!***************************************************************
! Program to do averaging over the domain and yield a single
! data point at each time step
!
! This program will read in GRADS data from a RAMS simulation, 
! perform an averaging function, and output the single point time
! series in GRADS format
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
!   3. selection of averaging function
!
! Output
!   The output will be two files which make a GRADS data set. One
!   file is the control file and the other is the binary data.
!

program main
  use gdata_utils
  use azavg_utils
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128

  character (len=LargeString) :: Infiles, OfileBase

  logical :: DoCylVol, DoCold, DoWarm, DoVertVel, DoAbsVertVel, DoLwp

  type (GradsControlFiles) :: GctlFiles
  integer, dimension(:), allocatable :: StmIx, StmIy

  real :: DeltaX, DeltaY, MinW, MaxW, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, LwpThresh
  real :: Xstart, Xinc, Ystart, Yinc
  real, dimension(:), allocatable :: MinP, Xcoords, Ycoords, Zcoords

  ! Data arrays: cloud, tempc, precipr
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of cloud, tempc, precipr in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  type (GradsVar) :: W, Press, TempC, CintLiq, OutFilter

  integer :: i
  integer :: ZcoordLoc

  integer :: ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt
  real :: MinLon, MaxLon, MinLat, MaxLat

  logical :: SelectThisPoint

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBase, LargeString, DoCylVol, DoCold, DoWarm, DoVertVel, DoAbsVertVel, DoLwp, MinW, MaxW, LwpThresh, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ)
  call String2List(Infiles, ':', GctlFiles%Fnames, MaxFiles, GctlFiles%Nfiles, 'input files')

  write (*,*) 'Creating GRADS data filter:'
  write (*,*) '  GRADS input control files:'
  do i = 1, GctlFiles%Nfiles
    write (*,*) '  ', i, ': ', trim(GctlFiles%Fnames(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  Data selection specs: '
  if (DoCylVol) then
    write (*,*) '    Cylindrical volume:'
    write (*,*) '      Minimum radius: ', MinR
    write (*,*) '      Maximum radius: ', MaxR
    write (*,*) '      Minimum angle: ', MinPhi
    write (*,*) '      Maximum angle: ', MaxPhi
    write (*,*) '      Minimum height: ', MinZ
    write (*,*) '      Maximum height: ', MaxZ
  endif
  if (DoCold) then
    write (*,*) '    Cold (T <= 0 deg C)'
  endif
  if (DoWarm) then
    write (*,*) '    Warm (T > 0 deg C)'
  endif
  if (DoVertVel) then
    write (*,*) '    Vertical velocity:'
    write (*,*) '      Minimum W: ', MinW
    write (*,*) '      Maximum W: ', MaxW
  endif
  if (DoAbsVertVel) then
    write (*,*) '    Absolute value of vertical velocity:'
    write (*,*) '      Minimum W: ', MinW
    write (*,*) '      Maximum W: ', MaxW
  endif
  if (DoLwp) then
    write (*,*) '    LWP:'
    write (*,*) '      LWP (column integrated liquid water) threshold: ', LwpThresh
  endif
  write (*,*) ''
  flush(6)

  ! Read the GRADS data description files and collect the information about the data
  call ReadGradsCtlFiles(GctlFiles)

  ! See if we have access to all of the required GRADS variables
  !
  ! Note that press, tempc, w are 3d variables and cint_liq is a 2d variable (only one
  ! z level).
  if (DoCylVol) then
    ! Need pressure (to find the storm center)
    call InitGvarFromGdescrip(GctlFiles, Press, 'press')
  endif
  if (DoWarm .or. DoCold) then
    ! Need temperature (in deg. C)
    call InitGvarFromGdescrip(GctlFiles, TempC, 'tempc')
  endif
  if (DoVertVel .or. DoAbsVertVel) then
    ! Need vertical velocity
    call InitGvarFromGdescrip(GctlFiles, W, 'w')
  endif
  if (DoLwp) then
    ! Need column integrated liquid water
    call InitGvarFromGdescrip(GctlFiles, CintLiq, 'cint_liq')
  endif

  ! If we got to here all the required variables are available. Check to see if they have consistent
  ! dimensions.
  ! Use 3D values for Nx, Ny, Nz, Nt if any 3D variable is being used. If column integreated
  ! liquid is the only variable being used, then set Nx, Ny, Nz, Nt to the 2D values (from col. int. liq.)
  if (DoCylVol) then
    ! check consistency of dimensions of the variables
    if (DoWarm .or. DoCold) then
      ! using Tempc, check that it is consistent with Press
      if (.not. (GvarDimsMatch(Press, TempC, .false.))) then
        write (*,*) 'ERROR: dimensions of press and tempc do not match'
        stop
      endif
    endif
    if (DoVertVel .or. DoAbsVertVel) then
      ! using W, check that it is consistent with Press
      if (.not. (GvarDimsMatch(Press, W, .false.))) then
        write (*,*) 'ERROR: dimensions of press and w do not match'
        stop
      endif
    endif
    if (DoLwp) then
      ! using CintLiq, check that it is consistent with Press
      if (.not. (GvarDimsMatch(Press, CintLiq, .true.))) then
        write (*,*) 'ERROR: dimensions of press and cint_liq do not match'
        stop
      endif
    endif

    ! Variables are consistent, record the dimensions and coordinates for later
    call SetOutDimsCoords(Press, Nx, Ny, Nz, Nt, Xcoords, Ycoords, Zcoords, MinLon, MaxLon, MinLat, MaxLat)
  else if (DoWarm .or. DoCold) then
    ! check consistency of dimensions of the variables
    if (DoVertVel .or. DoAbsVertVel) then
      ! using W, check that it is consistent with TempC
      if (.not. (GvarDimsMatch(TempC, W, .false.))) then
        write (*,*) 'ERROR: dimensions of tempc and w do not match'
        stop
      endif
    endif
    if (DoLwp) then
      ! using CintLiq, check that it is consistent with TempC
      if (.not. (GvarDimsMatch(TempC, CintLiq, .true.))) then
        write (*,*) 'ERROR: dimensions of tempc and cint_liq do not match'
        stop
      endif
    endif

    ! Variables are consistent, record the dimensions and coordinates for later
    call SetOutDimsCoords(TempC, Nx, Ny, Nz, Nt, Xcoords, Ycoords, Zcoords, MinLon, MaxLon, MinLat, MaxLat)
  else if (DoVertVel .or. DoAbsVertVel) then
    if (DoLwp) then
      ! using CintLiq, check that it is consistent with W
      if (.not. (GvarDimsMatch(W, CintLiq, .true.))) then
        write (*,*) 'ERROR: dimensions of w and cint_liq do not match'
        stop
      endif
    endif

    ! Variables are consistent, record the dimensions and coordinates for later
    call SetOutDimsCoords(W, Nx, Ny, Nz, Nt, Xcoords, Ycoords, Zcoords, MinLon, MaxLon, MinLat, MaxLat)
  else
    ! only using column integrated liquid, record the dimensions and coordinates for later
    call SetOutDimsCoords(CintLiq, Nx, Ny, Nz, Nt, Xcoords, Ycoords, Zcoords, MinLon, MaxLon, MinLat, MaxLat)
  endif

  DeltaX = (Xcoords(2) - Xcoords(1)) * 1000.0
  DeltaY = (Ycoords(2) - Ycoords(1)) * 1000.0

  write (*,*) 'Horizontal grid info:'
  write (*,*) '  X range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', MinLon, MaxLon, Xcoords(1), Xcoords(Nx)
  write (*,*) '  Y range (min lat, max lat) --> (min y, max y): '
  write (*,*) '    ', MinLat, MaxLat, Ycoords(1), Ycoords(Ny)
  write (*,*) ''
  write (*,*) 'Vertical grid info:'
  do iz = 1, Nz
    write (*,*) '  ', iz, ' --> ', Zcoords(iz)
  end do
  write (*,*) ''
  flush(6)

  ! Read in the GRADS variable data
  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', Nx
  write (*,*) '  Number of y (latitude) points:           ', Ny
  write (*,*) '  Number of z (vertical level) points:     ', Nz
  write (*,*) '  Number of t (time) points:               ', Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''
  write (*,*) '  Grid delta x: ', DeltaX
  write (*,*) '  Grid delta y: ', DeltaY
  write (*,*) ''
  flush(6)

  write (*,*) 'Locations of variables in GRADS data (file number, var number):'
  if (DoCylVol) then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
  end if
  if (DoWarm .or. DoCold) then
    write (*,'(a20,a,a2,i3,a1)') 'tempc: (', trim(TempC%DataFile), ', ', TempC%Vnum, ')'
  end if
  if (DoVertVel .or. DoAbsVertVel) then
    write (*,'(a20,a,a2,i3,a1)') 'w: (', trim(W%DataFile), ', ', W%Vnum, ')'
  end if
  if (DoLwp) then
    write (*,'(a20,a,a2,i3,a1)') 'cint_liq: (', trim(CintLiq%DataFile), ', ', CintLiq%Vnum, ')'
  end if
  write (*,*) ''

  if (DoCylVol) then
    call ReadGradsData(Press)
  end if
  if (DoWarm .or. DoCold) then
    call ReadGradsData(TempC)
  end if
  if (DoVertVel .or. DoAbsVertVel) then
    call ReadGradsData(W)
  end if
  if (DoLwp) then
    call ReadGradsData(CintLiq)
  end if
  write (*,*) ''
  flush(6)

  ! Allocate the output array and create the filter
  if (DoCylVol) then
    call InitGvarFromGvar(Press, OutFilter, 'gfilter')
  else if (DoWarm .or. DoCold) then
    call InitGvarFromGvar(TempC, OutFilter, 'gfilter')
  else if (DoVertVel .or. DoAbsVertVel) then
    call InitGvarFromGvar(W, OutFilter, 'gfilter')
  else if (DoLwp) then
    call InitGvarFromGvar(CintLiq, OutFilter, 'gfilter')
  endif

  if (DoCylVol) then
    ! Generate the storm center for all time steps
    allocate (StmIx(1:Nt), StmIy(1:Nt), MinP(1:Nt))
    call RecordStormCenter(Press, StmIx, StmIy, MinP)
  
    do it = 1, Nt
      write (*,*) 'Timestep: ', it
      write (*,'(a,i3,a,i3,a,g,a,g,a)') '  Storm Center: (', StmIx(it), ', ', StmIy(it), &
         ') --> (', Xcoords(StmIx(it)), ', ', Ycoords(StmIy(it)), ')'
      write (*,*) '  Minumum pressure: ', MinP(it)
    end do
    write (*,*) ''
    flush(6)
  end if

  ! Apply the filter. We are forcing the user to specify at least one filtering mechamism. If
  ! multiple filters are specified, then take the intersection of all filters.
  do it = 1, Nt
    do iz = 1, Nz
      do ix = 1, Nx
        do iy = 1, Ny
          SelectThisPoint = .true.
          if (DoCylVol) then
            SelectThisPoint = SelectThisPoint .and. InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, &
                  MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords) 
          end if
          if (DoCold) then
            SelectThisPoint = SelectThisPoint .and. (TempC%Vdata(it,iz,ix,iy) .le. 0.0)
          end if
          if (DoWarm) then
            SelectThisPoint = SelectThisPoint .and. (TempC%Vdata(it,iz,ix,iy) .gt. 0.0)
          end if
          if (DoVertVel) then
            SelectThisPoint = SelectThisPoint .and. ((W%Vdata(it,iz,ix,iy) .ge. MinW) .and. (W%Vdata(it,iz,ix,iy) .le. MaxW))
          end if
          if (DoAbsVertVel) then
            SelectThisPoint = SelectThisPoint .and. ((abs(W%Vdata(it,iz,ix,iy)) .ge. MinW) .and. (abs(W%Vdata(it,iz,ix,iy)) .le. MaxW))
          end if
          if (DoLwp) then
            SelectThisPoint = SelectThisPoint .and. (CintLiq%Vdata(ix,iy,1,it) .ge. LwpThresh)
          end if

          if (SelectThisPoint) then
            OutFilter%Vdata(it,iz,ix,iy) = 1.0
          else
            OutFilter%Vdata(it,iz,ix,iy) = 0.0
          end if
        end do
      end do
    end do
  end do

  ! Write out the output filter
  call WriteGrads(OutFilter, OfileBase, 'gfilter')

  stop

contains
!**********************************************************************
! Subroutines go here. Want these contained within the main program
! so that the interfaces to these subroutine are 'explicit' (ie, like
! ANSI C external declarations). This is needed to get the passing
! of allocatable arrays working properly.
!**********************************************************************

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   Outfiles - output GRADS file, base name for two files
!   AvgFunc - averaging function selection
!

subroutine GetMyArgs(Infiles, OfileBase, StringSize, DoCylVol, DoCold, DoWarm, DoVertVel, DoAbsVertVel, DoLwp, MinW, MaxW, LwpThresh, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ)
  implicit none

  character (len=*) :: Infiles, OfileBase
  integer :: StringSize
  real :: MinW, MaxW, LwpThresh, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  logical :: DoCylVol, DoCold, DoWarm, DoVertVel, DoAbsVertVel, DoLwp

  integer :: iargc, i, Nargs
  character (len=StringSize) :: Arg, FilterFunc

  logical :: BadArgs

  ! set defaults
  DoCylVol     = .false.
  DoCold       = .false.
  DoWarm       = .false.
  DoVertVel    = .false.
  DoAbsVertVel = .false.
  DoLwp        = .false.
  Infiles   = ''
  OfileBase = ''
  MinW      = -999.0
  MaxW      = -999.0
  LwpThresh = -1.0
  MinR      = -1.0
  MaxR      = -1.0
  MinPhi    = -1.0
  MaxPhi    = -1.0
  MinZ      = -1.0
  MaxZ      = -1.0

  ! walk through the list of arguments
  !   see the USAGE string at the end of this routine
  !
  !   first arg --> Infiles
  !   second arg --> OfileBase
  !
  !   remaining arguments are the filter specs
  BadArgs      = .false.
  i = 1
  Nargs = iargc()
  do while (i .le. Nargs)
    call getarg(i, Arg)

    if (i .eq. 1) then
      Infiles = Arg
      i = i + 1
    else if (i .eq. 2) then
      OfileBase = Arg
      i = i + 1
    else
      FilterFunc = Arg
      i = i + 1

      if (FilterFunc .eq. 'cylvol') then
        DoCylVol = .true.

        if ((i+5) .gt. Nargs) then
          write (*,*) 'ERROR: not enough arguments (6) left to fully specify the clyvol <filter_spec>'
          BadArgs = .true.
        else
          ! have enough args
          call getarg(i, Arg)
          read(Arg, '(f)') MinR
          call getarg(i+1, Arg)
          read(Arg, '(f)') MaxR
          call getarg(i+2, Arg)
          read(Arg, '(f)') MinPhi
          call getarg(i+3, Arg)
          read(Arg, '(f)') MaxPhi
          call getarg(i+4, Arg)
          read(Arg, '(f)') MinZ
          call getarg(i+5, Arg)
          read(Arg, '(f)') MaxZ
          i = i + 6
    
          if ((MinR .lt. 0.0) .or. (MaxR .lt. 0.0) .or. (MaxR .le. MinR)) then
            write (*,*) 'ERROR: <min_r> and <max_r> must be >= 0.0, and <max_r> must be > <min_r>'
            write (*,*) ''
            BadArgs = .true.
          end if
          if ((MinPhi .lt. 0.0) .or. (MaxPhi .lt. 0.0) .or. (MaxPhi .le. MinPhi)) then
            write (*,*) 'ERROR: <min_phi> and <max_phi> must be >= 0.0, and <max_phi> must be > <min_phi>'
            write (*,*) ''
            BadArgs = .true.
          end if
          if ((MinZ .lt. 0.0) .or. (MaxZ .lt. 0.0) .or. (MaxZ .le. MinZ)) then
            write (*,*) 'ERROR: <min_z> and <max_z> must be >= 0.0, and <max_z> must be > <min_z>'
            write (*,*) ''
            BadArgs = .true.
          end if
        end if
      else if (FilterFunc .eq. 'cold') then
        DoCold = .true.
      else if (FilterFunc .eq. 'warm') then
        DoWarm = .true.
      else if (FilterFunc .eq. 'vertvel') then
        DoVertVel = .true.
  
        if ((i+1) .gt. Nargs) then
          write (*,*) 'ERROR: not enough arguments (2) left to fully specify the vertvel <filter_spec>'
          BadArgs = .true.
        else
          call getarg(i, Arg)
          read(Arg, '(f)') MinW
          call getarg(i+1, Arg)
          read(Arg, '(f)') MaxW
          i = i + 2
    
          if (MaxW .le. MinW) then
            write (*,*) 'ERROR: <max_w> must be > <min_w>'
            write (*,*) ''
            BadArgs = .true.
          end if
        end if
      else if (FilterFunc .eq. 'absvertvel') then
        DoAbsVertVel = .true.
  
        if ((i+1) .gt. Nargs) then
          write (*,*) 'ERROR: not enough arguments (2) left to fully specify the absvertvel <filter_spec>'
          BadArgs = .true.
        else
          call getarg(i, Arg)
          read(Arg, '(f)') MinW
          call getarg(i+1, Arg)
          read(Arg, '(f)') MaxW
          i = i + 2
    
          if (MaxW .le. MinW) then
            write (*,*) 'ERROR: <max_w> must be > <min_w>'
            write (*,*) ''
            BadArgs = .true.
          end if
        end if
      else if (FilterFunc .eq. 'lwp') then
        DoLwp = .true.
  
        if (i .gt. Nargs) then
          write (*,*) 'ERROR: not enough arguments (1) left to fully specify the lwp <filter_spec>'
          BadArgs = .true.
        else
          call getarg(i, Arg)
          read(Arg, '(f)') LwpThresh
          i = i + 1
  
          if (LwpThresh .lt. 0.0) then
            write (*,*) 'ERROR: <lpw_> must be >= 0.0'
            write (*,*) ''
            BadArgs = .true.
          end if
        end if
      else
        write (*,*) 'ERROR: <filter>, ', trim(FilterFunc), ', must be one of:'
        write (*,*) '          cylvol'
        write (*,*) '          cold'
        write (*,*) '          warm'
        write (*,*) '          vertvel'
        write (*,*) '          absvertvel'
        write (*,*) '          lwp'
        write (*,*) ''
        BadArgs = .true.
      endif
    endif
  enddo

  if ((.not.DoCylVol) .and. (.not.DoCold) .and. (.not.DoWarm) &
       .and. (.not.DoVertVel) .and. (.not.DoAbsVertVel) .and. (.not.DoLwp)) then
    write (*,*) 'ERROR: must specify at least one <filter_spec>'
    write (*,*) ''
    BadArgs = .true.
  endif

  if (DoCold .and. DoWarm) then
    write (*,*) 'ERROR: cannot specify both cold and warm <filter_spec>'
    write (*,*) ''
    BadArgs = .true.
  end if

  if (DoVertVel .and. DoAbsVertVel) then
    write (*,*) 'ERROR: cannot specify both vertvel and absvertvel <filter_spec>'
    write (*,*) ''
    BadArgs = .true.
  end if

  if (BadArgs) then
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_files> <filter_spec> [<filter_spec>...]'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, colon seprated list, data + rates names, this program will tag on .ctl, .dat suffixes'
    write (*,*) '        <filter_spec>: <filter> <filter_args>'
    write (*,*) ''
    write (*,*) '        <filter> <filter_args>:'
    write (*,*) '            cylvol <min_r> <max_r> <min_phi> <max_phi> <min_z> <max_z>'
    write (*,*) '                Select data within a cylindrical (r,phi,z) volume centered on the storm center'
    write (*,*) '                  <min_r> <max_r> (in km)'
    write (*,*) '                  <max_phi> <max_phi> (in radians)'
    write (*,*) '                  <min_z> <max_z> (in m)'
    write (*,*) '            cold'
    write (*,*) '                Select data within regions where T <= 0 deg C'
    write (*,*) '            warm'
    write (*,*) '                Select data within regions where T > 0 deg C'
    write (*,*) '            vertvel <min_w> <max_w>'
    write (*,*) '                Select data within regions where <min_w> <= w <= <max_w>'
    write (*,*) '                  <min_w> and <max_w> in m/s'
    write (*,*) '            absvertvel <min_w> <max_w>'
    write (*,*) '                Select data within regions where <min_w> <= |w| <= <max_w>'
    write (*,*) '            lwp <lwp_thresh>'
    write (*,*) '                Select data within regions where column integrated liquid (lwp) >= <lwp_thresh>'
    write (*,*) ''
    write (*,*) '        Note that more than one <filter_spec> can be specified. When this is done the data selection'
    write (*,*) '        becomes the intersection of all the filter specs. Eg., specifying both cylvol and cold'
    write (*,*) '        will select only the data that is both inside the cylindrical volume and in regions <= 0 deg C'
    write (*,*) ''

    stop
  end if

  return
end subroutine


!**********************************************************************
! SetOutDimsCoords()
!
! This routine will set the coordinate and dimension data 

subroutine SetOutDimsCoords(Gvar, Nx, Ny, Nz, Nt, Xcoords, Ycoords, Zcoords, MinLon, MaxLon, MinLat, MaxLat)
  use gdata_utils
  use azavg_utils
  implicit none

  type (GradsVar) :: Gvar
  integer :: Nx, Ny, Nz, Nt
  real, dimension(:), allocatable :: Xcoords, Ycoords, Zcoords
  real ::  MinLon, MaxLon, MinLat, MaxLat

  integer :: iz

  Nx = Gvar%Nx
  Ny = Gvar%Ny
  Nz = Gvar%Nz
  Nt = Gvar%Nt

  MinLon = Gvar%Xcoords(1)
  MaxLon = Gvar%Xcoords(Gvar%Nx)
  MinLat = Gvar%Ycoords(1)
  MaxLat = Gvar%Ycoords(Gvar%Ny)

  ! Calculate the x,y coordinates (in km) for doing selection by radius from storm center.
  call ConvertGridCoords(Gvar, Xcoords, Ycoords)
  allocate (Zcoords(1:Nz))
  do iz = 1, Nz
      Zcoords(iz) = Gvar%Zcoords(iz)
  enddo

  return
end subroutine

end program main
