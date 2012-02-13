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
  use gdata_mod
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  character (len=LargeString) :: Infiles, OfileBase

  logical :: DoCylVol, DoCold, DoWarm, DoVertVel, DoAbsVertVel, DoLwp

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles
  integer :: Nfiles
  integer, dimension(:), allocatable :: StmIx, StmIy

  real :: DeltaX, DeltaY, MinW, MaxW, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, LwpThresh
  real :: Xstart, Xinc, Ystart, Yinc
  real, dimension(:), allocatable :: MinP, Xcoords, Ycoords, Zcoords

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars
  type (GradsOutDescription) :: GoutDescrip

  ! Data arrays: cloud, tempc, precipr
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of cloud, tempc, precipr in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: W, Press, TempC, CintLiq, OutFilter
  type (GradsVarLocation) :: WLoc, PressLoc, TempcLoc, CintLiqLoc

  integer :: i
  integer :: Ierror
  integer :: ZcoordLoc

  integer :: ix, iy, iz, it

  logical :: SelectThisPoint, InsideCylVol

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBase, LargeString, DoCylVol, DoCold, DoWarm, DoVertVel, DoAbsVertVel, DoLwp, MinW, MaxW, LwpThresh, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Creating GRADS data filter:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
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
  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''
  flush(6)

  ! See if we have access to all of the required GRADS variables
  ! These checks also set Nx, Ny, Nz, Nt, Nvars for the following code section
  !
  ! Note that press, tempc, w are 3d variables and cint_liq is a 2d variable (only one
  ! z level). CheckDataDescripOneVar() will make sure that Nz gets set to the 3d variable's
  ! number of z levels even though the call using cint_liq is last. If cint_liq is the only
  ! variable being read in, then Nz will be set to one.
  if (DoCylVol) then
    ! Need pressure (to find the storm center)
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
  end if
  if (DoWarm .or. DoCold) then
    ! Need temperature (in deg. C)
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, TempcLoc, 'tempc')
  end if
  if (DoVertVel .or. DoAbsVertVel) then
    ! Need vertical velocity
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, WLoc, 'w')
  end if
  if (DoLwp) then
    ! Need column integrated liquid water
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CintLiqLoc, 'cint_liq')
  end if

  ! Calculate the x,y coordinates (in km) for doing selection by radius from storm center.
  ! Also for setting DeltaX and DeltaY (in m)
  ! Note that we want to access a 3d variable in GdataDescrip() if one exists in order
  ! to get the correct list of Zcoords. Look at what's in GdataDescrip() and point ZcoordLoc
  ! to an appropriate variable.
  allocate (Xcoords(1:Nx), Ycoords(1:Ny), Zcoords(1:Nz), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  ZcoordLoc = -1
  if (DoCylVol .or. DoWarm .or. DoCold .or. DoVertVel .or. DoAbsVertVel) then
     ! just grab the first one that exists
     if (DoCylVol) then
       ZcoordLoc = PressLoc%Fnum
     else if (DoWarm .or. DoCold) then
       ZcoordLoc = TempcLoc%Fnum
     else if (DoVertVel .or. DoAbsVertVel) then
       ZcoordLoc = WLoc%Fnum
     end if
  else
    ZcoordLoc = CintLiqLoc%Fnum
  end if

  call ConvertGridCoords(Nx, Ny, GdataDescrip(1), Xcoords, Ycoords)
  do iz = 1, Nz
    Zcoords(iz) = GdataDescrip(ZcoordLoc)%Zcoords(iz)
  end do

  DeltaX = (Xcoords(2) - Xcoords(1)) * 1000.0
  DeltaY = (Ycoords(2) - Ycoords(1)) * 1000.0

  write (*,*) 'Horizontal grid info:'
  write (*,*) '  X range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', GdataDescrip(1)%Xcoords(1), GdataDescrip(1)%Xcoords(Nx), Xcoords(1), Xcoords(Nx)
  write (*,*) '  Y range (min lat, max lat) --> (min y, max y): '
  write (*,*) '    ', GdataDescrip(1)%Ycoords(1), GdataDescrip(1)%Ycoords(Ny), Ycoords(1), Ycoords(Ny)
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
  write (*,*) '  Total number of grid variables:          ', Nvars
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''
  write (*,*) '  Grid delta x: ', DeltaX
  write (*,*) '  Grid delta y: ', DeltaY
  write (*,*) ''
  flush(6)

  write (*,*) 'Locations of variables in GRADS data (file number, var number):'
  if (DoCylVol) then
    write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'
  end if
  if (DoWarm .or. DoCold) then
    write (*,'(a20,i3,a2,i3,a1)') 'tempc: (', TempcLoc%Fnum, ', ', TempcLoc%Vnum, ')'
  end if
  if (DoVertVel .or. DoAbsVertVel) then
    write (*,'(a20,i3,a2,i3,a1)') 'w: (', WLoc%Fnum, ', ', WLoc%Vnum, ')'
  end if
  if (DoLwp) then
    write (*,'(a20,i3,a2,i3,a1)') 'cint_liq: (', CintLiqLoc%Fnum, ', ', CintLiqLoc%Vnum, ')'
  end if
  write (*,*) ''

  if (DoCylVol) then
    allocate (Press(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory allocation failed'
      stop
    end if
    call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
  end if
  if (DoWarm .or. DoCold) then
    allocate (TempC(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory allocation failed'
      stop
    end if
    call ReadGradsData(GdataDescrip, 'tempc', TempcLoc, TempC, Nx, Ny, Nz, Nt)
  end if
  if (DoVertVel .or. DoAbsVertVel) then
    allocate (W(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory allocation failed'
      stop
    end if
    call ReadGradsData(GdataDescrip, 'w', WLoc, W, Nx, Ny, Nz, Nt)
  end if
  if (DoLwp) then
    allocate (CintLiq(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory allocation failed'
      stop
    end if
    call ReadGradsData(GdataDescrip, 'cint_liq', CintLiqLoc, CintLiq, Nx, Ny, 1, Nt)
  end if
  write (*,*) ''
  flush(6)

  ! Allocate the output array and create the filter
  allocate (OutFilter(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Ouput data array memory allocation failed'
    stop
  end if

  if (DoCylVol) then
    ! Generate the storm center for all time steps
    allocate (StmIx(1:Nt), StmIy(1:Nt), MinP(1:Nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Ouput data array memory allocation failed'
      stop
    end if
    call RecordStormCenter(Nx, Ny, Nz, Nt, Press, StmIx, StmIy, MinP)
  
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
            SelectThisPoint = SelectThisPoint .and. (TempC(ix,iy,iz,it) .le. 0.0)
          end if
          if (DoWarm) then
            SelectThisPoint = SelectThisPoint .and. (TempC(ix,iy,iz,it) .gt. 0.0)
          end if
          if (DoVertVel) then
            SelectThisPoint = SelectThisPoint .and. ((W(ix,iy,iz,it) .ge. MinW) .and. (W(ix,iy,iz,it) .le. MaxW))
          end if
          if (DoAbsVertVel) then
            SelectThisPoint = SelectThisPoint .and. ((abs(W(ix,iy,iz,it)) .ge. MinW) .and. (abs(W(ix,iy,iz,it)) .le. MaxW))
          end if
          if (DoLwp) then
            SelectThisPoint = SelectThisPoint .and. (CintLiq(ix,iy,1,it) .ge. LwpThresh)
          end if

          if (SelectThisPoint) then
            OutFilter(ix,iy,iz,it) = 1.0
          else
            OutFilter(ix,iy,iz,it) = 0.0
          end if
        end do
      end do
    end do
  end do

  ! Write out the output filter
  Xstart = GdataDescrip(1)%Xcoords(1)
  Xinc = GdataDescrip(1)%Xcoords(2)-GdataDescrip(1)%Xcoords(1)
  Ystart = GdataDescrip(1)%Ycoords(1)
  Yinc = GdataDescrip(1)%Ycoords(2)-GdataDescrip(1)%Ycoords(1)

  call BuildGoutDescrip(Nx, Ny, Nz, Nt, OutFilter, OfileBase, GdataDescrip(1)%UndefVal, 'gfilter', &
          Xstart, Xinc, Ystart, Yinc, GdataDescrip(ZcoordLoc)%Zcoords, GdataDescrip(1)%Tstart, &
          GdataDescrip(1)%Tinc, GoutDescrip, 'gfilter')

  call WriteGrads(GoutDescrip, OutFilter)

  ! Clean up
  deallocate (Xcoords, Ycoords, Zcoords, OutFilter, stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory de-allocation failed'
    stop
  end if
  if (DoCylVol) then
    deallocate (Press, StmIx, StmIy, MinP, stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory de-allocation failed'
      stop
    end if
  end if
  if (DoWarm .or. DoCold) then
    deallocate (TempC, stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory de-allocation failed'
      stop
    end if
  end if
  if (DoVertVel .or. DoAbsVertVel) then
    deallocate (W, stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory de-allocation failed'
      stop
    end if
  end if
  if (DoLwp) then
    deallocate (CintLiq, stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory de-allocation failed'
      stop
    end if
  end if

  stop
end

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
