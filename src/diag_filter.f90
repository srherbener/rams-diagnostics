!***************************************************************
! Program to apply filters to input data and output "mask" data
! that shows where points exist that were selected by the filtering.
!
! Args
!   1. directory containing input files
!   2. suffix to tag onto the end of input file names
!   3. output file name 
!   4. selection of averaging function
!
! Output
!   The output will be an hdf5 file containing '1's where data was
!   selected and '0's else where. The output will always be a 3D
!   field so that it can be applied downstream to either 2D or 3D fields.
!

program diag_filter
  use rhdf5_utils
  use diag_utils
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128

  character (len=LargeString) :: InDir, InSuffix, OutFile

  logical :: DoCylVol, DoCold, DoWarm, DoVertVel, DoAbsVertVel, DoCwp

  integer, dimension(:), allocatable :: StmIx, StmIy

  real :: DeltaX, DeltaY, MinW, MaxW, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CwpThresh
  real :: Xstart, Xinc, Ystart, Yinc
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm

  ! Data arrays: cloud, tempc, precipr
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of cloud, tempc, precipr in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  type (Rhdf5Var) :: W, Press, TempC, Cwp, OutFilter, Xcoords, Ycoords, Zcoords, Tcoords
  type (Rhdf5Var) :: MinP, StormX, StormY
  character (len=MediumString) :: Wfile, PressFile, TempFile, CwpFile

  integer :: i
  integer :: ZcoordLoc

  integer :: ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt

  logical :: SelectThisPoint

  ! Get the command line arguments
  call GetMyArgs(InDir, InSuffix, OutFile, LargeString, DoCylVol, DoCold, DoWarm, DoVertVel, DoAbsVertVel, DoCwp, MinW, MaxW, CwpThresh, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ)

  write (*,*) 'Creating HDF5 data filter:'
  write (*,*) '  Input directory: ', trim(InDir)
  write (*,*) '  Input file name suffix: ', trim(InSuffix)
  write (*,*) '  Output file name:  ', trim(OutFile)
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
  if (DoCwp) then
    write (*,*) '    CWP:'
    write (*,*) '      CWP (column integrated liquid water) threshold: ', CwpThresh
  endif
  write (*,*) ''
  flush(6)

  ! See if we have access to all of the required variables

  ! Always need vertical velocity, even if not doing W filtering, since this
  ! will serve as the model of the 3D field.
  Wfile = trim(InDir) // '/w' // trim(InSuffix)
  W%vname = 'w'
  call rhdf5_read_init(Wfile, W)

  if (DoCylVol) then
    ! Need pressure (to find the storm center)
    PressFile = trim(InDir) // '/sea_press' // trim(InSuffix)
    Press%vname = 'sea_press'
    call rhdf5_read_init(PressFile, Press)
    if (.not. (DimsMatch(W, Press))) then
      write (*,*) 'ERROR: horizontal dimensions of w and sea_press do not match'
      stop
    endif
  endif
  if (DoWarm .or. DoCold) then
    ! Need temperature (in deg. C)
    TempFile = trim(InDir) // '/tempc' // trim(InSuffix)
    TempC%vname = 'tempc'
    call rhdf5_read_init(TempFile, TempC)
    if (.not. (DimsMatch(W, TempC))) then
      write (*,*) 'ERROR: horizontal dimensions of w and tempc do not match'
      stop
    endif
  endif
  if (DoCwp) then
    ! Need column integrated liquid water
    CwpFile = trim(InDir) // '/vint_cloud' // trim(InSuffix)
    Cwp%vname = 'vertint_cloud'
    call rhdf5_read_init(CwpFile, Cwp)
    if (.not. (DimsMatch(W, Cwp))) then
      write (*,*) 'ERROR: horizontal dimensions of w and vint_cloud do not match'
      stop
    endif
  endif

  ! Set the output dimensions and coordinates to those of the W field
  call SetOutDimsCoords(Wfile, W, Nx, Ny, Nz, Nt, Xcoords, Ycoords, Zcoords, Tcoords)

  ! Convert lat (x coords) and lon (y coords) to distances in km
  allocate(XcoordsKm(Nx))
  allocate(YcoordsKm(Ny))
  call ConvertGridCoords(Nx, Ny, Nz, Nt, Xcoords%vdata, Ycoords%vdata, XcoordsKm, YcoordsKm)

  DeltaX = (XcoordsKm(2) - XcoordsKm(1)) * 1000.0
  DeltaY = (YcoordsKm(2) - YcoordsKm(1)) * 1000.0

  write (*,*) 'Horizontal grid info:'
  write (*,*) '  X range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', Xcoords%vdata(1), Xcoords%vdata(Nx), XcoordsKm(1), XcoordsKm(Nx)
  write (*,*) '  Y range (min lat, max lat) --> (min y, max y): '
  write (*,*) '    ', Ycoords%vdata(1), Ycoords%vdata(Ny), YcoordsKm(1), YcoordsKm(Ny)
  write (*,*) ''
  write (*,*) '  Delta x of domain:     ', DeltaX
  write (*,*) '  Delta y of domain:     ', DeltaY
  write (*,*) ''
  write (*,*) 'Vertical grid info:'
  do iz = 1, Nz
    write (*,*) '  ', iz, ' --> ', Zcoords%vdata(iz)
  end do
  write (*,*) ''
  flush(6)

  ! Read in the variable data
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

  if (DoCylVol) then
    write (*,*) 'Reading variable: sea_press'
    write (*,*) '  HDF5 file: ', trim(PressFile)
    write (*,*) ''
    call rhdf5_read(PressFile, Press)
  end if
  if (DoWarm .or. DoCold) then
    write (*,*) 'Reading variable: tempc'
    write (*,*) '  HDF5 file: ', trim(TempFile)
    write (*,*) ''
    call rhdf5_read(TempFile, TempC)
  end if
  if (DoVertVel .or. DoAbsVertVel) then
    write (*,*) 'Reading variable: w'
    write (*,*) '  HDF5 file: ', trim(Wfile)
    write (*,*) ''
    call rhdf5_read(Wfile, W)
  end if
  if (DoCwp) then
    write (*,*) 'Reading variable: vint_cloud'
    write (*,*) '  HDF5 file: ', trim(CwpFile)
    write (*,*) ''
    call rhdf5_read(CwpFile, Cwp)
  end if
  write (*,*) ''
  flush(6)

  ! Allocate the output array and create the filter
  OutFilter%vname = 'filter'
  OutFilter%ndims = 4
  OutFilter%dims(1) = Nx
  OutFilter%dims(2) = Ny
  OutFilter%dims(3) = Nz
  OutFilter%dims(4) = Nt
  OutFilter%dimnames(1) = 'x'
  OutFilter%dimnames(2) = 'y'
  OutFilter%dimnames(3) = 'z'
  OutFilter%dimnames(4) = 't'
  OutFilter%units = 'logical'
  OutFilter%descrip = 'data selection locations'
  allocate(OutFilter%vdata(Nx*Ny*Nz*Nt))

  if (DoCylVol) then
    ! Generate the storm center for all time steps
    MinP%vname = 'min_press'
    MinP%ndims = 1
    MinP%dims(1) = Nt
    MinP%dimnames(1) = 't'
    MinP%units = 'mb'
    MinP%descrip = 'minimum SLP of storm'
    allocate(MinP%vdata(Nt))
    allocate (StmIx(1:Nt), StmIy(1:Nt))
    call RecordStormCenter(Nx, Ny, Nz, Nt, Press%vdata, StmIx, StmIy, MinP%vdata)

    StormX%vname = 'min_press_xloc'
    StormX%ndims = 1
    StormX%dims(1) = Nt
    StormX%dimnames(1) = 't'
    StormX%units = 'deg lon'
    StormX%descrip = 'longitude location of minimum SLP of storm'
    allocate(StormX%vdata(Nt))
    
    StormY%vname = 'min_press_yloc'
    StormY%ndims = 1
    StormY%dims(1) = Nt
    StormY%dimnames(1) = 't'
    StormY%units = 'deg lat'
    StormY%descrip = 'latitude location of minimum SLP of storm'
    allocate(StormY%vdata(Nt))
  
    do it = 1, Nt
      write (*,*) 'Timestep: ', it
      write (*,'(a,i3,a,i3,a,g,a,g,a)') '  Storm Center: (', StmIx(it), ', ', StmIy(it), &
         ') --> (', XcoordsKm(StmIx(it)), ', ', YcoordsKm(StmIy(it)), ')'
      write (*,*) '  Minumum pressure: ', MinP%vdata(it)

      StormX%vdata(it) = Xcoords%vdata(StmIx(it))
      StormY%vdata(it) = Ycoords%vdata(StmIy(it))
    end do
    write (*,*) ''
    flush(6)
  end if

  ! Apply the filter. We are forcing the user to specify at least one filtering mechamism. If
  ! multiple filters are specified, then take the intersection of all filters.
  do it = 1, Nt
    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          SelectThisPoint = .true.
          if (DoCylVol) then
            SelectThisPoint = SelectThisPoint .and. InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, &
                  MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, XcoordsKm, YcoordsKm, Zcoords%vdata) 
          end if
          if (DoCold) then
            SelectThisPoint = SelectThisPoint .and. &
              (MultiDimLookup(Nx, Ny, Nz, Nt, ix, iy, iz, it, Var3d=TempC%vdata) .le. 0.0)
          end if
          if (DoWarm) then
            SelectThisPoint = SelectThisPoint .and. &
              (MultiDimLookup(Nx, Ny, Nz, Nt, ix, iy, iz, it, Var3d=TempC%vdata) .gt. 0.0)
          end if
          if (DoVertVel) then
            SelectThisPoint = SelectThisPoint .and. &
              ((MultiDimLookup(Nx, Ny, Nz, Nt, ix, iy, iz, it, Var3d=W%vdata) .ge. MinW) .and. &
               (MultiDimLookup(Nx, Ny, Nz, Nt, ix, iy, iz, it, Var3d=W%vdata) .le. MaxW))
          end if
          if (DoAbsVertVel) then
            SelectThisPoint = SelectThisPoint .and. &
              ((abs(MultiDimLookup(Nx, Ny, Nz, Nt, ix, iy, iz, it, Var3d=W%vdata)) .ge. MinW) .and. &
               (abs(MultiDimLookup(Nx, Ny, Nz, Nt, ix, iy, iz, it, Var3d=W%vdata)) .le. MaxW))
          end if
          if (DoCwp) then
            SelectThisPoint = SelectThisPoint .and. &
              (MultiDimLookup(Nx, Ny, Nz, Nt, ix, iy, iz, it, Var2d=Cwp%vdata) .gt. CwpThresh)
          end if

          if (SelectThisPoint) then
            call MultiDimAssign(Nx, Ny, Nz, Nt, ix, iy, iz, it, 1.0, Var3d=OutFilter%vdata)
          else
            call MultiDimAssign(Nx, Ny, Nz, Nt, ix, iy, iz, it, 0.0, Var3d=OutFilter%vdata)
          end if
        end do
      end do
    end do
  end do

  ! write out the filter data
  ! third arg to rhdf5_write is "append" flag:
  !   0 - create new file
  !   1 - append to existing file
  write (*,*) 'Writing HDF5 output: ', trim(OutFile)
  write (*,*) ''
  call rhdf5_write(OutFile, OutFilter, 0)

  ! write out the location of the minimum pressure and the minimum pressure values
  if (DoCylVol) then
    call rhdf5_write(OutFile, MinP, 1)
    call rhdf5_write(OutFile, StormX, 1)
    call rhdf5_write(OutFile, StormY, 1)
  endif

  ! write out the coordinate data
  call rhdf5_write(OutFile, Xcoords, 1)
  call rhdf5_write(OutFile, Ycoords, 1)
  call rhdf5_write(OutFile, Zcoords, 1)
  call rhdf5_write(OutFile, Tcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(OutFile, Xcoords, 'x')
  call rhdf5_set_dimension(OutFile, Ycoords, 'y')
  call rhdf5_set_dimension(OutFile, Zcoords, 'z')
  call rhdf5_set_dimension(OutFile, Tcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(OutFile, OutFilter)
  if (DoCylVol) then
    call rhdf5_attach_dimensions(OutFile, MinP)
    call rhdf5_attach_dimensions(OutFile, StormX)
    call rhdf5_attach_dimensions(OutFile, StormY)
  endif

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

subroutine GetMyArgs(InDir, InSuffix, OutFile, StringSize, DoCylVol, DoCold, DoWarm, DoVertVel, DoAbsVertVel, DoCwp, MinW, MaxW, CwpThresh, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ)
  implicit none

  character (len=*) :: InDir, InSuffix, OutFile
  integer :: StringSize
  real :: MinW, MaxW, CwpThresh, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  logical :: DoCylVol, DoCold, DoWarm, DoVertVel, DoAbsVertVel, DoCwp

  integer :: iargc, i, Nargs
  character (len=StringSize) :: Arg, FilterFunc

  logical :: BadArgs

  ! set defaults
  DoCylVol     = .false.
  DoCold       = .false.
  DoWarm       = .false.
  DoVertVel    = .false.
  DoAbsVertVel = .false.
  DoCwp        = .false.
  InDir   = ''
  InSuffix = ''
  OutFile = ''
  MinW      = -999.0
  MaxW      = -999.0
  CwpThresh = -1.0
  MinR      = -1.0
  MaxR      = -1.0
  MinPhi    = -1.0
  MaxPhi    = -1.0
  MinZ      = -1.0
  MaxZ      = -1.0

  ! walk through the list of arguments
  !   see the USAGE string at the end of this routine
  !
  !   first arg --> InDir
  !   second arg --> InSuffix
  !   third arg --> OutFile
  !
  !   remaining arguments are the filter specs
  BadArgs      = .false.
  i = 1
  Nargs = iargc()
  do while (i .le. Nargs)
    call getarg(i, Arg)

    if (i .eq. 1) then
      InDir = Arg
      i = i + 1
    else if (i .eq. 2) then
      InSuffix = Arg
      i = i + 1
    else if (i .eq. 3) then
      OutFile = Arg
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
      else if (FilterFunc .eq. 'cwp') then
        DoCwp = .true.
  
        if (i .gt. Nargs) then
          write (*,*) 'ERROR: not enough arguments (1) left to fully specify the cwp <filter_spec>'
          BadArgs = .true.
        else
          call getarg(i, Arg)
          read(Arg, '(f)') CwpThresh
          i = i + 1
  
          if (CwpThresh .lt. 0.0) then
            write (*,*) 'ERROR: <cpw_thresh> must be >= 0.0'
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
        write (*,*) '          cwp'
        write (*,*) ''
        BadArgs = .true.
      endif
    endif
  enddo

  if ((.not.DoCylVol) .and. (.not.DoCold) .and. (.not.DoWarm) &
       .and. (.not.DoVertVel) .and. (.not.DoAbsVertVel) .and. (.not.DoCwp)) then
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
    write (*,*) 'USAGE: diag_filter <in_dir> <in_suffix> <out_data_files> <filter_spec> [<filter_spec>...]'
    write (*,*) '        <in_dir>: directory containing input files'
    write (*,*) '        <in_suffix>: suffix to tag onto the end of input file names'
    write (*,*) '           Note: input file names become: <in_dir>/<var_name><in_suffix>'
    write (*,*) '        <out_file>: output hdf5 file name'
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
    write (*,*) '            cwp <cwp_thresh>'
    write (*,*) '                Select data within regions where column integrated cloud water (cwp) >= <cwp_thresh>'
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

subroutine SetOutDimsCoords(Hfile, Hvar, Nx, Ny, Nz, Nt, Xcoords, Ycoords, Zcoords, Tcoords)
  use rhdf5_utils
  use diag_utils
  implicit none

  character (len=*) :: Hfile
  type (Rhdf5Var) :: Hvar, Xcoords, Ycoords, Zcoords, Tcoords
  integer :: Nx, Ny, Nz, Nt

  if (Hvar%ndims .eq. 4) then
    Nx = Hvar%dims(1)
    Ny = Hvar%dims(2)
    Nz = Hvar%dims(3)
    Nt = Hvar%dims(4)
  else
    Nx = Hvar%dims(1)
    Ny = Hvar%dims(2)
    Nz = 1
    Nt = Hvar%dims(3)
  endif

  ! Read in longitude, latitude and height values
  Xcoords%vname = 'x_coords'
  call rhdf5_read_init(Hfile, Xcoords)
  call rhdf5_read(Hfile, Xcoords)
  
  Ycoords%vname = 'y_coords'
  call rhdf5_read_init(Hfile, Ycoords)
  call rhdf5_read(Hfile, Ycoords)

  Zcoords%vname = 'z_coords'
  call rhdf5_read_init(Hfile, Zcoords)
  call rhdf5_read(Hfile, Zcoords)

  Tcoords%vname = 't_coords'
  call rhdf5_read_init(Hfile, Tcoords)
  call rhdf5_read(Hfile, Tcoords)

  return
end subroutine

end program diag_filter
