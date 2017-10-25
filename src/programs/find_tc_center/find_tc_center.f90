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

program find_tc_center
  use rhdf5_utils
  use diag_utils
  implicit none

  integer, parameter :: LargeString  = 512
  integer, parameter :: MediumString = 256
  integer, parameter :: LittleString = 128

  character (len=LargeString) :: PressFile, PressVar, TopoFile, TopoVar, OutFile

  type (Rhdf5Var) :: Press, Topo, Xcoords, Ycoords, Zcoords, Tcoords
  type (Rhdf5Var) :: Radius, Phi, MinPress, MinPressX, MinPressY, MinPressXidx, MinPressYidx
  type (Rhdf5Var) :: PressCentX, PressCentY, PressCentXidx, PressCentYidx
  type (Rhdf5Var) :: SpeedX, SpeedY
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm, StormX, StormY
  real :: DeltaX, DeltaY, MaxElev
  real :: MinLat, MinLon, MaxLat, MaxLon
  real :: PressRad
  integer :: MinPressIx, MinPressIy
  integer :: PressCentIx, PressCentIy

  integer :: Nx, Ny, Nz, Nt
  integer :: iz, it
  character (len=RHDF5_MAX_STRING) :: FileAcc
  integer :: PfileId, TfileId, OfileId

  logical, dimension(:,:), allocatable :: SelectGrid

  ! Get the command line arguments
  call GetMyArgs(LargeString, PressFile, PressVar, TopoFile, TopoVar, MaxElev, MinLat, MaxLat, MinLon, MaxLon, PressRad, OutFile)

  print'(a)',       'Creating HDF5 TC center data:'
  print'(a)',       ''
  print'(a)',       '  Sea level pressure data: '
  print'(a,x,a)',   '     File:', trim(PressFile)
  print'(a,x,a)',   '     Dataset:', trim(PressVar)
  print'(a)',       ''
  print'(a)',       '  Topography file: '
  print'(a,x,a)',   '     File:', trim(TopoFile)
  print'(a,x,a)',   '     Variable:', trim(TopoVar)
  print'(a)',       ''
  print'(a)',       '  Thresholds: '
  print'(a,f8.2)',  '     Maximum Elevation (m):', MaxElev
  print'(a,2f8.2)', '     Latitude Range (degrees):', MinLat, MaxLat
  print'(a,2f8.2)', '     Longitude Range (degrees):', MinLon, MaxLon
  print'(a,f8.2)',  '     Radius for Environmental Pressure (km):', PressRad
  print'(a)',       ''
  print'(a,x,a)',   '  Output file name:  ', trim(OutFile)
  print'(a)',       ''
  flush(6)

  ! See if we have access to all of the required variables. Check to make sure the horizontal
  ! and time dimensions match between all vars. Compare against the model file.
  !
  ! Since we will be processing one time step at a time, the effective number of dimensions
  ! on each of the input is decreased by one. We want to eliminate the time dimension which
  ! is always the last one. This works out conveniently since all we need to do is decrement
  ! the number of dimensions by one.

  Press%vname = trim(PressVar)
  Topo%vname = trim(TopoVar)

  call rhdf5_read_init(PressFile, Press)
  call rhdf5_read_init(TopoFile, Topo)

  if (.not. DimsMatch(Press, Topo)) then
    print'(a)', 'ERROR: horizontal and time dimensions of pressure and topography variables do not match'
    stop
  endif

  ! 2D vars
  Nx = Press%dims(1)
  Ny = Press%dims(2)
  Nz = 1
  Nt = Press%dims(3)

  ! Coordinates
  call SetOutCoords(PressFile, Xcoords, Ycoords, Zcoords, Tcoords)

  ! Set the output dimensions and coordinates to those of the selected input var
  
  ! Convert lat (x coords) and lon (y coords) to distances in km
  allocate(XcoordsKm(Nx))
  allocate(YcoordsKm(Ny))
  call ConvertGridCoords(Nx, Ny, Xcoords%vdata, Ycoords%vdata, XcoordsKm, YcoordsKm)

  DeltaX = (XcoordsKm(2) - XcoordsKm(1)) * 1000.0
  DeltaY = (YcoordsKm(2) - YcoordsKm(1)) * 1000.0

  print*, 'Horizontal grid info:'
  print*, '  X range (min lon, max lon) --> (min x, max x): '
  print*, '    ', Xcoords%vdata(1), Xcoords%vdata(Nx), XcoordsKm(1), XcoordsKm(Nx)
  print*, '  Y range (min lat, max lat) --> (min y, max y): '
  print*, '    ', Ycoords%vdata(1), Ycoords%vdata(Ny), YcoordsKm(1), YcoordsKm(Ny)
  print*, ''
  print*, '  Delta x of domain:     ', DeltaX
  print*, '  Delta y of domain:     ', DeltaY
  print*, ''
  print*, 'Vertical grid info:'
  do iz = 1, Nz
    print*, '  ', iz, ' --> ', Zcoords%vdata(iz)
  end do
  print*, ''
  flush(6)

  ! Stats
  print*, 'Gridded data information:'
  print*, '  Number of x (longitude) points:          ', Nx
  print*, '  Number of y (latitude) points:           ', Ny
  print*, '  Number of z (vertical level) points:     ', Nz
  print*, '  Number of t (time) points:               ', Nt
  print*, ''
  print*, '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  print*, ''
  print*, '  Grid delta x: ', DeltaX
  print*, '  Grid delta y: ', DeltaY
  print*, ''
  flush(6)

  ! prepare input vars for reading --> decrement their ndims by one
  Press%ndims = Press%ndims - 1
  Topo%ndims = Topo%ndims - 1

  ! prepare for writing the radius and Phi values (polar coords) into the output file.
  Radius%vname = 'radius'
  Radius%ndims = 2
  Radius%dims(1) = Nx
  Radius%dims(2) = Ny
  Radius%dimnames(1) = 'x'
  Radius%dimnames(2) = 'y'
  Radius%units = 'km'
  Radius%descrip = 'radius from storm center'
  allocate(Radius%vdata(Nx*Ny))

  Phi%vname = 'phi'
  Phi%ndims = 2
  Phi%dims(1) = Nx
  Phi%dims(2) = Ny
  Phi%dimnames(1) = 'x'
  Phi%dimnames(2) = 'y'
  Phi%units = 'radians'
  Phi%descrip = 'angle from storm center'
  allocate(Phi%vdata(Nx*Ny))

  ! ndims == 0 means this var will be a time series, f(t)

  ! Minimum Pressure value and location
  MinPress%vname = 'min_press'
  MinPress%ndims = 0
  MinPress%units = 'mb'
  MinPress%descrip = 'minimum SLP of storm'
  allocate(MinPress%vdata(1))

  MinPressX%vname = 'min_press_xloc'
  MinPressX%ndims = 0
  MinPressX%units = 'deg lon'
  MinPressX%descrip = 'longitude location of minimum SLP of storm'
  allocate(MinPressX%vdata(1))
  
  MinPressY%vname = 'min_press_yloc'
  MinPressY%ndims = 0
  MinPressY%units = 'deg lat'
  MinPressY%descrip = 'latitude location of minimum SLP of storm'
  allocate(MinPressY%vdata(1))

  MinPressXidx%vname = 'min_press_x_index'
  MinPressXidx%ndims = 0
  MinPressXidx%units = 'lon index'
  MinPressXidx%descrip = 'longitude index of minimum SLP of storm'
  allocate(MinPressXidx%vdata(1))
  
  MinPressYidx%vname = 'min_press_y_index'
  MinPressYidx%ndims = 0
  MinPressYidx%units = 'lat index'
  MinPressYidx%descrip = 'latitude index of minimum SLP of storm'
  allocate(MinPressYidx%vdata(1))

  ! Pressure Centroid location
  PressCentX%vname = 'press_cent_xloc'
  PressCentX%ndims = 0
  PressCentX%units = 'deg lon'
  PressCentX%descrip = 'longitude location of minimum SLP of storm'
  allocate(PressCentX%vdata(1))
  
  PressCentY%vname = 'press_cent_yloc'
  PressCentY%ndims = 0
  PressCentY%units = 'deg lat'
  PressCentY%descrip = 'latitude location of minimum SLP of storm'
  allocate(PressCentY%vdata(1))

  PressCentXidx%vname = 'press_cent_x_index'
  PressCentXidx%ndims = 0
  PressCentXidx%units = 'lon index'
  PressCentXidx%descrip = 'longitude index of minimum SLP of storm'
  allocate(PressCentXidx%vdata(1))
  
  PressCentYidx%vname = 'press_cent_y_index'
  PressCentYidx%ndims = 0
  PressCentYidx%units = 'lat index'
  PressCentYidx%descrip = 'latitude index of minimum SLP of storm'
  allocate(PressCentYidx%vdata(1))

  ! Open the input files and the output file
  FileAcc = 'R'

  call rhdf5_open_file(PressFile, FileAcc, 0, PfileId)
  print*, 'Reading HDF5 file: ', trim(PressFile)
  call rhdf5_open_file(TopoFile, FileAcc, 0, TfileId)
  print*, 'Reading HDF5 file: ', trim(TopoFile)

  print*, ''

  FileAcc = 'W'
  call rhdf5_open_file(OutFile, FileAcc, 1, OfileId)
  print*, 'Writing HDF5 file: ', trim(OutFile)
  print*, ''

  ! Allocate the selection grid. This is a 2D logical array that will
  ! contain the results of applying the MaxElev, and the "inside the box"
  ! tests so that these tests don't have to be repeated multiple times
  ! in the flow.
  allocate(SelectGrid(Nx,Ny))

  ! Do the filtering one time step at a time.
  !
  ! The necessary input files have been opened, data buffers allocated,
  ! plus the time dimensions have been "removed" from the descriptions
  ! of the input variables.
  !

  allocate(StormX(Nt))  ! use these to record storm location in the following loop
  allocate(StormY(Nt))  ! these will be used after the loop to caclulate storm motion (speed)
  do it = 1, Nt
    ! read the input vars
    call rhdf5_read_variable(PfileId, Press%vname, Press%ndims, it, Press%dims, rdata=Press%vdata)
    call rhdf5_read_variable(TfileId, Topo%vname,  Topo%ndims,  it, Topo%dims,  rdata=Topo%vdata)

    ! Mark the points that will be sampled accroding to their elevation, and location within
    ! the horizontal domain.
    call SetSelectGrid(Nx, Ny, Topo%vdata, Xcoords%vdata, Ycoords%vdata, MaxElev, MinLat, MaxLat, MinLon, MaxLon, SelectGrid)

    ! Find the storm center using minumum pressure - this is the guess that seeds the
    ! pressure centroid calculation
    call FindMinPressure(Nx, Ny, Press%vdata, SelectGrid, MinPressIx, MinPressIy, MinPress%vdata(1))

    ! Record the lat,lon of the miniumum pressure
    MinPressX%vdata(1) = Xcoords%vdata(MinPressIx)
    MinPressY%vdata(1) = Ycoords%vdata(MinPressIy)
    MinPressXidx%vdata(1) = float(MinPressIx)
    MinPressYidx%vdata(1) = float(MinPressIy)

    ! Calculate the polar coordinates of all horizontal points from the minimum pressure
    call CalcPolarCoords(Nx, Ny, Radius%vdata, Phi%vdata, XcoordsKm, YcoordsKm, MinPressIx, MinPressIy)

    ! Find the pressure centroid
    call FindPressureCent(Nx, Ny, Press%vdata, SelectGrid, Radius%vdata, XcoordsKm, YcoordsKm, PressRad, PressCentIx, PressCentIy)

    ! Record the lat,lon of the pressure centroid
    PressCentX%vdata(1) = Xcoords%vdata(PressCentIx)
    PressCentY%vdata(1) = Ycoords%vdata(PressCentIy)
    PressCentXidx%vdata(1) = float(PressCentIx)
    PressCentYidx%vdata(1) = float(PressCentIy)

    ! Record storm location in km
    StormX(it) = XcoordsKm(PressCentIx)
    StormY(it) = YcoordsKm(PressCentIy)

    ! Calculate the radius of all horizontal points from the pressure centroid
    call CalcPolarCoords(Nx, Ny, Radius%vdata, Phi%vdata, XcoordsKm, YcoordsKm, PressCentIx, PressCentIy)

    ! Write out the results
    ! Polor coords (from pressure centroid)
    call rhdf5_write_variable(OfileId, Radius%vname, Radius%ndims, it, Radius%dims, &
      Radius%units, Radius%descrip, Radius%dimnames, rdata=Radius%vdata)
    call rhdf5_write_variable(OfileId, Phi%vname, Phi%ndims, it, Phi%dims, &
      Phi%units, Phi%descrip, Phi%dimnames, rdata=Phi%vdata)

    ! Minimum pressure
    call rhdf5_write_variable(OfileId, MinPress%vname, MinPress%ndims, it, MinPress%dims, &
      MinPress%units, MinPress%descrip, MinPress%dimnames, rdata=MinPress%vdata)
    call rhdf5_write_variable(OfileId, MinPressX%vname, MinPressX%ndims, it, MinPressX%dims, &
      MinPressX%units, MinPressX%descrip, MinPressX%dimnames, rdata=MinPressX%vdata)
    call rhdf5_write_variable(OfileId, MinPressY%vname, MinPressY%ndims, it, MinPressY%dims, &
      MinPressY%units, MinPressY%descrip, MinPressY%dimnames, rdata=MinPressY%vdata)
    call rhdf5_write_variable(OfileId, MinPressXidx%vname, MinPressXidx%ndims, it, MinPressXidx%dims, &
      MinPressXidx%units, MinPressXidx%descrip, MinPressXidx%dimnames, rdata=MinPressXidx%vdata)
    call rhdf5_write_variable(OfileId, MinPressYidx%vname, MinPressYidx%ndims, it, MinPressYidx%dims, &
      MinPressYidx%units, MinPressYidx%descrip, MinPressYidx%dimnames, rdata=MinPressYidx%vdata)

    ! Pressure centroid
    call rhdf5_write_variable(OfileId, PressCentX%vname, PressCentX%ndims, it, PressCentX%dims, &
      PressCentX%units, PressCentX%descrip, PressCentX%dimnames, rdata=PressCentX%vdata)
    call rhdf5_write_variable(OfileId, PressCentY%vname, PressCentY%ndims, it, PressCentY%dims, &
      PressCentY%units, PressCentY%descrip, PressCentY%dimnames, rdata=PressCentY%vdata)
    call rhdf5_write_variable(OfileId, PressCentXidx%vname, PressCentXidx%ndims, it, PressCentXidx%dims, &
      PressCentXidx%units, PressCentXidx%descrip, PressCentXidx%dimnames, rdata=PressCentXidx%vdata)
    call rhdf5_write_variable(OfileId, PressCentYidx%vname, PressCentYidx%ndims, it, PressCentYidx%dims, &
      PressCentYidx%units, PressCentYidx%descrip, PressCentYidx%dimnames, rdata=PressCentYidx%vdata)

    ! cleanup
    deallocate(Press%vdata)
    deallocate(Topo%vdata)


    ! Write out status to screen every 50 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,50) .eq. 0) then
      print*, 'Working: Timestep: ', it

      print'(a,i3,a,i3,a,g15.2,a,g15.2,a)', '    Pressure Centroid: (', PressCentIx, ', ', PressCentIy, &
       ') --> (', XcoordsKm(PressCentIx), ', ', YcoordsKm(PressCentIy), ')'
      print*, '   Minumum Pressure: ', MinPress%vdata(1)

      print*, ''
    endif
  enddo

  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''

  ! Now have the time series of the location of the storm center
  ! Use these locations to caclulate the storm motion (speed). Keep
  ! the x and y components separated so that downstream processes
  ! can have access to these components.
  !
  ! Calculate storm motion from pressure centroid locations.
  SpeedX%vname       = 'storm_speed_x'
  SpeedX%ndims       = 1
  SpeedX%dims(1)     = Nt
  SpeedX%dimnames(1) = 't'
  SpeedX%units       = 'm/s'
  SpeedX%descrip     = 'x component of storm speed'
  allocate(SpeedX%vdata(Nt))

  SpeedY%vname       = 'storm_speed_y'
  SpeedY%ndims       = 1
  SpeedY%dims(1)     = Nt
  SpeedY%dimnames(1) = 't'
  SpeedY%units       = 'm/s'
  SpeedY%descrip     = 'y component of storm speed'
  allocate(SpeedY%vdata(Nt))

  call CalcStormSpeed(Nt, StormX, StormY, Tcoords%vdata(1), SpeedX%vdata(1), SpeedY%vdata(1))

  call rhdf5_write(OutFile, SpeedX, 1)
  call rhdf5_write(OutFile, SpeedY, 1)

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
  call rhdf5_attach_dimensions(OutFile, Radius)
  call rhdf5_attach_dimensions(OutFile, Phi)

  call rhdf5_attach_dimensions(OutFile, MinPress)
  call rhdf5_attach_dimensions(OutFile, MinPressX)
  call rhdf5_attach_dimensions(OutFile, MinPressY)
  call rhdf5_attach_dimensions(OutFile, MinPressXidx)
  call rhdf5_attach_dimensions(OutFile, MinPressYidx)

  call rhdf5_attach_dimensions(OutFile, PressCentX)
  call rhdf5_attach_dimensions(OutFile, PressCentY)
  call rhdf5_attach_dimensions(OutFile, PressCentXidx)
  call rhdf5_attach_dimensions(OutFile, PressCentYidx)

  call rhdf5_attach_dimensions(OutFile, SpeedX)
  call rhdf5_attach_dimensions(OutFile, SpeedY)

  ! cleanup
  call rhdf5_close_file(OfileId)
  call rhdf5_close_file(PfileId)
  call rhdf5_close_file(TfileId)

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

 subroutine GetMyArgs(Nstr, PressFile, PressVar, TopoFile, TopoVar, MaxElev, MinLat, MaxLat, MinLon, MaxLon, PressRad, OutFile)
  implicit none

  integer :: Nstr
  character (len=Nstr) :: PressFile, PressVar, TopoFile, TopoVar, OutFile
  real :: MaxElev, MinLat, MaxLat, MinLon, MaxLon, PressRad

  integer :: iargc, i, Nargs
  character (len=Nstr) :: Arg
  logical :: BadArg

  PressFile = ''
  PressVar  = ''
  TopoFile  = ''
  TopoVar   = ''
  MaxElev   = 0.0
  MinLat    = 0.0
  MaxLat    = 0.0
  MinLon    = 0.0
  MaxLon    = 0.0
  OutFile   = ''

  Nargs = iargc()

  if (Nargs .ne. 11) then
    print'(a)', 'ERROR: must specify exactly 5 argments'
    print'(a)', ''

    print'(a)', 'USAGE: find_tc_center <PressFile> <PressVar> <TopoFile> <TopoVar> <MaxElev> <MinLat> <MaxLat> <MinLon> <MaxLon> <OutFile>'
    print'(a)', '        <PressFile>: sea level pressure input file'
    print'(a)', '        <PressVar>: sea level pressure input variable name'
    print'(a)', '        <TopoFile>: topography input file'
    print'(a)', '        <TopoVar>: topography input variable name'
    print'(a)', '        <MaxElev>: only select points with elevation < MaxElev (meters)'
    print'(a)', '        <MinLat> <MaxLat> <MinLon> <MaxLon>: only select points within the box defined by these arguments (degrees)'
    print'(a)', '        <PressRad>: radius from storm center for sampling environmental pressure (km)'
    print'(a)', '        <OutFile>: output file name'
    print'(a)', ''

    stop
  end if

  ! Walk through list of arguments
  BadArg = .false.
  i = 1
  do while (i .le. Nargs)
    call getarg(i, Arg)

    if (i .eq. 1) then
      PressFile = Arg
      i = i + 1
    else if (i .eq. 2) then
      PressVar = Arg
      i = i + 1
    else if (i .eq. 3) then
      TopoFile = Arg
      i = i + 1
    else if (i .eq. 4) then
      TopoVar = Arg
      i = i + 1
    else if (i .eq. 5) then
      read(Arg,'(f15.7)') MaxElev
      i = i + 1
    else if (i .eq. 6) then
      read(Arg,'(f15.7)') MinLat
      i = i + 1
    else if (i .eq. 7) then
      read(Arg,'(f15.7)') MaxLat
      i = i + 1
    else if (i .eq. 8) then
      read(Arg,'(f15.7)') MinLon
      i = i + 1
    else if (i .eq. 9) then
      read(Arg,'(f15.7)') MaxLon
      i = i + 1
    else if (i .eq. 10) then
      read(Arg,'(f15.7)') PressRad
      i = i + 1
    else if (i .eq. 11) then
      OutFile = Arg
      i = i + 1
    endif
  enddo

  if (MaxElev .lt. 0) then
    print'(a)', 'ERROR: <MaxElev> must be greater than or equal to zero'
    BadArg = .true.
  endif
  if (MinLat .ge. MaxLat) then
    print'(a)', 'ERROR: <MinLat> must be less than <MaxLat>'
    BadArg = .true.
  endif
  if (MinLon .ge. MaxLon) then
    print'(a)', 'ERROR: <MinLon> must be less than <MaxLon>'
    BadArg = .true.
  endif

  if (BadArg) then
    stop
  endif
  
  return
end subroutine GetMyArgs

end program find_tc_center
