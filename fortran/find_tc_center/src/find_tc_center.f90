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
  type (Rhdf5Var) :: Radius, MinP, StormX, StormY, StormXidx, StormYidx
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm
  real :: DeltaX, DeltaY, MaxElev
  real :: MinLat, MinLon, MaxLat, MaxLon
  integer :: StmIx, StmIy

  integer :: Nx, Ny, Nz, Nt
  integer :: iz, it
  character (len=RHDF5_MAX_STRING) :: FileAcc
  integer :: PfileId, TfileId, OfileId

  ! Get the command line arguments
  call GetMyArgs(LargeString, PressFile, PressVar, TopoFile, TopoVar, OutFile)
  ! Don't consider points that have elevation > MaxElev
  MaxElev = 1.0
  ! Define box to restrict search for min pressure
  !     Lat range: 10 to 23
  !     Lon range: -39 to -21
  MinLat = 10.0
  MaxLat = 23.0

  MinLon = -39.0
  MaxLon = -21.0

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
  print'(a,f8.2)',  '     Maximum Elevation:', MaxElev
  print'(a,2f8.2)', '     Latitude Range:', MinLat, MaxLat
  print'(a,2f8.2)', '     Longitude Range:', MinLon, MaxLon
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
  call ConvertGridCoords(Nx, Ny, Nz, Xcoords%vdata, Ycoords%vdata, XcoordsKm, YcoordsKm)

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

  ! prepare for writing the radius values into the output file.
  Radius%vname = 'radius'
  Radius%ndims = 2
  Radius%dims(1) = Nx
  Radius%dims(2) = Ny
  Radius%dimnames(1) = 'x'
  Radius%dimnames(2) = 'y'
  Radius%units = 'km'
  Radius%descrip = 'radius from storm center'
  allocate(Radius%vdata(Nx*Ny))

  ! ndims == 0 means this var will be a time series, f(t)
  MinP%vname = 'min_press'
  MinP%ndims = 0
  MinP%units = 'mb'
  MinP%descrip = 'minimum SLP of storm'
  allocate(MinP%vdata(1))

  StormX%vname = 'min_press_xloc'
  StormX%ndims = 0
  StormX%units = 'deg lon'
  StormX%descrip = 'longitude location of minimum SLP of storm'
  allocate(StormX%vdata(1))
  
  StormY%vname = 'min_press_yloc'
  StormY%ndims = 0
  StormY%units = 'deg lat'
  StormY%descrip = 'latitude location of minimum SLP of storm'
  allocate(StormY%vdata(1))

  StormXidx%vname = 'min_press_x_index'
  StormXidx%ndims = 0
  StormXidx%units = 'lon index'
  StormXidx%descrip = 'longitude index of minimum SLP of storm'
  allocate(StormXidx%vdata(1))
  
  StormYidx%vname = 'min_press_y_index'
  StormYidx%ndims = 0
  StormYidx%units = 'lat index'
  StormYidx%descrip = 'latitude index of minimum SLP of storm'
  allocate(StormYidx%vdata(1))

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

  ! Do the filtering one time step at a time.
  !
  ! The necessary input files have been opened, data buffers allocated,
  ! plus the time dimensions have been "removed" from the descriptions
  ! of the input variables.
  !

  do it = 1, Nt
    ! read the input vars
    call rhdf5_read_variable(PfileId, Press%vname, Press%ndims, it, Press%dims, rdata=Press%vdata)
    call rhdf5_read_variable(TfileId, Topo%vname,  Topo%ndims,  it, Topo%dims,  rdata=Topo%vdata)

    ! zero out radius
    Radius%vdata = 0.0

    ! Find the storm center
    call FindStormCenter(Nx, Ny, Press%vdata, Topo%vdata, Xcoords%vdata, Ycoords%vdata, &
            MaxElev, MinLat, MaxLat, MinLon, MaxLon, StmIx, StmIy, MinP%vdata(1))

    ! Record the lat,lon of the storm center
    StormX%vdata(1) = Xcoords%vdata(StmIx)
    StormY%vdata(1) = Ycoords%vdata(StmIy)
    StormXidx%vdata(1) = float(StmIx)
    StormYidx%vdata(1) = float(StmIy)

    ! Write out the results
    call rhdf5_write_variable(OfileId, Radius%vname, Radius%ndims, it, Radius%dims, &
      Radius%units, Radius%descrip, Radius%dimnames, rdata=Radius%vdata)
    call rhdf5_write_variable(OfileId, MinP%vname, MinP%ndims, it, MinP%dims, &
      MinP%units, MinP%descrip, MinP%dimnames, rdata=MinP%vdata)
    call rhdf5_write_variable(OfileId, StormX%vname, StormX%ndims, it, StormX%dims, &
      StormX%units, StormX%descrip, StormX%dimnames, rdata=StormX%vdata)
    call rhdf5_write_variable(OfileId, StormY%vname, StormY%ndims, it, StormY%dims, &
      StormY%units, StormY%descrip, StormY%dimnames, rdata=StormY%vdata)
    call rhdf5_write_variable(OfileId, StormXidx%vname, StormXidx%ndims, it, StormXidx%dims, &
      StormXidx%units, StormXidx%descrip, StormXidx%dimnames, rdata=StormXidx%vdata)
    call rhdf5_write_variable(OfileId, StormYidx%vname, StormYidx%ndims, it, StormYidx%dims, &
      StormYidx%units, StormYidx%descrip, StormYidx%dimnames, rdata=StormYidx%vdata)

    ! cleanup
    deallocate(Press%vdata)
    deallocate(Topo%vdata)


    ! Write out status to screen every 100 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,100) .eq. 0) then
      print*, 'Working: Timestep: ', it

      print'(a,i3,a,i3,a,g,a,g,a)', '    Storm Center: (', StmIx, ', ', StmIy, &
       ') --> (', XcoordsKm(StmIx), ', ', YcoordsKm(StmIy), ')'
      print*, '   Minumum pressure: ', MinP%vdata(1)

      print*, ''
    endif
  enddo

  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''

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
  call rhdf5_attach_dimensions(OutFile, MinP)
  call rhdf5_attach_dimensions(OutFile, StormX)
  call rhdf5_attach_dimensions(OutFile, StormY)
  call rhdf5_attach_dimensions(OutFile, StormXidx)
  call rhdf5_attach_dimensions(OutFile, StormYidx)
  call rhdf5_attach_dimensions(OutFile, Radius)

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

 subroutine GetMyArgs(Nstr, PressFile, PressVar, TopoFile, TopoVar, OutFile)
  implicit none

  integer :: Nstr
  character (len=Nstr) :: PressFile, PressVar, TopoFile, TopoVar, OutFile

  integer :: iargc, i, j, Nargs, Nfld
  character (len=Nstr) :: Arg


  PressFile = ''
  PressVar = ''
  TopoFile = ''
  TopoVar = ''
  OutFile = ''

  Nargs = iargc()

  if (Nargs .ne. 5) then
    print'(a)', 'ERROR: must specify exactly 5 argments'
    print'(a)', ''

    print'(a)', 'USAGE: find_tc_center <PressFile> <PressVar> <TopoFile> <TopoVar> <OutFile>'
    print'(a)', '        <PressFile>: sea level pressure input file'
    print'(a)', '        <PressVar>: sea level pressure input variable name'
    print'(a)', '        <TopoFile>: topography input file'
    print'(a)', '        <TopoVar>: topography input variable name'
    print'(a)', '        <OutFile>: output file name'
    print'(a)', ''

    stop
  end if

  ! Walk through list of arguments
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
      OutFile = Arg
      i = i + 1
    endif
  enddo
  
  return
end subroutine GetMyArgs

!**********************************************************************
! FindStormCenter
!
! This routine will locate the storm center using the simple hueristic
! of the center being where the minimum surface pressure exists.
!
! Argument it holds the time step that you want to analyze. (iStmCtr, jStmCtr)
! hold the grid position of the minumum pressure value on the first vertical
! level (iz = 1).
!

subroutine FindStormCenter(Nx, Ny, Press, Topo, Lon, Lat, MaxElev, &
                           MinLat, MaxLat, MinLon, MaxLon, StmCtrX, &
                           StmCtrY, MinP)
  implicit none

  integer :: Nx, Ny
  real, dimension(Nx,Ny) :: Press, Topo
  real, dimension(Nx) :: Lon
  real, dimension(Ny) :: Lat
  real :: MaxElev, MinLat, MaxLat, MinLon, MaxLon
  integer :: StmCtrX, StmCtrY
  real :: MinP

  integer :: ix, iy
  logical :: SelectPoint

  MinP = 1e10 ! ridiculously large pressure
  StmCtrX = 0 
  StmCtrY = 0 

  do ix = 1, Nx
    do iy = 1, Ny
      SelectPoint = Topo(ix,iy) .lt. MaxElev
      SelectPoint = SelectPoint .and. (Lon(ix) .ge. MinLon)
      SelectPoint = SelectPoint .and. (Lon(ix) .le. MaxLon)
      SelectPoint = SelectPoint .and. (Lat(iy) .ge. MinLat)
      SelectPoint = SelectPoint .and. (Lat(iy) .le. MaxLat)

      if (SelectPoint) then
        if (Press(ix,iy) .lt. MinP) then
          MinP = Press(ix,iy)
          StmCtrX = ix
          StmCtrY = iy
        endif
      endif
    enddo
  enddo
  
  return
end subroutine FindStormCenter

end program find_tc_center
