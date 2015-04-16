!***************************************************************
! Program to do azimuthial averaging
!
! This program will read in GRADS data from a RAMS simulation, find
! the storm center and perform azimuthial averaging on the given
! quantity.
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
!   3. number of radial bands to split data into
!   4. RAMS quantity to perform the averaging on
!
! Output
!   The output will be two files which make a GRADS data set. One
!   file is the control file and the other is the binary data.
!

program azavg
  use rhdf5_utils
  use diag_utils

  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10
  integer, parameter :: MaxCoords=1000
  real, parameter :: UndefVal=-999.0
  integer, parameter :: MaxArgFields = 20

  integer NumRbands
  character (len=MediumString) :: InDir
  character (len=MediumString) :: InSuffix
  character (len=MediumString) :: OutFile
  character (len=MediumString) :: FilterFile
  character (len=MediumString) :: StormCenterFile
  character (len=LittleString) :: VarToAvg, RevuVar
  logical :: DoHorizVel, DoTangential
  character (len=MediumString), dimension(MaxArgFields) :: ArgList
  integer :: Nfields
  logical :: DoHist
  character (len=LargeString) :: BinsFile
  integer :: NumBins
  real, dimension(:), allocatable :: HistBins

  ! Data arrays: need one for w (vertical velocity), press (pressure)
  ! and the var we are doing the averaging on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number

  type (Rhdf5Var) :: U, V, Avar, Aavg, Rcoords, Zcoords, Tcoords, Bcoords, Xcoords, Ycoords
  type (Rhdf5Var) :: Filter, Radius, StormXindex, StormYindex
  character (len=LittleString) :: rh5f_facc
  integer :: rh5f_filter, rh5f_center, rh5f_u, rh5f_v, rh5f_avar, rh5f_aavg
  character (len=MediumString) :: Ufile, Vfile, AvarFile
  integer :: i_storm, j_storm
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm
  real :: StormX, StormY, MaxRadius

  integer :: i
  integer :: AvarNelems

  integer :: ix, iy, it
  integer :: Nx, Ny, Nz, Nt, FilterNz

  real :: DeltaX, DeltaY, RbandInc

  ! Get the command line arguments
  call GetMyArgs(InDir, InSuffix, OutFile, FilterFile, StormCenterFile, NumRbands, MaxRadius, VarToAvg, &
    RevuVar)

  if (VarToAvg(1:5) .eq. 'hist:') then
    call String2List(VarToAvg, ':', ArgList, MaxArgFields, Nfields, 'hist spec')
    if (Nfields .eq. 3) then
      ! got the right amount of fields
      !   field    value
      !    1       'hist'
      !    2       variable name
      !    3       bins file
      DoHist = .true.
      VarToAvg = trim(ArgList(2))
      BinsFile = trim(ArgList(3))
    else
      write (*,*) 'ERROR: average function hist requires five fields: hist:<var>:<num_bins>:<bin_start>:<bin_size>'
      stop
    endif
  else
    ! Not doing histograms
    DoHist = .false.
  endif

  ! Do the check for horizontal velocity cases in order to support
  ! both averaging and creation of a histogram.
  if (VarToAvg == 'speed_t') then
    DoHorizVel = .true.
    DoTangential = .true.
  else
    if (VarToAvg == 'speed_r') then
      DoHorizVel = .true.
      DoTangential = .false.
    else
      DoHorizVel = .false.
      DoTangential = .false.
    end if
  end if

  write (*,*) 'Calculating azimuthal average for RAMS data:'
  write (*,*) '  Input directory: ', trim(InDir)
  write (*,*) '  Input file name suffix: ', trim(InSuffix)
  write (*,*) '  Output file name:  ', trim(OutFile)
  write (*,*) '  Filter file name: ', trim(FilterFile)
  write (*,*) '  Storm center file name: ', trim(StormCenterFile)
  write (*,*) '  Number of radial bands: ', NumRbands
  if (DoHorizVel) then
    if (DoTangential) then
      write (*,*) '  RAMS variable that is being averaged: Tangential Horizontal Velocity'
    else
      write (*,*) '  RAMS variable that is being averaged: Radial Horizontal Velocity'
    end if
  else
    write (*,*) '  RAMS variable that is being averaged: ', trim(VarToAvg)
    write (*,*) '  REVU variable name: ', trim(RevuVar)
  end if
  if (DoHist) then
    write (*,*) '  Creating histogram counts:'
    write (*,*) '    Bins File: ', trim(BinsFile)
  endif
  write (*,*) ''

  ! Read the variable information from the HDF5 files and check for consistency.
  !
  ! Always use filter file
  !
  ! If doing horizontal velocity, ie VarToAvg is 'speed_t' or 'speed_r', set up U and V
  ! otherwise set up the variable named in VarToAvg
  !

  ! Get variable and coordinate descriptions from FilterFile
  ! Get the z coordinate description from the variable file since this
  ! can change between 2D and 3D variables.
  Filter%vname = 'filter'
  Radius%vname = 'radius'
  StormXindex%vname = 'press_cent_x_index'
  StormYindex%vname = 'press_cent_y_index'
  Xcoords%vname = 'x_coords'
  Ycoords%vname = 'y_coords'
  Tcoords%vname = 't_coords'
  call rhdf5_read_init(FilterFile, Filter)
  call rhdf5_read_init(FilterFile, Xcoords)
  call rhdf5_read_init(FilterFile, Ycoords)
  call rhdf5_read_init(FilterFile, Tcoords)

  call rhdf5_read_init(StormCenterFile, Radius)
  call rhdf5_read_init(StormCenterFile, StormXindex)
  call rhdf5_read_init(StormCenterFile, StormYindex)

  ! Read in the 1D variable data
  call rhdf5_read(FilterFile, Xcoords)
  call rhdf5_read(FilterFile, Ycoords)
  call rhdf5_read(FilterFile, Tcoords)

  call rhdf5_read(StormCenterFile, StormXindex)
  call rhdf5_read(StormCenterFile, StormYindex)

  if (DoHorizVel) then
    Ufile = trim(InDir) // '/u' // trim(InSuffix)
    U%vname = 'u'
    call rhdf5_read_init(Ufile, U)

    Zcoords%vname = 'z_coords'
    call rhdf5_read_init(Ufile, Zcoords)
    call rhdf5_read(Ufile, Zcoords)

    Vfile = trim(InDir) // '/v' // trim(InSuffix)
    V%vname = 'v'
    call rhdf5_read_init(Vfile, V)

    ! Initialize the elements in Avar
    Avar%vname = trim(VarToAvg)
    Avar%ndims = U%ndims
    Avar%dims = U%dims
    Avar%dimnames = U%dimnames
    Avar%units = 'm/s'
    if (DoTangential) then
      Avar%descrip = 'tangential wind speed'
    else
      Avar%descrip = 'radial wind speed'
    endif
  else
    AvarFile = trim(InDir) // '/' // trim(VarToAvg) // trim(InSuffix)
    Avar%vname = trim(RevuVar)
    call rhdf5_read_init(AvarFile, Avar)

    Zcoords%vname = 'z_coords'
    call rhdf5_read_init(AvarFile, Zcoords)
    call rhdf5_read(AvarFile, Zcoords)
  endif

  ! check that the variable dimensions (size and coordinate values) match up, if this
  ! isn't true, then the subsequent anlysis will be bogus
  !
  ! third arg of GvarDimsMatch() is true if one of the vars is 2D, else false
  if (DoHorizVel) then
    if (.not. (DimsMatch(Filter, U) .and. DimsMatch(Filter, V) .and. &
               DimsMatch(Filter, Radius))) then
      write (*,*) 'ERROR: dimensions of u, v, filter, and radius do not match'
      stop
    endif
  else
    if (.not. (DimsMatch(Filter, Radius) .and. DimsMatch(Filter, Avar))) then
      write (*,*) 'ERROR: dimensions of filter, radius, and ', trim(VarToAvg), ' do not match'
      stop
    endif
  endif
 
  ! Avar is either 3D (x,y,z,t) or 2D (x,y,t)
  if (Avar%ndims .eq. 4) then
    ! (x,y,z,t)
    Nx = Avar%dims(1)
    Ny = Avar%dims(2)
    Nz = Avar%dims(3)
    Nt = Avar%dims(4)
  elseif (Avar%ndims .eq. 3) then
    ! (x,y,t)
    Nx = Avar%dims(1)
    Ny = Avar%dims(2)
    Nz = 1
    Nt = Avar%dims(3)
  else
    write (*,*) 'ERROR: <var_to_average> must possess either 2D or 3D spatial dimensions'
    stop
  endif

  FilterNz = Filter%dims(3)

  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:               ', Nx
  write (*,*) '  Number of y (latitude) points:                ', Ny
  write (*,*) '  Number of z (vertical level) points (var):    ', Nz
  write (*,*) '  Number of z (vertical level) points (filter): ', FilterNz
  write (*,*) '  Number of t (time) points:                    ', Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''

  write (*,*) ''

  ! If doing histogram, read in the bins and the number of bins. Do this before
  ! the next section since NumBins is being used to set up the output coordinate values.
  !
  ! If not doing a histogram, set NumBins to one, and place a zero in HistBins.
  ! This will allow the code that creates Bcoords (below) to just use NumBins
  ! and HistBins without checking DoHist.
  if (DoHist) then
    ! ReadBinsFile will allocate HistBins (it needs to read in NumBins before
    ! HistBins gets allocated).
    call ReadBinsFile(BinsFile, NumBins, HistBins)
    write (*,*) 'Histogram Bins:'
    do i = 1, NumBins
      write (*,*) '  ', HistBins(i)
    enddo
    write (*,*) ''
  else
    ! Need to set NumBins to one since NumBins is being used to set the size
    ! of the y dimension in the output variable, Aavg.
    NumBins = 1
    allocate(HistBins(1))
    HistBins(1) = 0.0
  endif

  ! Convert the GRADS grid coordinates from longitude, latitude to flat plane (x and y).
  ! XcoordsKm, YcoordsKm are in units of km
  allocate(XcoordsKm(Nx))
  allocate(YcoordsKm(Ny))
  call ConvertGridCoords(Nx, Ny, Xcoords%vdata, Ycoords%vdata, XcoordsKm, YcoordsKm)

  write (*,*) 'Horzontal Grid Coordinate Info:'
  write (*,*) '  X Range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', Xcoords%vdata(1), Xcoords%vdata(Nx), XcoordsKm(1), XcoordsKm(Nx)
  write (*,*) '  Y Range (min lat, max lat) --> (min y, max y): '
  write (*,*) '    ', Ycoords%vdata(1), Ycoords%vdata(Ny), YcoordsKm(1), YcoordsKm(Ny)
  write (*,*) ''

  MaxRadius = anint(MaxRadius)           ! round to nearest integer
  DeltaX = XcoordsKm(Nx) - XcoordsKm(1)
  DeltaY = YcoordsKm(Ny) - YcoordsKm(1)
  RbandInc = MaxRadius / real(NumRbands)

  Rcoords%vname = 'x_coords'
  Rcoords%ndims = 1
  Rcoords%dims(1) = NumRbands
  Rcoords%dimnames(1) = 'x'
  Rcoords%units = 'degrees_east'
  Rcoords%descrip = 'radius in meters'
  allocate (Rcoords%vdata(NumRbands))
  do ix = 1, NumRbands
    Rcoords%vdata(ix) = real(ix) * RbandInc * 1000.0
  enddo

  write (*,*) 'Radial band information:'
  write (*,*) '  Delta x of domain:     ', DeltaX
  write (*,*) '  Delta y of domain:     ', DeltaY
  write (*,*) '  Radial distance:       ', MaxRadius
  write (*,*) '  Radial band increment: ', RbandInc
  write (*,*) ''

  ! Open files for read and write. The variables in the HDF5 files include the
  ! time dimension which is always the last dimension. We want to read and write
  ! one time step at at time which requires the variables to not indlude the time
  ! dimension so we need to remove the time dimension from the variables we are
  ! using. This turns out to be easy - since time is always the last dimension, we
  ! simply decrement the number of dimensions by one to remove time.
  !
  ! Chop off the time dimensions of the variables that have time dimensions
  !
  rh5f_facc = 'R'
  call rhdf5_open_file(FilterFile, rh5f_facc, 0, rh5f_filter)
  call rhdf5_open_file(StormCenterFile, rh5f_facc, 0, rh5f_center)

  Filter%ndims = Filter%ndims - 1
  Radius%ndims = Radius%ndims - 1

  ! set up the input variable data
  rh5f_facc = 'R'
  if (DoHorizVel) then
    call rhdf5_open_file(Ufile, rh5f_facc, 0, rh5f_u)
    call rhdf5_open_file(Vfile, rh5f_facc, 0, rh5f_v)
    U%ndims = U%ndims - 1
    V%ndims = V%ndims - 1
  else
    call rhdf5_open_file(AvarFile, rh5f_facc, 0, rh5f_avar)
  endif

  Avar%ndims = Avar%ndims - 1
  AvarNelems = 1
  do i = 1, Avar%ndims
    AvarNelems = AvarNelems * Avar%dims(i)
  enddo

  ! set up the output variable
  rh5f_facc = 'W'
  call rhdf5_open_file(OutFile, rh5f_facc, 1, rh5f_aavg)
  write (*,*) 'Writing HDF5 output: ', trim(OutFile)
  write (*,*) ''

  ! NumBins will be one when not doing histograms. Otherwise it will be
  ! the number of bins the user specified.

  Aavg%vname = trim(VarToAvg)
  Aavg%ndims = 3
  Aavg%dims(1) = NumRbands
  Aavg%dims(2) = NumBins
  Aavg%dims(3) = Nz
  Aavg%dimnames(1) = 'x'
  Aavg%dimnames(2) = 'y'
  Aavg%dimnames(3) = 'z'
  Aavg%units = Avar%units 
  Aavg%descrip = 'azimuthally averaged ' // trim(VarToAvg) 
  allocate(Aavg%vdata(Aavg%dims(1)*Aavg%dims(2)*Aavg%dims(3)))

  ! Do the averaging - one time step at a time
  do it = 1, Nt
    ! find storm center in km
    i_storm = nint(StormXindex%vdata(it))
    j_storm = nint(StormYindex%vdata(it))
    StormX = XcoordsKm(i_storm)
    StormY = YcoordsKm(j_storm)

    ! Read in the filter, radius and variable data
    call rhdf5_read_variable(rh5f_filter, Filter%vname, Filter%ndims, it, Filter%dims, rdata=Filter%vdata)
    call rhdf5_read_variable(rh5f_center, Radius%vname, Radius%ndims, it, Radius%dims, rdata=Radius%vdata)

    if (DoHorizVel) then
      call rhdf5_read_variable(rh5f_u, U%vname, U%ndims, it, U%dims, rdata=U%vdata)
      call rhdf5_read_variable(rh5f_v, V%vname, V%ndims, it, V%dims, rdata=V%vdata)

      ! convert u,v to tangential or radial
      allocate(Avar%vdata(AvarNelems))
      call ConvertHorizVelocity(Nx, Ny, Nz, U%vdata, V%vdata, StormX, StormY, &
        XcoordsKm, YcoordsKm, Avar%vdata, DoTangential)

      ! Free up variable memory
      deallocate(U%vdata)
      deallocate(V%vdata)
    else
      call rhdf5_read_variable(rh5f_avar, Avar%vname, Avar%ndims, it, Avar%dims, rdata=Avar%vdata)
    endif

    ! do the averaging and write out the results to the output file
    call AzimuthalAverage(Nx, Ny, Nz, FilterNz, NumRbands, NumBins, Avar%vdata, Aavg%vdata, &
      RbandInc, Filter%vdata, Radius%vdata, HistBins, DoHist, UndefVal)

    call rhdf5_write_variable(rh5f_aavg, Aavg%vname, Aavg%ndims, it, Aavg%dims, &
      Aavg%units, Aavg%descrip, Aavg%dimnames, rdata=Aavg%vdata)

    ! Free up variable memory
    deallocate(Filter%vdata)
    deallocate(Radius%vdata)
    deallocate(Avar%vdata)
    
    ! Write out status to screen every 100 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,50) .eq. 0) then
      write (*,*) 'Working: Timestep, Time: ', it, Tcoords%vdata(it)

      write (*,'(a,4f15.4)') '   Storm center: lon, lat, x, y: ', Xcoords%vdata(i_storm), Ycoords%vdata(j_storm), StormX, StormY
      write (*,*) ''
      flush(6)
    endif
  enddo


  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''
  flush(6)

  ! close the files
  call rhdf5_close_file(rh5f_filter)
  call rhdf5_close_file(rh5f_center)
  if (DoHorizVel) then
    call rhdf5_close_file(rh5f_u)
    call rhdf5_close_file(rh5f_v)
  else
    call rhdf5_close_file(rh5f_avar)
  endif
  call rhdf5_close_file(rh5f_aavg)

  ! write out the coordinate data
  ! need to create a dummy y coordinate to keep grads happy
  ! NumBins will be one and HistBins(1) will be zero when we are
  ! not doing histograms.
  Bcoords%vname = 'y_coords'
  Bcoords%ndims = 1
  Bcoords%dims(1) = NumBins
  Bcoords%dimnames(1) = 'y'
  Bcoords%units = 'degrees_north'
  Bcoords%descrip = 'dummy coordinates'
  allocate (Bcoords%vdata(NumBins))
  do iy = 1, NumBins
    Bcoords%vdata(iy) = HistBins(iy)
  enddo
  
  call rhdf5_write(OutFile, Rcoords, 1)
  call rhdf5_write(OutFile, Bcoords, 1)
  call rhdf5_write(OutFile, Zcoords, 1)
  call rhdf5_write(OutFile, Tcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(OutFile, Rcoords, 'x')
  call rhdf5_set_dimension(OutFile, Bcoords, 'y')
  call rhdf5_set_dimension(OutFile, Zcoords, 'z')
  call rhdf5_set_dimension(OutFile, Tcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(OutFile, Aavg)

  stop

Contains

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the command line arguments
!

subroutine GetMyArgs(InDir, InSuffix, OutFile, FilterFile, StormCenterFile, NumRbands, &
                     MaxRadius, VarToAvg, RevuVar)

  implicit none

  integer :: NumRbands
  real :: MaxRadius
  character (len=*) :: InDir, InSuffix, OutFile, FilterFile, StormCenterFile, VarToAvg, RevuVar

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 9) then
    write (*,*) 'ERROR: must supply exactly 9 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_dir> <in_suffix> <out_file> <filter_file> <storm_center_file> <num_radial_bands> <max_radius> <var_to_average> <revu_var_name>'
    write (*,*) '        <in_dir>: directory where input files live'
    write (*,*) '        <in_suffix>: suffix on input file names'
    write (*,*) '        <out_file>: name of output file, HDF5 format'
    write (*,*) '        <filter_file>: file containing the filter mask data'
    write (*,*) '        <storm_center_file>: file containing the radius, storm center, min pressure data'
    write (*,*) '        <num_radial_bands>: number of bands to split data into'
    write (*,*) '        <max_radius>: maximum radius for averaging (determines size of radial bands)'
    write (*,*) '        <var_to_average>: name of variable to do the averaging on'
    write (*,*) '             special var names:'
    write (*,*) '                speed_t: horizontal tangential speed'
    write (*,*) '                speed_r: horizontal radial speed'
    write (*,*) '                hist: create histogram counts instead of averages'
    write (*,*) '                  hist:<var>:<bins_file>'
    write (*,*) '                     <var>: same as <var_to_average>'
    write (*,*) '                     <bins_file>: file containing the edge values of the bins'
    write (*,*) '        <revu_var_name>: name of variable in the REVU file'
    write (*,*) ''
    stop
  end if

  call getarg(1, InDir)
  call getarg(2, InSuffix)
  call getarg(3, OutFile)
  call getarg(4, FilterFile)
  call getarg(5, StormCenterFile)

  call getarg(6, arg)       !this is a string
  read (arg, '(i15)') NumRbands !convert to integer

  call getarg(7, arg)
  read (arg, '(f15.7)') MaxRadius ! convert to real

  call getarg(8, VarToAvg)

  call getarg(9, RevuVar)

  return
end subroutine GetMyArgs

end program azavg
