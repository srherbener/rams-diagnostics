!***************************************************************
! Program to do averaging over the domain and yield a single
! data point at each time step
!
! This program will read in HDF5 data from a RAMS simulation, 
! perform an averaging function, and output the single point time
! series in HDF5 format
!

program tsavg
  use rhdf5_utils
  use diag_utils
  implicit none

  integer, parameter :: LargeString  = 512
  integer, parameter :: MediumString = 256
  integer, parameter :: LittleString = 128
  integer, parameter :: MaxArgFields = 20

  real, parameter :: UndefVal = -999.0

  character (len=MediumString) :: InDir
  character (len=MediumString) :: InSuffix
  character (len=MediumString) :: OutFile
  character (len=MediumString) :: FilterFile
  character (len=LittleString) :: AvgFunc
  character (len=MediumString), dimension(MaxArgFields) :: ArgList
  integer :: Nfields
  character (len=LittleString) :: VarDim
  logical :: UseFilter

  ! Data arrays
  ! Dims: x, y, z, t
  type (Rhdf5Var) :: InXcoords, InYcoords, InZcoords, InTcoords
  type (Rhdf5Var) :: OutXcoords, OutYcoords, OutZcoords, OrigDimSize
  type (Rhdf5Var) :: U, V, AzWind, Speed10m, Dens, Var, Filter, TserAvg, RadMaxWind
  character (len=MediumString) :: Ufile, Vfile, AzWindFile, Speed10mFile, DensFile, VarFile, InCoordFile
  character (len=LittleString) :: VarFprefix
  character (len=LittleString) :: rh5f_facc
  
  integer :: NumBins
  real :: BinStart, BinInc
  real, dimension(:), allocatable :: Bins

  integer :: rh5f_azwind, rh5f_u, rh5f_v, rh5f_speed10m, rh5f_dens, rh5f_var, rh5f_filter, rh5f_out

  integer :: id, ib, ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY
  real, dimension(:), allocatable :: InXcoordsKm, InYcoordsKm
  logical :: BadDims

  integer :: i_rmw

  ! Get the command line arguments
  call GetMyArgs(InDir, InSuffix, OutFile, AvgFunc, FilterFile, UseFilter)
  if ((AvgFunc(1:4) .eq. 'min:') .or. (AvgFunc(1:4) .eq. 'max:') .or. (AvgFunc(1:4) .eq. 'hda:')) then
    call String2List(AvgFunc, ':', ArgList, MaxArgFields, Nfields, 'avg spec')
    if (Nfields .eq. 4) then
      ! got the right amount of fields
      !   field    value
      !    1       'min', 'max', or 'hda' 
      !    2       name of variable inside the REVU file
      !    3       prefix for the REVU file name
      !    4       dimensionality of variable
      AvgFunc    = trim(ArgList(1))
      Var%vname  = trim(ArgList(2))
      VarFprefix = trim(ArgList(3))
      VarFile    = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      VarDim     = trim(ArgList(4))
    else
      write (*,*) 'ERROR: average function hda requires four fields: hda:<var>:<file>:<dim>'
      stop
    endif
  endif

  if (AvgFunc(1:5) .eq. 'hist:') then
    call String2List(AvgFunc, ':', ArgList, MaxArgFields, Nfields, 'hist spec')
    if (Nfields .eq. 7) then
      ! got the right amount of fields
      !   field    value
      !    1       'hist'
      !    2       name of variable inside the REVU file
      !    3       prefix for the REVU file name
      !    4       dimensionality of variable
      !    5       number of bins
      !    6       bin start
      !    7       bin increment
      AvgFunc    = trim(ArgList(1))
      Var%vname  = trim(ArgList(2))
      VarFprefix = trim(ArgList(3))
      VarFile    = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      VarDim     = trim(ArgList(4))
      read(ArgList(5), '(i)') NumBins
      read(ArgList(6), '(f)') BinStart
      read(ArgList(7), '(f)') BinInc
    else
      write (*,*) 'ERROR: average function hist requires seven fields: hist:<var>:<file>:<dim>:<num_bins>:<bin_start>:<bin_inc>'
      stop
    endif
  endif

  write (*,*) 'Time seris of average for RAMS data:'
  write (*,*) '  Input directory: ', trim(InDir)
  write (*,*) '  Input file suffix: ', trim(InSuffix)
  write (*,*) '  Output file:  ', trim(OutFile)
  write (*,*) '  Averaging function: ', trim(AvgFunc)
  if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. (AvgFunc .eq. 'hda')) then
    write (*,*) '    Variable name: ', trim(Var%vname)
    write (*,*) '    File name: ', trim(VarFile)
    write (*,*) '    Dimensionality: ', trim(VarDim)
  else if (AvgFunc .eq. 'hist') then
    write (*,*) '    Variable name: ', trim(Var%vname)
    write (*,*) '    File name: ', trim(VarFile)
    write (*,*) '    Dimensionality: ', trim(VarDim)
    write (*,*) '    Binning specs:'
    write (*,*) '      Number of bins: ', NumBins
    write (*,*) '      Bins start at: ', BinStart
    write (*,*) '      Delta between bins: ', BinInc
  endif
  write (*,*) '  Filter file: ', trim(FilterFile)
  write (*,*) '    Using filter: ', UseFilter
  write (*,*) ''
  flush(6)

  ! set up file and variable names
  AzWindFile = trim(InDir) // '/speed_t' // trim(InSuffix)
  AzWind%vname = 'speed_t'

  ! FilterFile is set by command line arguments
  Filter%vname = 'filter'

  DensFile = trim(InDir) // '/dn0' // trim(InSuffix)
  Dens%vname = 'dn0'

  Ufile = trim(InDir) // '/u' // trim(InSuffix)
  U%vname = 'u'

  Vfile = trim(InDir) // '/v' // trim(InSuffix)
  V%vname = 'v'

  Speed10mFile = trim(InDir) // '/speed10m' // trim(InSuffix)
  Speed10m%vname = 'speed10m'

  ! Check that the dimensions are consistent between the variables needed for
  ! the selected averaging function.
  !
  ! There is no associated filter with the max_azwind since a filter has already been
  ! applied by the azavg program (which created the azwind data). All other functions
  ! need the filter data.
  !
  ! Expect 3D vars to be: (x,y,z,t)
  !        2D vars to be: (x,y,t)
  !

  ! Set Filter%ndims to one here. This will result with Filter%ndims being equal to zero (since
  ! we need to chop of the time dimension to read the filter time step by time step) when we
  ! are not using a filter. When using a filter, Filter%ndims will get set by what's contained
  ! in FilterFile.
  if (AvgFunc .eq. 'max_azwind') then
    ! max_azwind: must not use filter (azavg already applied a filter)
    if (UseFilter) then
      write (*,*) 'ERROR: cannot use a filter with function: max_azwind'
      stop
    endif
    
    ! Read in the filter and use it to check against all the other variables
    call rhdf5_read_init(AzWindFile, AzWind)

    Nx = AzWind%dims(1)
    Ny = AzWind%dims(2)
    Nz = AzWind%dims(3)
    Nt = AzWind%dims(4)
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
           (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist')) then
    call rhdf5_read_init(VarFile, Var)

    if (VarDim .eq. '2d') then
      Nx = Var%dims(1)
      Ny = Var%dims(2)
      Nz = 1
      Nt = Var%dims(3)
    else
      Nx = Var%dims(1)
      Ny = Var%dims(2)
      Nz = Var%dims(3)
      Nt = Var%dims(4)
    endif

    ! filter is optional
    if (UseFilter) then
      call rhdf5_read_init(FilterFile, Filter)
      if (.not.(DimsMatch(Filter, Var))) then
        write (*,*) 'ERROR: dimensions of filter do not match dimensions of input variable: ', trim(Var%vname)
        stop
      endif
    endif
  else if ((AvgFunc .eq. 'horiz_ke') .or. (AvgFunc .eq. 'storm_int')) then
    ! horiz_ke and storm_int require a filter
    if (.not. UseFilter) then
      write (*,*) 'ERROR: must use a filter with functions: horiz_ke and storm_int'
      stop
    endif
    
    ! Read in the filter and use it to check against all the other variables
    call rhdf5_read_init(FilterFile, Filter)
  
    Nx = Filter%dims(1)
    Ny = Filter%dims(2)
    Nz = Filter%dims(3)
    Nt = Filter%dims(4)

    BadDims = .false.
    if (AvgFunc .eq. 'horiz_ke') then
      call rhdf5_read_init(DensFile, Dens)
      BadDims = BadDims .or. (.not.(DimsMatch(Filter, Dens)))
  
      call rhdf5_read_init(Ufile, U)
      BadDims = BadDims .or. (.not.(DimsMatch(Filter, U)))
  
      call rhdf5_read_init(Vfile, V)
      BadDims = BadDims .or. (.not.(DimsMatch(Filter, V)))

      if (BadDims) then
        write (*,*) 'ERROR: dimensions of filter, dn0, u and v do not match'
        stop
      endif
    else if (AvgFunc .eq. 'storm_int') then
      ! speed10m is a 2D variable
      call rhdf5_read_init(Speed10mFile, Speed10m)
      BadDims = BadDims .or. (.not.(DimsMatch(Filter, Speed10m)))
      Nz = 1
  
      if (BadDims) then
        write (*,*) 'ERROR: dimensions of filter and speed10m do not match'
        stop
      endif
    endif
  endif

  ! Set up the dimensions for reading in the input field data, one time step per read. 
  ! In other words, remove the time dimension from the input dimensions since we will 
  ! be incrementing through every time step in a loop. The time dimension is always the
  ! last dimension so what this boils down to is to decrement number of dimensions by one.
  write (*,*) 'Input variable information:'
  if (AvgFunc .eq. 'max_azwind') then
    AzWind%ndims = AzWind%ndims - 1
    write (*,*) '  Number of dimensions: ', AzWind%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, AzWind%ndims
      write (*,*), '    ', trim(AzWind%dimnames(id)), ': ', AzWind%dims(id)
    enddo
  else if (AvgFunc .eq. 'horiz_ke') then
    Dens%ndims = Dens%ndims - 1
    U%ndims = U%ndims - 1
    V%ndims = V%ndims - 1
    write (*,*) '  Number of dimensions: ', Dens%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Dens%ndims
      write (*,*), '    ', trim(Dens%dimnames(id)), ': ', Dens%dims(id)
    enddo
  else if (AvgFunc .eq. 'storm_int') then
    Speed10m%ndims = Speed10m%ndims - 1
    write (*,*) '  Number of dimensions: ', Speed10m%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Speed10m%ndims
      write (*,*), '    ', trim(Speed10m%dimnames(id)), ': ', Speed10m%dims(id)
    enddo
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
           (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist')) then
    Var%ndims = Var%ndims - 1
    write (*,*) '  Number of dimensions: ', Var%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Var%ndims
      write (*,*), '    ', trim(Var%dimnames(id)), ': ', Var%dims(id)
    enddo
  endif
  write (*,*) ''

  if (UseFilter) then
    Filter%ndims = Filter%ndims - 1

    write (*,*) 'Filter variable information:'
    write (*,*) '  Number of dimensions: ', Filter%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Filter%ndims
      write (*,*), '    ', trim(Filter%dimnames(id)), ': ', Filter%dims(id)
    enddo
    write (*,*) ''
  endif

  ! Set up the dimensions for the output and allocate the output data array. Always
  ! set up as if the output were 3D. This is done so that the output file can
  ! be read into GRADS which expects 3D variables. Always have (x,y,z) for the
  ! dimension names, but set the sizes of the dimensions according to the averaging
  ! function asked for.
  TserAvg%vname = trim(AvgFunc)
  if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist')) then
    TserAvg%vname = trim(TserAvg%vname) // '_' // trim(VarFprefix)
  endif
  TserAvg%descrip = 'time series averaged ' // trim(AvgFunc) 
  TserAvg%ndims = 3 
  TserAvg%dimnames(1) = 'x' 
  TserAvg%dimnames(2) = 'y' 
  TserAvg%dimnames(3) = 'z' 

  if (AvgFunc .eq. 'max_azwind') then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = 'm/s'
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max')) then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = Var%units
  else if (AvgFunc .eq. 'hda') then
    ! all z points
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = Nz
    TserAvg%units = Var%units
  else if (AvgFunc .eq. 'hist') then
    ! put bin values in the x dimension
    TserAvg%dims(1) = NumBins
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = Var%units
  else if (AvgFunc .eq. 'horiz_ke') then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = 'Joules'
  else if (AvgFunc .eq. 'storm_int') then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = 'int'
  endif

  allocate(TserAvg%vdata(TserAvg%dims(1)*TserAvg%dims(2)*TserAvg%dims(3)))

  ! If doing max_azwind, set up the output var for the radius of max wind
  if (AvgFunc .eq. 'max_azwind') then
    RadMaxWind%vname = 'rmw'
    RadMaxWind%descrip = 'time series averaged ' // trim(AvgFunc) 
    RadMaxWind%units = 'km'
    RadMaxWind%ndims = 3 
    RadMaxWind%dimnames(1) = 'x' 
    RadMaxWind%dimnames(2) = 'y' 
    RadMaxWind%dimnames(3) = 'z' 
    RadMaxWind%dims(1) = 1
    RadMaxWind%dims(2) = 1
    RadMaxWind%dims(3) = 1
    allocate(RadMaxWind%vdata(RadMaxWind%dims(1)*RadMaxWind%dims(2)*RadMaxWind%dims(3)))
  endif

  ! Report the dimensions
  write (*,*) 'Output variable information:'
  write (*,*) '  Name: ', trim(TserAvg%vname)
  write (*,*) '  Units: ', trim(TserAvg%units)
  write (*,*) '  Description: ', trim(TserAvg%descrip)
  write (*,*) '  Number of dimensions: ', TserAvg%ndims
  write (*,*) '  Dimension sizes:'
  do id = 1, TserAvg%ndims
    write (*,*), '    ', trim(TserAvg%dimnames(id)), ': ', TserAvg%dims(id)
  enddo
  write (*,*) ''
  write (*,*) '  Number of time steps: ', Nt
  write (*,*) ''
  flush(6)

  ! Read in the input coordinates
  if (AvgFunc .eq. 'max_azwind') then
    InCoordFile = trim(AzWindFile)
  else if (AvgFunc .eq. 'horiz_ke') then
    InCoordFile = trim(DensFile)
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
           (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist')) then
    InCoordFile = trim(VarFile)
  else if (AvgFunc .eq. 'storm_int') then
    InCoordFile = trim(Speed10mFile)
  endif

  InXcoords%vname = 'x_coords'
  call rhdf5_read_init(InCoordFile, InXcoords)
  call rhdf5_read(InCoordFile, InXcoords)

  InYcoords%vname = 'y_coords'
  call rhdf5_read_init(InCoordFile, InYcoords)
  call rhdf5_read(InCoordFile, InYcoords)

  if (AvgFunc .eq. 'storm_int') then
    ! fake it for the storm intensity metric
    InZcoords%vname = 'z_coords'
    InZcoords%ndims = 1
    InZcoords%dims(1) = 1
    InZcoords%dimnames(1) = 'z'
    InZcoords%units = 'meter'
    InZcoords%descrip = 'sigma-z'
    allocate(InZcoords%vdata(1))
    InZcoords%vdata(1) = 10.0
  else
    InZcoords%vname = 'z_coords'
    call rhdf5_read_init(InCoordFile, InZcoords)
    call rhdf5_read(InCoordFile, InZcoords)
  endif

  InTcoords%vname = 't_coords'
  call rhdf5_read_init(InCoordFile, InTcoords)
  call rhdf5_read(InCoordFile, InTcoords)

  ! Need to get delta for x and y (in meters) and the z heights (in meters)
  ! for DoHorizKe to be able to calculate volume * density -> mass
  allocate(InXcoordsKm(Nx))
  allocate(InYcoordsKm(Ny))
  if (AvgFunc .eq. 'max_azwind') then
    ! x,y coords are in meters, convert to km
    do ix = 1, Nx
      InXcoordsKm(ix) = InXcoords%vdata(ix) / 1000.0
    enddo
    do iy = 1, Ny
      InYcoordsKm(iy) = InYcoords%vdata(iy) / 1000.0
    enddo
  else
    ! x,y coords are in degrees lon,lat respectively
    call ConvertGridCoords(Nx, Ny, Nz, InXcoords%vdata, InYcoords%vdata, InXcoordsKm, InYcoordsKm)
  endif

  DeltaX = InXcoordsKm(2) - InXcoordsKm(1)
  DeltaY = InYcoordsKm(2) - InYcoordsKm(1)

  write (*,*) 'Horizontal grid info:'
  write (*,*) '  X range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', InXcoords%vdata(1), InXcoords%vdata(Nx), InXcoordsKm(1), InXcoordsKm(Nx)
  write (*,*) '  Y range (min lat, max lat) --> (min y, max y): '
  write (*,*) '    ', InYcoords%vdata(1), InYcoords%vdata(Ny), InYcoordsKm(1), InYcoordsKm(Ny)
  write (*,*) ''
  write (*,*) 'Vertical grid info:'
  do iz = 1, Nz
    write (*,*) '  ', iz, ' --> ', InZcoords%vdata(iz)
  end do
  write (*,*) ''
  flush(6)

  ! If doing histogram, calculate the bin values
  if (AvgFunc .eq. 'hist') then
    allocate(Bins(NumBins))
    Bins(1) = BinStart
    do ib = 2, NumBins
      Bins(ib) = Bins(ib-1) + BinInc
    enddo
  endif

  ! Prepare the output coordinates
  !
  ! Create dummy coordinates (for GRADS sake) and write out the
  ! time series as a 4D var, (x,y,z,t) regardless of how many
  ! true dimensions exist in the output.
  OutXcoords%vname = 'x_coords'
  OutXcoords%ndims = 1
  OutXcoords%dimnames(1) = 'x'
  OutXcoords%units = 'degrees_east'
  OutXcoords%descrip = 'longitude'
  if (AvgFunc .eq. 'hist') then
    OutXcoords%dims(1) = NumBins
    allocate(OutXcoords%vdata(NumBins))
    do ib = 1, NumBins
      OutXcoords%vdata(ib) = Bins(ib)
    enddo
  else
    OutXcoords%dims(1) = 1
    allocate(OutXcoords%vdata(1))
    OutXcoords%vdata(1) = 1.0
  endif
  
  OutYcoords%vname = 'y_coords'
  OutYcoords%ndims = 1
  OutYcoords%dimnames(1) = 'y'
  OutYcoords%units = 'degrees_north'
  OutYcoords%descrip = 'latitude'
  OutYcoords%dims(1) = 1
  allocate(OutYcoords%vdata(1))
  OutYcoords%vdata(1) = 1.0
  
  OutZcoords%vname = 'z_coords'
  OutZcoords%ndims = 1
  OutZcoords%dimnames(1) = 'z'
  OutZcoords%units = 'meter'
  OutZcoords%descrip = 'sigma-z'
  if (AvgFunc .eq. 'hda') then
    OutZcoords%dims(1) = Nz
    allocate(OutZcoords%vdata(Nz))
    do iz = 1, Nz
      OutZcoords%vdata(iz) = InZcoords%vdata(iz)
    enddo
  else
    OutZcoords%dims(1) = 1
    allocate(OutZcoords%vdata(1))
    OutZcoords%vdata(1) = 1.0
  endif

  ! Perform the averaging function.
  rh5f_facc = 'W'
  call rhdf5_open_file(OutFile, rh5f_facc, 1, rh5f_out)

  rh5f_facc = 'R'
  if (UseFilter) then
    call rhdf5_open_file(FilterFile, rh5f_facc, 1, rh5f_filter)
  endif
  if (AvgFunc .eq. 'max_azwind') then
    call rhdf5_open_file(AzWindFile, rh5f_facc, 1, rh5f_azwind)
  else if (AvgFunc .eq. 'horiz_ke') then
    call rhdf5_open_file(DensFile, rh5f_facc, 1, rh5f_dens)
    call rhdf5_open_file(Ufile, rh5f_facc, 1, rh5f_u)
    call rhdf5_open_file(Vfile, rh5f_facc, 1, rh5f_v)
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
           (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist'))  then
    call rhdf5_open_file(VarFile, rh5f_facc, 1, rh5f_var)
  else if (AvgFunc .eq. 'storm_int') then
    call rhdf5_open_file(Speed10mFile, rh5f_facc, 1, rh5f_speed10m)
  endif

  do it = 1, Nt
    ! if using a filter read in the data
    if (UseFilter) then
      call rhdf5_read_variable(rh5f_filter, Filter%vname, Filter%ndims, it, Filter%dims, rdata=Filter%vdata)
    endif

    ! do the averaging function
    if (AvgFunc .eq. 'max_azwind') then
      call rhdf5_read_variable(rh5f_azwind, AzWind%vname, AzWind%ndims, it, AzWind%dims, rdata=AzWind%vdata)
      call DoMaxAzWind(Nx, Ny, Nz, AzWind%vdata, UndefVal, TserAvg%vdata(1), i_rmw)
      RadMaxWind%vdata(1) = InXcoordsKm(i_rmw)
      deallocate(AzWind%vdata)
    else if (AvgFunc .eq. 'horiz_ke') then
      call rhdf5_read_variable(rh5f_dens, Dens%vname, Dens%ndims, it, Dens%dims, rdata=Dens%vdata)
      call rhdf5_read_variable(rh5f_u, U%vname, U%ndims, it, U%dims, rdata=U%vdata)
      call rhdf5_read_variable(rh5f_v, V%vname, V%ndims, it, V%dims, rdata=V%vdata)

      call DoHorizKe(Nx, Ny, Nz, Dens%vdata, U%vdata, V%vdata, Filter%vdata, DeltaX, DeltaY, InZcoords%vdata, TserAvg%vdata(1))

      deallocate(Dens%vdata)
      deallocate(U%vdata)
      deallocate(V%vdata)
    else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
             (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist'))  then
      call rhdf5_read_variable(rh5f_var, Var%vname, Var%ndims, it, Var%dims, rdata=Var%vdata)

      if (AvgFunc .eq. 'min') then
        call DoMin(Nx, Ny, Nz, Filter%dims(3), Var%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata(1))
      else if (AvgFunc .eq. 'max') then
        call DoMax(Nx, Ny, Nz, Filter%dims(3), Var%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata(1))
      else if (AvgFunc .eq. 'hda') then
        call DoHda(Nx, Ny, Nz, Filter%dims(3), Var%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata)
      else if (AvgFunc .eq. 'hist') then
        call DoHist(Nx, Ny, Nz, Filter%dims(3), NumBins, Var%vdata, Filter%vdata, UseFilter, UndefVal, Bins, TserAvg%vdata)
      endif

      deallocate(Var%vdata)
    else if (AvgFunc .eq. 'storm_int') then
      call rhdf5_read_variable(rh5f_speed10m, Speed10m%vname, Speed10m%ndims, it, Speed10m%dims, rdata=Speed10m%vdata)
      call DoStormInt(Nx, Ny, Nz, Speed10m%vdata, Filter%vdata, TserAvg%vdata(1))
      deallocate(Speed10m%vdata)
    endif

    ! if using a filter, deallocate the space for the next time around
    if (UseFilter) then
      deallocate(Filter%vdata)
    endif

    ! write out the averaged data
    call rhdf5_write_variable(rh5f_out, TserAvg%vname, TserAvg%ndims, it, TserAvg%dims, &
       TserAvg%units, TserAvg%descrip, TserAvg%dimnames, rdata=TserAvg%vdata)
    if (AvgFunc .eq. 'max_azwind') then
      call rhdf5_write_variable(rh5f_out, RadMaxWind%vname, RadMaxWind%ndims, it, RadMaxWind%dims, &
         RadMaxWind%units, RadMaxWind%descrip, RadMaxWind%dimnames, rdata=RadMaxWind%vdata)
    endif

    ! print a message for the user on longer jobs so that it can be
    ! seen that progress is being made
    if (modulo(it,100) .eq. 0) then
      write (*,*) 'Working: Number of time steps processed so far: ', it
    endif
  enddo

  rh5f_facc = 'W'
  call rhdf5_close_file(rh5f_out)

  rh5f_facc = 'R'
  if (UseFilter) then
    call rhdf5_close_file(rh5f_filter)
  endif
  if (AvgFunc .eq. 'max_azwind') then
    call rhdf5_close_file(rh5f_azwind)
  else if (AvgFunc .eq. 'horiz_ke') then
    call rhdf5_close_file(rh5f_dens)
    call rhdf5_close_file(rh5f_u)
    call rhdf5_close_file(rh5f_v)
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
           (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist'))  then
    call rhdf5_close_file(rh5f_var)
  else if (AvgFunc .eq. 'storm_int') then
    call rhdf5_close_file(rh5f_speed10m)
  endif
 
  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''

  ! Finish off output file
  ! write out the coordinate data
  call rhdf5_write(OutFile, OutXcoords, 1)
  call rhdf5_write(OutFile, OutYcoords, 1)
  call rhdf5_write(OutFile, OutZcoords, 1)
  call rhdf5_write(OutFile, InTcoords, 1)

  ! If doing hist function, write out the input dimension sizes
  ! for downstream analyses. Eg. if you want to do fractional
  ! area calculations then the counts in the histogram do not
  ! tell you how many total points are in the domain.
  if (AvgFunc .eq. 'hist') then
    ! common settings
    OrigDimSize%ndims = 1
    OrigDimSize%dims(1) = 1
    OrigDimSize%units = 'number'
    allocate(OrigDimSize%vdata(1))

    ! X
    OrigDimSize%vname = 'Nx'
    OrigDimSize%dimnames(1) = 'x'
    OrigDimSize%descrip = 'number of domain x points'
    OrigDimSize%vdata(1) = float(Nx)
    call rhdf5_write(OutFile, OrigDimSize, 1)

    ! Y
    OrigDimSize%vname = 'Ny'
    OrigDimSize%dimnames(1) = 'y'
    OrigDimSize%descrip = 'number of domain y points'
    OrigDimSize%vdata(1) = float(Ny)
    call rhdf5_write(OutFile, OrigDimSize, 1)

    ! Z
    OrigDimSize%vname = 'Nz'
    OrigDimSize%dimnames(1) = 'z'
    OrigDimSize%descrip = 'number of domain z points'
    OrigDimSize%vdata(1) = float(Nz)
    call rhdf5_write(OutFile, OrigDimSize, 1)

    ! T
    OrigDimSize%vname = 'Nt'
    OrigDimSize%dimnames(1) = 't'
    OrigDimSize%descrip = 'number of domain t points'
    OrigDimSize%vdata(1) = float(Nt)
    call rhdf5_write(OutFile, OrigDimSize, 1)

    deallocate(OrigDimSize%vdata)
  endif

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(OutFile, OutXcoords, 'x')
  call rhdf5_set_dimension(OutFile, OutYcoords, 'y')
  call rhdf5_set_dimension(OutFile, OutZcoords, 'z')
  call rhdf5_set_dimension(OutFile, InTcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(OutFile, TserAvg)
  if (AvgFunc .eq. 'max_azwind') then
    call rhdf5_attach_dimensions(OutFile, RadMaxWind)
  endif
  
  stop

contains
!**********************************************************************
! SUBROUTINES
!**********************************************************************

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the command line arguments
!
subroutine GetMyArgs(InDir, InSuffix, OutFile, AvgFunc, FilterFile, UseFilter)
  implicit none

  character (len=*) :: InDir, InSuffix, OutFile, AvgFunc, FilterFile
  logical :: UseFilter

  integer :: iargc
  character (len=128) :: arg
  integer :: Nitems

  logical :: BadArgs

  if (iargc() .ne. 5) then
    write (*,*) 'ERROR: must supply exactly 5 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: tsavg <in_dir> <in_suffix> <out_file> <avg_function> <filter_file>'
    write (*,*) '        <in_dir>: directory where input files live'
    write (*,*) '        <in_suffix>: suffix on input file names'
    write (*,*) '        <out_file>: name of output file, HDF5 format'
    write (*,*) '        <avg_function>: averaging function to use on input data'
    write (*,*) '            horiz_ke -> total kinetic energy form horizontal winds'
    write (*,*) '            storm_int -> storm intensity metric from horizontal wind speeds'
    write (*,*) '            max_azwind -> max value of azimuthially averaged wind'
    write (*,*) '            min:<var>:<file>:<dim> -> domain minimum for <var>'
    write (*,*) '              <var>: revu var name inside the file'
    write (*,*) '              <file>: prefix for the revu file'
    write (*,*) '              <dim>: dimensionality of variable, either "2d" or "3d"'
    write (*,*) '            max:<var>:<file>:<dim> -> domain maximum for <var>'
    write (*,*) '              <var>: revu var name inside the file'
    write (*,*) '              <file>: prefix for the revu file'
    write (*,*) '              <dim>: dimensionality of variable, either "2d" or "3d"'
    write (*,*) '            hda:<var>:<file>:<dim> -> horizontal domain average at each z level for <var>'
    write (*,*) '              <var>: revu var name inside the file'
    write (*,*) '              <file>: prefix for the revu file'
    write (*,*) '              <dim>: dimensionality of variable, either "2d" or "3d"'
    write (*,*) '            hist:<var>:<file>:<dim>:<num_bins>:<bin_start>:<bin_inc>'
    write (*,*) '              <var>,<file>,<dim> same as for hda'
    write (*,*) '              <num_bins>: number of bins'
    write (*,*) '              <bin_start>: starting value for bins'
    write (*,*) '              <bin_inc>: delta between bins'
    write (*,*) '        <filter_file>: file containing the filter mask'
    stop
  end if

  call getarg(1, InDir)
  call getarg(2, InSuffix)
  call getarg(3, OutFile)
  call getarg(4, AvgFunc)
  call getarg(5, FilterFile)

  if ((FilterFile .eq. 'none') .or. (FilterFile .eq. 'NONE')) then
    UseFilter = .false.
  else
    UseFilter = .true.
  endif

  BadArgs = .false.

  if ((AvgFunc .ne. 'horiz_ke')       .and. &
      (AvgFunc(1:4) .ne. 'min:')  .and. &
      (AvgFunc(1:4) .ne. 'max:')  .and. &
      (AvgFunc(1:4) .ne. 'hda:')  .and. &
      (AvgFunc(1:5) .ne. 'hist:')  .and. &
      (AvgFunc .ne. 'storm_int')  .and. &
      (AvgFunc .ne. 'max_azwind')) then
    write (*,*) 'ERROR: <avg_function> must be one of:'
    write (*,*) '          horiz_ke'
    write (*,*) '          min:<var>:<file>:<dim>'
    write (*,*) '          max:<var>:<file>:<dim>'
    write (*,*) '          hda:<var>:<file>:<dim>'
    write (*,*) '          hist:<var>:<file>:<dim>:<num_bins>:<bin_start>:<bin_inc>'
    write (*,*) '          storm_int'
    write (*,*) '          max_azwind'
    write (*,*) ''
    BadArgs = .true.
  end if

  if (BadArgs) then
    stop
  end if

  return
end subroutine GetMyArgs

!************************************************************************************
! DoMaxAzWind()
!
! This subroutine will find the maximum wind speed in AzWind and copy that
! to TserAvg. The intent of this metric is to mimic the Saffir-Simpson scale
! which uses the average 10m tangential wind speed that has been sustained
! over 10 minutes.
!
! Typically, don't have enough RAMS files to do 10m time averaging so as a
! proxy use the azimuthally averaged tangetial wind speed.
!
! Keep track of the index into the radius values that produced the maximum
! wind speed. Pass this back to the caller so that the caller can look
! up the radius value. Radius is the first dimension, x.
!
subroutine DoMaxAzWind(Nx, Ny, Nz, AzWind, UndefVal, AzWindMax, Irmw)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: AzWind
  real :: UndefVal
  real :: AzWindMax
  integer :: Irmw

  integer :: ix,iy,iz

  ! dimension order is: x,y,z

  ! RAMS places the first z level below the surface so use the second level
  ! to approximate the 10m winds.

  AzWindMax = 0.0
  Irmw = 0
  iz = 2
  do iy = 1, Ny
    do ix = 1, Nx
      if ((AzWind(ix,iy,iz) .gt. AzWindMax) .and. (anint(AzWind(ix,iy,iz)) .ne. UndefVal)) then
        AzWindMax = AzWind(ix,iy,iz)
        Irmw = ix
      endif
    end do
  end do

  return
end subroutine DoMaxAzWind

!**************************************************************************************
! DoHorizKe()
!
! This routine will calculate the total kinetic energy over the given cylindrical
! volume. Do not want average since we want the size of the storm reflected in
! this diagnostic.

subroutine DoHorizKe(Nx, Ny, Nz, Dens, U, V, Filter, DeltaX, DeltaY, Zcoords, TserAvg)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: Dens, U, V, Filter
  real :: DeltaX, DeltaY
  real, dimension(Nz) :: Zcoords
  real :: TserAvg

  integer ix,iy,iz
  integer NumPoints
  real SumKe, CurrKe, LevThickness

  ! KE is 1/2 * m * v^2
  !   - calculate this a every point inside the defined cylindrical volume
  !   - get m by density * volume
  !   - v^2 is based on horiz velocity so is equal to U^2 + V^2
  !
  ! Zcoords are technically the center points of the levels in the RAMS simulation. Since we don't have
  ! the level definition from the RAMS runs here, just use the difference from the i+1st z coord minus the
  ! ith z coord to approzimate the ith level thickness. This will be close enough for the measurement.
  !
  ! The integrated kinetic engery measurement done on hurricanes is taken over the domain on the 10m
  ! level. Just use the first model level above the surface (z == 2) for this measurement which will
  ! be close enough.

  SumKe = 0.0
  NumPoints = 0

  iz = 2

  do iy = 1, Ny
    do ix = 1, Nx
      if (anint(Filter(ix,iy,iz)) .eq. 1.0) then
        if (iz .eq. Nz) then
          ! Use the level below for this case (since no level above)
          LevThickness = Zcoords(iz) - Zcoords(iz-1)
        else
          LevThickness = Zcoords(iz+1) - Zcoords(iz)
        end if
        CurrKe = 0.5 * DeltaX * DeltaY * LevThickness * Dens(ix,iy,iz) * (U(ix,iy,iz)**2 + V(ix,iy,iz)**2)
        SumKe = SumKe + CurrKe
        NumPoints = NumPoints + 1
      endif
    enddo
  enddo

  TserAvg = SumKe;
  
  return
end subroutine DoHorizKe

!**************************************************************************************
! DoHda()
!
! This routine will do horizontal domain average for all z levels.
!
subroutine DoHda(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomAvg)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(Nz) :: DomAvg

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint
  integer :: NumPoints

  do iz = 1, Nz
    DomAvg(iz) = 0.0
    NumPoints = 0

    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          DomAvg(iz) = DomAvg(iz) + Var(ix,iy,iz)
          NumPoints = NumPoints + 1
        endif
      enddo
    enddo

    if (NumPoints .eq. 0) then
      DomAvg(iz) = UndefVal
    else
      DomAvg(iz) = DomAvg(iz) / float(NumPoints)
    endif
  enddo

  return
end subroutine DoHda

!**************************************************************************************
! DoMin()
!
! This routine will return the domain minimum value. Values equal to UndefVal will
! be ignored.
!
subroutine DoMin(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomMin)
  implicit none

  real, parameter :: BigPosNum = 10.0e50

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real :: DomMin

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint

  DomMin = BigPosNum
  do iz = 1, Nz
    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          if (Var(ix,iy,iz) .lt. DomMin) then
            DomMin = Var(ix,iy,iz)
          endif
        endif
      enddo
    enddo
  enddo

  if (DomMin .eq. BigPosNum) then
    ! all entries were UndefVal
    DomMin = UndefVal
  endif

  return
end subroutine DoMin

!**************************************************************************************
! DoMax()
!
! This routine will return the domain maximum value. Values equal to UndefVal will
! be ignored.
!
subroutine DoMax(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomMax)
  implicit none

  real, parameter :: BigNegNum = -10.0e50

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real :: DomMax

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint

  DomMax = BigNegNum
  do iz = 1, Nz
    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          if (Var(ix,iy,iz) .gt. DomMax) then
            DomMax = Var(ix,iy,iz)
          endif
        endif
      enddo
    enddo
  enddo

  if (DomMax .eq. BigNegNum) then
    ! all entries were UndefVal
    DomMax = UndefVal
  endif

  return
end subroutine DoMax

!**************************************************************************************
! DoHist()
!
! This routine will do histogram binning over all of the domain.
!

subroutine DoHist(Nx, Ny, Nz, FilterNz, Nb, Var, Filter, UseFilter, UndefVal, Bins, Counts)
  implicit none

  integer :: Nx, Ny, Nz, Nb, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(Nb) :: Bins, Counts

  integer :: ib, ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint

  ! Emulate matlab histc command, ie the values in Bins are treated as the edges
  ! of the bin ranges. Do not count any values from Var that fall outside the
  ! bin range (< Bins(1) or > Bins(Nb)). The count is incremented for bin ib when:
  !
  !   Bins(ib) <= Var(ix,iy,iz) < Bins(ib+1)
  !
  ! The count in Bin(Nb), last bin, is incremented when:
  !
  !   Var(ix,iy,iz) == Bins(ib)
  !
  do ib = 1, Nb
    Counts(ib) = 0.0
  enddo

  ! Build the histogram (counts)
  do iz = 1, Nz
    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 1, Ny
      do ix = 1, Nx
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          ! Check all bins except the last.
          !
          ! Exiting out of the loop when finding the bin will help a lot when the
          ! distribution of values is biased toward smaller values. After exiting
          ! out of the loop you can either just check the last bin (which will be wasted)
          ! or put in a logical variable and check that variable saying you can skip the
          ! check of the last bin. Since you would have to check the logical variable and
          ! the last bin every time you might as well just check the last bin instead.
          do ib = 1, Nb-1
             if ((Bins(ib) .le. Var(ix,iy,iz)) .and. (Var(ix,iy,iz) .lt. Bins(ib+1))) then
                Counts(ib) = Counts(ib) + 1.0
                exit
             endif
          enddo
          ! check the last bin
          if (Bins(Nb) .eq. Var(ix,iy,iz)) then
            Counts(Nb) = Counts(Nb) + 1.0
          endif
        endif
      enddo
    enddo
  enddo

  return
end subroutine DoHist

!**************************************************************************************
! DoStormInt()
!
! This routine will calculate a storm intensity metric based on the horizontal wind
! speeds. The metric is based on the hurricane categories from the Saffir-Simpson
! and their associated wind speeds.
!
!  convert the surface wind data into weighted coverage data
!  run through each time point and each x,y location of the sfc wind data
!  weight each count by the Saffir-Simpson category number
!
!   Wind speed   Weight
!    <  33 m/s -> 0  (not hurricane force)
!    33-42 m/s -> 1  (Category 1 hurricane)
!    43-49 m/s -> 2  (etc.)
!    50-58 m/s -> 3
!    59-69 m/s -> 4
!    >= 70 m/s -> 5
!

subroutine DoStormInt(Nx, Ny, Nz, Speed10m, Filter, TserAvg)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: Filter
  real, dimension(Nx,Ny) :: Speed10m
  real :: TserAvg

  integer ix,iy,iz
  integer nCat0, nCat1, nCat2, nCat3, nCat4, nCat5, NumPoints
  real Wspeed, SiMetric

  nCat0 = 0
  nCat1 = 0
  nCat2 = 0
  nCat3 = 0
  nCat4 = 0
  nCat5 = 0

  iz = 2 ! use next to bottom layer in the filter

  do iy = 1, Ny
    do ix = 1, Nx
      if (anint(Filter(ix,iy,iz)) .eq. 1.0) then
        ! Count up the number of grid points with wind speeds fitting each of the
        ! Saffir-Simpson categories. Then form the metric by weighting each category
        ! count.
        Wspeed = Speed10m(ix,iy)
        if (Wspeed .ge. 70.0) then
          nCat5 = nCat5 + 1
        else if (Wspeed .ge. 59.0) then
          nCat4 = nCat4 + 1
        else if (Wspeed .ge. 50.0) then
          nCat3 = nCat3 + 1
        else if (Wspeed .ge. 43.0) then
          nCat2 = nCat2 + 1
        else if (Wspeed .ge. 33.0) then
          nCat1 = nCat1 + 1
        else
          nCat0 = nCat0 + 1
        endif
      endif
    enddo
  enddo

  !Linear weighting
  NumPoints = nCat0 + nCat1 + nCat2 + nCat3 + nCat4 + nCat5
  SiMetric = float(nCat1) + (float(nCat2)*2.0) + (float(nCat3)*3.0) + (float(nCat4)*4.0) + (float(nCat5)*5.0)

  if (NumPoints .eq. 0) then
    TserAvg = 0.0
  else
    TserAvg = SiMetric / float(NumPoints)
  end if
  
  return
end subroutine DoStormInt


!!! !*****************************************************************************
!!! ! DoCloud()
!!! !
!!! ! This subroutine will perform the cloud droplet total mass time series
!!! ! averaging. Can select between supercooled or warm rain droplets.
!!! !
!!! 
!!! subroutine DoCloud(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, Dens, CintLiq, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   logical :: DoSc
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
!!!   type (GradsVar) :: Cloud, TempC, Dens, CintLiq, TserAvg
!!!   integer, dimension(1:Cloud%Nt) :: StmIx, StmIy
!!!   real, dimension(1:Cloud%Nx) :: Xcoords
!!!   real, dimension(1:Cloud%Ny) :: Ycoords
!!!   real, dimension(1:Cloud%Nz) :: Zcoords
!!! 
!!!   real, dimension(0:Cloud%Nz) :: ZmHeights
!!!   integer :: ix, iy, iz, it, NumPoints
!!! 
!!!   call SetZmHeights (Cloud%Nz, ZmHeights)
!!! 
!!!   ! Convert the cloud mixing ratio to mass for each grid point. 
!!!   !
!!!   ! Mixing ratios in GRADS files are g/kg
!!!   ! Density is in kg/m**3
!!!   ! Heights are in m
!!!   ! So, express the mass value in g using the formula
!!!   !   (mix ratio) * (density) * (layer thickness) * (layer horiz area)
!!!   ! The layer thickness for layer k is: ZmHeights(k) - ZmHeights(k-1)
!!! 
!!!   do it = 1, Cloud%Nt
!!!     ! Sum up the cloud droplet mass over the specified radial band. Only include the
!!!     ! grid points where tempc is 0 or less (supercooled)
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, Cloud%Nz
!!!       do ix = 1, Cloud%Nx
!!!         do iy = 1, Cloud%Ny
!!!           if ((InsideCylVol(Cloud%Nx, Cloud%Ny, Cloud%Nz, Cloud%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords))  .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
!!!             if (DoSc) then
!!!               if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
!!!                 TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + (Cloud%Vdata(it,iz,ix,iy) * Dens%Vdata(it,iz,ix,iy) * (ZmHeights(iz)-ZmHeights(iz-1)))
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             else
!!!               if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
!!!                 TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + (Cloud%Vdata(it,iz,ix,iy) * Dens%Vdata(it,iz,ix,iy) * (ZmHeights(iz)-ZmHeights(iz-1)))
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     ! At this point TserAvg%Vdata(it,1,1,1) holds g/m**2, multiply by grid cell horizontal area. Note this assumes
!!!     ! each grid cell has the same horizontal area.
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) * DeltaX * DeltaY
!!!     if (NumPoints .eq. 0) then
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       write (*,*) 'ScCloud: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     endif
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoWup()
!!! !
!!! ! This subroutine will do the average vertical velocity in regions of significant
!!! ! updrafts.
!!! !
!!! 
!!! subroutine DoWup(DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, W, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: W, TserAvg
!!!   integer, dimension(1:W%Nt) :: StmIx, StmIy
!!!   real, dimension(1:W%Nx) :: Xcoords
!!!   real, dimension(1:W%Ny) :: Ycoords
!!!   real, dimension(1:W%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it, NumPoints
!!! 
!!!   do it = 1, W%Nt
!!!     ! Average w over regions where significant updrafts occur
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, W%Nz
!!!       do ix = 1, W%Nx
!!!         do iy = 1, W%Ny
!!!           if (InsideCylVol(W%Nx, W%Ny, W%Nz, W%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             if (W%Vdata(it,iz,ix,iy) .ge. Wthreshold) then
!!!               TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + W%Vdata(it,iz,ix,iy)
!!!               NumPoints = NumPoints + 1
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) / float(NumPoints)
!!!       write (*,*) 'Wup: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!   end do
!!! 
!!! end subroutine
!!! 
!!! 
!!! !************************************************************************************
!!! ! DoCloudDiam()
!!! !
!!! ! This routine will calculate an averaged cloud droplet diameter. Can select between
!!! ! warm rain droplets or supercooled droplets.
!!! 
!!! subroutine DoCloudDiam(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, CloudDiam, CintLiq, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   logical :: DoSc
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
!!!   type (GradsVar) :: Cloud, TempC, CloudDiam, CintLiq, TserAvg
!!!   integer, dimension(1:CloudDiam%Nt) :: StmIx, StmIy
!!!   real, dimension(1:CloudDiam%Nx) :: Xcoords
!!!   real, dimension(1:CloudDiam%Ny) :: Ycoords
!!!   real, dimension(1:CloudDiam%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it, NumPoints
!!!   real SumQ, SumQD
!!!   real MaxQ, Climit, SumD
!!! 
!!!   do it = 1, CloudDiam%Nt
!!!     ! Calculate a mass-weighted mean diameter for supercooled cloud droplets.
!!!     ! 
!!!     !    Mean diameter (TserAvg value) = Sum(cloud * cloud_d) / Sum(cloud)
!!!     !    where cloud and cloud_d are only included in the sum when tempc is <= 0
!!!     !
!!!     ! Find the max Q (mass) value and use it to filter data - select data points if
!!!     ! the Q is within 20% of the max Q (.2 to 1.0). Then form the average of the diameters
!!!     ! of the selected points.
!!!     SumQD = 0.0
!!!     SumQ = 0.0
!!!     SumD = 0.0
!!!     NumPoints = 0
!!!     do iz = 1, CloudDiam%Nz
!!!       do ix = 1, CloudDiam%Nx
!!!         do iy = 1, CloudDiam%Ny
!!!           if ((InsideCylVol(CloudDiam%Nx, CloudDiam%Ny, CloudDiam%Nz, CloudDiam%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
!!!             if (DoSc) then
!!!               if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
!!!                  SumQD = SumQD + (Cloud%Vdata(it,iz,ix,iy) * CloudDiam%Vdata(it,iz,ix,iy))
!!!                  SumQ = SumQ + Cloud%Vdata(it,iz,ix,iy)
!!!                  SumD = SumD + CloudDiam%Vdata(it,iz,ix,iy)
!!!                  NumPoints = NumPoints + 1
!!!               end if
!!!             else
!!!               if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
!!!                  SumQD = SumQD + (Cloud%Vdata(it,iz,ix,iy) * CloudDiam%Vdata(it,iz,ix,iy))
!!!                  SumQ = SumQ + Cloud%Vdata(it,iz,ix,iy)
!!!                  SumD = SumD + CloudDiam%Vdata(it,iz,ix,iy)
!!!                  NumPoints = NumPoints + 1
!!!               end if
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (SumQ .eq. 0.0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = SumQD / SumQ
!!!       write (*,*) 'ScCloudDiam: Ts:', it, ', NumPoints: ', NumPoints
!!!     end if
!!! 
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoCloudConc()
!!! !
!!! ! This subroutine will do the average cloud droplet concentration
!!! 
!!! subroutine DoCloudConc(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, DoSc, TempC, CloudConc, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: CloudConc, TempC, TserAvg
!!!   integer, dimension(1:CloudConc%Nt) :: StmIx, StmIy
!!!   real, dimension(1:CloudConc%Nx) :: Xcoords
!!!   real, dimension(1:CloudConc%Ny) :: Ycoords
!!!   real, dimension(1:CloudConc%Nz) :: Zcoords
!!!   logical :: DoSc
!!! 
!!!   integer ix,iy,iz,it
!!!   integer NumPoints
!!!   real SumCloudConc
!!! 
!!!   ! Calculate the average cloud droplet concentration near the eyewall region.
!!!   ! 
!!!   ! Call the cloud base to be between 1000m and 1200m. The z level corresponding to
!!!   ! that is z = 4 which is 1138m. Cover from surface to the cloud base level. This
!!!   ! results in using z = 1 to z = 4.
!!! 
!!!   do it = 1, CloudConc%Nt
!!!     SumCloudConc = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, CloudConc%Nz
!!!       do ix = 1, CloudConc%Nx
!!!         do iy = 1, CloudConc%Ny
!!!           if (InsideCylVol(CloudConc%Nx, CloudConc%Ny, CloudConc%Nz, CloudConc%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             if (DoSc) then
!!!               if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
!!!                 SumCloudConc = SumCloudConc + CloudConc%Vdata(it,iz,ix,iy)
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             else
!!!               if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
!!!                 SumCloudConc = SumCloudConc + CloudConc%Vdata(it,iz,ix,iy)
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = SumCloudConc / float(NumPoints)
!!!       write (*,*) 'CloudConc: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoPrecipR()
!!! !
!!! ! This subroutine will do the average cloud droplet concentration near the eyewall
!!! ! region.
!!! 
!!! subroutine DoPrecipR(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, PrecipR, CintLiq, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
!!!   ! PrecipR is 2D var
!!!   type (GradsVar) :: PrecipR, CintLiq, TserAvg
!!!   integer, dimension(1:PrecipR%Nt) :: StmIx, StmIy
!!!   real, dimension(1:PrecipR%Nx) :: Xcoords
!!!   real, dimension(1:PrecipR%Ny) :: Ycoords
!!!   real, dimension(1:PrecipR%Nz) :: Zcoords
!!! 
!!!   integer :: ix,iy,iz,it
!!!   integer :: NumPoints
!!!   real :: SumPrecip
!!! 
!!!   do it = 1, PrecipR%Nt
!!!     SumPrecip = 0.0
!!!     NumPoints = 0
!!! 
!!!     do ix = 1, PrecipR%Nx
!!!       do iy = 1, PrecipR%Ny
!!!         if (InsideCylVol(PrecipR%Nx, PrecipR%Ny, PrecipR%Nz, PrecipR%Nt, ix, iy, 1, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords) .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
!!!           SumPrecip = SumPrecip + PrecipR%Vdata(it,1,ix,iy)
!!!           NumPoints = NumPoints + 1
!!!         end if
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) '  WARNING: no data points selected for this time step'
!!!     else
!!!       ! At this point TserAvg%Vdata(it,1,1,1) holds mm/hr, multiply by grid cell horizontal area.
!!!       ! Note this assumes each grid cell has the same horizontal area. What this does
!!!       ! is convert mm/hr to kg/hr when assuming that the density of water is 1000kg/m**3.
!!!       !   mm/hr * m**2 * 1000 kg/m**3 * 0.001 m/mm -> kg/hr
!!!       TserAvg%Vdata(it,1,1,1) = (SumPrecip / float(NumPoints)) * DeltaX * DeltaY
!!!       write (*,*) 'Precip: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!     write (*,*) ''
!!!     flush(6)
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoCcnConc()
!!! !
!!! ! This subroutine will do the average CCN concentration.
!!! !
!!! 
!!! subroutine DoCcnConc(DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CcnConc, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: CcnConc, TserAvg
!!!   integer, dimension(1:CcnConc%Nt) :: StmIx, StmIy
!!!   real, dimension(1:CcnConc%Nx) :: Xcoords
!!!   real, dimension(1:CcnConc%Ny) :: Ycoords
!!!   real, dimension(1:CcnConc%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it, NumPoints
!!! 
!!!   do it = 1, CcnConc%Nt
!!!     ! Average w over regions where significant updrafts occur
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, CcnConc%Nz
!!!       do ix = 1, CcnConc%Nx
!!!         do iy = 1, CcnConc%Ny
!!!           if (InsideCylVol(CcnConc%Nx, CcnConc%Ny, CcnConc%Nz, CcnConc%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + CcnConc%Vdata(it,iz,ix,iy)
!!!             NumPoints = NumPoints + 1
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) / float(NumPoints)
!!!       write (*,*) 'Wup: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!   end do
!!! 
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoTestCvs()
!!! !
!!! ! This subroutine will perform a test on the cylindrical volume selection routine.
!!! ! Just runs through the entire grid and outputs a '1' when selection occurs otherwise
!!! ! outputs a '0'. Then view the result in grads and see if selection is correct.
!!! 
!!! subroutine DoTestCvs(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, TestSelect)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: TestSelect
!!!   integer, dimension(1:TestSelect%Nt) :: StmIx, StmIy
!!!   real, dimension(1:TestSelect%Nx) :: Xcoords
!!!   real, dimension(1:TestSelect%Ny) :: Ycoords
!!!   real, dimension(1:TestSelect%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it
!!!   integer NumPoints
!!! 
!!!   write (*,*) 'Testing cylindrical volume selection:'
!!!   do it = 1, TestSelect%Nt
!!!     NumPoints = 0
!!!     do iz = 1, TestSelect%Nz
!!!       do ix = 1, TestSelect%Nx
!!!         do iy = 1, TestSelect%Ny
!!!           if (InsideCylVol(TestSelect%Nx, TestSelect%Ny, TestSelect%Nz, TestSelect%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             TestSelect%Vdata(it,iz,ix,iy) = 1.0
!!!             NumPoints = NumPoints + 1
!!!           else
!!!             TestSelect%Vdata(it,iz,ix,iy) = 0.0
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!!     ! mark the storm center
!!!     TestSelect%Vdata(it,iz,StmIx(it),StmIy(it)) = 2.0
!!!     write (*,*) '  Timestep, Number of points selected: ', it, NumPoints
!!!   end do
!!! end subroutine
!!! 
!!! !******************************************************************************
!!! ! ConvertTinc()
!!! !
!!! ! This function will convert the time increment spec'd in the GRAD control
!!! ! file into a number of seconds.
!!! !
!!! 
!!! real function ConvertTinc(Tinc)
!!!   implicit none
!!! 
!!!   character (len=*) :: Tinc
!!! 
!!!   character (len=128) :: Tval, Tunits, Tfmt
!!!   integer :: i, Uloc, Tlen, Itval
!!!   
!!!   ! Walk through the Tinc string. Concatenate the numbers onto Tval and
!!!   ! the alaph characters onto Tunits. This algorithm assumes that GRADS
!!!   ! will use a format like <numeric_value><units> for Tinc where <units>
!!!   ! is a string like 'hr' or 'mn'.
!!!   Uloc = 0
!!!   Tlen = len_trim(Tinc)
!!!   do i = 1, Tlen
!!!     if ((Tinc(i:i) .ge. 'A') .and. (Tinc(i:i) .le. 'z')) then
!!!       if (Uloc .eq. 0) then
!!!         ! the first alpha character is the beginning of the spec for units
!!!         Uloc = i
!!!       end if
!!!     end if
!!!   end do
!!! 
!!!   write (Tfmt, '(a,i2.2,a,i2.2,a)') '(a', (Uloc-1), 'a' , ((Tlen-Uloc)+1), ')'
!!!   read (Tinc, Tfmt) Tval, Tunits
!!! 
!!!   read(Tval, '(i)') Itval
!!!   if (Tunits .eq. 'hr') then
!!!     ConvertTinc = float(Itval) * 3600.0
!!!   else
!!!     if (Tunits .eq. 'mn') then
!!!       ConvertTinc = float(Itval) * 60.0
!!!     else
!!!       ConvertTinc = float(Itval)
!!!     end if
!!!   end if
!!!   return
!!! end function
!!! 
!!! !******************************************************************************
!!! ! CalcRates()
!!! !
!!! ! This routine will calculate time derivatives of the input data using
!!! ! a centered difference method.
!!! !
!!! 
!!! subroutine CalcRates(DeltaT, TserAvg, Rates)
!!!   use gdata_utils
!!!   implicit none
!!! 
!!!   type (GradsVar) :: TserAvg, Rates
!!!   real :: DeltaT
!!! 
!!!   real :: f1, f2
!!!   integer :: ix, iy, iz, it
!!! 
!!!   ! use a centered difference, uses points at t-1, t and t+1
!!!   do ix = 1, TserAvg%Nx
!!!     do iy = 1, TserAvg%Ny
!!!       do iz = 1, TserAvg%Nz
!!!         do it = 2, TserAvg%Nt-1
!!!           f1 = (TserAvg%Vdata(it,1,1,1) + TserAvg%Vdata(it-1,1,1,1)) / 2.0
!!!           f2 = (TserAvg%Vdata(it+1,1,1,1) + TserAvg%Vdata(it,1,1,1)) / 2.0
!!! 
!!!           Rates%Vdata(it-1,1,1,1) = (f2 - f1) / DeltaT
!!!         end do
!!!       end do
!!!     end do
!!!   end do
!!! end subroutine

end program tsavg
