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

  type FileSpec
    character (len=RHDF5_MAX_STRING) :: fname
    character (len=RHDF5_MAX_STRING) :: vname
  endtype

  type ArgList
    character (len=LittleString) :: AvgFunc
    character (len=LittleString) :: VelInType
    real :: HkeZthick
    real :: HfracLimit
    real :: Zbot
    real :: Ztop
    real :: PcprateLimit
    logical :: SubtractSmotion
    type (FileSpec) :: Output
    type (FileSpec) :: Filter
    type (FileSpec) :: StormCenter
    type (FileSpec) :: V1bins
    type (FileSpec) :: V2bins
    type (FileSpec) :: Var1
    type (FileSpec) :: Var2
    type (FileSpec) :: Var3
    type (FileSpec) :: Dens
  endtype

  type (ArgList) :: Args

  ! Data arrays
  ! Dims: x, y, z, t
  type (Rhdf5Var) :: InXcoords, InYcoords, InZcoords, InTcoords
  type (Rhdf5Var) :: OutXcoords, OutYcoords, OutZcoords, OrigDimSize
  type (Rhdf5Var) :: Var1, Var2, Var3, Filter, StormCenter, Dens, TserAvg, RadMaxWind, Rad34ktWind, Rad50ktWind, Rad64ktWind
  type (Rhdf5Var) :: Ustorm, Vstorm
  character (len=RHDF5_MAX_STRING) :: InCoordFile
  character (len=RHDF5_MAX_STRING) :: rh5f_facc
  integer :: Kbot, Ktop

  logical :: UseFilter, UseStormCenter, UseVar2, UseVar3, UseDens
  logical :: Var1Is3d, Var2Is3d
  
  integer :: V1nbins, V2nbins
  real, dimension(:), allocatable :: V1bins, V2bins

  real, dimension(:,:), allocatable :: HorizSpeed

  integer :: rh5f_var1, rh5f_var2, rh5f_var3, rh5f_dens, rh5f_filter, rh5f_out

  integer :: id, ib, ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt, FilterNz
  integer :: Var1Nz, Var2Nz
  real :: DeltaX, DeltaY
  real, dimension(:), allocatable :: InXcoordsKm, InYcoordsKm
  logical :: BadDims

  ! Get the command line arguments
  call GetMyArgs(Args)

  ! For convenience
  UseFilter      = (trim(Args%Filter%fname) .ne. 'none')
  UseStormCenter = (trim(Args%StormCenter%fname) .ne. 'none')
  UseVar2        = (trim(Args%Var2%fname) .ne. 'none')
  UseVar3        = (trim(Args%Var3%fname) .ne. 'none')
  UseDens        = (trim(Args%Dens%fname) .ne. 'none')

  write (*,*) 'Time seris of average for RAMS data:'
  write (*,*) '  Output file:  ', trim(Args%Output%fname)
  write (*,*) '  Averaging function: ', trim(Args%AvgFunc)
  if (trim(Args%AvgFunc) .eq. 'horiz_ke') then
    write (*,*) '    Input Type: ', trim(Args%VelInType)
    write (*,*) '    Z thickness: ', Args%HkeZthick
  else if ((trim(Args%AvgFunc) .eq. 'min') .or. (trim(Args%AvgFunc) .eq. 'max') .or. (trim(Args%AvgFunc) .eq. 'hda') .or. (trim(Args%AvgFunc) .eq. 'turb_mmts')) then
    write (*,*) '    Variable name: ', trim(Args%Var1%vname)
    write (*,*) '    File name: ', trim(Args%Var1%fname)
  else if (trim(Args%AvgFunc) .eq. 'hist') then
    write (*,*) '    Variable name: ', trim(Args%Var1%vname)
    write (*,*) '    File name: ', trim(Args%Var1%fname)
    write (*,*) '    Bins File: ', trim(Args%V1bins%fname)
  else if (trim(Args%AvgFunc) .eq. 'hfrac') then
    write (*,*) '    Variable name: ', trim(Args%Var1%vname)
    write (*,*) '    File name: ', trim(Args%Var1%fname)
    write (*,*) '    Limit for fraction calculation: ', Args%HfracLimit 
  else if (trim(Args%AvgFunc) .eq. 'pop') then
    write (*,*) '    Precip rate variable name: ', trim(Args%Var1%vname)
    write (*,*) '    Precip rate File name: ', trim(Args%Var1%fname)
    write (*,*) '    Liquid water path variable name: ', trim(Args%Var2%vname)
    write (*,*) '    Liquid water path file name: ', trim(Args%Var2%fname)
    write (*,*) '    Lower troposhperic static stability variable name: ', trim(Args%Var3%vname)
    write (*,*) '    Lower troposhperic static stability file name: ', trim(Args%Var3%fname)
    write (*,*) '    Precip rate threshold: ', Args%PcprateLimit
    write (*,*) '    Liquid water path bins file: ', trim(Args%V1bins%fname)
    write (*,*) '    Lower tropospheric static stability bins file: ', trim(Args%V2bins%fname)
  else if (trim(Args%AvgFunc) .eq. 'hist2d') then
    write (*,*) '    X variable name: ', trim(Args%Var1%vname)
    write (*,*) '    X file name: ', trim(Args%Var1%fname)
    write (*,*) '    X variable bins file: ', trim(Args%V1bins%fname)
    write (*,*) '    Y variable name: ', trim(Args%Var2%vname)
    write (*,*) '    Y file name: ', trim(Args%Var2%vname)
    write (*,*) '    Y variable bins file: ', trim(Args%V2bins%fname)
  else if (trim(Args%AvgFunc) .eq. 'turb_cov') then
    write (*,*) '    X variable name: ', trim(Args%Var1%vname)
    write (*,*) '    X file name: ', trim(Args%Var1%fname)
    write (*,*) '    Y variable name: ', trim(Args%Var2%vname)
    write (*,*) '    Y file name: ', trim(Args%V2bins%fname)
  else if (trim(Args%AvgFunc) .eq. 'ltss') then
    write (*,*) '    Theta variable name: ', trim(Args%Var1%vname)
    write (*,*) '    Theta file name: ', trim(Args%Var1%fname)
    write (*,*) '    Z at bottom: ', Args%Zbot
    write (*,*) '    Z at top: ', Args%Ztop
  endif
  if (UseFilter) then
    write (*,*) '  Filter file: ', trim(Args%Filter%fname)
  endif
  if (UseStormCenter) then
    write (*,*) '  Storm center file: ', trim(Args%StormCenter%fname)
  endif
  write (*,*) ''
  flush(6)

  ! set up input variables
  Filter%vname = Args%Filter%vname
  StormCenter%vname = Args%StormCenter%vname

  Var1%vname = Args%Var1%vname
  Var2%vname = Args%Var2%vname
  Var3%vname = Args%Var3%vname

  Dens%vname = Args%Dens%vname

  ! Check that the dimensions are consistent between the variables needed for
  ! the selected averaging function.
  !
  ! Expect 3D vars to be: (x,y,z,t)
  !        2D vars to be: (x,y,t)
  BadDims = .false.
  
  ! always have var1, use var1 to set Nx,Ny,Nz,Nt
  call rhdf5_read_init(Args%Var1%fname, Var1)
  if (Var1%ndims .eq. 4) then
    Nx = Var1%dims(1)
    Ny = Var1%dims(2)
    Nz = Var1%dims(3)
    Nt = Var1%dims(4)
    Var1Is3d = .true.
  else
    Nx = Var1%dims(1)
    Ny = Var1%dims(2)
    Nz = 1
    Nt = Var1%dims(3)
    Var1Is3d = .false.
  endif
  Var1Nz = Nz

  ! if using a filter, check it against Var1
  if (UseFilter) then
    call rhdf5_read_init(Args%Filter%fname, Filter)
    FilterNz = Filter%dims(3)

    if (.not.(DimsMatch(Var1, Filter))) then
      write (*,*) 'ERROR: dimensions do not match between var1 and filter: ', trim(Var1%vname), trim(Filter%vname)
      BadDims = .true.
    endif
  else
    FilterNz = 1
  endif

  select case (trim(Args%AvgFunc))
    case ('horiz_ke')
      call rhdf5_read_init(Args%Dens%fname, Dens)
      if (.not.(DimsMatch(Var1, Dens))) then
        write (*,*) 'ERROR: dimensions do not match between var1 and dens: ', trim(Var1%vname), trim(Dens%vname)
        BadDims = .true.
      endif

      if (trim(Args%VelInType) .eq. 'uv') then
        call rhdf5_read_init(Args%Var2%fname, Var2)
        if (.not.(DimsMatch(Var1, Var2))) then
          write (*,*) 'ERROR: dimensions do not match between var1 and var2: ', trim(Var1%vname), trim(Var2%vname)
          BadDims = .true.
        endif
      endif

    case ('turb_cov')
      call rhdf5_read_init(Args%Var2%fname, Var2)
      if (.not.(DimsMatch(Var1, Var2))) then
        write (*,*) 'ERROR: dimensions do not match between var1 and var2: ', trim(Var1%vname), trim(Var2%vname)
        BadDims = .true.
      endif

    case ('hist2d')
      call rhdf5_read_init(Args%Var2%fname, Var2)
      if (.not.(DimsMatch(Var1, Var2))) then
        write (*,*) 'ERROR: dimensions do not match between var1 and var2: ', trim(Var1%vname), trim(Var2%vname)
        BadDims = .true.
      endif

    case ('pop')
      call rhdf5_read_init(Args%Var2%fname, Var2) ! lwp, 2D (var1 is pcprate, 2D), so check these two variables
      call rhdf5_read_init(Args%Var3%fname, Var3) ! ltss, 1D so don't need to check dims

      if (.not.(DimsMatch(Var1, Var2))) then
        write (*,*) 'ERROR: dimensions do not match between var1 and var2: ', trim(Var1%vname), trim(Var2%vname)
        BadDims = .true.
      endif

  endselect

  ! check the z dimensions of Var1 vs Var2 when Var2 is being used
  if (UseVar2) then
    if (Var2%ndims .eq. 4) then
      Var2Nz = Var2%dims(3)
      Var2Is3d = .true.
    else
      Var2Nz = 1
      Var2Is3d = .false.
    endif

    ! If both Var1 and Var2 are 3D, then make sure the z dims match
    if (Var1Is3d .and. Var2Is3d) then
      if (Var1Nz .ne. Var2Nz) then
        write (*,*) 'ERROR: var1 and var2 are both 3d, but their z dimensions do not match'
        BadDims = .true.
      endif
    endif

    ! Make sure if doing turb_cov, that var1 and var2 z dimensions match
    if (trim(Args%AvgFunc) .eq. 'turb_cov') then
      if (Var1Nz .ne. Var2Nz) then
        write (*,*) 'ERROR: var1 and var2 need to have their z dimensions match for "turb_cov" function'
        BadDims = .true.
      endif
    endif

    ! At this point Nz is set to Var1Nz. Make sure that Nz is set to the greater
    ! of Var1Nz and Var2Nz. This will allow for a mix of 2d and 3d variables for
    ! the hist2d function.
    if (Var2Nz .gt. Var1Nz) then
      Nz = Var2Nz
    endif
  endif

  if (BadDims) then
    stop
  endif

  ! If subtracting the storm motion, read in the entire time series for the storm
  ! motion (from the storm center file). The storm motion is stored in the file
  ! as one dimenion, (t), so it's more efficient to read in the entire time
  ! series in one shot.
  if (Args%SubtractSmotion) then
    Ustorm%vname = 'storm_speed_x'
    call rhdf5_read_init(Args%StormCenter%fname, Ustorm)
    call rhdf5_read(Args%StormCenter%fname, Ustorm)

    Vstorm%vname = 'storm_speed_y'
    call rhdf5_read_init(Args%StormCenter%fname, Vstorm)
    call rhdf5_read(Args%StormCenter%fname, Vstorm)
  endif

  ! Set up the dimensions for reading in the input field data, one time step per read. 
  ! In other words, remove the time dimension from the input dimensions since we will 
  ! be incrementing through every time step in a loop. The time dimension is always the
  ! last dimension so what this boils down to is to decrement number of dimensions by one.
  
  ! Always have var1
  Var1%ndims = Var1%ndims - 1
  write (*,*) 'Var1 information'
  write (*,*) '  Dataset name: ', trim(Var1%vname)
  write (*,*) '  Number of dimensions: ', Var1%ndims
  write (*,*) '  Dimension sizes:'
  do id = 1, Var1%ndims
    write (*,*), '    ', trim(Var1%dimnames(id)), ': ', Var1%dims(id)
  enddo
  write (*,*) ''

  ! Check for existence of other vars
  if (UseVar2) then
    Var2%ndims = Var2%ndims - 1
    write (*,*) 'Var2 information'
    write (*,*) '  Dataset name: ', trim(Var2%vname)
    write (*,*) '  Number of dimensions: ', Var2%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Var2%ndims
      write (*,*), '    ', trim(Var2%dimnames(id)), ': ', Var2%dims(id)
    enddo
    write (*,*) ''
  endif

  if (UseVar3) then
    Var3%ndims = Var3%ndims - 1
    write (*,*) 'Var3 information'
    write (*,*) '  Dataset name: ', trim(Var3%vname)
    write (*,*) '  Number of dimensions: ', Var3%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Var3%ndims
      write (*,*), '    ', trim(Var3%dimnames(id)), ': ', Var3%dims(id)
    enddo
    write (*,*) ''
  endif

  if (UseFilter) then
    Filter%ndims = Filter%ndims - 1
    write (*,*) 'Filter information'
    write (*,*) '  Dataset name: ', trim(Filter%vname)
    write (*,*) '  Number of dimensions: ', Filter%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Filter%ndims
      write (*,*), '    ', trim(Filter%dimnames(id)), ': ', Filter%dims(id)
    enddo
    write (*,*) ''
  endif

  if (UseDens) then
    Dens%ndims = Dens%ndims - 1
    write (*,*) 'Density information'
    write (*,*) '  Dataset name: ', trim(Dens%vname)
    write (*,*) '  Number of dimensions: ', Dens%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Dens%ndims
      write (*,*), '    ', trim(Dens%dimnames(id)), ': ', Dens%dims(id)
    enddo
    write (*,*) ''
  endif


  ! If doing 'hist', read in the bins. Do it before the next section since V1nbins is being
  ! use to set up the output variable coordinates.
  if ((Args%AvgFunc .eq. 'hist') .or. (Args%AvgFunc .eq. 'pop') .or. (Args%AvgFunc .eq. 'hist2d')) then
    call ReadBinsFile(Args%V1bins%fname, V1nbins, V1bins)
    if ((Args%AvgFunc .eq. 'pop') .or. (Args%AvgFunc .eq. 'hist2d')) then
      call ReadBinsFile(Args%V2bins%fname, V2nbins, V2bins)
    endif
  endif

  ! Set up the dimensions for the output and allocate the output data array. Always
  ! set up as if the output were 3D. This is done so that the output file can
  ! be read into GRADS which expects 3D variables. Always have (x,y,z) for the
  ! dimension names, but set the sizes of the dimensions according to the averaging
  ! function asked for.
  TserAvg%vname = trim(Args%Output%vname)
  TserAvg%descrip = 'time series averaged ' // trim(TserAvg%vname) 
  TserAvg%ndims = 3 
  TserAvg%dimnames(1) = 'x' 
  TserAvg%dimnames(2) = 'y' 
  TserAvg%dimnames(3) = 'z' 

  if ((Args%AvgFunc .eq. 'min') .or. (Args%AvgFunc .eq. 'max') .or. &
      (Args%AvgFunc .eq. 'max_azwind') .or. (Args%AvgFunc .eq. 'min_azslp')) then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = Var1%units
  else if ((Args%AvgFunc .eq. 'hda') .or. (Args%AvgFunc .eq. 'hfrac')) then
    ! y has size 2, one for the sum and the other for the count
    ! all z levels
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 2
    TserAvg%dims(3) = Nz
    TserAvg%units = Var1%units
  else if ((Args%AvgFunc .eq. 'turb_cov') .or. (Args%AvgFunc .eq. 'turb_mmts')) then
    ! y has size 4
    !   For turb_cov
    !     y(1) - sum of mean values of x
    !     y(2) - sum of mean values of y
    !     y(3) - sum of x'y' products
    !     y(4) - number of points
    !   For turb_mmts
    !     y(1) - sum of mean values of x
    !     y(2) - sum of x'x' products
    !     y(3) - sum of x'x'x' products
    !     y(4) - number of points
    ! all z levels
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 4
    TserAvg%dims(3) = Nz
    if (Args%AvgFunc .eq. 'turb_cov') then
      TserAvg%units = trim(Var1%units) // ' ' // trim(Var2%units)
    else if (Args%AvgFunc .eq. 'turb_mmts') then
      TserAvg%units = Var1%units
    endif
  else if (Args%AvgFunc .eq. 'hist') then
    ! put bin values in the x dimension
    TserAvg%dims(1) = V1nbins
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = Nz
    TserAvg%units = Var1%units
  else if (Args%AvgFunc .eq. 'pop') then
    ! put LWP bin values in the x dimension
    ! put LTSS bin values in the y dimension
    ! put Nr and Nt results in z dimension
    TserAvg%dims(1) = V1nbins
    TserAvg%dims(2) = V2nbins
    TserAvg%dims(3) = 2
    TserAvg%units = trim(Var1%units) // ':' // trim(Var2%units)
  else if (Args%AvgFunc .eq. 'hist2d') then
    TserAvg%dims(1) = V1nbins
    TserAvg%dims(2) = V2nbins
    TserAvg%dims(3) = Nz
    TserAvg%units = trim(Var1%units) // ':' // trim(Var2%units)
  else if (Args%AvgFunc .eq. 'ltss') then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = Var1%units
  else if (Args%AvgFunc .eq. 'horiz_ke') then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = 'Joules'
  else if (Args%AvgFunc .eq. 'storm_int') then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = 'int'
  endif

  allocate(TserAvg%vdata(TserAvg%dims(1)*TserAvg%dims(2)*TserAvg%dims(3)))

  ! If doing max_azwind, set up the output var for the radius of max wind
  if (Args%AvgFunc .eq. 'max_azwind') then
    RadMaxWind%vname = 'rmw'
    RadMaxWind%descrip = 'time series averaged ' // trim(Args%AvgFunc) 
    RadMaxWind%units = 'km'
    RadMaxWind%ndims = 3 
    RadMaxWind%dimnames(1) = 'x' 
    RadMaxWind%dimnames(2) = 'y' 
    RadMaxWind%dimnames(3) = 'z' 
    RadMaxWind%dims(1) = 1
    RadMaxWind%dims(2) = 1
    RadMaxWind%dims(3) = 1
    allocate(RadMaxWind%vdata(RadMaxWind%dims(1)*RadMaxWind%dims(2)*RadMaxWind%dims(3)))

    Rad34ktWind%vname = 'r34kt'
    Rad34ktWind%descrip = 'time series averaged ' // trim(Args%AvgFunc) 
    Rad34ktWind%units = 'km'
    Rad34ktWind%ndims = 3 
    Rad34ktWind%dimnames(1) = 'x' 
    Rad34ktWind%dimnames(2) = 'y' 
    Rad34ktWind%dimnames(3) = 'z' 
    Rad34ktWind%dims(1) = 1
    Rad34ktWind%dims(2) = 1
    Rad34ktWind%dims(3) = 1
    allocate(Rad34ktWind%vdata(Rad34ktWind%dims(1)*Rad34ktWind%dims(2)*Rad34ktWind%dims(3)))

    Rad50ktWind%vname = 'r50kt'
    Rad50ktWind%descrip = 'time series averaged ' // trim(Args%AvgFunc) 
    Rad50ktWind%units = 'km'
    Rad50ktWind%ndims = 3 
    Rad50ktWind%dimnames(1) = 'x' 
    Rad50ktWind%dimnames(2) = 'y' 
    Rad50ktWind%dimnames(3) = 'z' 
    Rad50ktWind%dims(1) = 1
    Rad50ktWind%dims(2) = 1
    Rad50ktWind%dims(3) = 1
    allocate(Rad50ktWind%vdata(Rad50ktWind%dims(1)*Rad50ktWind%dims(2)*Rad50ktWind%dims(3)))

    Rad64ktWind%vname = 'r64kt'
    Rad64ktWind%descrip = 'time series averaged ' // trim(Args%AvgFunc) 
    Rad64ktWind%units = 'km'
    Rad64ktWind%ndims = 3 
    Rad64ktWind%dimnames(1) = 'x' 
    Rad64ktWind%dimnames(2) = 'y' 
    Rad64ktWind%dimnames(3) = 'z' 
    Rad64ktWind%dims(1) = 1
    Rad64ktWind%dims(2) = 1
    Rad64ktWind%dims(3) = 1
    allocate(Rad64ktWind%vdata(Rad64ktWind%dims(1)*Rad64ktWind%dims(2)*Rad64ktWind%dims(3)))
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

  ! Read in the input coordinates. Make sure to use a 3D variable if
  ! one is available.
  InCoordFile = trim(Args%Var1%fname)
  if ((.not. Var1Is3d) .and. Var2Is3d) then
    InCoordFile = trim(Args%Var2%fname)
  endif

  InXcoords%vname = 'x_coords'
  call rhdf5_read_init(InCoordFile, InXcoords)
  call rhdf5_read(InCoordFile, InXcoords)

  InYcoords%vname = 'y_coords'
  call rhdf5_read_init(InCoordFile, InYcoords)
  call rhdf5_read(InCoordFile, InYcoords)

  if (Args%AvgFunc .eq. 'storm_int') then
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
  if ((Args%AvgFunc .eq. 'max_azwind') .or. (Args%AvgFunc .eq. 'min_azslp')) then
    ! x,y coords are in meters, convert to km
    do ix = 1, Nx
      InXcoordsKm(ix) = InXcoords%vdata(ix) / 1000.0
    enddo
    do iy = 1, Ny
      InYcoordsKm(iy) = InYcoords%vdata(iy) / 1000.0
    enddo
  else
    ! x,y coords are in degrees lon,lat respectively
    call ConvertGridCoords(Nx, Ny, InXcoords%vdata, InYcoords%vdata, InXcoordsKm, InYcoordsKm)
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

  ! if doing ltss, find the indices associated with Zbot and Ztop
  if (Args%AvgFunc .eq. 'ltss') then
    Kbot = Nt
    do iz = Nz,1,-1
      if (InZcoords%vdata(iz) .ge. Args%Zbot) then
        Kbot = iz
      endif
    enddo

    Ktop = 1
    do iz = 1,Nz
      if (InZcoords%vdata(iz) .le. Args%Ztop) then
        Ktop = iz
      endif
    enddo

    write (*,*) 'Height indices for LTSS calculation:'
    write (*,*) '   Zbot: ', Args%Zbot
    write (*,*) '     Index for Zbot: ', Kbot
    write (*,*) '   Ztop: ', Args%Ztop
    write (*,*) '     Index for Ztop: ', Ktop
    write (*,*) ''
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
  if ((Args%AvgFunc .eq. 'hist') .or. (Args%AvgFunc .eq. 'pop') .or. (Args%AvgFunc .eq. 'hist2d')) then
    OutXcoords%dims(1) = V1nbins
    allocate(OutXcoords%vdata(V1nbins))
    do ib = 1, V1nbins
      OutXcoords%vdata(ib) = V1bins(ib)
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
  if ((Args%AvgFunc .eq. 'pop') .or. (Args%AvgFunc .eq. 'hist2d')) then
    OutYcoords%dims(1) = V2nbins
    allocate(OutYcoords%vdata(V2nbins))
    do ib = 1, V2nbins
      OutYcoords%vdata(ib) = V2bins(ib)
    enddo
  elseif ((Args%AvgFunc .eq. 'hfrac') .or. (Args%AvgFunc .eq. 'hda')) then
    OutYcoords%dims(1) = 2
    allocate(OutYcoords%vdata(2))
    OutYcoords%vdata(1) = 1.0
    OutYcoords%vdata(2) = 2.0
  elseif ((Args%AvgFunc .eq. 'turb_cov') .or. (Args%AvgFunc .eq. 'turb_mmts')) then
    OutYcoords%dims(1) = 4
    allocate(OutYcoords%vdata(4))
    OutYcoords%vdata(1) = 1.0
    OutYcoords%vdata(2) = 2.0
    OutYcoords%vdata(3) = 3.0
    OutYcoords%vdata(4) = 4.0
  else
    OutYcoords%dims(1) = 1
    allocate(OutYcoords%vdata(1))
    OutYcoords%vdata(1) = 1.0
  endif
  
  OutZcoords%vname = 'z_coords'
  OutZcoords%ndims = 1
  OutZcoords%dimnames(1) = 'z'
  OutZcoords%units = 'meter'
  OutZcoords%descrip = 'sigma-z'
  if ((Args%AvgFunc .eq. 'hda') .or. (Args%AvgFunc .eq. 'hist') .or. (Args%AvgFunc .eq. 'hist2d') .or.  &
      (Args%AvgFunc .eq. 'hfrac') .or. (Args%AvgFunc .eq. 'turb_cov') .or. (Args%AvgFunc .eq. 'turb_mmts')) then
    OutZcoords%dims(1) = Nz
    allocate(OutZcoords%vdata(Nz))
    do iz = 1, Nz
      OutZcoords%vdata(iz) = InZcoords%vdata(iz)
    enddo
  elseif (Args%AvgFunc .eq. 'pop') then
    ! need a size of 2 for the z dimension
    !   index 1 --> number of grids that are raining
    !   index 2 --> number of total grids (selected per bin on lwp)
    OutZcoords%dims(1) = 2
    allocate(OutZcoords%vdata(2))
    ! dummy values
    OutZcoords%vdata(1) = 1.0
    OutZcoords%vdata(2) = 2.0
  else
    OutZcoords%dims(1) = 1
    allocate(OutZcoords%vdata(1))
    OutZcoords%vdata(1) = 1.0
  endif

  ! Perform the averaging function.
  rh5f_facc = 'W'
  call rhdf5_open_file(Args%Output%fname, rh5f_facc, 1, rh5f_out)

  rh5f_facc = 'R'

  call rhdf5_open_file(Args%Var1%fname, rh5f_facc, 1, rh5f_var1)
  if (UseVar2) then
    call rhdf5_open_file(Args%Var2%fname, rh5f_facc, 1, rh5f_var2)
  endif
  if (UseVar3) then
    call rhdf5_open_file(Args%Var3%fname, rh5f_facc, 1, rh5f_var3)
  endif
  if (UseFilter) then
    call rhdf5_open_file(Args%Filter%fname, rh5f_facc, 1, rh5f_filter)
  endif
  if (UseDens) then
    call rhdf5_open_file(Args%Dens%fname, rh5f_facc, 1, rh5f_dens)
  endif

  do it = 1, Nt
    ! read the input files
    call rhdf5_read_variable(rh5f_var1, Var1%vname, Var1%ndims, it, Var1%dims, rdata=Var1%vdata)
    if (UseVar2) then
      call rhdf5_read_variable(rh5f_var2, Var2%vname, Var2%ndims, it, Var2%dims, rdata=Var2%vdata)
    endif
    if (UseVar3) then
      call rhdf5_read_variable(rh5f_var3, Var3%vname, Var3%ndims, it, Var3%dims, rdata=Var3%vdata)
    endif
    if (UseFilter) then
      call rhdf5_read_variable(rh5f_filter, Filter%vname, Filter%ndims, it, Filter%dims, rdata=Filter%vdata)
    endif
    if (UseDens) then
      call rhdf5_read_variable(rh5f_dens, Dens%vname, Dens%ndims, it, Dens%dims, rdata=Dens%vdata)
    endif

    ! If subtracting storm motion, GetMyArgs has already checked that Var1 and Var2 are
    ! restricted to combinations of "u" and "v". Allow for Var1 and Var2 to be either
    ! of "u" or "v" here. This will allow for changes in the restrictions in GetMyArgs
    ! (as long as they stick with combinations of "u" and "v") without updating this code.
    if (Args%SubtractSmotion) then
      if (trim(Var1%vname) .eq. 'u') then
        call TranslateField(Nx, Ny, Nz, Var1%vdata, Ustorm%vdata(it))
      elseif (trim(Var1%vname) .eq. 'v') then
        call TranslateField(Nx, Ny, Nz, Var1%vdata, Vstorm%vdata(it))
      endif

      if (trim(Var2%vname) .eq. 'u') then
        call TranslateField(Nx, Ny, Nz, Var2%vdata, Ustorm%vdata(it))
      elseif (trim(Var2%vname) .eq. 'v') then
        call TranslateField(Nx, Ny, Nz, Var2%vdata, Vstorm%vdata(it))
      endif
    endif

    ! do the averaging function
    select case (trim(Args%AvgFunc))

      case ('max_azwind')
        call DoMaxAzWind(Nx, Ny, Nz, Var1%vdata, InXcoordsKm, UndefVal, TserAvg%vdata(1), RadMaxWind%vdata(1), &
                         Rad34ktWind%vdata(1), Rad50ktWind%vdata(1), Rad64ktWind%vdata(1))

      case ('min_azslp')
        call DoMinAzSlp(Nx, Ny, Nz, Var1%vdata, UndefVal, TserAvg%vdata(1))

      case ('horiz_ke')
        if (trim(Args%VelInType) .eq. 'uv') then
          ! Convert U and V to velocity magnitude, and then calculate horizontal KE
          allocate(HorizSpeed(Nx,Ny))
          call UvToSpeed(Nx, Ny, Nz, Var1%vdata, Var2%vdata, HorizSpeed)
          call DoHorizKe(Nx, Ny, Nz, FilterNz, Dens%vdata, HorizSpeed, Filter%vdata, DeltaX, DeltaY, Args%HkeZthick, TserAvg%vdata(1))
          deallocate(HorizSpeed)
        elseif (trim(Args%VelInType) .eq. 's10') then
          call DoHorizKe(Nx, Ny, Nz, FilterNz, Dens%vdata, Var1%vdata, Filter%vdata, DeltaX, DeltaY, Args%HkeZthick, TserAvg%vdata(1))
        endif

      case ('min')
        call DoMin(Nx, Ny, Nz, FilterNz, Var1%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata(1))

      case ('max')
        call DoMax(Nx, Ny, Nz, FilterNz, Var1%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata(1))

      case ('hda')
        call DoHda(Nx, Ny, Nz, FilterNz, Var1%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata)

      case ('hist')
        call DoHist(Nx, Ny, Nz, FilterNz, V1nbins, Var1%vdata, Filter%vdata, UseFilter, UndefVal, V1bins, TserAvg%vdata)

      case ('hfrac')
        call DoHfrac(Nx, Ny, Nz, FilterNz, Args%HfracLimit, Var1%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata)

      case ('turb_mmts')
        call DoTurbMmts(Nx, Ny, Nz, FilterNz, Var1%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata)

      case ('pop')
        call DoPop(Nx, Ny, Nz, FilterNz, V1nbins, V2nbins, Var1%vdata, Var2%vdata, Var3%vdata(1), Filter%vdata, UseFilter, UndefVal, Args%PcprateLimit, V1bins, V2bins, TserAvg%vdata)

      case ('hist2d')
        call DoHist2d(Nx, Ny, Nz, Var1Nz, Var2Nz, FilterNz, V1nbins, V2nbins, Var1%vdata, Var2%vdata, Filter%vdata, UseFilter, UndefVal, V1bins, V2bins, TserAvg%vdata)

      case ('turb_cov')
        call DoTurbCov(Nx, Ny, Nz, FilterNz, Var1%vdata, Var2%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata)

      case ('ltss')
        call DoLtss(Nx, Ny, Nz, FilterNz, Var1%vdata, Filter%vdata, UseFilter, UndefVal, Kbot, Ktop, TserAvg%vdata(1))

      case ('storm_int')
        call DoStormInt(Nx, Ny, FilterNz, Var1%vdata, Filter%vdata, TserAvg%vdata(1))

    endselect

    ! Deallocate the variable data memory
    deallocate(Var1%vdata)
    if (UseVar2) then
      deallocate(Var2%vdata)
    endif
    if (UseVar3) then
      deallocate(Var3%vdata)
    endif
    if (UseFilter) then
      deallocate(Filter%vdata)
    endif
    if (UseDens) then
      deallocate(Dens%vdata)
    endif

    ! write out the averaged data
    call rhdf5_write_variable(rh5f_out, TserAvg%vname, TserAvg%ndims, it, TserAvg%dims, &
       TserAvg%units, TserAvg%descrip, TserAvg%dimnames, rdata=TserAvg%vdata)
    if (Args%AvgFunc .eq. 'max_azwind') then
      call rhdf5_write_variable(rh5f_out, RadMaxWind%vname, RadMaxWind%ndims, it, RadMaxWind%dims, &
         RadMaxWind%units, RadMaxWind%descrip, RadMaxWind%dimnames, rdata=RadMaxWind%vdata)
      call rhdf5_write_variable(rh5f_out, Rad34ktWind%vname, Rad34ktWind%ndims, it, Rad34ktWind%dims, &
         Rad34ktWind%units, Rad34ktWind%descrip, Rad34ktWind%dimnames, rdata=Rad34ktWind%vdata)
      call rhdf5_write_variable(rh5f_out, Rad50ktWind%vname, Rad50ktWind%ndims, it, Rad50ktWind%dims, &
         Rad50ktWind%units, Rad50ktWind%descrip, Rad50ktWind%dimnames, rdata=Rad50ktWind%vdata)
      call rhdf5_write_variable(rh5f_out, Rad64ktWind%vname, Rad64ktWind%ndims, it, Rad64ktWind%dims, &
         Rad64ktWind%units, Rad64ktWind%descrip, Rad64ktWind%dimnames, rdata=Rad64ktWind%vdata)
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
  call rhdf5_close_file(rh5f_var1)
  if (UseVar2) then
    call rhdf5_close_file(rh5f_var2)
  endif
  if (UseVar3) then
    call rhdf5_close_file(rh5f_var3)
  endif
  if (UseFilter) then
    call rhdf5_close_file(rh5f_filter)
  endif
  if (UseDens) then
    call rhdf5_close_file(rh5f_dens)
  endif
 
  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''

  ! Finish off output file
  ! write out the coordinate data
  call rhdf5_write(Args%Output%fname, OutXcoords, 1)
  call rhdf5_write(Args%Output%fname, OutYcoords, 1)
  call rhdf5_write(Args%Output%fname, OutZcoords, 1)
  call rhdf5_write(Args%Output%fname, InTcoords, 1)

  ! If doing hist function, write out the input dimension sizes
  ! for downstream analyses. Eg. if you want to do fractional
  ! area calculations then the counts in the histogram do not
  ! tell you how many total points are in the domain.
  if ((Args%AvgFunc .eq. 'hist') .or. (Args%AvgFunc .eq. 'hist2d')) then
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
    call rhdf5_write(Args%Output%fname, OrigDimSize, 1)

    ! Y
    OrigDimSize%vname = 'Ny'
    OrigDimSize%dimnames(1) = 'y'
    OrigDimSize%descrip = 'number of domain y points'
    OrigDimSize%vdata(1) = float(Ny)
    call rhdf5_write(Args%Output%fname, OrigDimSize, 1)

    ! Z
    OrigDimSize%vname = 'Nz'
    OrigDimSize%dimnames(1) = 'z'
    OrigDimSize%descrip = 'number of domain z points'
    OrigDimSize%vdata(1) = float(Nz)
    call rhdf5_write(Args%Output%fname, OrigDimSize, 1)

    ! T
    OrigDimSize%vname = 'Nt'
    OrigDimSize%dimnames(1) = 't'
    OrigDimSize%descrip = 'number of domain t points'
    OrigDimSize%vdata(1) = float(Nt)
    call rhdf5_write(Args%Output%fname, OrigDimSize, 1)

    deallocate(OrigDimSize%vdata)
  endif

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(Args%Output%fname, OutXcoords, 'x')
  call rhdf5_set_dimension(Args%Output%fname, OutYcoords, 'y')
  call rhdf5_set_dimension(Args%Output%fname, OutZcoords, 'z')
  call rhdf5_set_dimension(Args%Output%fname, InTcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(Args%Output%fname, TserAvg)
  if (Args%AvgFunc .eq. 'max_azwind') then
    call rhdf5_attach_dimensions(Args%Output%fname, RadMaxWind)
    call rhdf5_attach_dimensions(Args%Output%fname, Rad34ktWind)
    call rhdf5_attach_dimensions(Args%Output%fname, Rad50ktWind)
    call rhdf5_attach_dimensions(Args%Output%fname, Rad64ktWind)
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
subroutine GetMyArgs(Args)
  use getoptions

  implicit none

  type (ArgList) :: Args

  character :: optval
  logical :: BadArgs
  logical :: FirstParam
  character (len=MediumString), dimension(MaxArgFields) :: ArgList
  integer :: Nfields

  ! Default
  Args%SubtractSmotion = .false.

  ! Initialization
  Args%AvgFunc = 'none'

  Args%HkeZthick = 0.0
  Args%HfracLimit = 0.0
  Args%Zbot = 0.0
  Args%Ztop = 0.0
  Args%PcprateLimit = 0.0
  Args%VelInType = 'none'

  Args%Output%fname = 'none'
  Args%Output%vname = 'none'

  Args%Filter%fname = 'none'
  Args%Filter%vname = 'none'

  Args%StormCenter%fname = 'none'
  Args%StormCenter%vname = 'none'

  Args%V1bins%fname = 'none'
  Args%V1bins%vname = 'none'

  Args%V2bins%fname = 'none'
  Args%V2bins%vname = 'none'

  Args%Var1%fname = 'none'
  Args%Var1%vname = 'none'
  
  Args%Var2%fname = 'none'
  Args%Var2%vname = 'none'
  
  Args%Var3%fname = 'none'
  Args%Var3%vname = 'none'

  Args%Dens%fname = 'none'
  Args%Dens%vname = 'none'

  ! loop through all of the command line tokens
  ! optarg is a character string variable that the getoptions module supplies
  !   optarg gets set to the argument value for an option that uses an argument
  ! getopt returns single character:
  !    '>' finished
  !    '!' unknown option
  !    '.' command line argument (no option)
  BadArgs = .false.
  FirstParam = .true.
  do
    optval = getopt('m')

    select case (optval)
      case ('>')  ! finished
        exit

      case ('!')  ! unrecognized argument
        write(*,*) 'ERROR: unknown option: ', trim(optarg)
        write(*,*) ''
        BadArgs = .true.

      case ('m')
        Args%SubtractSmotion = .true.

      case ('.')
        if (FirstParam) then
          ! first parameter -> <arg_func>[:<args>]
          call String2List(optarg, ':', ArgList, MaxArgFields, Nfields, 'avgerage function spec')
          Args%AvgFunc = trim(ArgList(1))

          ! Grab args for those functions that use them
          select case (trim(Args%AvgFunc))
            case ('horiz_ke')
              if (Nfields .eq. 3) then
                Args%VelInType = trim(ArgList(2))
                read(ArgList(3), '(f)') Args%HkeZthick
              else
                write (*,*) 'ERROR: average function horiz_ke requires three fields: horiz_ke:<in_type>:<thickness>'
                write(*,*) ''
                BadArgs = .true.
              endif

            case ('hfrac')
              if (Nfields .eq. 2) then
                read(ArgList(2), '(f)') Args%HfracLimit
              else
                write (*,*) 'ERROR: average function hfrac requires two fields: hfrac:<limit>'
                write(*,*) ''
                BadArgs = .true.
              endif

            case ('ltss')
              if (Nfields .eq. 3) then
                read(ArgList(2), '(f)') Args%Zbot
                read(ArgList(3), '(f)') Args%Ztop
              else
                write (*,*) 'ERROR: average function ltss requires three fields: hfrac:<z_bot>:<z_top>'
                write(*,*) ''
                BadArgs = .true.
              endif
 
            case ('pop')
              if (Nfields .eq. 2) then
                read(ArgList(2), '(f)') Args%PcprateLimit
              else
                write (*,*) 'ERROR: average function pop requires two fields: pop:<pcp_limit>'
                write(*,*) ''
                BadArgs = .true.
              endif

          endselect

          FirstParam = .false.
        else
          ! file spec ->   <file_type>:<file_name>:<variable_name>
          call String2List(optarg, ':', ArgList, MaxArgFields, Nfields, 'file spec') 
          if (Nfields .eq. 3) then
            select case (trim(ArgList(1)))
              case ('filter')
                Args%Filter%fname = trim(ArgList(2))
                Args%Filter%vname = trim(ArgList(3))
  
              case ('storm')
                Args%StormCenter%fname = trim(ArgList(2))
                Args%StormCenter%vname = trim(ArgList(3))
  
              case ('v1bins')
                Args%V1bins%fname = trim(ArgList(2))
                Args%V1bins%vname = trim(ArgList(3))
  
              case ('v2bins')
                Args%V2bins%fname = trim(ArgList(2))
                Args%V2bins%vname = trim(ArgList(3))
  
              case ('var1')
                Args%Var1%fname = trim(ArgList(2))
                Args%Var1%vname = trim(ArgList(3))
  
              case ('var2')
                Args%Var2%fname = trim(ArgList(2))
                Args%Var2%vname = trim(ArgList(3))
  
              case ('var3')
                Args%Var3%fname = trim(ArgList(2))
                Args%Var3%vname = trim(ArgList(3))
  
              case ('dens')
                Args%Dens%fname = trim(ArgList(2))
                Args%Dens%vname = trim(ArgList(3))
  
              case ('out')
                Args%Output%fname = trim(ArgList(2))
                Args%Output%vname = trim(ArgList(3))
  
              case default
                write(*,*) 'ERROR: unknown <file_type>: ', trim(ArgList(1))
                write(*,*) ''
                BadArgs = .true.

            endselect
          else
            write(*,*) 'ERROR: must use <file_type>:<file_name>:<variable_name> for file specs: ', trim(optarg)
            write(*,*) ''
            BadArgs = .true.
          endif
        endif

    endselect
  enddo

  ! Check for required files
  ! Leave empty cases so that default can throw error if an unknown avg func is used
  if ((trim(Args%Var1%fname) .eq. 'none') .or. (trim(Args%Output%fname) .eq. 'none')) then
    write (*,*) 'ERROR: must always specify var1 and out files'
    write (*,*) ''
    BadArgs = .true.
  endif

  if (Args%SubtractSmotion) then
    if (trim(Args%StormCenter%fname) .eq. 'none') then
      write (*,*) 'ERROR: must specify storm center file when using -m option'
      write (*,*) ''
      BadArgs = .true.
    endif
  endif

  select case (trim(Args%AvgFunc))
    case ('horiz_ke')
      if (trim(Args%Dens%fname) .eq. 'none') then
        write (*,*) 'ERROR: must specify dens file with horiz_ke'
        write (*,*) ''
        BadArgs = .true.
      endif

      if (trim(Args%VelInType) .eq. 'uv') then
        if (trim(Args%Var2%fname) .eq. 'none') then
          write (*,*) 'ERROR: must specify var2 file with horiz_ke and <in_type> == "uv"'
          write (*,*) ''
          BadArgs = .true.
        endif
      endif

    case ('min')

    case ('max')

    case ('hda')

    case ('hist')
      if (trim(Args%V1bins%fname) .eq. 'none') then
        write (*,*) 'ERROR: must specify v1bins file with hist'
        write (*,*) ''
        BadArgs = .true.
      endif

    case ('turb_cov')
      if (trim(Args%Var2%fname) .eq. 'none') then
        write (*,*) 'ERROR: must specify var2 file with turb_cov'
        write (*,*) ''
        BadArgs = .true.
      endif

    case ('turb_mmts')

    case ('hfrac')

    case ('hist2d')
      if ((trim(Args%Var2%fname) .eq. 'none') .or. &
          (trim(Args%V1bins%fname) .eq. 'none') .or. &
          (trim(Args%V2bins%fname) .eq. 'none')) then
        write (*,*) 'ERROR: must always specify var2, v1bins and v2bins files with hist2d'
        write (*,*) ''
        BadArgs = .true.
      endif

    case ('ltss')

    case ('min_azslp')
      if (UseFilter) then
        write (*,*) 'ERROR: cannot use a filter with function: min_azslp'
        write (*,*) ''
        BadArgs = .true.
      endif

    case ('max_azwind')
      if (UseFilter) then
        write (*,*) 'ERROR: cannot use a filter with function: max_azwind'
        write (*,*) ''
        BadArgs = .true.
      endif

    case ('pop')
      if ((trim(Args%Var2%fname) .eq. 'none') .or. &
          (trim(Args%Var3%fname) .eq. 'none') .or. &
          (trim(Args%V1bins%fname) .eq. 'none') .or. &
          (trim(Args%V2bins%fname) .eq. 'none')) then
        write (*,*) 'ERROR: must specify var2, var3, v1bins and v2bins files with pop'
        write (*,*) ''
        BadArgs = .true.
      endif

    case default
      write (*,*) 'ERROR: <avg_function> must be one of:'
      write (*,*) '          horiz_ke'
      write (*,*) '          min'
      write (*,*) '          max'
      write (*,*) '          hda'
      write (*,*) '          hist'
      write (*,*) '          turb_cov'
      write (*,*) '          turb_mmts'
      write (*,*) '          hfrac'
      write (*,*) '          hist2d'
      write (*,*) '          ltss'
      write (*,*) '          storm_int'
      write (*,*) '          min_azslp'
      write (*,*) '          max_azwind'
      write (*,*) '          pop'
      write (*,*) ''
      BadArgs = .true.

  endselect

  ! Check to make sure u and v are being used with the -m option.
  ! For single var averaging functions, Var1 can be be either u or v.
  ! For two var averaging functions, Var1 must be u and Var2 must be v.
  if (Args%SubtractSmotion) then
    if (trim(Args%Var2%vname) .eq. 'none') then
      ! Single var averaging function
      if ((trim(Args%Var1%vname) .ne. 'u') .and. (trim(Args%Var1%vname) .ne. 'v')) then
        write (*,*) 'ERROR: must specify var1 as "u" or "v" when using -m option'
        write (*,*) ''
        BadArgs = .true.
      endif
    else
      ! Two var averaging functions
      if ((trim(Args%Var1%vname) .ne. 'u') .or. (trim(Args%Var2%vname) .ne. 'v')) then
        write (*,*) 'ERROR: must specify var1 as "u" and var2 as "v" when using -m option'
        write (*,*) ''
        BadArgs = .true.
      endif
    endif
  endif

  ! check other specs for valid values
  if (Args%VelInType .ne. 'none') then
    if ((Args%VelInType .ne. 'uv') .and. (Args%VelInType .ne. 's10')) then
      write (*,*) 'ERROR: <in_type> for average functions must be one of: "uv" or "s10"'
      write (*,*) ''
      BadArgs = .true.
    endif
  endif

  if (BadArgs) then
    write (*,*) 'USAGE: tsavg [-m] <avg_function> <file_list>'
    write (*,*) ''
    write (*,*) '   -m: subtract storm motion from input winds'
    write (*,*) ''
    write (*,*) '   <avg_function>: averaging function to use on input data'
    write (*,*) '            horiz_ke:<in_type>:<thickness> -> total kinetic energy form horizontal winds'
    write (*,*) '              <in_type>: "uv" calculate from lowest level of u and v fields, set var1 to u, var2 to v'
    write (*,*) '                         "s10" calculate from 10m wind speed field, set var1 to speed10m'
    write (*,*) '              <thickness>: depth (m) of slab that KE is being calculated in' 
    write (*,*) '                          Powell and Reinhold, 2007 use 1.0 m'
    write (*,*) ''
    write (*,*) '            min -> domain minimum'
    write (*,*) ''
    write (*,*) '            max -> domain maximum'
    write (*,*) ''
    write (*,*) '            hda -> horizontal domain average at each z level'
    write (*,*) ''
    write (*,*) '            hist -> histogram across each z level'
    write (*,*) ''
    write (*,*) '            turb_cov -> calculate turbulence covariance for each z level'
    write (*,*) ''
    write (*,*) '            turb_mmts -> calculate turbulence moments for each z level'
    write (*,*) ''
    write (*,*) '            hfrac:<limit> -> fractional occurrence across each z level'
    write (*,*) '              <limit>: if var > threshold, that hoizontal grid cell gets a 1, otherwise a zero.'
    write (*,*) '                       Then the fraction for that level = Number of cells with 1s divided by total Number of cells'
    write (*,*) ''
    write (*,*) '            hist2d -> construct 2d histogram'
    write (*,*) ''
    write (*,*) '            ltss::<z_bot>:<z_top> -> calculate lower tropospheric static stability'
    write (*,*) '              var1 must be set to theta variable (3d)'
    write (*,*) '              <z_bot>: height (Z) for bottom'
    write (*,*) '              <z_top>: height (Z) for top'
    write (*,*) ''
    write (*,*) '            storm_int -> storm intensity metric from horizontal wind speeds'
    write (*,*) ''
    write (*,*) '            max_azwind -> max value of azimuthially averaged wind'
    write (*,*) '              <in_type>: "uv" calculate from lowest level of u and v fields,'
    write (*,*) '                         "s10" calculate from 10m wind speed field,'
    write (*,*) ''
    write (*,*) '            min_azslp -> min value of azimuthially averaged sea-level pressure'
    write (*,*) ''
    write (*,*) '            pop:<pcp_limit>'
    write (*,*) '              var1 must be set to precip rate variable (2d):'
    write (*,*) '                <pcp_limit>: threshold to determine if raining'
    write (*,*) '                  precip rate >= <pcp_limit> --> raining'
    write (*,*) '                  precip rate <  <pcp_limit> --> not raining'
    write (*,*) '              var2 must be set to liquid water path variable (2d):'
    write (*,*) '              var3 must be set to lower tropospheric static stabilty variable (1d):'
    write (*,*) ''
    write (*,*) '   <file_list>: supply specs for input and output files'
    write (*,*) '     <file_list> := <file_spec> ...'
    write (*,*) '       <file_spec> := <file_type>:<file_name>:<variable_name>'
    write (*,*) '          <file_type> is one of:'
    write (*,*) '             "filter": filter mask, optional'
    write (*,*) '             "storm": storm track and motion, required only if -m used'
    write (*,*) '             "v1bins": histogram bin edges, required for <avg_function>: hist pop hist2d'
    write (*,*) '             "v2bins": histogram bin edges, required for <avg_function>: pop hist2d'
    write (*,*) '             "var1": variable to average'
    write (*,*) '             "var2": variable to average'
    write (*,*) '             "var3": variable to average'
    write (*,*) '             "u": zonal component of momentum, required only if -t used'
    write (*,*) '             "v": meridional component of momentum, required only if -t used'
    write (*,*) '             "dens": density'
    write (*,*) '             "out": output file, required'
    write (*,*) '          <file_name> is the complete path to the file'
    write (*,*) '          <variable_name> is the name of the variable inside the file'
    write (*,*) ''

    stop
  endif

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
! Also keep track of the greatest radial index where the wind speed was
! >= 34 kt, 50 kt and 64 kt. AzWind is in m/s so the corresponding values
! are:
!   34 kt -> 17.5 m/s
!   50 kt -> 25.7 m/s
!   64 kt -> 32.9 m/s
!
subroutine DoMaxAzWind(Nx, Ny, Nz, AzWind, Radius, UndefVal, AzWindMax, RadMax, Rad34kt, Rad50kt, Rad64kt)
  implicit none

  ! 1 kt = 0.514444 m/s
  real, parameter :: w_34kt_in_ms = 34.0 * 0.514444
  real, parameter :: w_50kt_in_ms = 50.0 * 0.514444
  real, parameter :: w_64kt_in_ms = 64.0 * 0.514444

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: AzWind
  real, dimension(Nx) :: Radius
  real :: UndefVal
  real :: AzWindMax
  real :: RadMax, Rad34kt, Rad50kt, Rad64kt

  integer :: ix,iy,iz
  integer :: i_rmw, i_34kt, i_50kt, i_64kt

  ! dimension order is: x,y,z

  ! RAMS places the first z level below the surface so use the second level
  ! to approximate the 10m winds.

  AzWindMax = 0.0
  i_rmw = 0
  i_34kt = 0
  i_50kt = 0
  i_64kt = 0
  iy = 1 ! dummy dimension in azmuthally averaged wind speed
  if (Nz .eq. 1) then
    iz = 1
  else
    iz = 2
  endif
  do ix = 1, Nx
    ! Max wind
    if ((AzWind(ix,iy,iz) .gt. AzWindMax) .and. (anint(AzWind(ix,iy,iz)) .ne. UndefVal)) then
      AzWindMax = AzWind(ix,iy,iz)
      i_rmw = ix
    endif
    ! Last index where wind >= 34 kt
    if ((AzWind(ix,iy,iz) .ge. w_34kt_in_ms) .and. (anint(AzWind(ix,iy,iz)) .ne. UndefVal)) then
      i_34kt = ix
    endif
    ! Last index where wind >= 50 kt
    if ((AzWind(ix,iy,iz) .ge. w_50kt_in_ms) .and. (anint(AzWind(ix,iy,iz)) .ne. UndefVal)) then
      i_50kt = ix
    endif
    ! Last index where wind >= 64 kt
    if ((AzWind(ix,iy,iz) .ge. w_64kt_in_ms) .and. (anint(AzWind(ix,iy,iz)) .ne. UndefVal)) then
      i_64kt = ix
    endif
  end do

  RadMax = Radius(i_rmw)
  Rad34kt = InterpRadius(Nx, Ny, Nz, AzWind, Radius, w_34kt_in_ms, i_34kt, UndefVal)
  Rad50kt = InterpRadius(Nx, Ny, Nz, AzWind, Radius, w_50kt_in_ms, i_50kt, UndefVal)
  Rad64kt = InterpRadius(Nx, Ny, Nz, AzWind, Radius, w_64kt_in_ms, i_64kt, UndefVal)

  return
end subroutine DoMaxAzWind

!************************************************************************************
! DoMinAzSlp()
!
! This subroutine will find the min sea-level pressure in AzSlp and copy that
! to TserAvg. The intent of this metric is to mimic the usage of mininum central
! pressure for an intensity measurement.
!
!
subroutine DoMinAzSlp(Nx, Ny, Nz, AzSlp, UndefVal, AzSlpMin)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: AzSlp
  real :: UndefVal
  real :: AzSlpMin

  integer :: ix,iy,iz

  ! dimension order is: x,y,z

  ! RAMS places the first z level below the surface so use the second level
  ! to approximate the 10m winds.

  AzSlpMin = 1000000
  iy = 1 ! dummy dimension in azmuthally averaged sea level pressure
  iz = 1 ! dummy dimension in azmuthally averaged sea level pressure
  do ix = 1, Nx
    ! Min SLP
    if ((AzSlp(ix,iy,iz) .lt. AzSlpMin) .and. (anint(AzSlp(ix,iy,iz)) .ne. UndefVal)) then
      AzSlpMin = AzSlp(ix,iy,iz)
    endif
  end do

  return
end subroutine DoMinAzSlp

!**************************************************************************************
! InterpRadius()
!
! This function will interpolate a radial value (x) given the cooresponding wind speed.
!
real function InterpRadius(Nx, Ny, Nz, AzWind, Radius, W, Ileft, UndefVal)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: AzWind
  real, dimension(Nx) :: Radius
  real :: W, UndefVal
  integer :: Ileft

  integer :: ix1, ix2, iy, iz
  real :: W1, W2, R1, R2

  iy = 1 ! dummy dimension in azmuthally averaged wind speed
  if (Nz .eq. 1) then
    iz = 1
  else
    iz = 2
  endif

  if (Ileft .eq. 0) then
    ! Didn't find any entry in AzWind >= W -> return UndefVal
    InterpRadius = UndefVal
  else if (Ileft .eq. Nx) then
    ! The last entry in AzWind is >= W -> can't interpolate so just return the last entry
    InterpRadius = AzWind(Nx,iy,iz)
  else
    ! W has fallen between two entries in AzWind -> linearly interpolate radius value
    ! Ileft points to the last entry in AzWind that was >= W (left side of interval)
    ix1 = Ileft
    ix2 = Ileft + 1

    R1 = Radius(ix1)
    R2 = Radius(ix2)
    W1 = AzWind(ix1,iy,iz)
    W2 = AzWind(ix2,iy,iz)

    InterpRadius = R1 + (((W - W1)/(W2 - W1)) * (R2 - R1))
  endif

end function InterpRadius

!**************************************************************************************
! DoHorizKe()
!
! This routine will calculate the total kinetic energy over the given cylindrical
! volume. Do not want average since we want the size of the storm reflected in
! this diagnostic.

subroutine DoHorizKe(Nx, Ny, Nz, FilterNz, Dens, Speed, Filter, DeltaX, DeltaY, Zthick, TserAvg)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Dens
  real, dimension(Nx,Ny,FilterNz) :: Filter
  real, dimension(Nx,Ny) :: Speed
  real :: DeltaX, DeltaY
  real :: Zthick
  real :: TserAvg

  integer :: ix,iy,iz
  integer :: filter_z
  integer :: NumPoints
  real :: SumKe, CurrKe

  ! KE is 1/2 * m * v^2
  !   - calculate this a every point inside the defined cylindrical volume
  !   - get m by density * volume
  !   - v^2 is based on horiz velocity so is equal to Speed^2
  !
  ! The integrated kinetic engery measurement done on hurricanes is taken over the domain on the 10m
  ! level. Just use the first model level above the surface (z == 2) for this measurement which will
  ! be close enough.
  iz = 2

  SumKe = 0.0
  NumPoints = 0

  if (FilterNz .eq. 1) then
    filter_z = 1
  else
    filter_z = iz
  endif

  do iy = 2, Ny-1
    do ix = 2, Nx-1
      if (anint(Filter(ix,iy,filter_z)) .eq. 1.0) then
        CurrKe = 0.5 * DeltaX * DeltaY * Zthick * Dens(ix,iy,iz) * Speed(ix,iy)**2.
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
subroutine DoHda(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomStats)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(2,Nz) :: DomStats ! (1,z) <-- domain sum, (2,z) <-- domain count

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint
  integer :: NumPoints
  real    :: DomSum

  do iz = 1, Nz
    DomSum = 0.0
    NumPoints = 0

    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (FilterNz .eq. 1) then
      filter_z = 1
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
          DomSum = DomSum + Var(ix,iy,iz)
          NumPoints = NumPoints + 1
        endif
      enddo
    enddo

    DomStats(1,iz) = DomSum
    DomStats(2,iz) = float(NumPoints)
  enddo

  return
end subroutine DoHda

!**************************************************************************************
! DoHfrac()
!
! This routine will do horizontal domain fraction calculations
!
subroutine DoHfrac(Nx, Ny, Nz, FilterNz, Threshold, Var, Filter, UseFilter, UndefVal, DomStats)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal, Threshold
  real, dimension(2,Nz) :: DomStats ! (1,x) <-- fraction count, (2,x) <-- total count

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint
  integer :: NumPoints, NumFracPoints

  do iz = 1, Nz
    NumPoints = 0
    NumFracPoints = 0

    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (FilterNz .eq. 1) then
      filter_z = 1
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
          if (Var(ix,iy,iz) .gt. Threshold) then
            NumFracPoints = NumFracPoints + 1
          endif
          NumPoints = NumPoints + 1
        endif
      enddo
    enddo

    DomStats(1,iz) = float(NumFracPoints)
    DomStats(2,iz) = float(NumPoints)
  enddo

  return
end subroutine DoHfrac

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

    if (FilterNz .eq. 1) then
      filter_z = 1
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

    if (FilterNz .eq. 1) then
      filter_z = 1
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
  real, dimension(Nb,Nz) :: Counts
  real, dimension(Nb) :: Bins

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
  ! Build the histogram (counts)
  do iz = 1, Nz
    do ib = 1, Nb
      Counts(ib,iz) = 0.0
    enddo

    if (FilterNz .eq. 1) then
      filter_z = 1
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal

        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        
        ib = FindBin(Nb, Bins, Var(ix,iy,iz))
        SelectPoint = SelectPoint .and. (ib .ne. -1)
        
        if (SelectPoint) then
          Counts(ib,iz) = Counts(ib,iz) + 1.0
        endif
      enddo
    enddo
  enddo

  return
end subroutine DoHist

!**************************************************************************************
! DoPop()
!
! This routine will create counts for generating a probability of precip statistic.
! The data is binned on LWP, and each bin gets two counts: one for the number of grid cells that
! it is raining in that bin and the other for the total number of cells in that bin. A cell counts
! as raining when its precip rate is >= to the precip limit
!

subroutine DoPop(Nx, Ny, Nz, FilterNz, Nb_lwp, Nb_ltss, PrecipRate, Lwp, Ltss, Filter, UseFilter, UndefVal, PrecipRateLimit, Bins_lwp, Bins_ltss, Counts)
  implicit none

  integer :: Nx, Ny, Nz, Nb_lwp, Nb_ltss, FilterNz
  real, dimension(Nx,Ny,Nz) :: PrecipRate, Lwp
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal, PrecipRateLimit, Ltss
  real, dimension(Nb_lwp) :: Bins_lwp
  real, dimension(Nb_ltss) :: Bins_ltss
  real, dimension(Nb_lwp,Nb_ltss,2) :: Counts

  integer :: ib_lwp, ib_ltss, ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint

  ! Emulate matlab histc command, ie the values in Bins are treated as the edges
  ! of the bin ranges. Do not count any values from Lwp that fall outside the
  ! bin range (< Bins(1) or > Bins(Nb)). The count is incremented for bin ib when:
  !
  !   Bins(ib) <= Lwp(ix,iy,iz) < Bins(ib+1)
  !
  ! The count in Bin(Nb), last bin, is incremented when:
  !
  !   Lwp(ix,iy,iz) == Bins(ib)
  !
  ! Same for the Ltts bins
  !
  do ib_lwp = 1, Nb_lwp
    do ib_ltss = 1, Nb_ltss
      Counts(ib_lwp,ib_ltss,1) = 0.0
      Counts(ib_lwp,ib_ltss,2) = 0.0
    enddo
  enddo

  ! Figure out the LTSS bin
  ! Then loop through the LWP and figure out the counts
  ib_ltss = FindBin(Nb_ltss, Bins_ltss, Ltss)

  ! Only if an LTSS bin has been selected
  if (ib_ltss .ne. -1) then
    do iz = 1, Nz
      if (FilterNz .eq. 1) then
        filter_z = 1
      else
        filter_z = iz
      endif
  
      do iy = 2, Ny-1
        do ix = 2, Nx-1
          SelectPoint = anint(Lwp(ix,iy,iz)) .ne. UndefVal
          if (UseFilter) then
            SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
          endif
          if (SelectPoint) then
            ib_lwp = FindBin(Nb_lwp, Bins_lwp, Lwp(ix,iy,iz))
            if (ib_lwp .ne. -1) then
              Counts(ib_lwp,ib_ltss,1) = Counts(ib_lwp,ib_ltss,1) + 1.0
              if (PrecipRate(ix,iy,iz) .ge. PrecipRateLimit) then
                Counts(ib_lwp,ib_ltss,2) = Counts(ib_lwp,ib_ltss,2) + 1.0
              endif
            endif
          endif
        enddo
      enddo
    enddo
  endif

  return
end subroutine DoPop

!**************************************************************************************
! DoHist2d()
!
! This routine will create a 2D histogram.
!
subroutine DoHist2d(Nx, Ny, Nz, Var1Nz, Var2Nz, FilterNz, V1nbins, V2nbins, Xvar, Yvar, Filter, UseFilter, UndefVal, V1bins, V2bins, Counts)
  implicit none

  integer :: Nx, Ny, Nz, Var1Nz, Var2Nz, FilterNz, V1nbins, V2nbins
  real, dimension(Nx,Ny,Var1Nz) :: Xvar
  real, dimension(Nx,Ny,Var2Nz) :: Yvar
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(V1nbins) :: V1bins
  real, dimension(V2nbins) :: V2bins
  real, dimension(V1nbins,V2nbins,Nz) :: Counts

  integer :: ix, iy, iz, iz_x, iz_y, iz_filter
  integer :: ib_x, ib_y

  logical :: SelectPoint

  ! go level by level
  do iz = 1, Nz
    ! zero out the counts
    do ib_x = 1, V1nbins
      do ib_y = 1, V2nbins
        Counts(ib_x,ib_y,iz) = 0.0
      enddo
    enddo

    ! Since Xvar and Yvar can have differing z dimensions, check to make sure we are not
    ! going past the max z coordinate.
    if (iz .gt. Var1Nz) then
      iz_x = Var1Nz
    else
      iz_x = iz
    endif
    if (iz .gt. Var2Nz) then
      iz_y = Var2Nz
    else
      iz_y = iz
    endif
   
    if (FilterNz .eq. 1) then
      iz_filter = 1
    else
      iz_filter = iz
    endif

    ! Walk through each grid cell on this level. Simultaneously bin the Xvar and Yvar values.
    ! Increment the count for this grid cell only if the Xvar and Yvar values got placed into
    ! an Xbin and Ybin.
    do iy = 2, Ny-1
      do ix = 2, Nx-1
        ! don't select if either x or y is undefined
        SelectPoint = ((anint(Xvar(ix,iy,iz_x)) .ne. UndefVal) .and. &
                       (anint(Yvar(ix,iy,iz_y)) .ne. UndefVal))

        ! if using a filter, check the filter value
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,iz_filter)) .eq. 1.0)
        endif

        ! get the bins and only count if both x and y bins were found
        ib_x = FindBin(V1nbins, V1bins, Xvar(ix,iy,iz_x))
        ib_y = FindBin(V2nbins, V2bins, Yvar(ix,iy,iz_y))
        SelectPoint = SelectPoint .and. (ib_x .ne. -1) .and. (ib_y .ne. -1)
        
        if (SelectPoint) then
          Counts(ib_x,ib_y,iz) = Counts(ib_x,ib_y,iz) + 1.0
        endif
      enddo
    enddo
  enddo
end subroutine DoHist2d

!**************************************************************************************
! DoLtss()
!
! This routine will calculate the lower tropsheric static stability (LTSS) from the
! theta field. 
!
! Use an average of differences. Tried difference of averages and the traces are noisy.

subroutine DoLtss(Nx, Ny, Nz, FilterNz, Theta, Filter, UseFilter, UndefVal, Kbot, Ktop, Ltss)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Theta
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal, Ltss
  integer :: Kbot, Ktop

  integer :: ix, iy
  integer :: filter_zbot, filter_ztop
  logical :: SelectPoint
  integer :: Npts

  ! Figure out which levels to use in the filter data
  if (FilterNz .eq. 1) then
    filter_zbot = 1
    filter_ztop = 1
  else
    filter_zbot = Kbot
    filter_ztop = Ktop
  endif

  ! Calculate theta difference at each point and then take the average of
  ! these differences.
  Ltss = 0.0
  Npts = 0
  do iy = 2, Ny-1
    do ix = 2, Nx-1
      SelectPoint = (anint(Theta(ix,iy,Kbot)) .ne. UndefVal) .and. (anint(Theta(ix,iy,Ktop)) .ne. UndefVal)
      if (UseFilter) then
        SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_zbot)) .eq. 1.0) &
                                  .and. (anint(Filter(ix,iy,filter_ztop)) .eq. 1.0)
      endif

      ! Sum up the theta differences
      if (SelectPoint) then
        Ltss = Ltss + (Theta(ix,iy,Ktop) - Theta(ix,iy,Kbot))
        Npts = Npts + 1
      endif
    enddo
  enddo

  if (Npts .eq. 0) then
    Ltss = UndefVal
  else
    Ltss = Ltss / float(Npts)
  endif

  return
end subroutine DoLtss

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

subroutine DoStormInt(Nx, Ny, FilterNz, Speed10m, Filter, TserAvg)
  implicit none

  integer :: Nx, Ny, FilterNz
  real, dimension(Nx,Ny,FilterNz) :: Filter
  real, dimension(Nx,Ny) :: Speed10m
  real :: TserAvg

  integer ix,iy
  integer nCat0, nCat1, nCat2, nCat3, nCat4, nCat5, NumPoints
  real Wspeed, SiMetric

  nCat0 = 0
  nCat1 = 0
  nCat2 = 0
  nCat3 = 0
  nCat4 = 0
  nCat5 = 0

  do iy = 2, Ny-1
    do ix = 2, Nx-1
      if (anint(Filter(ix,iy,1)) .eq. 1.0) then
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

!**************************************************************************************
! DoTurbCov()
!
! This routine will do summing for the calculation of the mean, variance and skew
! (first three moments) of the input variable.
!
subroutine DoTurbCov(Nx, Ny, Nz, FilterNz, Xvar, Yvar, Filter, UseFilter, UndefVal, DomStats)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Xvar, Yvar
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(4,Nz) :: DomStats ! (1,Nz) <-- sum of mean values from running mean of Xvar
                                    ! (2,Nz) <-- sum of mean values from running mean of Yvar
                                    ! (3,Nz) <-- sum of (x-xmean)*(y-ymean)
                                    ! (4,Nz) <-- number of points

  integer :: ix, iy, iz
  integer :: Xstart, Xend, Ystart, Yend
  integer :: filter_z
  logical :: SelectPoint
  integer :: NumPoints
  double precision    :: DomS1
  double precision    :: DomS2
  double precision    :: DomS3
  real    :: Xmean, Ymean

  ! RAMS reserves the first and last x and y values for lateral
  ! boundaries. These only contain valid field data under certain
  ! circumstances such as cyclic boundary cases. Most of the time
  ! we want these to be excluded so for now always exclude them
  ! (shouldn't hurt results with cyclic boundaries where the
  ! boundary values could have been included).
  Xstart = 2
  Xend = Nx - 1
  Ystart = 2
  Yend = Ny - 1

  do iz = 1, Nz
    DomS1 = 0.0d+0
    DomS2 = 0.0d+0
    DomS3 = 0.0d+0
    NumPoints = 0

    if (FilterNz .eq. 1) then
      filter_z = 1
    else
      filter_z = iz
    endif

    do iy = Ystart, Yend
      do ix = Xstart, Xend
        SelectPoint = ((anint(Xvar(ix,iy,iz)) .ne. UndefVal) .and. (anint(Yvar(ix,iy,iz)) .ne. UndefVal))
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          Xmean = Calc2dMean(Nx, Ny, Nz, Xvar, ix, iy, iz, Xstart, Xend, Ystart, Yend)
          Ymean = Calc2dMean(Nx, Ny, Nz, Yvar, ix, iy, iz, Xstart, Xend, Ystart, Yend)

          DomS1 = DomS1 + Xmean
          DomS2 = DomS2 + Ymean
          DomS3 = DomS3 + (dble(Xvar(ix,iy,iz)-Xmean)*dble(Yvar(ix,iy,iz)-Ymean))

          NumPoints = NumPoints + 1
        endif
      enddo
    enddo

    DomStats(1,iz) = real(DomS1)
    DomStats(2,iz) = real(DomS2)
    DomStats(3,iz) = real(DomS3)
    DomStats(4,iz) = float(NumPoints)
  enddo

  return
end subroutine DoTurbCov

!**************************************************************************************
! DoTurbMmts()
!
! This routine will do summing for the calculation of the mean, variance and skew
! (first three moments) of the input variable.
!
subroutine DoTurbMmts(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomStats)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(4,Nz) :: DomStats ! (1,Nz) <-- sum of mean values from running mean
                                    ! (2,Nz) <-- sum of (x-xmean)**2
                                    ! (3,Nz) <-- sum of (x-xmean)**3
                                    ! (4,Nz) <-- number of points

  integer :: ix, iy, iz
  integer :: Xstart, Xend, Ystart, Yend
  integer :: filter_z
  logical :: SelectPoint
  integer :: NumPoints
  double precision    :: DomS1
  double precision    :: DomS2
  double precision    :: DomS3
  real    :: Vmean

  ! RAMS reserves the first and last x and y values for lateral
  ! boundaries. These only contain valid field data under certain
  ! circumstances such as cyclic boundary cases. Most of the time
  ! we want these to be excluded so for now always exclude them
  ! (shouldn't hurt results with cyclic boundaries where the
  ! boundary values could have been included).
  Xstart = 2
  Xend = Nx - 1
  Ystart = 2
  Yend = Ny - 1

  do iz = 1, Nz
    DomS1 = 0.0d+0
    DomS2 = 0.0d+0
    DomS3 = 0.0d+0
    NumPoints = 0

    if (FilterNz .eq. 1) then
      filter_z = 1
    else
      filter_z = iz
    endif

    do iy = Ystart, Yend
      do ix = Xstart, Xend
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          Vmean = Calc2dMean(Nx, Ny, Nz, Var, ix, iy, iz, Xstart, Xend, Ystart, Yend)
          DomS1 = DomS1 + Vmean
          DomS2 = DomS2 + (dble(Var(ix,iy,iz)-Vmean)**2)
          DomS3 = DomS3 + (dble(Var(ix,iy,iz)-Vmean)**3)
          NumPoints = NumPoints + 1
        endif
      enddo
    enddo

    DomStats(1,iz) = real(DomS1)
    DomStats(2,iz) = real(DomS2)
    DomStats(3,iz) = real(DomS3)
    DomStats(4,iz) = float(NumPoints)
  enddo

  return
end subroutine DoTurbMmts

!*****************************************************************************
! Calc2dMean()
!
! This function will find the mean for the given ix, iy, iz point. The mean
! is taken from the surrounding N points about (ix,iy) on the level given
! by iz. When called for each (x,y) point on a level, this effective does
! a 2D running average over that level.
!
! 
!
real function Calc2dMean(Nx, Ny, Nz, Var, ix, iy, iz, xs, xe, ys, ye)
  integer :: Nx, Ny, Nz, ix, iy, iz, xs, xe, ys, ye
  real, dimension(Nx,Ny,Nz) :: Var

  integer :: i, j
  integer :: Fsize
  real :: VarVal
  real :: Npts

  ! use a 5x5 box for computing the average
  Fsize = 2
  Npts = 25.0

  Calc2dMean = 0.0
  do j = iy-Fsize, iy+Fsize
    do i = ix-Fsize, ix+Fsize
      ! If the box extends beyond the borders of Var, then fill with zeros.
      if ((i .lt. xs) .or. (j .lt. ys) .or. (i .gt. xe) .or. (j .gt. ye)) then
         VarVal = 0.0
      else
         VarVal = Var(i,j,iz)
      endif

      Calc2dMean = Calc2dMean + VarVal
    enddo
  enddo
  
  Calc2dMean = Calc2dMean / Npts

end function Calc2dMean

!************************************************************************
! UvToSpeed()
!
! This function will convert the U, V fields into speed (vector magnitude)
! on the lowest above ground model level.

subroutine UvToSpeed(Nx, Ny, Nz, U, V, HorizSpeed)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: U, V
  real, dimension(Nx,Ny) :: HorizSpeed

  integer :: ix, iy, iz

  ! RAMS has first model level underground, so use z = 2 (1st model level above ground)
  iz = 2
  do iy = 1, Ny
    do ix = 1, Nx
      HorizSpeed(ix,iy) = sqrt(U(ix,iy,iz)**2. + V(ix,iy,iz)**2.)
    enddo
  enddo

  return
end subroutine UvToSpeed

end program tsavg
