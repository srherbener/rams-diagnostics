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

  type FileSpec
     character (len=RHDF5_MAX_STRING) :: fname
     character (len=RHDF5_MAX_STRING) :: vname
  endtype

  type ArgList
    logical :: DoHist
    logical :: DoHorizVel
    logical :: DoTangential
    logical :: SubtractSmotion
    integer :: NumRbands
    real :: MaxRadius
    type (FileSpec) :: Output
    type (FileSpec) :: Filter
    type (FileSpec) :: StormCenter
    type (FileSpec) :: Bins
    type (FileSpec) :: U
    type (FileSpec) :: V
    type (FileSpec) :: Avar
  endtype

  type (ArgList) :: Args
  integer :: NumBins
  real, dimension(:), allocatable :: HistBins

  ! Data arrays: need one for w (vertical velocity), press (pressure)
  ! and the var we are doing the averaging on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number

  type (Rhdf5Var) :: U, V, Ustorm, Vstorm, Avar, Aavg, Rcoords, Zcoords, Tcoords, Bcoords, Xcoords, Ycoords
  type (Rhdf5Var) :: Filter, Radius, StormXindex, StormYindex
  character (len=RHDF5_MAX_STRING) :: rh5f_facc
  integer :: rh5f_filter, rh5f_center, rh5f_u, rh5f_v, rh5f_avar, rh5f_aavg
  integer :: i_storm, j_storm
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm
  real :: StormX, StormY

  integer :: i
  integer :: AvarNelems

  integer :: ix, iy, it
  integer :: Nx, Ny, Nz, Nt, FilterNz

  real :: DeltaX, DeltaY, RbandInc

  ! Get the command line arguments
  call GetMyArgs(Args)

  write (*,'("Calculating azimuthal average for RAMS data:")')
  write (*,'("  Number of radial bands: ",i4)') Args%NumRbands
  write (*,'("  Maximum radius: ",f10.2)') Args%MaxRadius
  write (*,'("  Subtracting storm motion from input winds: ",l)') Args%SubtractSmotion
  write (*,'("  Input file specs:")')
  if (Args%DoHorizVel) then
    write (*,'("    U: ",a," (",a,")")') trim(Args%U%fname), trim(Args%U%vname)
    write (*,'("    V: ",a," (",a,")")') trim(Args%V%fname), trim(Args%V%vname)
  else
    write (*,'("    Avar: ",a," (",a,")")') trim(Args%Avar%fname), trim(Args%Avar%vname)
  endif
  write (*,'("    Output: ",a," (",a,")")') trim(Args%Output%fname), trim(Args%Output%vname)
  write (*,'("    Filter: ",a," (",a,")")') trim(Args%Filter%fname), trim(Args%Filter%vname)
  write (*,'("    Storm center: ",a," (",a,")")') trim(Args%StormCenter%fname), trim(Args%StormCenter%vname)

  if (Args%DoHist) then
    write (*,'("  Creating histogram counts:")')
    write (*,'("    Bins file: ",a)') trim(Args%Bins%fname)
  endif

  write (*,*) ''

  ! Read the variable information from the HDF5 files and check for consistency.
  !
  ! Always use filter file
  !
  ! Get variable and coordinate descriptions from the filter file
  ! Get the z coordinate description from the variable file since this
  ! can change between 2D and 3D variables.
  Filter%vname = trim(Args%Filter%vname)
  Radius%vname = 'radius'
  StormXindex%vname = 'press_cent_x_index'
  StormYindex%vname = 'press_cent_y_index'
  Xcoords%vname = 'x_coords'
  Ycoords%vname = 'y_coords'
  Tcoords%vname = 't_coords'
  call rhdf5_read_init(Args%Filter%fname, Filter)
  call rhdf5_read_init(Args%Filter%fname, Xcoords)
  call rhdf5_read_init(Args%Filter%fname, Ycoords)
  call rhdf5_read_init(Args%Filter%fname, Tcoords)

  call rhdf5_read_init(Args%StormCenter%fname, Radius)
  call rhdf5_read_init(Args%StormCenter%fname, StormXindex)
  call rhdf5_read_init(Args%StormCenter%fname, StormYindex)

  ! Read in the 1D variable data
  call rhdf5_read(Args%Filter%fname, Xcoords)
  call rhdf5_read(Args%Filter%fname, Ycoords)
  call rhdf5_read(Args%Filter%fname, Tcoords)

  call rhdf5_read(Args%StormCenter%fname, StormXindex)
  call rhdf5_read(Args%StormCenter%fname, StormYindex)

  if (Args%DoHorizVel) then
    U%vname = trim(Args%U%vname)
    call rhdf5_read_init(Args%U%fname, U)

    Zcoords%vname = 'z_coords'
    call rhdf5_read_init(Args%U%fname, Zcoords)
    call rhdf5_read(Args%U%fname, Zcoords)

    V%vname = trim(Args%V%vname)
    call rhdf5_read_init(Args%V%fname, V)

    ! Initialize the elements in Avar
    Avar%ndims = U%ndims
    Avar%dims = U%dims
    Avar%dimnames = U%dimnames
    Avar%units = 'm/s'
    if (Args%DoTangential) then
      Avar%vname = 'speed_t'
      Avar%descrip = 'tangential wind speed'
    else
      Avar%vname = 'speed_r'
      Avar%descrip = 'radial wind speed'
    endif
  else
    Avar%vname = trim(Args%Avar%vname)
    call rhdf5_read_init(Args%Avar%fname, Avar)

    Zcoords%vname = 'z_coords'
    call rhdf5_read_init(Args%Avar%fname, Zcoords)
    call rhdf5_read(Args%Avar%fname, Zcoords)
  endif

  ! If subtracting storm motion, initialize the storm u and v vars
  if (Args%SubtractSmotion) then
    Ustorm%vname = 'storm_speed_x'
    call rhdf5_read_init(Args%StormCenter%fname, Ustorm)

    Vstorm%vname = 'storm_speed_y'
    call rhdf5_read_init(Args%StormCenter%fname, Vstorm)
  endif

  ! check that the variable dimensions (size and coordinate values) match up, if this
  ! isn't true, then the subsequent anlysis will be bogus
  !
  ! third arg of GvarDimsMatch() is true if one of the vars is 2D, else false
  if (Args%DoHorizVel) then
    if (.not. (DimsMatch(Filter, U) .and. DimsMatch(Filter, V) .and. &
               DimsMatch(Filter, Radius))) then
      write (*,*) 'ERROR: dimensions of u, v, filter, and radius do not match'
      stop
    endif
  else
    if (.not. (DimsMatch(Filter, Radius) .and. DimsMatch(Filter, Avar))) then
      write (*,*) 'ERROR: dimensions of filter, radius, and ', trim(Avar%vname), ' do not match'
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
  if (Args%DoHist) then
    ! ReadBinsFile will allocate HistBins (it needs to read in NumBins before
    ! HistBins gets allocated).
    call ReadBinsFile(Args%Bins%fname, NumBins, HistBins)
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

  Args%MaxRadius = anint(Args%MaxRadius)           ! round to nearest integer
  DeltaX = XcoordsKm(Nx) - XcoordsKm(1)
  DeltaY = YcoordsKm(Ny) - YcoordsKm(1)
  RbandInc = Args%MaxRadius / real(Args%NumRbands)

  Rcoords%vname = 'x_coords'
  Rcoords%ndims = 1
  Rcoords%dims(1) = Args%NumRbands
  Rcoords%dimnames(1) = 'x'
  Rcoords%units = 'degrees_east'
  Rcoords%descrip = 'radius in meters'
  allocate (Rcoords%vdata(Args%NumRbands))
  do ix = 1, Args%NumRbands
    Rcoords%vdata(ix) = real(ix) * RbandInc * 1000.0
  enddo

  write (*,*) 'Radial band information:'
  write (*,*) '  Delta x of domain:     ', DeltaX
  write (*,*) '  Delta y of domain:     ', DeltaY
  write (*,*) '  Radial distance:       ', Args%MaxRadius
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
  call rhdf5_open_file(Args%Filter%fname, rh5f_facc, 0, rh5f_filter)
  call rhdf5_open_file(Args%StormCenter%fname, rh5f_facc, 0, rh5f_center)

  Filter%ndims = Filter%ndims - 1
  Radius%ndims = Radius%ndims - 1

  ! set up the input variable data
  rh5f_facc = 'R'
  if (Args%DoHorizVel) then
    call rhdf5_open_file(Args%U%fname, rh5f_facc, 0, rh5f_u)
    call rhdf5_open_file(Args%V%fname, rh5f_facc, 0, rh5f_v)
    U%ndims = U%ndims - 1
    V%ndims = V%ndims - 1
  else
    call rhdf5_open_file(Args%Avar%fname, rh5f_facc, 0, rh5f_avar)
  endif

  Avar%ndims = Avar%ndims - 1
  AvarNelems = 1
  do i = 1, Avar%ndims
    AvarNelems = AvarNelems * Avar%dims(i)
  enddo

  ! set up the output variable
  rh5f_facc = 'W'
  call rhdf5_open_file(Args%Output%fname, rh5f_facc, 1, rh5f_aavg)
  write (*,*) 'Writing HDF5 output: ', trim(Args%Output%fname)
  write (*,*) ''

  ! NumBins will be one when not doing histograms. Otherwise it will be
  ! the number of bins the user specified.

  Aavg%vname = trim(Args%Output%vname)
  Aavg%ndims = 3
  Aavg%dims(1) = Args%NumRbands
  Aavg%dims(2) = NumBins
  Aavg%dims(3) = Nz
  Aavg%dimnames(1) = 'x'
  Aavg%dimnames(2) = 'y'
  Aavg%dimnames(3) = 'z'
  Aavg%units = Avar%units 
  Aavg%descrip = 'azimuthally averaged ' // trim(Args%Output%vname) 
  allocate(Aavg%vdata(Aavg%dims(1)*Aavg%dims(2)*Aavg%dims(3)))

  ! If subtracting out storm motion, read in the storm u and v
  ! Read in entire dataset (time step number == 0, Argument 4 to rhdf5_read_variable). Because
  ! we are doing this outside the time step loop, don't alter the number of dimensions on [UV]storm.
  if (Args%SubtractSmotion) then
    call rhdf5_read_variable(rh5f_center, Ustorm%vname, Ustorm%ndims, 0, Ustorm%dims, rdata=Ustorm%vdata)
    call rhdf5_read_variable(rh5f_center, Vstorm%vname, Vstorm%ndims, 0, Vstorm%dims, rdata=Vstorm%vdata)
  endif

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

    if (Args%DoHorizVel) then
      call rhdf5_read_variable(rh5f_u, U%vname, U%ndims, it, U%dims, rdata=U%vdata)
      call rhdf5_read_variable(rh5f_v, V%vname, V%ndims, it, V%dims, rdata=V%vdata)

      ! If subtracting storm motion, adjust u and v
      if (Args%SubtractSmotion) then
        call TranslateField(Nx, Ny, Nz, U%vdata, Ustorm%vdata(it))
        call TranslateField(Nx, Ny, Nz, V%vdata, Vstorm%vdata(it))
      endif

      ! convert u,v to tangential or radial
      allocate(Avar%vdata(AvarNelems))
      call ConvertHorizVelocity(Nx, Ny, Nz, U%vdata, V%vdata, StormX, StormY, &
        XcoordsKm, YcoordsKm, Avar%vdata, Args%DoTangential)

      ! Free up variable memory
      deallocate(U%vdata)
      deallocate(V%vdata)
    else
      call rhdf5_read_variable(rh5f_avar, Avar%vname, Avar%ndims, it, Avar%dims, rdata=Avar%vdata)

      ! If subtracting storm motion, adjust u and v
      if (Args%SubtractSmotion) then
        if (trim(Avar%vname) .eq. 'u') then
          call TranslateField(Nx, Ny, Nz, Avar%vdata, Ustorm%vdata(it))
        elseif (trim(Avar%vname) .eq. 'v') then
          call TranslateField(Nx, Ny, Nz, Avar%vdata, Vstorm%vdata(it))
        endif
      endif
    endif

    ! do the averaging and write out the results to the output file
    call AzimuthalAverage(Nx, Ny, Nz, FilterNz, Args%NumRbands, NumBins, Avar%vdata, Aavg%vdata, &
      RbandInc, Filter%vdata, Radius%vdata, HistBins, Args%DoHist, UndefVal)

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
      if (Args%SubtractSmotion) then
        write (*,'(a,4f15.4)') '   Storm motion: u, v: ', Ustorm%vdata(it), Vstorm%vdata(it)
      endif
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
  if (Args%DoHorizVel) then
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
  
  call rhdf5_write(Args%Output%fname, Rcoords, 1)
  call rhdf5_write(Args%Output%fname, Bcoords, 1)
  call rhdf5_write(Args%Output%fname, Zcoords, 1)
  call rhdf5_write(Args%Output%fname, Tcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(Args%Output%fname, Rcoords, 'x')
  call rhdf5_set_dimension(Args%Output%fname, Bcoords, 'y')
  call rhdf5_set_dimension(Args%Output%fname, Zcoords, 'z')
  call rhdf5_set_dimension(Args%Output%fname, Tcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(Args%Output%fname, Aavg)

  stop

Contains

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
  character (len=MediumString), dimension(MaxArgFields) :: ArgList
  integer :: Nfields

  ! default values
  Args%DoHist = .false.
  Args%DoHorizVel = .false.
  Args%DoTangential = .false.
  Args%NumRbands = 50
  Args%MaxRadius = 500.0
  Args%SubtractSmotion = .false.

  ! initialization
  Args%Output%fname = 'none'
  Args%Output%vname = 'none'

  Args%Filter%fname = 'none'
  Args%Filter%vname = 'none'

  Args%StormCenter%fname = 'none'
  Args%StormCenter%vname = 'none'

  Args%Bins%fname = 'none'
  Args%Bins%vname = 'none'

  Args%U%fname = 'none'
  Args%U%vname = 'none'

  Args%V%fname = 'none'
  Args%V%vname = 'none'

  Args%Avar%fname = 'none'
  Args%Avar%vname = 'none'

  ! loop through all of the command line tokens
  ! optarg is a character string variable that the getoptions module supplies
  !   optarg gets set to the argument value for an option that uses an argument
  ! getopt returns single character:
  !    '>' finished
  !    '!' unknown option
  !    '.' command line argument (no option)
  BadArgs = .false.
  do
    optval = getopt('b:hmr:t:')

    select case (optval)
      case ('>')  ! finished
        exit

      case ('!')  ! unrecognized argument
        write(*,*) 'ERROR: unknown option: ', trim(optarg)
        write(*,*) ''
        BadArgs = .true.

      case ('b')
        read(optarg, '(i)') Args%NumRbands

      case ('h')
        Args%DoHist = .true.

      case ('m')
        Args%SubtractSmotion = .true.

      case ('r')
        read(optarg, '(f)') Args%MaxRadius

      case ('t')
        Args%DoHorizVel = .true.
        select case (optarg)
          case ('h_tan')
            Args%DoTangential = .true.

          case ('h_rad')
            Args%DoTangential = .false.

          case default
            write(*,*) 'ERROR: must use one of "h_tan" or "h_rad" for argement to -t option: ', trim(optarg)
            write(*,*) ''
            BadArgs = .true.
        endselect

      case ('.')
        ! file spec ->   <file_type>:<file_name>:<variable_name>
        call String2List(optarg, ':', ArgList, MaxArgFields, Nfields, 'file spec') 
        if (Nfields .eq. 3) then
          select case (ArgList(1))
            case ('filter')
              Args%Filter%fname = trim(ArgList(2))
              Args%Filter%vname = trim(ArgList(3))

            case ('storm')
              Args%StormCenter%fname = trim(ArgList(2))
              Args%StormCenter%vname = trim(ArgList(3))

            case ('bins')
              Args%Bins%fname = trim(ArgList(2))
              Args%Bins%vname = trim(ArgList(3))

            case ('avar')
              Args%Avar%fname = trim(ArgList(2))
              Args%Avar%vname = trim(ArgList(3))

            case ('u')
              Args%U%fname = trim(ArgList(2))
              Args%U%vname = trim(ArgList(3))

            case ('v')
              Args%V%fname = trim(ArgList(2))
              Args%V%vname = trim(ArgList(3))

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
    endselect
  enddo

  ! check for required files
  if ((trim(Args%Filter%fname) .eq. 'none') .or. (trim(Args%Output%fname) .eq. 'none') .or. &
      (trim(Args%StormCenter%fname) .eq. 'none')) then
    write(*,*) 'ERROR: must always specify filter, out, and storm files'
    write(*,*) ''
    BadArgs = .true.
  endif

  if (Args%DoHist) then
    if (trim(Args%Bins%fname) .eq. 'none') then
      write(*,*) 'ERROR: must specify bins file when using -h option'
      write(*,*) ''
      BadArgs = .true.
    endif
  endif

  if (Args%DoHorizVel) then
    if ((trim(Args%U%fname) .eq. 'none') .or. (trim(Args%V%fname) .eq. 'none')) then
      write(*,*) 'ERROR: must specify u and v files when using -t option'
      write(*,*) ''
      BadArgs = .true.
    endif
  else
    if (trim(Args%Avar%fname) .eq. 'none') then
      write(*,*) 'ERROR: must specify avar file when not using -t option'
      write(*,*) ''
      BadArgs = .true.
    endif
  endif

  if (BadArgs) then
    write (*,*) 'USAGE: azavg [-b <num_radial_bands>] [-r <max_radius>] [-h] [-m] [-t <input_type>] <list_of_files>'
    write (*,*) ''
    write (*,*) '   -b: split radial distance into <num_radial_bands> bands'
    write (*,*) '   -r: select data within <max_radius> for averaging (determines size of radial bands)'
    write (*,*) '   -h: create histogram instead of average (requires a bins file)'
    write (*,*) '   -m: subtract storm motion from input winds'
    write (*,*) '   -t: for input data processing, valid values of <input_type>:'
    write (*,*) '     "h_tan": do tangential speed of horizontal winds'
    write (*,*) '     "h_rad": do radial speed of horizontal winds'
    write (*,*) '     both of these values require u and v files'
    write (*,*) ''
    write (*,*) '   <list_of_files>: supply specs for input and output files'
    write (*,*) '     <list_of_files> := <file_spec> ...'
    write (*,*) '       <file_spec> := <file_type>:<file_name>:<variable_name>'
    write (*,*) '          <file_type> is one of:'
    write (*,*) '             "filter": filter mask, required'
    write (*,*) '             "storm": storm track and motion, required'
    write (*,*) '             "bins": histogram bin edges, required only if -h used'
    write (*,*) '             "avar": variable to average, required except when -t is used'
    write (*,*) '             "u": zonal component of momentum, required only if -t used'
    write (*,*) '             "v": meridional component of momentum, required only if -t used'
    write (*,*) '             "out": output file, required'
    write (*,*) '          <file_name> is the complete path to the file'
    write (*,*) '          <variable_name> is the name of the variable inside the file'
    write (*,*) ''
    stop
  endif

  return
end subroutine GetMyArgs

end program azavg
