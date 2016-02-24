!***************************************************************
! Program to create CFADs
!
! This program will read in HDF5 data from a RAMS simulation, use a filter
! to select data, and create a CFAD.
!

program cfad
  use rhdf5_utils
  use diag_utils

  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10
  integer, parameter :: MaxCoords=1000
  real, parameter :: UndefVal=-999.0

  integer NumRbands
  real :: MaxRadius, WfilterMin, WfilterMax
  character (len=MediumString) :: InDir
  character (len=MediumString) :: InSuffix
  character (len=MediumString) :: OutFile
  character (len=MediumString) :: FilterFile
  character (len=MediumString) :: RadiusFile
  character (len=LittleString) :: VarName, RevuVar
  logical :: DoHorizVel, DoTangential

  integer :: NumBins
  real :: BinStart, BinInc

  type (Rhdf5Var) :: U, V, Avar, CfadOut, Rcoords, Bcoords, Zcoords, Tcoords, Lon, Lat
  type (Rhdf5Var) :: Filter, Radius, StormLon, StormLat
  character (len=MediumString) :: Ufile, Vfile, AvarFile
  integer, dimension(:), allocatable :: StmIx, StmIy
  real, dimension(:), allocatable :: MinP, Xcoords, Ycoords, StormX, StormY

  integer :: i

  integer :: ix, iy, iz, it, ilin, ib
  integer :: Nx, Ny, Nz, Nt

  real :: DeltaX, DeltaY, RbandInc
  real :: TestData, TestGridX, TestGridY

  ! Get the command line arguments
  call GetMyArgs(InDir, InSuffix, OutFile, FilterFile, RadiusFile, NumRbands, &
    NumBins, BinStart, BinInc, VarName, RevuVar, DoHorizVel, DoTangential)

  write (*,*) 'Calculating CFAD for RAMS data:'
  write (*,*) '  Input directory: ', trim(InDir)
  write (*,*) '  Input file name suffix: ', trim(InSuffix)
  write (*,*) '  Output file name:  ', trim(OutFile)
  write (*,*) '  Filter file name: ', trim(FilterFile)
  write (*,*) '  Radius file name:  ', trim(RadiusFile)
  write (*,*) '  Distance between cfad bins: ', BinInc
  if (DoHorizVel) then
    if (DoTangential) then
      write (*,*) '  RAMS variable that is being processed: Tangential Horizontal Velocity'
    else
      write (*,*) '  RAMS variable that is being processed: Radial Horizontal Velocity'
    end if
  else
    write (*,*) '  RAMS variable that is being processed: ', trim(VarName)
    write (*,*) '  REVU variable name: ', trim(RevuVar)
  end if
  write (*,*) ''

  if (trim(VarName) /= 'test') then
    ! Not running a test so grab the data from HDF5 input files and perform the averaging.

    ! Read the variable information from the HDF5 files and check for consistency.
    !
    ! Always use filter and radius files
    !
    ! If doing horizontal velocity, ie VarName is 'speed_t' or 'speed_r', set up U and V
    ! otherwise set up the variable named in VarName
    !

    Filter%vname = 'filter'
    call rhdf5_read_init(FilterFile, Filter)

    Radius%vname = 'radius'
    call rhdf5_read_init(RadiusFile, Radius)

    if (DoHorizVel) then
      Ufile = trim(InDir) // '/u' // trim(InSuffix)
      U%vname = 'u'
      call rhdf5_read_init(Ufile, U)

      Zcoords%vname = 'z_coords'
      call rhdf5_read_init(Ufile, Zcoords)

      Vfile = trim(InDir) // '/v' // trim(InSuffix)
      V%vname = 'v'
      call rhdf5_read_init(Vfile, V)

      ! Initialize the elements in Avar
      Avar%vname = trim(VarName)
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
      AvarFile = trim(InDir) // '/' // trim(VarName) // trim(InSuffix)
      Avar%vname = trim(RevuVar)
      call rhdf5_read_init(AvarFile, Avar)

      Zcoords%vname = 'z_coords'
      call rhdf5_read_init(AvarFile, Zcoords)
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
        write (*,*) 'ERROR: dimensions of filter, radius, and ', trim(VarName), ' do not match'
        stop
      endif
    endif
 
    ! Always read in Filter so use it to record Nx, Ny, Nz, Nt. At this point we have verified
    ! that Nx, Ny, Nz, Nt are consitent for all the variables.
    Nx = Filter%dims(1)
    Ny = Filter%dims(2)
    Nz = Filter%dims(3)
    Nt = Filter%dims(4)

    write (*,*) 'Gridded data information:'
    write (*,*) '  Number of x (longitude) points:          ', Nx
    write (*,*) '  Number of y (latitude) points:           ', Ny
    write (*,*) '  Number of z (vertical level) points:     ', Nz
    write (*,*) '  Number of t (time) points:               ', Nt
    write (*,*) ''
    write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
    write (*,*) ''
  
    write (*,*) ''

    ! Read in the longitude and latitude values, use the Filter file
    Lon%vname = 'x_coords'
    call rhdf5_read_init(FilterFile, Lon)
    call rhdf5_read(FilterFile, Lon)
    
    Lat%vname = 'y_coords'
    call rhdf5_read_init(FilterFile, Lat)
    call rhdf5_read(FilterFile, Lat)

    ! Get the time values
    Tcoords%vname = 't_coords'
    call rhdf5_read_init(FilterFile, Tcoords)
    call rhdf5_read(FilterFile, Tcoords)

    ! read in the radius data
    write (*,*) 'Reading variable: radius'
    write (*,*) '  HDF5 file: ', trim(RadiusFile)
    write (*,*) ''
    call rhdf5_read(RadiusFile, Radius)

    ! Convert the GRADS grid coordinates from longitude, latitude to flat plane (x and y).
    ! Xcoords, Ycoords are in units of km
    allocate(Xcoords(Nx))
    allocate(Ycoords(Ny))
    call ConvertGridCoords(Nx, Ny, Nz, Nt, Lon%vdata, Lat%vdata, Xcoords, Ycoords)

    write (*,*) 'Horzontal Grid Coordinate Info:'
    write (*,*) '  X Range (min lon, max lon) --> (min x, max x): '
    write (*,*) '    ', Lon%vdata(1), Lon%vdata(Nx), Xcoords(1), Xcoords(Nx)
    write (*,*) '  Y Range (min lat, max lat) --> (min y, max y): '
    write (*,*) '    ', Lat%vdata(1), Lat%vdata(Ny), Ycoords(1), Ycoords(Ny)
    write (*,*) ''
  
    ! Figure out how big to make each radial band. Assume that the storm center stays near
    ! the center of the grid domain.
    MaxRadius = 0.0
    do i = 1, Nt*Ny*Nx
      if (Radius%vdata(i) .gt. MaxRadius) then
        MaxRadius = Radius%vdata(i)
      endif
    enddo
    MaxRadius = anint(MaxRadius)  ! round to nearest integer
    DeltaX = Xcoords(Nx) - Xcoords(1)
    DeltaY = Ycoords(Ny) - Ycoords(1)
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

    ! Fill in the Bcoords array (which are the center
    ! values for the CFAD bins)
    !
    Bcoords%vname = 'y_coords'
    Bcoords%ndims = 1
    Bcoords%dims(1) = NumBins
    Bcoords%dimnames(1) = 'y'
    Bcoords%units = 'degrees_north'
    Bcoords%descrip = 'CFAD bins'
    allocate (Bcoords%vdata(NumBins))
  
    Bcoords%vdata(1) = BinStart
    do ib = 2, NumBins
      Bcoords%vdata(ib) = Bcoords%vdata(ib-1) + BinInc
    enddo

    write (*,*) 'Bin information: '
    write (*,*) '  Number of bins: ', NumBins
    write (*,*) '  Bin start value: ', BinStart
    write (*,*) '  Bin increment value: ', BinInc
    write (*,*) '  Bin values:'
    do ib = 1, NumBins
      write (*,*) '    ', ib, ' --> ', Bcoords%vdata(ib)
    end do
    write (*,*) ''
  
    ! StormLon and StormLat are in the RadiusFile
    StormLon%vname = 'min_press_xloc'
    call rhdf5_read_init(RadiusFile, StormLon)
    call rhdf5_read(RadiusFile, StormLon)

    StormLat%vname = 'min_press_yloc'
    call rhdf5_read_init(RadiusFile, StormLat)
    call rhdf5_read(RadiusFile, StormLat)

    ! Convert the lat/lon values to x/y values
    allocate(StormX(Nt))
    allocate(StormY(Nt))
    call ConvertStormCenter(Nx, Ny, Nt, Lon%vdata, Xcoords, StormLon%vdata, StormX, &
      Lat%vdata, Ycoords, StormLat%vdata, StormY)

    ! Dump out the storm center locations
    write (*,*) 'Storm center:'
    do it = 1, Nt
      write (*,'(a,5f15.4)') '  time, lat, lon, x, y: ', Tcoords%vdata(it), StormLon%vdata(it), StormLat%vdata(it), StormX(it), StormY(it)
    enddo
    write (*,*) ''

    ! Read in the data for the vars
    write (*,*) 'Reading variable: filter'
    write (*,*) '  HDF5 file: ', trim(FilterFile)
    write (*,*) ''
    call rhdf5_read(FilterFile, Filter)

    if (DoHorizVel) then
      write (*,*) 'Reading variable: u'
      write (*,*) '  HDF5 file: ', trim(Ufile)
      write (*,*) ''
      call rhdf5_read(Ufile, U)

      call rhdf5_read(Ufile, Zcoords)

      write (*,*) 'Reading variable: v'
      write (*,*) '  HDF5 file: ', trim(Vfile)
      write (*,*) ''
      call rhdf5_read(Vfile, V)

      allocate(Avar%vdata(Nx*Ny*Nz*Nt))
      call ConvertHorizVelocity(Nx, Ny, Nz, Nt, U%vdata, V%vdata, StormX, StormY, &
        Xcoords, Ycoords, Avar%vdata, DoTangential)
    else
      write (*,*) 'Reading variable: ', trim(VarName)
      write (*,*) '  HDF5 file: ', trim(AvarFile)
      write (*,*) ''
      call rhdf5_read(AvarFile, Avar)

      call rhdf5_read(AvarFile, Zcoords)
    end if

  else
    ! Performing a test, load up the variables with data that will produce a known result and
    ! finish out the program (averaging and output).

    ! Create one horizontal slice, 100 x 100 tiles, one z level, 5 time points.
    !   Make the x,y increments 1km
    ! Move the storm center around a little bit but keep it near the center.
    ! Make 10 radial bands, 

    Nx = 101
    Ny = 101
    Nz = 1
    Nt = 5

!    allocate (W%vdata(Nx*Ny*Nz*Nt))

    Zcoords%ndims = 1
    Zcoords%dims(1) = Nz
    Zcoords%dimnames(1) = 'z'
    Zcoords%units = 'meter'
    Zcoords%descrip = 'sigma-z'
    allocate (Zcoords%vdata(Nz))
    do iz = 1, Nz
      Zcoords%vdata(iz) = real(iz)
    enddo

    Tcoords%ndims = 1
    Tcoords%dims(1) = Nt
    Tcoords%dimnames(1) = 't'
    Tcoords%units = 'second'
    Tcoords%descrip = 'seconds since 1970-01-01 00:00:00 00:00'
    allocate (Tcoords%vdata(Nt))
    do it = 1, Nt
      Tcoords%vdata(it) = real(it)
    enddo

    allocate (Avar%vdata(Nx*Ny*Nz*Nt))
    Avar%units = 'test'
    allocate (StmIx(Nt), StmIy(Nt), MinP(Nt))
    allocate (Xcoords(Nx), Ycoords(Ny))

    NumRbands = 10
    TestGridX = 100.0
    TestGridY = 100.0

    MaxRadius = sqrt(TestGridX*TestGridX + TestGridY*TestGridY) / 2.0
    RbandInc = MaxRadius / real(NumRbands)

    Rcoords%vname = 'x_coords'
    Rcoords%ndims = 1
    Rcoords%dims(1) = NumRbands
    Rcoords%dimnames(1) = 'x'
    Rcoords%units = 'degrees_east'
    Rcoords%descrip = 'radius in meters'
    allocate (Rcoords%vdata(NumRbands))
    do ix = 1, NumRbands
      Rcoords%vdata(ix) = real(ix) * 5.0 * 1000.0
    enddo

    ! load up the coordinates
    DeltaX = TestGridX / real(Nx - 1)
    DeltaY = TestGridY / real(Ny - 1)
    do ix = 1, Nx
      Xcoords(ix) = real(ix-1) * DeltaX
    end do
    do iy = 1, Ny
      Ycoords(iy) = real(iy-1) * DeltaY
    end do

    do it = 1, Nt
      ! Storm center drifts from roughly the center of the grid
      StmIx(it) = int(Nx/2) + it
      StmIy(it) = int(Ny/2) - it
      MinP(it) = 980.0 - real(it)  ! Just make it up, it's only used in diagnostic msg

      do iz = 1, Nz
        do iy = 1, Ny
          do ix = 1, Nx
            ! The vdata elements in W and Avar are linear arrays, so need to translate
            ! the multi-dimension indecies to the equivalent linear index
            ilin = ix + ((iy-1)*(Nx)) + ((iz-1)*(Nx*Ny)) + ((it-1)*(Nx*Ny*Nz))

            ! Fill up W with the it value so you can see if screening works with different
            ! WfilterMin, WfilterMax values
!            W%vdata(ilin) = real(it)

            ! Make the data match the radial band number after the averaging takes place.
            !   int(sqrt(DeltaX*DeltaX + DeltaY*DeltaY) / RbandInc) + 1) gives you the radial band number
            !   just multiply it by the variable "it" so you see increasing slope for
            ! successive time steps
            DeltaX = Xcoords(StmIx(it)) - Xcoords(ix)
            DeltaY = Ycoords(StmIy(it)) - Ycoords(iy)
            Avar%vdata(ilin) = real((int(sqrt(DeltaX*DeltaX + DeltaY*DeltaY) / RbandInc) + 1) * it)

          enddo 
        enddo 
      enddo 
    enddo 
  endif
  
  ! Do the averaging and write the result into an HDF5 file
  CfadOut%vname = trim(VarName)
  CfadOut%ndims = 4
  CfadOut%dims(1) = NumRbands
  CfadOut%dims(2) = NumBins
  CfadOut%dims(3) = Nz
  CfadOut%dims(4) = Nt
  CfadOut%dimnames(1) = 'x'
  CfadOut%dimnames(2) = 'y'
  CfadOut%dimnames(3) = 'z'
  CfadOut%dimnames(4) = 't'
  CfadOut%units = Avar%units 
  CfadOut%descrip = 'azimuthally averaged ' // trim(VarName) 
  allocate(CfadOut%vdata(NumRbands*NumBins*Nz*Nt))

!  call AzimuthalAverage(Nx, Ny, Nz, Nt, Nz, NumRbands, Avar%vdata, CfadOut%vdata, &
!    RbandInc, Filter%vdata, Radius%vdata, UndefVal)
  call BuildCfad(Nx, Ny, Nz, Nt, NumRbands, NumBins, Avar%vdata, CfadOut%vdata, &
    RbandInc, Bcoords%vdata, Filter%vdata, Radius%vdata, UndefVal)

  ! third arg to rhdf5_write is "append" flag:
  !   0 - create new file
  !   1 - append to existing file
  write (*,*) 'Writing HDF5 output: ', trim(OutFile)
  write (*,*) ''
  call rhdf5_write(OutFile, CfadOut, 0)

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
  call rhdf5_attach_dimensions(OutFile, CfadOut)

  stop

Contains

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the command line arguments
!

subroutine GetMyArgs(InDir, InSuffix, OutFile, FilterFile, RadiusFile, NumRbands, &
    NumBins, BinStart, BinInc, VarName, RevuVar, DoHorizVel, DoTangential)

  implicit none

  integer :: NumRbands, NumBins
  character (len=*) :: InDir, InSuffix, OutFile, FilterFile, RadiusFile, VarName, RevuVar
  logical :: DoHorizVel, DoTangential
  real :: BinStart, BinInc

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 11) then
    write (*,*) 'ERROR: must supply exactly 11 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: cfad <in_dir> <in_suffix> <out_data_file> <filter_file> <radius_file> <num_radial_bands> <num_cfad_bins> <bin_start> <bin_inc> <var_to_average> <revu_var_name>'
    write (*,*) '        <in_dir>: directory where input files live'
    write (*,*) '        <in_suffix>: suffix on input file names'
    write (*,*) '        <out_file>: name of output file, HDF5 format'
    write (*,*) '        <filter_file_list>: file containing the filter mask'
    write (*,*) '        <radius_file>: name of file containing radius information, HDF5 format'
    write (*,*) '        <num_radial_bands>: number of bands to split data into'
    write (*,*) '        <num_cfad_bins>: number of bins in the cfad'
    write (*,*) '        <bin_start>: value associated with first bin in cfad'
    write (*,*) '        <bin_inc>: distance between bins in cfad'
    write (*,*) '        <var_to_average>: name of variable to do the averaging on'
    write (*,*) '             special var names:'
    write (*,*) '                speed_t: horizontal tangential speed'
    write (*,*) '                speed_r: horizontal radial speed'
    write (*,*) '        <revu_var_name>: name of variable in the REVU file'
    write (*,*) ''
    stop
  end if

  call getarg(1, InDir)
  call getarg(2, InSuffix)
  call getarg(3, OutFile)
  call getarg(4, FilterFile)
  call getarg(5, RadiusFile)

  call getarg(6, arg)       !this is a string
  read (arg, '(i)') NumRbands !convert to integer

  call getarg(7, arg)
  read (arg, '(i)') NumBins

  call getarg(8, arg)
  read (arg, '(f)') BinStart

  call getarg(9, arg)
  read (arg, '(f)') BinInc

  call getarg(10, VarName)
  if (VarName == 'speed_t') then
    DoHorizVel = .true.
    DoTangential = .true.
  else
    if (VarName == 'speed_r') then
      DoHorizVel = .true.
      DoTangential = .false.
    else
      DoHorizVel = .false.
      DoTangential = .false.
    end if
  end if

  call getarg(11, RevuVar)

  return
end subroutine GetMyArgs

end program cfad
