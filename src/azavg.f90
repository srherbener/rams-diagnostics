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
  use azavg_utils
  use revu_utils

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
  character (len=LittleString) :: VarToAvg, RevuVar
  logical :: DoHorizVel, DoTangential

  ! Data arrays: need one for w (vertical velocity), press (pressure)
  ! and the var we are doing the averaging on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number

  type (Rhdf5Var) :: U, V, W, Press, Avar, Aavg
  character (len=MediumString) :: Ufile, Vfile, Wfile, PressFile, AvarFile, AavgFile
  integer, dimension(:), allocatable :: StmIx, StmIy
  real, dimension(:), allocatable :: MinP, Lon, Lat, Xcoords, Ycoords

  integer :: i
  logical :: VarIs2d

  integer :: ix, iy, iz, it, ilin
  integer :: Nx, Ny, Nz, Nt, VarNz

  real :: DeltaX, DeltaY, RbandInc
  real :: TestData, TestGridX, TestGridY

  ! Get the command line arguments
  call GetMyArgs(InDir, InSuffix, OutFile, NumRbands, MaxRadius, WfilterMin, WfilterMax, VarToAvg, &
    RevuVar, DoHorizVel, DoTangential, VarIs2d)

  write (*,*) 'Calculating azimuthal average for RAMS data:'
  write (*,*) '  Input directory: ', trim(InDir)
  write (*,*) '  Input file name suffix: ', trim(InSuffix)
  write (*,*) '  Output file name:  ', trim(OutFile)
  write (*,*) '  Number of radial bands: ', NumRbands
  write (*,*) '  Maximum radius: ', MaxRadius
  write (*,*) '  W filter interval :     ', WfilterMin, WfilterMax
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
  if (VarIs2d) then
    write (*,*) '    Variable contains 2D data'
  else
    write (*,*) '    Variable contains 3D data'
  end if
  write (*,*) ''


  if (trim(VarToAvg) /= 'test') then
    ! Not running a test so grab the data from HDF5 input files and perform the averaging.

    ! Read the variable information from the HDF5 files and check for consistency.
    !
    ! Always use w (filtering) and press (find storm center).
    ! If doing horizontal velocity, ie VarToAvg is 'speed_t' or 'speed_r', set up U and V
    ! otherwise set up the variable named in VarToAvg
    !
    Wfile = trim(InDir) // '/w' // trim(InSuffix)
    W%vname = 'w'
    call rhdf5_read_init(Wfile, W)

    PressFile = trim(InDir) // '/press' // trim(InSuffix)
    Press%vname = 'press'
    call rhdf5_read_init(PressFile, Press)

    if (DoHorizVel) then
      Ufile = trim(InDir) // '/u' // trim(InSuffix)
      U%vname = 'u'
      call rhdf5_read_init(Ufile, U)

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
    endif

    ! check that the variable dimensions (size and coordinate values) match up, if this
    ! isn't true, then the subsequent anlysis will be bogus
    !
    ! third arg of GvarDimsMatch() is true if one of the vars is 2D, else false
    if (DoHorizVel) then
      if (.not. (DimsMatch(W, U) .and. DimsMatch(W, V) .and. &
                 DimsMatch(W, Press))) then
        write (*,*) 'ERROR: dimensions of u, v, w, and press do not match'
        stop
      endif
    else
      if (.not. (DimsMatch(W, Press) .and. DimsMatch(W, Avar))) then
        write (*,*) 'ERROR: dimensions of w, press, and ', trim(VarToAvg), ' do not match'
        stop
      endif
    endif
 
    ! Always read in W so use it to record Nx, Ny, Nz, Nt. At this point we have verified
    ! that Nx, Ny, Nz, Nt are consitent for all the variables.
    Nx = W%dims(1)
    Ny = W%dims(2)
    Nz = W%dims(3)
    Nt = W%dims(4)
    if (Avar%ndims .eq. 4) then
      VarNz = Avar%dims(3)
    else
      VarNz = 1
    endif

    write (*,*) 'Gridded data information:'
    write (*,*) '  Number of x (longitude) points:            ', Nx
    write (*,*) '  Number of y (latitude) points:             ', Ny
    write (*,*) '  Number of z (vertical level) points (3D):  ', Nz
    write (*,*) '  Number of z (vertical level) points (var): ', VarNz
    write (*,*) '  Number of t (time) points:                 ', Nt
    write (*,*) ''
    write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
    write (*,*) ''
  
    write (*,*) ''

    ! Read in the longitude and latitude values, use the W file
    allocate(Lon(Nx))
    allocate(Lat(Ny))
    call GetLonLat(Nx, Ny, Lon, Lat, Wfile)

    ! Convert the GRADS grid coordinates from longitude, latitude to flat plane (x and y).
    ! Always have W which has been verified to be consistent with all other vars
    ! so use W for this call
    ! Xcoords, Ycoords are in units of km
    allocate(Xcoords(Nx))
    allocate(Ycoords(Ny))
    call ConvertGridCoords(Nx, Ny, Nz, Nt, Lon, Lat, Xcoords, Ycoords)

    write (*,*) 'Horzontal Grid Coordinate Info:'
    write (*,*) '  X Range (min lon, max lon) --> (min x, max x): '
    write (*,*) '    ', Lon(1), Lon(Nx), Xcoords(1), Xcoords(Nx)
    write (*,*) '  Y Range (min lat, max lat) --> (min y, max y): '
    write (*,*) '    ', Lat(1), Lat(Ny), Ycoords(1), Ycoords(Ny)
    write (*,*) ''
  
    ! Figure out how big to make each radial band. Assume that the storm center stays near
    ! the center of the grid domain.
    DeltaX = Xcoords(Nx) - Xcoords(1)
    DeltaY = Ycoords(Ny) - Ycoords(1)
    RbandInc = MaxRadius / real(NumRbands)
  
    write (*,*) 'Radial band information:'
    write (*,*) '  Delta x of domain:     ', DeltaX
    write (*,*) '  Delta y of domain:     ', DeltaY
    write (*,*) '  Radial distance:       ', MaxRadius
    write (*,*) '  Radial band increment: ', RbandInc
    write (*,*) ''
  
    ! Read in the data for the vars
    write (*,*) 'Reading variable: w'
    write (*,*) '  HDF5 file: ', trim(Wfile)
    write (*,*) ''
    allocate(W%vdata(Nx*Ny*Nz*Nt))
    call rhdf5_read(Wfile, W)

    write (*,*) 'Reading variable: press'
    write (*,*) '  HDF5 file: ', trim(PressFile)
    write (*,*) ''
    allocate(Press%vdata(Nx*Ny*Nz*Nt))
    call rhdf5_read(PressFile, Press)

    allocate(StmIx(Nt))
    allocate(StmIy(Nt))
    allocate(MinP(Nt))
    call RecordStormCenter(Nx, Ny, Nz, Nt, Press%vdata, StmIx, StmIy, MinP)

    if (DoHorizVel) then
      write (*,*) 'Reading variable: u'
      write (*,*) '  HDF5 file: ', trim(Ufile)
      write (*,*) ''
      allocate(U%vdata(Nx*Ny*Nz*Nt))
      call rhdf5_read(Ufile, U)

      write (*,*) 'Reading variable: v'
      write (*,*) '  HDF5 file: ', trim(Vfile)
      write (*,*) ''
      allocate(V%vdata(Nx*Ny*Nz*Nt))
      call rhdf5_read(Vfile, V)

      allocate(Avar%vdata(Nx*Ny*VarNz*Nt))
      call ConvertHorizVelocity(Nx, Ny, VarNz, Nt, U%vdata, V%vdata, StmIx, StmIy, &
        Xcoords, Ycoords, Avar%vdata, DoTangential)
    else
      write (*,*) 'Reading variable: ', trim(VarToAvg)
      write (*,*) '  HDF5 file: ', trim(AvarFile)
      write (*,*) ''
      allocate(Avar%vdata(Nx*Ny*VarNz*Nt))
      call rhdf5_read(AvarFile, Avar)
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
    VarNz = Nz
    Nt = 5

    allocate (W%vdata(Nx*Ny*Nz*Nt))
    allocate (Avar%vdata(Nx*Ny*VarNz*Nt))
    Avar%units = 'test'
    allocate (StmIx(Nt), StmIy(Nt), MinP(Nt))
    allocate (Xcoords(Nx), Ycoords(Ny))

    NumRbands = 10
    TestGridX = 100.0
    TestGridY = 100.0

    MaxRadius = sqrt(TestGridX*TestGridX + TestGridY*TestGridY) / 2.0
    RbandInc = MaxRadius / real(NumRbands)

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
            W%vdata(ilin) = real(it)

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
  
  ! Do the averaging and write the result into GRADS files
  Aavg%vname = trim(VarToAvg)
  Aavg%ndims = 3
  Aavg%dims(1) = NumRbands
  Aavg%dims(2) = VarNz
  Aavg%dims(3) = Nt
  Aavg%dimnames(1) = 'r'
  Aavg%dimnames(2) = 'z'
  Aavg%dimnames(3) = 't'
  Aavg%units = Avar%units 
  Aavg%descrip = 'azimuthally averaged ' // trim(VarToAvg) 
  allocate(Aavg%vdata(NumRbands*VarNz*Nt))

  call AzimuthalAverage(Nx, Ny, Nz, Nt, VarNz, NumRbands, W%vdata, StmIx, StmIy, MinP, &
    Avar%vdata, Aavg%vdata, Xcoords, Ycoords, MaxRadius, RbandInc, WfilterMin, WfilterMax, UndefVal)

  write (*,*) 'Writing HDF5 output: ', trim(OutFile)
  write (*,*) ''
  call rhdf5_write(OutFile, Aavg)

  stop

Contains

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the command line arguments
!

subroutine GetMyArgs(InDir, InSuffix, OutFile, NumRbands, MaxRadius, WfilterMin, WfilterMax, &
  VarToAvg, RevuVar, DoHorizVel, DoTangential, VarIs2d)

  implicit none

  integer :: NumRbands
  real :: MaxRadius, WfilterMin, WfilterMax
  character (len=*) :: InDir, InSuffix, OutFile, VarToAvg, RevuVar
  logical :: DoHorizVel, DoTangential, VarIs2d

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 10) then
    write (*,*) 'ERROR: must supply exactly 10 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_dir> <in_suffix> <out_data_file> <num_radial_bands> <max_radius> <w_fliter_min> <w_filter_max> <var_to_average> <dim_of_var>'
    write (*,*) '        <in_dir>: directory where input files live'
    write (*,*) '        <in_suffix>: suffix on input file names'
    write (*,*) '        <out_file>: name of output file, HDF5 format'
    write (*,*) '        <num_radial_bands>: number of bands to split data into'
    write (*,*) '        <max_radius>: maximum radius from storm center'
    write (*,*) '        <w_filter_min>: min (left end) of interval used to filter out data'
    write (*,*) '        <w_filter_max>: max (right end) of interval used to filter out data'
    write (*,*) '        <var_to_average>: name of variable to do the averaging on'
    write (*,*) '        <revu_var_name>: name of variable in the REVU file'
    write (*,*) '        <dim_of_var>: indicates if <var_to_average> is 2d or 3d'
    write (*,*) '           <dim_of_var> must be either "2d" or "3d"'
    write (*,*) ''
    stop
  end if

  call getarg(1, InDir)
  call getarg(2, InSuffix)
  call getarg(3, OutFile)

  call getarg(4, arg)       !this is a string
  read (arg, '(i)') NumRbands !convert to integer

  call getarg(5, arg)
  read (arg, '(f)') MaxRadius

  call getarg(6, arg)
  read (arg, '(f)') WfilterMin

  call getarg(7, arg)
  read (arg, '(f)') WfilterMax
  if (WfilterMin .ge. WfilterMax) then
    write (*,*) 'ERROR: must specify <w_filter_min> to be less than <w_filter_max>'
    stop
  endif

  call getarg(8, VarToAvg)
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

  call getarg(9, RevuVar)

  call getarg(10, arg)
  if (arg .eq. '2d') then
    VarIs2d = .true.
  else if (arg .eq. '3d') then
    VarIs2d = .false.
  else
    write (*,*) 'ERROR: must use either "2d" or "3d" for <dim_of_var> argument'
    stop
  end if

  return
end subroutine GetMyArgs

end program azavg
