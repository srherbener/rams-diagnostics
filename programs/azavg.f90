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

program main
  use gdata_utils
  use azavg_utils
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128

  integer NumRbands
  real :: WfilterMin, WfilterMax
  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarToAvg
  logical :: DoHorizVel, DoTangential

  type (GradsControlFiles) :: GctlFiles

  ! Data arrays: need one for w (vertical velocity), press (pressure)
  ! and the var we are doing the averaging on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  type (GradsVar) :: U, V, W, Press, Avar, AzAvg
  integer, dimension(:), allocatable :: StmIx, StmIy
  real, dimension(:), allocatable :: MinP, Xcoords, Ycoords

  integer :: i
  logical :: VarIs2d

  integer :: ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt, VarNz

  real :: DeltaX, DeltaY, RadialDist, RbandInc
  real :: TestData, TestGridX, TestGridY
  real, dimension(1:MaxCoords) :: TestZcoords

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBase, NumRbands, WfilterMin, WfilterMax, VarToAvg, DoHorizVel, DoTangential, VarIs2d)
  call String2List(Infiles, ':', GctlFiles%Fnames, MaxFiles, GctlFiles%Nfiles, 'input files')

  write (*,*) 'Calculating azimuthal average for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, GctlFiles%Nfiles
    write (*,*) '  ', i, ': ', trim(GctlFiles%Fnames(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  Number of radial bands: ', NumRbands
  write (*,*) '  W filter interval :     ', WfilterMin, WfilterMax
  if (DoHorizVel) then
    if (DoTangential) then
      write (*,*) '  RAMS variable that is being averaged: Tangential Horizontal Velocity'
    else
      write (*,*) '  RAMS variable that is being averaged: Radial Horizontal Velocity'
    end if
  else
    write (*,*) '  RAMS variable that is being averaged: ', trim(VarToAvg)
  end if
  if (VarIs2d) then
    write (*,*) '    Variable contains 2D data'
  else
    write (*,*) '    Variable contains 3D data'
  end if
  write (*,*) ''

  if (trim(VarToAvg) /= 'test') then
    ! Not running a test so grab the data from GRADS input files and perform the averaging.

    ! Read the GRADS data description files and collect the information about the data
    call ReadGradsCtlFiles(GctlFiles)

    ! Initialize the GRADS vars - this loads up the GradsVar structure except for the actual
    ! data. Reading in the data is deferred until later after we know everything is okay with
    ! the GRADS vars since it takes a while to read in the actual data.
    !
    ! Always use w (filtering) and press (find storm center).
    ! If doing horizontal velocity, ie VarToAvg is 'speed_t' or 'speed_r', set up U and V
    ! otherwise set up the variable named in VarToAvg
    call InitGvarFromGdescrip(GctlFiles, W, 'w')
    call InitGvarFromGdescrip(GctlFiles, Press, 'press')
    if (DoHorizVel) then
      call InitGvarFromGdescrip(GctlFiles, U, 'u')
      call InitGvarFromGdescrip(GctlFiles, V, 'v')
      ! In this case Avar gets built from U and V instead of getting
      ! read in directly from a GRADS data file. In this case initialize
      ! Avar from data in U.
      call InitGvarFromGvar(U, Avar, VarToAvg)
    else
      call InitGvarFromGdescrip(GctlFiles, Avar, VarToAvg)
    endif

    ! check that the variable dimensions (size and coordinate values) match up, if this
    ! isn't true, then the subsequent anlysis will be bogus
    !
    ! third arg of GvarDimsMatch() is true if one of the vars is 2D, else false
    if (DoHorizVel) then
      if (.not. (GvarDimsMatch(W, U, .false.) .and. GvarDimsMatch(W, V, .false.) .and. &
                 GvarDimsMatch(W, Press, .false.))) then
        write (*,*) 'ERROR: dimensions of u, v, w, and press do not match'
        stop
      endif
    else
      if (.not. (GvarDimsMatch(W, Press, .false.) .and. GvarDimsMatch(W, Avar, VarIs2d))) then
        write (*,*) 'ERROR: dimensions of w, press, and ', trim(VarToAvg), ' do not match'
        stop
      endif
    endif
 
    ! Always read in W so use it to record Nx, Ny, Nz, Nt. At this point we have verified
    ! that Nx, Ny, Nz, Nt are consitent for all the variables.
    Nx = W%Nx
    Ny = W%Ny
    Nz = W%Nz
    Nt = W%Nt
    VarNz = Avar%Nz

    write (*,*) 'Gridded data information:'
    write (*,*) '  Number of x (longitude) points:            ', Nx
    write (*,*) '  Number of y (latitude) points:             ', Ny
    write (*,*) '  Number of z (vertical level) points (3D):  ', Nz
    write (*,*) '  Number of z (vertical level) points (var): ', VarNz
    write (*,*) '  Number of t (time) points:                 ', Nt
    write (*,*) ''
    write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
    write (*,*) ''
  
    write (*,*) 'Locations of variables in GRADS data (file, var number):'
    write (*,'(a20,a,a2,i3,a1)') 'w: (', trim(W%DataFile), ', ', W%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    if (DoHorizVel) then
      write (*,'(a20,a,a2,i3,a1)') 'speed - u: (', trim(U%DataFile), ', ', U%Vnum, ')'
      write (*,'(a20,a,a2,i3,a1)') 'speed - v: (', trim(V%DataFile), ', ', V%Vnum, ')'
    else
      write (*,'(a17,a3,a,a2,i3,a1)') trim(VarToAvg), ': (', trim(Avar%DataFile), ', ', Avar%Vnum, ')'
    end if
    write (*,*) ''

    ! Convert the GRADS grid coordinates from longitude, latitude to flat plane (x and y).
    ! Always have W which has been verified to be consistent with all other vars
    ! so use W for this call
    ! Xcoords, Ycoords are in units of km
    call ConvertGridCoords(W, Xcoords, Ycoords)
  
    write (*,*) 'Horzontal Grid Coordinate Info:'
    write (*,*) '  X Range (min lon, max lon) --> (min x, max x): '
    write (*,*) '    ', W%Xcoords(1), W%Xcoords(W%Nx), Xcoords(1), Xcoords(Nx)
    write (*,*) '  Y Range: '
    write (*,*) '    ', W%Ycoords(1), W%Ycoords(W%Ny), Ycoords(1), Ycoords(Ny)
    write (*,*) ''
  
    ! Figure out how big to make each radial band. Assume that the storm center stays near
    ! the center of the grid domain. Take the diagonal distance of the domain, cut it in
    ! half and then break this up into NumRbands sections of equal length.
  
    DeltaX = Xcoords(Nx) - Xcoords(1)
    DeltaY = Ycoords(Ny) - Ycoords(1)
    RadialDist = sqrt(DeltaX*DeltaX + DeltaY*DeltaY) / 2.0
    RbandInc = RadialDist / real(NumRbands)
  
    write (*,*) 'Radial band information:'
    write (*,*) '  Delta x of domain:     ', DeltaX
    write (*,*) '  Delta y of domain:     ', DeltaY
    write (*,*) '  Radial distance:       ', RadialDist
    write (*,*) '  Radial band increment: ', RbandInc
    write (*,*) ''
  
    ! Read in the data for the vars using the description and location information
    call ReadGradsData(W)
    call ReadGradsData(Press)
    call RecordStormCenter(Press, StmIx, StmIy, MinP)
    if (DoHorizVel) then
      call ReadGradsData(U)
      call ReadGradsData(V)
      call ConvertHorizVelocity(U, V, StmIx, StmIy, Xcoords, Ycoords, Avar, DoTangential)
    else
      call ReadGradsData(Avar)
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
    VarNz = 1
    Nt = 5

    NumRbands = 10
    TestGridX = 100.0
    TestGridY = 100.0

    RadialDist = sqrt(TestGridX*TestGridX + TestGridY*TestGridY) / 2.0
    RbandInc = RadialDist / real(NumRbands)

    ! load up the coordinates
    allocate (Xcoords(1:Nx), Ycoords(1:Ny))
    DeltaX = TestGridX / real(Nx - 1)
    DeltaY = TestGridY / real(Ny - 1)
    do ix = 1, Nx
      Xcoords(ix) = real(ix-1) * DeltaX
    end do
    do iy = 1, Ny
      Ycoords(iy) = real(iy-1) * DeltaY
    end do
    do iz = 1, Nz
      TestZcoords(iz) = real(iz)
    end do

    call InitGradsVar(W, 'w', Nx, Ny, Nz, Nt, 0.0, DeltaX, 0.0, DeltaY, TestZcoords, &
                      '00:00Z24aug1998', '01mn', 10.0e30, '<NONE>', 0, 0)
    call InitGradsVar(Avar, VarToAvg, Nx, Ny, Nz, Nt, 0.0, DeltaX, 0.0, DeltaY, TestZcoords, &
                      '00:00Z24aug1998', '01mn', 10.0e30, '<NONE>', 0, 0)

    allocate (StmIx(1:Nt), StmIy(1:Nt), MinP(1:Nt))

    do it = 1, Nt
      do iz = 1, Nz
        do iy = 1, Ny
          do ix = 1, Nx

            ! Storm center drifts from roughly the center of the grid
            StmIx(it) = int(Nx/2) + it
            StmIy(it) = int(Ny/2) - it
            MinP(it) = 980.0 - real(it)  ! Just make it up, it's only used in diagnostic msg
                                         ! in AzimuthalAverage().

            ! Fill up W with the it value so you can see if screening works with different
            ! WfilterMin, WfilterMax values
            W%Vdata(iz,it,ix,iy) = real(it)

            ! Make the data match the radial band number after the averaging takes place.
            !   int(sqrt(DeltaX*DeltaX + DeltaY*DeltaY) / RbandInc) + 1) gives you the radial band number
            !   just multiply it by the variable "it" so you see increasing slope for
            ! successive time steps
            DeltaX = Xcoords(StmIx(it)) - Xcoords(ix)
            DeltaY = Ycoords(StmIy(it)) - Ycoords(iy)
            Avar%Vdata(it,iz,ix,iy) = real((int(sqrt(DeltaX*DeltaX + DeltaY*DeltaY) / RbandInc) + 1) * it)

          end do 
        end do 
      end do 
    end do 

  end if
  
  ! Initialize the output GRADS var: AzAvg
  call InitGradsVar(AzAvg, VarToAvg, NumRbands, 1, Avar%Nz, Avar%Nt, &
                    0.0, RbandInc, 0.0, 1.0, Avar%Zcoords, Avar%Tstart, Avar%Tinc, &
                    Avar%UndefVal, '<NONE>', 0, 0)

  ! Do the averageing and write the result into GRADS files
  call AzimuthalAverage(NumRbands, W, StmIx, StmIy, MinP, Avar, AzAvg, &
          Xcoords, Ycoords, RadialDist, RbandInc, WfilterMin, WfilterMax, Avar%UndefVal)

  call WriteGrads(AzAvg, OfileBase, 'azavg')

  stop
end

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   OfileBase - output GRADS file, base name for two files
!   NumRbands - number of radial bands to split data into
!   VarToAvg - RAMS variable to do the averaging on
!

subroutine GetMyArgs(Infiles, OfileBase, NumRbands, WfilterMin, WfilterMax, VarToAvg, DoHorizVel, DoTangential, VarIs2d)
  implicit none

  integer :: NumRbands
  real :: WfilterMin, WfilterMax
  character (len=*) :: Infiles, OfileBase, VarToAvg
  logical :: DoHorizVel, DoTangential, VarIs2d

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 7) then
    write (*,*) 'ERROR: must supply exactly 7 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_file> <num_radial_bands> <w_fliter_min> <w_filter_max> <var_to_average> <dim_of_var>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <num_radial_bands>: number of bands to split data into'
    write (*,*) '        <w_filter_min>: min (left end) of interval used to filter out data'
    write (*,*) '        <w_filter_max>: max (right end) of interval used to filter out data'
    write (*,*) '        <var_to_average>: name of RAMS variable to do the averaging on'
    write (*,*) '        <dim_of_var>: indicates if <var_to_average> is 2d or 3d'
    write (*,*) '           <dim_of_var> must be either "2d" or "3d"'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)

  call getarg(3, arg)       !this is a string
  read (arg, '(i)') NumRbands !convert to integer

  call getarg(4, arg)
  read (arg, '(f)') WfilterMin

  call getarg(5, arg)
  read (arg, '(f)') WfilterMax
  if (WfilterMin .ge. WfilterMax) then
    write (*,*) 'ERROR: must specify <w_filter_min> to be less than <w_filter_max>'
    stop
  endif

  call getarg(6, VarToAvg)
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

  call getarg(7, arg)
  if (arg .eq. '2d') then
    VarIs2d = .true.
  else if (arg .eq. '3d') then
    VarIs2d = .false.
  else
    write (*,*) 'ERROR: must use either "2d" or "3d" for <dim_of_var> argument'
    stop
  end if

  return
end subroutine
