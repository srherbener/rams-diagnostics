!***************************************************************
! Program to create CFAD data (Contoured Frequency by Altitude Diagrams)
!
! This program will read in GRADS data from a RAMS simulation, find
! the storm center and create normalized histogram type data which
! can be plotted into a CFAD. The histogram data is constructed for
! radial bands around the storm center.
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
  real :: Wthreshold
  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarName
  logical :: DoHorizVel, DoTangential

  type (GradsControlFiles) :: GctlFiles

  integer :: Nx, Ny, Nz, Nt

  ! Data arrays: need one for w (vertical velocity), press (pressure)
  ! and the var we are doing the averaging on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  type (GradsVar) :: U, V, Press, Avar, Cfad
  integer, dimension(:), allocatable :: StmIx, StmIy
  real, dimension(:), allocatable :: MinP, Xcoords, Ycoords
  integer :: NumBins, ib, iBinCenter
  real :: BinInc, BinStart, BinFilterMin, BinFilterMax
  real, dimension(:), allocatable :: BinVals

  integer :: i

  integer :: ix, iy, iz, it

  real :: DeltaX, DeltaY, RadialDist, RbandInc
  real :: TestData, TestGridX, TestGridY
  integer :: iTestRband, iTestCount
  real, dimension(1:MaxCoords) :: TestZcoords

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBase, NumRbands, BinFilterMin, BinFilterMax, VarName, &
                 DoHorizVel, DoTangential)
  call String2List(Infiles, ':', GctlFiles%Fnames, MaxFiles, GctlFiles%Nfiles, 'input files')

  write (*,*) 'Calculating azimuthal average for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, GctlFiles%Nfiles
    write (*,*) '  ', i, ': ', trim(GctlFiles%Fnames(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  Number of radial bands: ', NumRbands
  if (DoHorizVel) then
    if (DoTangential) then
      write (*,*) '  RAMS variable that is being analyzed: Tangential Horizontal Velocity'
    else
      write (*,*) '  RAMS variable that is being analyzed: Radial Horizontal Velocity'
    end if
  else
    write (*,*) '  RAMS variable that is being analyzed: ', trim(VarName)
  end if
  write (*,*) ''

  if (trim(VarName) /= 'test') then
    ! Not running a test so grab the data from GRADS input files and perform the averaging.

    ! Read the GRADS data description files and collect the information about the data
    call ReadGradsCtlFiles(GctlFiles)

    ! Initialize the GRADS vars
    call InitGvarFromGdescrip(GctlFiles, Press, 'press')
    if (DoHorizVel) then
      call InitGvarFromGdescrip(GctlFiles, U, 'u')
      call InitGvarFromGdescrip(GctlFiles, V, 'v')
      ! In this case, Avar gets built from U and V instead of getting
      ! read in directly from a GRADS data file. In this case initialize
      ! Avar from data in U
      call InitGvarFromGvar(U, Avar, VarName)
    else
      call InitGvarFromGdescrip(GctlFiles, Avar, VarName)
    endif
  
    ! check that the variable dimensions (size and coordinate values) match up, if this
    ! isn't true, then the subsequent anlysis will be bogus
    !
    ! third arg of GvarDimsMatch() is true if one of the vars is 2D, else false
    if (DoHorizVel) then
      if (.not. (GvarDimsMatch(Press, U, .false.) .and. GvarDimsMatch(Press, V, .false.))) then
        write (*,*) 'ERROR: dimensions of u, v, and press do not match'
        stop
      endif
    else
      if (.not. (GvarDimsMatch(Press, Avar, .false.))) then
        write (*,*) 'ERROR: dimensions of press, and ', trim(VarName), ' do not match'
        stop
      endif
    endif
 
    ! Always read in Press so use it to record Nx, Ny, Nz, Nt. At this point we have verified
    ! that Nx, Ny, Nz, Nt are consitent for all the variables.
    Nx = Press%Nx
    Ny = Press%Ny
    Nz = Press%Nz
    Nt = Press%Nt
  
    write (*,*) 'Gridded data information:'
    write (*,*) '  Number of x (longitude) points:          ', Nx
    write (*,*) '  Number of y (latitude) points:           ', Ny
    write (*,*) '  Number of z (vertical level) points:     ', Nz
    write (*,*) '  Number of t (time) points:               ', Nt
    write (*,*) ''
    write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
    write (*,*) ''
  
    write (*,*) 'Locations of variables in GRADS data (file, var number):'
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    if (DoHorizVel) then
      write (*,'(a20,a,a2,i3,a1)') 'speed - u: (', trim(U%DataFile), ', ', U%Vnum, ')'
      write (*,'(a20,a,a2,i3,a1)') 'speed - v: (', trim(V%DataFile), ', ', V%Vnum, ')'
    else
      write (*,'(a17,a3,a,a2,i3,a1)') trim(VarName), ': (', trim(Avar%DataFile), ', ', Avar%Vnum, ')'
    end if
    write (*,*) ''
  
    ! Convert the GRADS grid coordinates from longitude, latitude to flat plane (x and y).
    ! Always have Press which has been verified to be consistent with all other vars
    ! so use Press for this call
    ! Xcoords, Ycoords are in units of km
    call ConvertGridCoords(Press, Xcoords, Ycoords)
  
    write (*,*) 'Horzontal Grid Coordinate Info:'
    write (*,*) '  X Range (min lon, max lon) --> (min x, max x): '
    write (*,*) '    ', Press%Xcoords(1), Press%Xcoords(Nx), Xcoords(1), Xcoords(Nx)
    write (*,*) '  Y Range: '
    write (*,*) '    ', Press%Ycoords(1), Press%Ycoords(Ny), Ycoords(1), Ycoords(Ny)
    write (*,*) ''
  
    ! Figure out how big to make each radial band. Assume that the storm center stays near the center of
    ! the grid domain. Take the diagonal distance of the domain, cut it in half and then break this up
    ! into NumRbands sections of equal length.
  
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

    ! Set up the CFAD bins
    ! only know how to do 'w' for now -> set up 11 bins: -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5
    ! data goes in a bin if it's abs value is >= to abs value of bin, neg numbers in the
    ! negative bins, positive numbers in the positive bins, 0 bin is for numbers between -1, and 1.

    if (VarName .eq. 'w') then
      NumBins = 11
      BinStart = -5.0
      BinInc = 1.0
      iBinCenter = 6 ! index of middle entry in array of bin values 
    else
      write (*,*) 'Sorry: only know how to do "w" for now'
      stop
    endif
  
    allocate (BinVals(1:NumBins))

    BinVals(1) = BinStart
    do ib = 2, NumBins
      BinVals(ib) = BinVals(ib-1) + BinInc
    end do

    write (*,*) 'Bin information: '
    write (*,*) '  Number of bins: ', NumBins
    write (*,*) '  Bin start value: ', BinStart
    write (*,*) '  Bin increment value: ', BinInc
    write (*,*) '  Bin filter min: ', BinFilterMin
    write (*,*) '  Bin filter max: ', BinFilterMax
    do ib = 1, NumBins
      write (*,*) '    ', ib, ' --> ', BinVals(ib)
    end do
    write (*,*) ''

    ! Read in the data for the vars using the description and location information
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

    ! Create one horizontal slice, 10 x 10 tiles, 3 z levels, 5 time points.
    !   Make the x,y increments 1km
    ! Move the storm center around a little bit but keep it near the center.
    ! Make 10 radial bands, 

    Nx = 5
    Ny = 5
    Nz = 3
    Nt = 5

    NumRbands = 4
    TestGridX = 100.0
    TestGridY = 100.0

    NumBins = 5
    BinStart = 1.0
    BinInc = 1.0
    iBinCenter = 2

    allocate (BinVals(1:NumBins))

    BinVals(1) = BinStart
    do ib = 2, NumBins
      BinVals(ib) = BinVals(ib-1) + BinInc
    end do

    BinFilterMin = -10.0 ! pick range where all data will fall outside
    BinFilterMax = -10.0 ! of - ie, shut off the filtering for the test case

    RadialDist = sqrt(TestGridX*TestGridX + TestGridY*TestGridY) / 2.0
    RbandInc = RadialDist / real(NumRbands)

    ! load up the coordinates
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

    call InitGradsVar(Avar, VarName, Nx, Ny, Nz, Nt, 0.0, DeltaX, 0.0, DeltaY, TestZcoords, &
                      '00:00Z24aug1998', '01mn', 10.0e30, '<NONE>', 0, 0)

    allocate(StmIx(1:Nt), StmIy(1:Nt), MinP(1:Nt))

    iTestCount = 1
    do it = 1, Nt
      do iz = 1, Nz
        do iy = 1, Ny
          do ix = 1, Nx

            ! Storm center drifts from roughly the center of the grid
            StmIx(it) = int(Nx/2) + it
            StmIy(it) = int(Ny/2) - it
            MinP(it) = 980.0 - real(it)  ! Just make it up, it's only used in diagnostic msg
                                         ! in AzimuthalAverage().

            ! Make the data match the radial band number after the averaging takes place.
            !   int(sqrt(DeltaX*DeltaX + DeltaY*DeltaY) / RbandInc) + 1) gives you the radial band number
            !   just multiply it by it so you see increasing slope for successive time steps
            DeltaX = Xcoords(StmIx(it)) - Xcoords(ix)
            DeltaY = Ycoords(StmIy(it)) - Ycoords(iy)
            iTestRband = int(sqrt(DeltaX*DeltaX + DeltaY*DeltaY) / RbandInc) + 1
            !TestData = real(iTestRband)
            TestData = real(mod(iTestCount,5))
            Avar%Vdata(it,iz,ix,iy) = TestData
            iTestCount = iTestCount + 1
          end do 
        end do 
      end do 
    end do 


  end if
  
  ! Allocate the output array and compute the azimuthal averaging
  call InitGradsVar(Cfad, VarName, NumRbands, NumBins, Avar%Nz, Avar%Nt, &
                    0.0, RbandInc, BinStart, BinInc, Avar%Zcoords, Avar%Tstart, Avar%Tinc, &
                    Avar%UndefVal, '<NONE>', 0, 0)

  ! Do the binning.
  ! The output array, Cfad has the dimensions:
  !    t -> time
  !    z -> levels
  !    x -> radial bands
  !    y -> histogram bins

  call BuildCfad(NumRbands, NumBins, iBinCenter, BinFilterMin, BinFilterMax, &
          StmIx, StmIy, MinP, Avar, Cfad, Xcoords, Ycoords, RadialDist, RbandInc, BinVals)

  call WriteGrads(Cfad, OfileBase, 'cfad')

  stop
end

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   OfileBase - output GRADS file, base name for two files
!   NumRbands - number of radial bands to split data into
!   VarName - RAMS variable to do the averaging on
!

subroutine GetMyArgs(Infiles, OfileBase, NumRbands, BinFilterMin, BinFilterMax, &
                     VarName, DoHorizVel, DoTangential)
  implicit none

  integer :: NumRbands
  real :: BinFilterMin, BinFilterMax
  character (len=*) :: Infiles, OfileBase, VarName
  logical :: DoHorizVel, DoTangential

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 6) then
    write (*,*) 'ERROR: must supply exactly 6 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_file> <num_radial_bands> <bin_filter_min> <bin_filter_max> <var_to_average>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <num_radial_bands>: number of bands to split data into'
    write (*,*) '        <bin_filter_min>: min value (left end) of interval for filtering data'
    write (*,*) '        <bin_filter_max>: max value (right end) of interval for filtering data'
    write (*,*) '        <var_to_analyze>: name of RAMS variable to do the analysis on'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)

  call getarg(3, arg)       !this is a string
  read (arg, '(i)') NumRbands !convert to integer

  call getarg(4, arg)       !this is a string
  read (arg, '(f)') BinFilterMin !convert to real

  call getarg(5, arg)       !this is a string
  read (arg, '(f)') BinFilterMax !convert to real

  if (BinFilterMin .ge. BinFilterMax) then
    write (*,*) 'ERROR: must specify <bin_filter_min> to be less than <bin_filter_max>'
    stop
  endif

  call getarg(6, VarName)
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

  return
end subroutine

!*****************************************************************************
! BuildCfad()
!
! This routine will sort data points out into radial bands and then do 
! the histogram binning within each radial band.
!

subroutine BuildCfad(NumRbands, NumBins, iBinCenter, BinFilterMin, &
          BinFilterMax, StmIx, StmIy, MinP, Avar, Cfad, Xcoords, Ycoords, RadialDist, &
          RbandInc, BinVals)

  use gdata_utils
  implicit none

  integer :: NumRbands, NumBins, iBinCenter
  type (GradsVar) :: Avar, Cfad
  integer, dimension(1:Avar%Nt) :: StmIx, StmIy
  real, dimension(1:Avar%Nt) :: MinP
  real, dimension(1:Avar%Nx) :: Xcoords
  real, dimension(1:Avar%Ny) :: Ycoords
  real, dimension(1:NumBins) :: BinVals
  real :: BinFilterMin, BinFilterMax, RadialDist, RbandInc

  integer :: ix, iy, iz, it, ib, ir
  real :: DeltaX, DeltaY, Radius, BinTotal
  integer :: iRband, iBin

  ! initialize the counts
  do ir = 1, NumRbands
    do ib = 1, NumBins
      do iz = 1, Avar%Nz
        do it = 1, Avar%Nt
          Cfad%Vdata(it, iz, ir, ib) = 0.0;
        end do 
      end do 
    end do 
  end do 

  ! create the histogram data
  write (*,*) 'Binning Data:'

  do it = 1, Avar%Nt
    if (modulo(it,10) .eq. 0) then
      write (*,*) '  Timestep: ', it
      write (*,'(a,i3,a,i3,a,g,a,g,a)') '    Storm Center: (', StmIx(it), ', ', StmIy(it), ') --> (', &
            Xcoords(StmIx(it)), ', ', Ycoords(StmIy(it)), ')'
      write (*,*) '    Minimum Pressue: ', MinP(it)
    end if

    do iz = 1, Avar%Nz
      do iy = 1, Avar%Ny
        do ix = 1, Avar%Nx
          ! Throw out undefined data points, and those falling inbetween BinFilterMin
          ! and BinFilterMax (which is specific to processing w - intereseted in the
          ! greater abs values)
          if ((Avar%Vdata(it,iz,ix,iy) .ne. Avar%UndefVal) .and. &
              ((Avar%Vdata(it,iz,ix,iy) .lt. BinFilterMin) .or. &
               (Avar%Vdata(it,iz,ix,iy) .gt. BinFilterMax))) then
             ! find the radial band for this data point
             DeltaX = Xcoords(ix)-Xcoords(StmIx(it))
             DeltaY = Ycoords(iy)-Ycoords(StmIy(it))
             Radius = sqrt(DeltaX*DeltaX + DeltaY*DeltaY)
             iRband = int(Radius / RbandInc) + 1
             ! iRband will go from 1 to n, but n might extend beyond the last radial
             ! band due to approximations made in deriving RbandInc
             if (iRband .gt. NumRbands) then
               iRband = NumRbands
             end if

             ! Find the bin for this data point, iBinCenter has the index (for BinVals) that
             ! is the "center point" of this bin. This means that numbers between
             ! BinVals(iBinCenter-1) and BinVals(iBinCenter+1) go in this bin; indexes
             ! less than iBinCenter hold numbers between BinVals(iBinCenter) and BinVals(iBinCenter-1)
             ! and indecies greater than iBinCenter hold numbers between BinVals(iBinCenter) and
             ! BinVals(iBinCenter+1).
             ! This code assumes that the bin center is not on the ends of the BinVals array.
             iBin = 0
             do ib = 1, NumBins
               ! on the bin center
               if (ib .eq. iBinCenter) then
                 if ((Avar%Vdata(it,iz,ix,iy) .gt. BinVals(ib-1)) .and. (Avar%Vdata(it,iz,ix,iy) .lt. BinVals(ib+1))) then
                   iBin = ib
                 endif
               endif 

               ! to the left of the center bin
               if (ib .lt. iBinCenter) then
                 if (Avar%Vdata(it,iz,ix,iy) .le. BinVals(ib)) then
                   if (ib .eq. 1) then
                     iBin = ib
                   else
                     if (Avar%Vdata(it,iz,ix,iy) .gt. BinVals(ib-1)) then
                       iBin = ib
                     endif
                   endif
                 endif
               endif 

               ! to the right of the center bin
               if (ib .gt. iBinCenter) then
                 if (Avar%Vdata(it,iz,ix,iy) .ge. BinVals(ib)) then
                   if (ib .eq. NumBins) then
                     iBin = ib
                   else
                     if (Avar%Vdata(it,iz,ix,iy) .lt. BinVals(ib+1)) then
                       iBin = ib
                     endif
                   endif
                 endif
               endif 
             end do

             if (iBin .eq. 0) then
               write (*,*) 'ERROR: did not find a bin for data at location:'
               write (*,*) 'ERROR:   it, iz, ix, iy: ', it, iz, ix, iy
               stop
             endif
             
             Cfad%Vdata(it, iz, iRband, iBin) = Cfad%Vdata(it, iz, iRband, iBin) + 1.0
          endif
        end do
      end do
    end do
  end do

  ! Have counts in Cfad now. Convert these to % numbers.

  do it = 1, Avar%Nt
    do iz = 1, Avar%Nz
      do ir = 1, NumRbands
        ! sum up counts in all the bins
        BinTotal = 0.0
        do ib = 1, NumBins
          BinTotal = BinTotal + Cfad%Vdata(it, iz, ir, ib)
        end do
        ! convert entries to per cent values
        do ib = 1, NumBins
          ! If we selected out all data points, then the bin total will
          ! be zero. In this case, leave the zeros in the Cfad array.
          if (BinTotal .ne. 0.0) then
            Cfad%Vdata(it, iz, ir, ib) = (Cfad%Vdata(it, iz, ir, ib) / BinTotal) * 100.0
          endif
        end do
      end do
    end do
  end do

  write (*,*) ''
  return
end subroutine
