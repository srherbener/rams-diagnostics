!***************************************************************
! Program to do averaging over the domain and yield a single
! data point at each time step
!
! This program will read in GRADS data from a RAMS simulation, 
! perform an averaging function, and output the single point time
! series in GRADS format
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
!   3. selection of averaging function
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

  character (len=LargeString) :: Infiles, Outfiles
  character (len=LittleString) :: AvgFunc

  logical :: DoRates, DoSc

  character (len=MediumString), dimension(1:MaxFiles) :: OfileBases
  type (GradsControlFiles) :: GctlFiles
  integer :: Nofiles
  integer, dimension(:), allocatable :: StmIx, StmIy

  real :: DeltaX, DeltaY, DeltaT, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
  real :: Xstart, Xinc, Ystart, Yinc
  real :: ConvertTinc
  real, dimension(:), allocatable :: ZmHeights
  real, dimension(1:1) :: DummyZcoords
  real, dimension(:), allocatable :: MinP, Xcoords, Ycoords, Zcoords

  ! Data arrays: cloud, tempc, precipr
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of cloud, tempc, precipr in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  type (GradsVar) :: Cloud, TempC, Dens, PrecipR, U, V, W, CcnConc, CloudDiam, CloudConc, Press, CintLiq, AzWind, TsAvg, TestSelect, Rates 

  integer :: i

  integer :: ix, iy, iz, it

  ! Get the command line arguments
  call GetMyArgs(Infiles, Outfiles, AvgFunc, DoRates, Wthreshold, CilThresh, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ)
  call String2List(Infiles, ':', GctlFiles%Fnames, MaxFiles, GctlFiles%Nfiles, 'input files')
  call String2List(Outfiles, ':', OfileBases, MaxFiles, Nofiles, 'output files')

  DoSc = .true.
  if ((AvgFunc .eq. 'wr_cloud') .or. &
      (AvgFunc .eq. 'wr_cloud_diam') .or. &
      (AvgFunc .eq. 'wr_cloud_conc') .or. &
      (AvgFunc .eq. 'wr_cloud2') .or. &
      (AvgFunc .eq. 'wr_cloud2_diam') .or. &
      (AvgFunc .eq. 'wr_cloud2_conc')) then
    DoSc = .false.
  end if

  write (*,*) 'Time seris of average for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, GctlFiles%Nfiles
    write (*,*) '  ', i, ': ', trim(GctlFiles%Fnames(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBases(1))
  write (*,*) '  Output rates file base name:  ', trim(OfileBases(2))
  write (*,*) '  Averaging function: ', trim(AvgFunc)
  write (*,*) '  Do rates?: ', DoRates
  write (*,*) '  Do supercoooled?: ', DoSc
  write (*,*) '  Data selection specs: '
  write (*,*) '    Minimum radius: ', MinR
  write (*,*) '    Maximum radius: ', MaxR
  write (*,*) '    Minimum angle: ', MinPhi
  write (*,*) '    Maximum angle: ', MaxPhi
  write (*,*) '    Minimum height: ', MinZ
  write (*,*) '    Maximum height: ', MaxZ
  write (*,*) '    W threshold: ', Wthreshold
  write (*,*) '    Column integrated liquid threshold: ', CilThresh
  write (*,*) ''
  flush(6)

  ! Read the GRADS data description files and collect the information about the data
  call ReadGradsCtlFiles(GctlFiles)

  ! Check the data description for consistency and locate the variables in the GRADS control files
  if ((AvgFunc .eq. 'sc_cloud') .or. (AvgFunc .eq. 'wr_cloud')) then
    call InitGvarFromGdescrip(GctlFiles, Press,   'press')
    call InitGvarFromGdescrip(GctlFiles, TempC,   'tempc')
    call InitGvarFromGdescrip(GctlFiles, Dens,    'dn0')
    call InitGvarFromGdescrip(GctlFiles, Cloud,   'cloud')
    call InitGvarFromGdescrip(GctlFiles, CintLiq, 'cint_liq')

    if (.not. (GvarDimsMatch(Press, TempC, .false.) .and. GvarDimsMatch(Press, Dens, .false.) .and. &
               GvarDimsMatch(Press, Cloud, .false.) .and. GvarDimsMatch(Press, CintLiq, .true.))) then
      write (*,*) 'ERROR: dimensions of press, tempc, dn0, cloud and cint_liq do not match'
      stop
    endif
  else if ((AvgFunc .eq. 'sc_cloud2') .or. (AvgFunc .eq. 'wr_cloud2')) then
    call InitGvarFromGdescrip(GctlFiles, Press,   'press')
    call InitGvarFromGdescrip(GctlFiles, TempC,   'tempc')
    call InitGvarFromGdescrip(GctlFiles, Dens,    'dn0')
    call InitGvarFromGdescrip(GctlFiles, Cloud,   'cloud2')
    call InitGvarFromGdescrip(GctlFiles, CintLiq, 'cint_liq')

    if (.not. (GvarDimsMatch(Press, TempC, .false.) .and. GvarDimsMatch(Press, Dens, .false.) .and. &
               GvarDimsMatch(Press, Cloud, .false.) .and. GvarDimsMatch(Press, CintLiq, .true.))) then
      write (*,*) 'ERROR: dimensions of press, tempc, dn0, cloud2 and cint_liq do not match'
      stop
    endif
  else if (AvgFunc .eq. 'precipr') then
    call InitGvarFromGdescrip(GctlFiles, Press,   'press')
    call InitGvarFromGdescrip(GctlFiles, PrecipR, 'precipr')
    call InitGvarFromGdescrip(GctlFiles, CintLiq, 'cint_liq')

    if (.not. (GvarDimsMatch(Press, PrecipR, .true.) .and. GvarDimsMatch(Press, CintLiq, .true.))) then
      write (*,*) 'ERROR: dimensions of press, precipr, and cint_liq do not match'
      stop
    endif
  else if (AvgFunc .eq. 'w_up') then
    call InitGvarFromGdescrip(GctlFiles, Press,   'press')
    call InitGvarFromGdescrip(GctlFiles, W,       'w')

    if (.not. (GvarDimsMatch(Press, W, .false.))) then
      write (*,*) 'ERROR: dimensions of press, and w do not match'
      stop
    endif
  else if ((AvgFunc .eq. 'sc_cloud_diam') .or. (AvgFunc .eq. 'wr_cloud_diam')) then
    call InitGvarFromGdescrip(GctlFiles, Press,     'press')
    call InitGvarFromGdescrip(GctlFiles, TempC,     'tempc')
    call InitGvarFromGdescrip(GctlFiles, Cloud,     'cloud')
    call InitGvarFromGdescrip(GctlFiles, CloudDiam, 'cloud_d')
    call InitGvarFromGdescrip(GctlFiles, CintLiq,   'cint_liq')

    if (.not. (GvarDimsMatch(Press, TempC, .false.) .and. GvarDimsMatch(Press, Cloud, .false.) .and. &
               GvarDimsMatch(Press, CloudDiam, .false.) .and. GvarDimsMatch(Press, CintLiq, .true.))) then
      write (*,*) 'ERROR: dimensions of press, tempc, dn0, cloud_d and cint_liq do not match'
      stop
    endif
  else if ((AvgFunc .eq. 'sc_cloud2_diam') .or. (AvgFunc .eq. 'wr_cloud2_diam')) then
    call InitGvarFromGdescrip(GctlFiles, Press,     'press')
    call InitGvarFromGdescrip(GctlFiles, TempC,     'tempc')
    call InitGvarFromGdescrip(GctlFiles, Cloud,     'cloud2')
    call InitGvarFromGdescrip(GctlFiles, CloudDiam, 'cloud2_d')
    call InitGvarFromGdescrip(GctlFiles, CintLiq,   'cint_liq')

    if (.not. (GvarDimsMatch(Press, TempC, .false.) .and. GvarDimsMatch(Press, Cloud, .false.) .and. &
               GvarDimsMatch(Press, CloudDiam, .false.) .and. GvarDimsMatch(Press, CintLiq, .true.))) then
      write (*,*) 'ERROR: dimensions of press, tempc, dn0, cloud2_d and cint_liq do not match'
      stop
    endif
  else if ((AvgFunc .eq. 'sc_cloud_conc') .or. (AvgFunc .eq. 'wr_cloud_conc')) then
    call InitGvarFromGdescrip(GctlFiles, Press,     'press')
    call InitGvarFromGdescrip(GctlFiles, TempC,     'tempc')
    call InitGvarFromGdescrip(GctlFiles, CloudConc, 'cloud_cm3')

    if (.not. (GvarDimsMatch(Press, TempC, .false.) .and. GvarDimsMatch(Press, CloudConc, .false.))) then
      write (*,*) 'ERROR: dimensions of press, tempc, cloud_cm3 do not match'
      stop
    endif
  else if ((AvgFunc .eq. 'sc_cloud2_conc') .or. (AvgFunc .eq. 'wr_cloud2_conc')) then
    call InitGvarFromGdescrip(GctlFiles, Press,     'press')
    call InitGvarFromGdescrip(GctlFiles, TempC,     'tempc')
    call InitGvarFromGdescrip(GctlFiles, CloudConc, 'cloud2_cm3')

    if (.not. (GvarDimsMatch(Press, TempC, .false.) .and. GvarDimsMatch(Press, CloudConc, .false.))) then
      write (*,*) 'ERROR: dimensions of press, tempc, cloud2_cm3 do not match'
      stop
    endif
  else if (AvgFunc .eq. 'horiz_ke') then
    call InitGvarFromGdescrip(GctlFiles, Press, 'press')
    call InitGvarFromGdescrip(GctlFiles, Dens,  'dn0')
    call InitGvarFromGdescrip(GctlFiles, U,     'u')
    call InitGvarFromGdescrip(GctlFiles, V,     'v')

    if (.not. (GvarDimsMatch(Press, Dens, .false.) .and. &
               GvarDimsMatch(Press, U, .false.) .and. GvarDimsMatch(Press, V, .false.))) then
      write (*,*) 'ERROR: dimensions of press, dn0, u and v do not match'
      stop
    endif
  else if (AvgFunc .eq. 'storm_int') then
    call InitGvarFromGdescrip(GctlFiles, Press, 'press')
    call InitGvarFromGdescrip(GctlFiles, U,     'u')
    call InitGvarFromGdescrip(GctlFiles, V,     'v')

    if (.not. (GvarDimsMatch(Press, U, .false.) .and. GvarDimsMatch(Press, V, .false.))) then
      write (*,*) 'ERROR: dimensions of press, u and v do not match'
      stop
    endif
  else if (AvgFunc .eq. 'max_azwind') then
    call InitGvarFromGdescrip(GctlFiles, AzWind, 'speed_t')
  else if (AvgFunc .eq. 'ccnconc') then
    call InitGvarFromGdescrip(GctlFiles, Press,   'press')
    call InitGvarFromGdescrip(GctlFiles, CcnConc, 'ccn_conc')

    if (.not. (GvarDimsMatch(Press, CcnConc, .false.))) then
      write (*,*) 'ERROR: dimensions of press and ccn_conc do not match'
      stop
    endif
  else if (AvgFunc .eq. 'test_cvs') then
    call InitGvarFromGdescrip(GctlFiles, Press,   'press')
  end if

  ! Calculate the x,y coordinates (in km) for doing selection by radius from storm center.
  ! Always have Press, so use it for the coordinate calculations
  if (AvgFunc .eq. 'max_azwind') then
    ! Convert string time increment to seconds
    DeltaT = ConvertTinc(AzWind%Tinc)
  else
    call ConvertGridCoords(Press, Xcoords, Ycoords)
    allocate (Zcoords(1:Press%Nz))
    do iz = 1, Press%Nz
        Zcoords(iz) = Press%Zcoords(iz)
    enddo
  
    DeltaX = (Xcoords(2) - Xcoords(1)) * 1000.0
    DeltaY = (Ycoords(2) - Ycoords(1)) * 1000.0
  
    ! Convert string time increment to seconds
    DeltaT = ConvertTinc(Press%Tinc)
  
    write (*,*) 'Horizontal grid info:'
    write (*,*) '  X range (min lon, max lon) --> (min x, max x): '
    write (*,*) '    ', Press%Xcoords(1), Press%Xcoords(Press%Nx), Xcoords(1), Xcoords(Press%Nx)
    write (*,*) '  Y range (min lat, max lat) --> (min y, max y): '
    write (*,*) '    ', Press%Ycoords(1), Press%Ycoords(Press%Ny), Ycoords(1), Ycoords(Press%Ny)
    write (*,*) ''
    write (*,*) 'Vertical grid info:'
    do iz = 1, Press%Nz
      write (*,*) '  ', iz, ' --> ', Zcoords(iz)
    end do
    write (*,*) ''
    flush(6)
  endif

  if (AvgFunc .eq. 'max_azwind') then
    write (*,*) 'Gridded data information:'
    write (*,*) '  Number of x (longitude) points:          ', AzWind%Nx
    write (*,*) '  Number of y (latitude) points:           ', AzWind%Ny
    write (*,*) '  Number of z (vertical level) points:     ', AzWind%Nz
    write (*,*) '  Number of t (time) points:               ', AzWind%Nt
    write (*,*) ''
    write (*,*) '  Number of data values per grid variable: ', AzWind%Nx*AzWind%Ny*AzWind%Nz*AzWind%Nt
    write (*,*) ''
    write (*,*) '  Time increment: ', trim(AzWind%Tinc), ' --> ', DeltaT
    write (*,*) ''
    flush(6)
  else
    ! Read in the GRADS variable data
    write (*,*) 'Gridded data information:'
    write (*,*) '  Number of x (longitude) points:          ', Press%Nx
    write (*,*) '  Number of y (latitude) points:           ', Press%Ny
    write (*,*) '  Number of z (vertical level) points:     ', Press%Nz
    write (*,*) '  Number of t (time) points:               ', Press%Nt
    write (*,*) ''
    write (*,*) '  Number of data values per grid variable: ', Press%Nx*Press%Ny*Press%Nz*Press%Nt
    write (*,*) ''
    write (*,*) '  Grid delta x: ', DeltaX
    write (*,*) '  Grid delta y: ', DeltaY
    write (*,*) '  Time increment: ', trim(Press%Tinc), ' --> ', DeltaT
    write (*,*) ''
    flush(6)
  endif

  write (*,*) 'Locations of variables in GRADS data (file, var number):'
  if ((AvgFunc .eq. 'sc_cloud') .or. (AvgFunc .eq. 'wr_cloud')) then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'tempc: (', trim(TempC%DataFile), ', ', TempC%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'dn0: (', trim(Dens%DataFile), ', ', Dens%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cloud: (', trim(Cloud%DataFile), ', ', Cloud%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cint_liq: (', trim(CintLiq%DataFile), ', ', CintLiq%Vnum, ')'

    call ReadGradsData(Press)
    call ReadGradsData(TempC)
    call ReadGradsData(Dens)
    call ReadGradsData(Cloud)
    call ReadGradsData(CintLiq)
  else if ((AvgFunc .eq. 'sc_cloud2') .or. (AvgFunc .eq. 'wr_cloud2')) then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'tempc: (', trim(TempC%DataFile), ', ', TempC%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'dn0: (', trim(Dens%DataFile), ', ', Dens%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cloud2: (', trim(Cloud%DataFile), ', ', Cloud%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cint_liq: (', trim(CintLiq%DataFile), ', ', CintLiq%Vnum, ')'

    call ReadGradsData(Press)
    call ReadGradsData(TempC)
    call ReadGradsData(Dens)
    call ReadGradsData(Cloud)
    call ReadGradsData(CintLiq)
  else if (AvgFunc .eq. 'precipr') then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'precipr: (', trim(PrecipR%DataFile), ', ', PrecipR%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cint_liq: (', trim(CintLiq%DataFile), ', ', CintLiq%Vnum, ')'

    call ReadGradsData(Press)
    call ReadGradsData(PrecipR)
    call ReadGradsData(CintLiq)
  else if (AvgFunc .eq. 'w_up') then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'w: (', trim(W%DataFile), ', ', W%Vnum, ')'

    call ReadGradsData(Press)
    call ReadGradsData(W)
  else if ((AvgFunc .eq. 'sc_cloud_diam') .or. (AvgFunc .eq. 'wr_cloud_diam')) then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'tempc: (', trim(TempC%DataFile), ', ', TempC%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cloud: (', trim(Cloud%DataFile), ', ', Cloud%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cloud_d: (', trim(CloudDiam%DataFile), ', ', CloudDiam%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cint_liq: (', trim(CintLiq%DataFile), ', ', CintLiq%Vnum, ')'

    call ReadGradsData(Press)
    call ReadGradsData(TempC)
    call ReadGradsData(Cloud)
    call ReadGradsData(CloudDiam)
    call ReadGradsData(CintLiq)
  else if ((AvgFunc .eq. 'sc_cloud2_diam') .or. (AvgFunc .eq. 'wr_cloud2_diam')) then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'tempc: (', trim(TempC%DataFile), ', ', TempC%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cloud2: (', trim(Cloud%DataFile), ', ', Cloud%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cloud2_d: (', trim(CloudDiam%DataFile), ', ', CloudDiam%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cint_liq: (', trim(CintLiq%DataFile), ', ', CintLiq%Vnum, ')'

    call ReadGradsData(Press)
    call ReadGradsData(TempC)
    call ReadGradsData(Cloud)
    call ReadGradsData(CloudDiam)
    call ReadGradsData(CintLiq)
  else if ((AvgFunc .eq. 'sc_cloud_conc') .or. (AvgFunc .eq. 'wr_cloud_conc')) then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'tempc: (', trim(TempC%DataFile), ', ', TempC%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cloud_cm3: (', trim(CloudConc%DataFile), ', ', CloudConc%Vnum, ')'

    call ReadGradsData(Press)
    call ReadGradsData(TempC)
    call ReadGradsData(CloudConc)
  else if ((AvgFunc .eq. 'sc_cloud2_conc') .or. (AvgFunc .eq. 'wr_cloud2_conc')) then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'tempc: (', trim(TempC%DataFile), ', ', TempC%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'cloud2_cm3: (', trim(CloudConc%DataFile), ', ', CloudConc%Vnum, ')'

    call ReadGradsData(Press)
    call ReadGradsData(TempC)
    call ReadGradsData(CloudConc)
  else if (AvgFunc .eq. 'horiz_ke') then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'dn0: (', trim(Dens%DataFile), ', ', Dens%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'u: (', trim(U%DataFile), ', ', U%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'v: (', trim(V%DataFile), ', ', V%Vnum, ')'

    call ReadGradsData(Press)
    call ReadGradsData(Dens)
    call ReadGradsData(U)
    call ReadGradsData(V)
  else if (AvgFunc .eq. 'storm_int') then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'u: (', trim(U%DataFile), ', ', U%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'v: (', trim(V%DataFile), ', ', V%Vnum, ')'

    call ReadGradsData(U)
    call ReadGradsData(V)
  else if (AvgFunc .eq. 'max_azwind') then
    write (*,'(a20,a,a2,i3,a1)') 'speed_t: (', trim(AzWind%DataFile), ', ', AzWind%Vnum, ')'

    call ReadGradsData(AzWind)
  else if (AvgFunc .eq. 'ccnconc') then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
    write (*,'(a20,a,a2,i3,a1)') 'ccn_conc: (', trim(CcnConc%DataFile), ', ', CcnConc%Vnum, ')'

    call ReadGradsData(Press)
    call ReadGradsData(CcnConc)
  else if (AvgFunc .eq. 'test_cvs') then
    write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'

    call ReadGradsData(Press)
  end if

  write (*,*) ''
  flush(6)

  ! Allocate the output array and do the averaging
  if (AvgFunc .eq. 'test_cvs') then
    call InitGradsVar(TestSelect, AvgFunc, Press%Nx, Press%Ny, Press%Nz, Press%Nt, &
                      Press%Xstart, Press%Xinc, Press%Ystart, Press%Yinc, &
                      Press%Zcoords, Press%Tstart, Press%Tinc, &
                      Press%UndefVal, '<NONE>', 0, 0)
  else
    DummyZcoords(1) = 0.0
    if (AvgFunc .eq. 'max_azwind') then
      call InitGradsVar(TsAvg, AvgFunc, 1, 1, 1, AzWind%Nt, &
                      0.0, 1.0, 0.0, 1.0, DummyZcoords, AzWind%Tstart, AzWind%Tinc, &
                      AzWind%UndefVal, '<NONE>', 0, 0)
    else
      call InitGradsVar(TsAvg, AvgFunc, 1, 1, 1, Press%Nt, &
                      0.0, 1.0, 0.0, 1.0, DummyZcoords, Press%Tstart, Press%Tinc, &
                      Press%UndefVal, '<NONE>', 0, 0)
    endif
    if (DoRates .eq. .true.) then
      if (AvgFunc .eq. 'max_azwind') then
        call InitGradsVar(Rates, AvgFunc, 1, 1, 1, AzWind%Nt-2, &
                        0.0, 1.0, 0.0, 1.0, DummyZcoords, AzWind%Tstart, AzWind%Tinc, &
                        AzWind%UndefVal, '<NONE>', 0, 0)
      else
        call InitGradsVar(Rates, AvgFunc, 1, 1, 1, Press%Nt-2, &
                        0.0, 1.0, 0.0, 1.0, DummyZcoords, Press%Tstart, Press%Tinc, &
                        Press%UndefVal, '<NONE>', 0, 0)
      endif
    endif
  end if

  ! Generate the storm center for all time steps
  if (AvgFunc .ne. 'max_azwind') then
    call RecordStormCenter(Press, StmIx, StmIy, MinP)

    do it = 1, Press%Nt
      write (*,*) 'Timestep: ', it
      write (*,'(a,i3,a,i3,a,g,a,g,a)') '  Storm Center: (', StmIx(it), ', ', StmIy(it), &
         ') --> (', Xcoords(StmIx(it)), ', ', Ycoords(StmIy(it)), ')'
      write (*,*) '  Minumum pressure: ', MinP(it)
    end do
    write (*,*) ''
    flush(6)
  endif

  ! call the averaging function

  if ((AvgFunc .eq. 'sc_cloud') .or. (AvgFunc .eq. 'wr_cloud')) then
    call DoCloud(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, Dens, CintLiq, TsAvg)
  else if ((AvgFunc .eq. 'sc_cloud2') .or. (AvgFunc .eq. 'wr_cloud2')) then
    call DoCloud(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, Dens, CintLiq, TsAvg)
  else if (AvgFunc .eq. 'precipr') then
    call DoPrecipR(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, PrecipR, CintLiq, TsAvg)
  else if (AvgFunc .eq. 'w_up') then
    call DoWup(DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, W, TsAvg)
  else if ((AvgFunc .eq. 'sc_cloud_diam') .or. (AvgFunc .eq. 'wr_cloud_diam')) then
    call DoCloudDiam(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, CloudDiam, CintLiq, TsAvg)
  else if ((AvgFunc .eq. 'sc_cloud2_diam') .or. (AvgFunc .eq. 'wr_cloud2_diam')) then
    call DoCloudDiam(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, CloudDiam, CintLiq, TsAvg)
  else if ((AvgFunc .eq. 'sc_cloud_conc') .or. (AvgFunc .eq. 'wr_cloud_conc')) then
    call DoCloudConc(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, DoSc, TempC, CloudConc, TsAvg)
  else if ((AvgFunc .eq. 'sc_cloud2_conc') .or. (AvgFunc .eq. 'wr_cloud2_conc')) then
    call DoCloudConc(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, DoSc, TempC, CloudConc, TsAvg)
  else if (AvgFunc .eq. 'horiz_ke') then
    call DoHorizKe(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, U, V, Dens, TsAvg)
  else if (AvgFunc .eq. 'storm_int') then
    call DoStormInt(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, U, V, TsAvg)
  else if (AvgFunc .eq. 'max_azwind') then
    call DoMaxAzWind(AzWind, TsAvg)
  else if (AvgFunc .eq. 'ccnconc') then
    call DoCcnConc(DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CcnConc, TsAvg)
  else if (AvgFunc .eq. 'test_cvs') then
    call DoTestCvs(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, TestSelect)
  end if

  ! Convert to rates (derivative wrt time) if requested
  if (DoRates .eq. .true.) then
    write (*,*) 'Converting average data to rates'
    write (*,*) ''

    if (AvgFunc .eq. 'test_cvs') then
      write (*,*) 'WARNING: Skipping conversion due to selection of test function'
      write (*,*) ''
    else
      ! Doing the derivative with a difference method which will reduce
      ! the number of time points by two.
      call CalcRates(DeltaT, TsAvg, Rates)
    end if
  end if

  if (AvgFunc .eq. 'test_cvs') then
    call WriteGrads(TestSelect, OfileBases(1), 'ts')
  else
    call WriteGrads(TsAvg, OfileBases(1), 'ts')
    if (DoRates .eq. .true.) then
      call WriteGrads(Rates, OfileBases(2), 'dt')
    end if
  end if

  stop
end

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   Outfiles - output GRADS file, base name for two files
!   AvgFunc - averaging function selection
!

subroutine GetMyArgs(Infiles, Outfiles, AvgFunc, DoRates, Wthreshold, CilThresh, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ)
  implicit none

  character (len=*) :: Infiles, Outfiles, AvgFunc
  real :: Wthreshold, CilThresh, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  logical :: DoRates

  integer :: iargc
  character (len=128) :: arg
  integer :: Nitems

  logical :: BadArgs

  if (iargc() .ne. 12) then
    write (*,*) 'ERROR: must supply exactly 12 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_files> <averaging_function> <do_rates_flag> <w_threshold> <cli_threshold> <min_r> <max_r> <min_phi> <max_phi> <min_z> <max_z>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, colon seprated list, data + rates names, this program will tag on .ctl, .dat suffixes'
    write (*,*) '        <averaging_function>: averaging function to use on input data'
    write (*,*) '            sc_cloud -> total supercooled cloud mass'
    write (*,*) '            sc_cloud_diam -> total supercooled cloud droplet mean diameter'
    write (*,*) '            sc_cloud_conc -> average supercooled cloud droplet concentration'
    write (*,*) '            wr_cloud -> total warm rain cloud mass'
    write (*,*) '            wr_cloud_diam -> total warm rain cloud droplet mean diameter'
    write (*,*) '            wr_cloud_conc -> average warm rain cloud droplet concentration'
    write (*,*) '            sc_cloud2 -> total supercooled cloud2 mass'
    write (*,*) '            sc_cloud2_diam -> total supercooled cloud2 droplet mean diameter'
    write (*,*) '            sc_cloud2_conc -> average supercooled cloud2 droplet concentration'
    write (*,*) '            wr_cloud2 -> total warm rain cloud2 mass'
    write (*,*) '            wr_cloud2_diam -> total warm rain cloud2 droplet mean diameter'
    write (*,*) '            wr_cloud2_conc -> average warm rain cloud2 droplet concentration'
    write (*,*) '            precipr -> total precipitation rate'
    write (*,*) '            w_up -> average w in regions of significant updrafts'
    write (*,*) '            ccnconc -> average ccn concentration'
    write (*,*) '            horiz_ke -> total kinetic energy form horizontal winds'
    write (*,*) '            storm_int -> storm intensity metric from horizontal wind speeds'
    write (*,*) '            max_azwind -> max value of azimuthially averaged wind'
    write (*,*) '            test_cvs -> test the cylindrical volume selection scheme'
    write (*,*) ''
    write (*,*) '        <do_rates_flag>: use "rates" or "no_rates"'
    write (*,*) '            "rates" -> output averaged data plus deriviative wrt time data'
    write (*,*) '            "no_rates" -> output averaged data only'
    write (*,*) ''
    write (*,*) '        The following args are used to select data'
    write (*,*) '          <w_threshold>: select data where w is >= <w_threshold>'
    write (*,*) '          <cil_threshold>: select data where column integrated liquid is >= <cil_thresh>'
    write (*,*) '          Select data inside cylindrical (r,phi,z) volume:'
    write (*,*) '            <min_r> <max_r> (in km)'
    write (*,*) '            <max_phi> <max_phi> (in radians)'
    write (*,*) '            <min_z> <max_z> (in m)'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, Outfiles)

  call getarg(3, AvgFunc)

  call getarg(4, arg)
  if (arg .eq. 'rates') then
    DoRates = .true.
  else
    DoRates = .false.
  end if

  call getarg(5, arg)
  read(arg, '(f)') Wthreshold

  call getarg(6, arg)
  read(arg, '(f)') CilThresh

  call getarg(7, arg)
  read(arg, '(f)') MinR

  call getarg(8, arg)
  read(arg, '(f)') MaxR

  call getarg(9, arg)
  read(arg, '(f)') MinPhi

  call getarg(10, arg)
  read(arg, '(f)') MaxPhi

  call getarg(11, arg)
  read(arg, '(f)') MinZ

  call getarg(12, arg)
  read(arg, '(f)') MaxZ

  BadArgs = .false.

  if ((AvgFunc .ne. 'sc_cloud')       .and. &
      (AvgFunc .ne. 'sc_cloud_diam')  .and. &
      (AvgFunc .ne. 'sc_cloud_conc')  .and. &
      (AvgFunc .ne. 'wr_cloud')       .and. &
      (AvgFunc .ne. 'wr_cloud_diam')  .and. &
      (AvgFunc .ne. 'wr_cloud_conc')  .and. &
      (AvgFunc .ne. 'sc_cloud2')      .and. &
      (AvgFunc .ne. 'sc_cloud2_diam') .and. &
      (AvgFunc .ne. 'sc_cloud2_conc') .and. &
      (AvgFunc .ne. 'wr_cloud2')      .and. &
      (AvgFunc .ne. 'wr_cloud2_diam') .and. &
      (AvgFunc .ne. 'wr_cloud2_conc') .and. &
      (AvgFunc .ne. 'precipr')        .and. &
      (AvgFunc .ne. 'w_up')           .and. &
      (AvgFunc .ne. 'horiz_ke')       .and. &
      (AvgFunc .ne. 'storm_int')      .and. &
      (AvgFunc .ne. 'max_azwind')     .and. &
      (AvgFunc .ne. 'ccnconc')        .and. &
      (AvgFunc .ne. 'test_cvs'))       then
    write (*,*) 'ERROR: <averaging_function> must be one of:'
    write (*,*) '          sc_cloud'
    write (*,*) '          sc_cloud_diam'
    write (*,*) '          sc_cloud_conc'
    write (*,*) '          wr_cloud'
    write (*,*) '          wr_cloud_diam'
    write (*,*) '          wr_cloud_conc'
    write (*,*) '          sc_cloud2'
    write (*,*) '          sc_cloud2_diam'
    write (*,*) '          sc_cloud2_conc'
    write (*,*) '          wr_cloud2'
    write (*,*) '          wr_cloud2_diam'
    write (*,*) '          wr_cloud2_conc'
    write (*,*) '          precipr'
    write (*,*) '          w_up'
    write (*,*) '          ccnconc'
    write (*,*) '          horiz_ke'
    write (*,*) '          storm_int'
    write (*,*) '          max_azwind'
    write (*,*) '          test_cvs'
    write (*,*) ''
    BadArgs = .true.
  end if

  if ((MinR .lt. 0.0) .or. (MaxR .lt. 0.0) .or. (MaxR .le. MinR)) then
    write (*,*) 'ERROR: <min_r> and <max_r> must be >= 0.0, and <max_r> must be > <min_r>'
    write (*,*) ''
    BadArgs = .true.
  end if

  if ((MinPhi .lt. 0.0) .or. (MaxPhi .lt. 0.0) .or. (MaxPhi .le. MinPhi)) then
    write (*,*) 'ERROR: <min_phi> and <max_phi> must be >= 0.0, and <max_phi> must be > <min_phi>'
    write (*,*) ''
    BadArgs = .true.
  end if

  if ((MinZ .lt. 0.0) .or. (MaxZ .lt. 0.0) .or. (MaxZ .le. MinZ)) then
    write (*,*) 'ERROR: <min_z> and <max_z> must be >= 0.0, and <max_z> must be > <min_z>'
    write (*,*) ''
    BadArgs = .true.
  end if

  if (BadArgs) then
    stop
  end if

  return
end subroutine

!*****************************************************************************
! DoCloud()
!
! This subroutine will perform the cloud droplet total mass time series
! averaging. Can select between supercooled or warm rain droplets.
!

subroutine DoCloud(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, Dens, CintLiq, TsAvg)
  use gdata_utils
  use azavg_utils
  implicit none

  logical :: DoSc
  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
  type (GradsVar) :: Cloud, TempC, Dens, CintLiq, TsAvg
  integer, dimension(1:Cloud%Nt) :: StmIx, StmIy
  real, dimension(1:Cloud%Nx) :: Xcoords
  real, dimension(1:Cloud%Ny) :: Ycoords
  real, dimension(1:Cloud%Nz) :: Zcoords

  real, dimension(0:Cloud%Nz) :: ZmHeights
  integer :: ix, iy, iz, it, NumPoints

  call SetZmHeights (Cloud%Nz, ZmHeights)

  ! Convert the cloud mixing ratio to mass for each grid point. 
  !
  ! Mixing ratios in GRADS files are g/kg
  ! Density is in kg/m**3
  ! Heights are in m
  ! So, express the mass value in g using the formula
  !   (mix ratio) * (density) * (layer thickness) * (layer horiz area)
  ! The layer thickness for layer k is: ZmHeights(k) - ZmHeights(k-1)

  do it = 1, Cloud%Nt
    ! Sum up the cloud droplet mass over the specified radial band. Only include the
    ! grid points where tempc is 0 or less (supercooled)

    TsAvg%Vdata(it,1,1,1) = 0.0
    NumPoints = 0

    do iz = 1, Cloud%Nz
      do ix = 1, Cloud%Nx
        do iy = 1, Cloud%Ny
          if ((InsideCylVol(Cloud%Nx, Cloud%Ny, Cloud%Nz, Cloud%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords))  .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
            if (DoSc) then
              if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
                TsAvg%Vdata(it,1,1,1) = TsAvg%Vdata(it,1,1,1) + (Cloud%Vdata(it,iz,ix,iy) * Dens%Vdata(it,iz,ix,iy) * (ZmHeights(iz)-ZmHeights(iz-1)))
                NumPoints = NumPoints + 1
              end if
            else
              if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
                TsAvg%Vdata(it,1,1,1) = TsAvg%Vdata(it,1,1,1) + (Cloud%Vdata(it,iz,ix,iy) * Dens%Vdata(it,iz,ix,iy) * (ZmHeights(iz)-ZmHeights(iz-1)))
                NumPoints = NumPoints + 1
              end if
            end if
          end if
        end do
      end do
    end do

    ! At this point TsAvg%Vdata(it,1,1,1) holds g/m**2, multiply by grid cell horizontal area. Note this assumes
    ! each grid cell has the same horizontal area.

    TsAvg%Vdata(it,1,1,1) = TsAvg%Vdata(it,1,1,1) * DeltaX * DeltaY
    if (NumPoints .eq. 0) then
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      write (*,*) 'ScCloud: Timestep:', it, ', Number of points selected: ', NumPoints
    endif
  end do
end subroutine

!************************************************************************************
! DoWup()
!
! This subroutine will do the average vertical velocity in regions of significant
! updrafts.
!

subroutine DoWup(DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, W, TsAvg)
  use gdata_utils
  use azavg_utils
  implicit none

  real :: DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  type (GradsVar) :: W, TsAvg
  integer, dimension(1:W%Nt) :: StmIx, StmIy
  real, dimension(1:W%Nx) :: Xcoords
  real, dimension(1:W%Ny) :: Ycoords
  real, dimension(1:W%Nz) :: Zcoords

  integer ix,iy,iz,it, NumPoints

  do it = 1, W%Nt
    ! Average w over regions where significant updrafts occur

    TsAvg%Vdata(it,1,1,1) = 0.0
    NumPoints = 0

    do iz = 1, W%Nz
      do ix = 1, W%Nx
        do iy = 1, W%Ny
          if (InsideCylVol(W%Nx, W%Ny, W%Nz, W%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            if (W%Vdata(it,iz,ix,iy) .ge. Wthreshold) then
              TsAvg%Vdata(it,1,1,1) = TsAvg%Vdata(it,1,1,1) + W%Vdata(it,iz,ix,iy)
              NumPoints = NumPoints + 1
            end if
          end if
        end do
      end do
    end do

    if (NumPoints .eq. 0) then
      TsAvg%Vdata(it,1,1,1) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TsAvg%Vdata(it,1,1,1) = TsAvg%Vdata(it,1,1,1) / float(NumPoints)
      write (*,*) 'Wup: Timestep:', it, ', Number of points selected: ', NumPoints
    end if
  end do

end subroutine


!************************************************************************************
! DoCloudDiam()
!
! This routine will calculate an averaged cloud droplet diameter. Can select between
! warm rain droplets or supercooled droplets.

subroutine DoCloudDiam(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, CloudDiam, CintLiq, TsAvg)
  use gdata_utils
  use azavg_utils
  implicit none

  logical :: DoSc
  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
  type (GradsVar) :: Cloud, TempC, CloudDiam, CintLiq, TsAvg
  integer, dimension(1:CloudDiam%Nt) :: StmIx, StmIy
  real, dimension(1:CloudDiam%Nx) :: Xcoords
  real, dimension(1:CloudDiam%Ny) :: Ycoords
  real, dimension(1:CloudDiam%Nz) :: Zcoords

  integer ix,iy,iz,it, NumPoints
  real SumQ, SumQD
  real MaxQ, Climit, SumD

  do it = 1, CloudDiam%Nt
    ! Calculate a mass-weighted mean diameter for supercooled cloud droplets.
    ! 
    !    Mean diameter (TsAvg value) = Sum(cloud * cloud_d) / Sum(cloud)
    !    where cloud and cloud_d are only included in the sum when tempc is <= 0
    !
    ! Find the max Q (mass) value and use it to filter data - select data points if
    ! the Q is within 20% of the max Q (.2 to 1.0). Then form the average of the diameters
    ! of the selected points.
    SumQD = 0.0
    SumQ = 0.0
    SumD = 0.0
    NumPoints = 0
    do iz = 1, CloudDiam%Nz
      do ix = 1, CloudDiam%Nx
        do iy = 1, CloudDiam%Ny
          if ((InsideCylVol(CloudDiam%Nx, CloudDiam%Ny, CloudDiam%Nz, CloudDiam%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
            if (DoSc) then
              if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
                 SumQD = SumQD + (Cloud%Vdata(it,iz,ix,iy) * CloudDiam%Vdata(it,iz,ix,iy))
                 SumQ = SumQ + Cloud%Vdata(it,iz,ix,iy)
                 SumD = SumD + CloudDiam%Vdata(it,iz,ix,iy)
                 NumPoints = NumPoints + 1
              end if
            else
              if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
                 SumQD = SumQD + (Cloud%Vdata(it,iz,ix,iy) * CloudDiam%Vdata(it,iz,ix,iy))
                 SumQ = SumQ + Cloud%Vdata(it,iz,ix,iy)
                 SumD = SumD + CloudDiam%Vdata(it,iz,ix,iy)
                 NumPoints = NumPoints + 1
              end if
            end if
          end if
        end do
      end do
    end do

    if (SumQ .eq. 0.0) then
      TsAvg%Vdata(it,1,1,1) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TsAvg%Vdata(it,1,1,1) = SumQD / SumQ
      write (*,*) 'ScCloudDiam: Ts:', it, ', NumPoints: ', NumPoints
    end if

  end do
end subroutine

!************************************************************************************
! DoCloudConc()
!
! This subroutine will do the average cloud droplet concentration

subroutine DoCloudConc(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, DoSc, TempC, CloudConc, TsAvg)
  use gdata_utils
  use azavg_utils
  implicit none

  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  type (GradsVar) :: CloudConc, TempC, TsAvg
  integer, dimension(1:CloudConc%Nt) :: StmIx, StmIy
  real, dimension(1:CloudConc%Nx) :: Xcoords
  real, dimension(1:CloudConc%Ny) :: Ycoords
  real, dimension(1:CloudConc%Nz) :: Zcoords
  logical :: DoSc

  integer ix,iy,iz,it
  integer NumPoints
  real SumCloudConc

  ! Calculate the average cloud droplet concentration near the eyewall region.
  ! 
  ! Call the cloud base to be between 1000m and 1200m. The z level corresponding to
  ! that is z = 4 which is 1138m. Cover from surface to the cloud base level. This
  ! results in using z = 1 to z = 4.

  do it = 1, CloudConc%Nt
    SumCloudConc = 0.0
    NumPoints = 0

    do iz = 1, CloudConc%Nz
      do ix = 1, CloudConc%Nx
        do iy = 1, CloudConc%Ny
          if (InsideCylVol(CloudConc%Nx, CloudConc%Ny, CloudConc%Nz, CloudConc%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            if (DoSc) then
              if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
                SumCloudConc = SumCloudConc + CloudConc%Vdata(it,iz,ix,iy)
                NumPoints = NumPoints + 1
              end if
            else
              if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
                SumCloudConc = SumCloudConc + CloudConc%Vdata(it,iz,ix,iy)
                NumPoints = NumPoints + 1
              end if
            end if
          end if
        end do
      end do
    end do

    if (NumPoints .eq. 0) then
      TsAvg%Vdata(it,1,1,1) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TsAvg%Vdata(it,1,1,1) = SumCloudConc / float(NumPoints)
      write (*,*) 'CloudConc: Timestep:', it, ', Number of points selected: ', NumPoints
    end if
  end do
end subroutine

!************************************************************************************
! DoPrecipR()
!
! This subroutine will do the average cloud droplet concentration near the eyewall
! region.

subroutine DoPrecipR(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, PrecipR, CintLiq, TsAvg)
  use gdata_utils
  use azavg_utils
  implicit none

  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
  ! PrecipR is 2D var
  type (GradsVar) :: PrecipR, CintLiq, TsAvg
  integer, dimension(1:PrecipR%Nt) :: StmIx, StmIy
  real, dimension(1:PrecipR%Nx) :: Xcoords
  real, dimension(1:PrecipR%Ny) :: Ycoords
  real, dimension(1:PrecipR%Nz) :: Zcoords

  integer :: ix,iy,iz,it
  integer :: NumPoints
  real :: SumPrecip

  do it = 1, PrecipR%Nt
    SumPrecip = 0.0
    NumPoints = 0

    do ix = 1, PrecipR%Nx
      do iy = 1, PrecipR%Ny
        if (InsideCylVol(PrecipR%Nx, PrecipR%Ny, PrecipR%Nz, PrecipR%Nt, ix, iy, 1, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords) .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
          SumPrecip = SumPrecip + PrecipR%Vdata(it,1,ix,iy)
          NumPoints = NumPoints + 1
        end if
      end do
    end do

    if (NumPoints .eq. 0) then
      TsAvg%Vdata(it,1,1,1) = 0.0
      write (*,*) '  WARNING: no data points selected for this time step'
    else
      ! At this point TsAvg%Vdata(it,1,1,1) holds mm/hr, multiply by grid cell horizontal area.
      ! Note this assumes each grid cell has the same horizontal area. What this does
      ! is convert mm/hr to kg/hr when assuming that the density of water is 1000kg/m**3.
      !   mm/hr * m**2 * 1000 kg/m**3 * 0.001 m/mm -> kg/hr
      TsAvg%Vdata(it,1,1,1) = (SumPrecip / float(NumPoints)) * DeltaX * DeltaY
      write (*,*) 'Precip: Timestep:', it, ', Number of points selected: ', NumPoints
    end if
    write (*,*) ''
    flush(6)
  end do
end subroutine

!**************************************************************************************
! DoHorizKe()
!
! This routine will calculate the total kinetic energy over the given cylindrical
! volume. Do not want average since we want the size of the storm reflected in
! this diagnostic.

subroutine DoHorizKe(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, U, V, Dens, TsAvg)
  use gdata_utils
  use azavg_utils
  implicit none

  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  type (GradsVar) :: U, V, Dens, TsAvg
  integer, dimension(1:U%Nt) :: StmIx, StmIy
  real, dimension(1:U%Nx) :: Xcoords
  real, dimension(1:U%Ny) :: Ycoords
  real, dimension(1:U%Nz) :: Zcoords

  integer ix,iy,iz,it
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

  do it = 1, U%Nt
    SumKe = 0.0
    NumPoints = 0

    do iz = 1, U%Nz
      do ix = 1, U%Nx
        do iy = 1, U%Ny
          if (InsideCylVol(U%Nx, U%Ny, U%Nz, U%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            if (iz .eq. U%Nz) then
              ! Use the level below for this case (since no level above)
              LevThickness = Zcoords(iz) - Zcoords(iz-1)
            else
              LevThickness = Zcoords(iz+1) - Zcoords(iz)
            end if
            CurrKe = 0.5 * DeltaX * DeltaY * LevThickness * Dens%Vdata(it,iz,ix,iy) * (U%Vdata(it,iz,ix,iy)**2 + V%Vdata(it,iz,ix,iy)**2)
            SumKe = SumKe + CurrKe
            NumPoints = NumPoints + 1
          end if
        end do
      end do
    end do

    TsAvg%Vdata(it,1,1,1) = SumKe;
    write (*,*) 'HorizKe: Timestep:', it, ', Number of points selected: ', NumPoints
  end do
  
end subroutine

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

subroutine DoStormInt(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, U, V, TsAvg)
  use gdata_utils
  use azavg_utils
  implicit none

  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  type (GradsVar) :: U, V, TsAvg
  integer, dimension(1:U%Nt) :: StmIx, StmIy
  real, dimension(1:U%Nx) :: Xcoords
  real, dimension(1:U%Ny) :: Ycoords
  real, dimension(1:U%Nz) :: Zcoords

  integer ix,iy,iz,it
  integer nCat0, nCat1, nCat2, nCat3, nCat4, nCat5, NumPoints
  real Wspeed, SiMetric

  do it = 1, U%Nt
    nCat0 = 0
    nCat1 = 0
    nCat2 = 0
    nCat3 = 0
    nCat4 = 0
    nCat5 = 0

    do iz = 1, U%Nz
      do ix = 1, U%Nx
        do iy = 1, U%Ny
          if (InsideCylVol(U%Nx, U%Ny, U%Nz, U%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            ! Count up the number of grid points with wind speeds fitting each of the
            ! Saffir-Simpson categories. Then form the metric by weighting each category
            ! count.
            Wspeed = sqrt(U%Vdata(it,iz,ix,iy)**2 + V%Vdata(it,iz,ix,iy)**2)
            if (Wspeed .ge. 70.0) then
              nCat5 = nCat5 + 1
            else
              if (Wspeed .ge. 59.0) then
                nCat4 = nCat4 + 1
              else
                if (Wspeed .ge. 50.0) then
                  nCat3 = nCat3 + 1
                else
                  if (Wspeed .ge. 43.0) then
                    nCat2 = nCat2 + 1
                  else
                    if (Wspeed .ge. 33.0) then
                      nCat1 = nCat1 + 1
                    else
                      nCat0 = nCat0 + 1
                    end if
                  end if
                end if
              end if
            end if
          end if
        end do
      end do
    end do

    !Linear weighting
    NumPoints = nCat0 + nCat1 + nCat2 + nCat3 + nCat4 + nCat5
    SiMetric = float(nCat1) + (float(nCat2)*2.0) + (float(nCat3)*3.0) + (float(nCat4)*4.0) + (float(nCat5)*5.0)

    if (NumPoints .eq. 0) then
      TsAvg%Vdata(it,1,1,1) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TsAvg%Vdata(it,1,1,1) = SiMetric / float(NumPoints)
      write (*,*) 'StormInt: Timestep:', it, ', Number of points selected: ', NumPoints
    end if
  end do
  
end subroutine

!************************************************************************************
! DoCcnConc()
!
! This subroutine will do the average CCN concentration.
!

subroutine DoCcnConc(DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CcnConc, TsAvg)
  use gdata_utils
  use azavg_utils
  implicit none

  real :: DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  type (GradsVar) :: CcnConc, TsAvg
  integer, dimension(1:CcnConc%Nt) :: StmIx, StmIy
  real, dimension(1:CcnConc%Nx) :: Xcoords
  real, dimension(1:CcnConc%Ny) :: Ycoords
  real, dimension(1:CcnConc%Nz) :: Zcoords

  integer ix,iy,iz,it, NumPoints

  do it = 1, CcnConc%Nt
    ! Average w over regions where significant updrafts occur

    TsAvg%Vdata(it,1,1,1) = 0.0
    NumPoints = 0

    do iz = 1, CcnConc%Nz
      do ix = 1, CcnConc%Nx
        do iy = 1, CcnConc%Ny
          if (InsideCylVol(CcnConc%Nx, CcnConc%Ny, CcnConc%Nz, CcnConc%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            TsAvg%Vdata(it,1,1,1) = TsAvg%Vdata(it,1,1,1) + CcnConc%Vdata(it,iz,ix,iy)
            NumPoints = NumPoints + 1
          end if
        end do
      end do
    end do

    if (NumPoints .eq. 0) then
      TsAvg%Vdata(it,1,1,1) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TsAvg%Vdata(it,1,1,1) = TsAvg%Vdata(it,1,1,1) / float(NumPoints)
      write (*,*) 'Wup: Timestep:', it, ', Number of points selected: ', NumPoints
    end if
  end do

end subroutine

!************************************************************************************
! DoMaxAzWind()
!
! This subroutine will simply find the maximum wind speed in AzWind and copy that
! to TsAvg.
subroutine DoMaxAzWind(AzWind, TsAvg)
  use gdata_utils
  implicit none

  type (GradsVar) :: AzWind, TsAvg

  integer :: ix,iy,iz,it

  do it = 1, AzWind%Nt
    TsAvg%Vdata(it,1,1,1) = 0.0
    do iz = 1, AzWind%Nz
      do ix = 1, AzWind%Nx
        do iy = 1, AzWind%Ny
          if (AzWind%Vdata(it,iz,ix,iy) .gt. TsAvg%Vdata(it,1,1,1)) then
            TsAvg%Vdata(it,1,1,1) = AzWind%Vdata(it,iz,ix,iy)
          endif
        end do
      end do
    end do
  end do

  return
end subroutine

!************************************************************************************
! DoTestCvs()
!
! This subroutine will perform a test on the cylindrical volume selection routine.
! Just runs through the entire grid and outputs a '1' when selection occurs otherwise
! outputs a '0'. Then view the result in grads and see if selection is correct.

subroutine DoTestCvs(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, TestSelect)
  use gdata_utils
  use azavg_utils
  implicit none

  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  type (GradsVar) :: TestSelect
  integer, dimension(1:TestSelect%Nt) :: StmIx, StmIy
  real, dimension(1:TestSelect%Nx) :: Xcoords
  real, dimension(1:TestSelect%Ny) :: Ycoords
  real, dimension(1:TestSelect%Nz) :: Zcoords

  integer ix,iy,iz,it
  integer NumPoints

  write (*,*) 'Testing cylindrical volume selection:'
  do it = 1, TestSelect%Nt
    NumPoints = 0
    do iz = 1, TestSelect%Nz
      do ix = 1, TestSelect%Nx
        do iy = 1, TestSelect%Ny
          if (InsideCylVol(TestSelect%Nx, TestSelect%Ny, TestSelect%Nz, TestSelect%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            TestSelect%Vdata(it,iz,ix,iy) = 1.0
            NumPoints = NumPoints + 1
          else
            TestSelect%Vdata(it,iz,ix,iy) = 0.0
          end if
        end do
      end do
    end do
    ! mark the storm center
    TestSelect%Vdata(it,iz,StmIx(it),StmIy(it)) = 2.0
    write (*,*) '  Timestep, Number of points selected: ', it, NumPoints
  end do
end subroutine

!******************************************************************************
! ConvertTinc()
!
! This function will convert the time increment spec'd in the GRAD control
! file into a number of seconds.
!

real function ConvertTinc(Tinc)
  implicit none

  character (len=*) :: Tinc

  character (len=128) :: Tval, Tunits, Tfmt
  integer :: i, Uloc, Tlen, Itval
  
  ! Walk through the Tinc string. Concatenate the numbers onto Tval and
  ! the alaph characters onto Tunits. This algorithm assumes that GRADS
  ! will use a format like <numeric_value><units> for Tinc where <units>
  ! is a string like 'hr' or 'mn'.
  Uloc = 0
  Tlen = len_trim(Tinc)
  do i = 1, Tlen
    if ((Tinc(i:i) .ge. 'A') .and. (Tinc(i:i) .le. 'z')) then
      if (Uloc .eq. 0) then
        ! the first alpha character is the beginning of the spec for units
        Uloc = i
      end if
    end if
  end do

  write (Tfmt, '(a,i2.2,a,i2.2,a)') '(a', (Uloc-1), 'a' , ((Tlen-Uloc)+1), ')'
  read (Tinc, Tfmt) Tval, Tunits

  read(Tval, '(i)') Itval
  if (Tunits .eq. 'hr') then
    ConvertTinc = float(Itval) * 3600.0
  else
    if (Tunits .eq. 'mn') then
      ConvertTinc = float(Itval) * 60.0
    else
      ConvertTinc = float(Itval)
    end if
  end if
  return
end function

!******************************************************************************
! CalcRates()
!
! This routine will calculate time derivatives of the input data using
! a centered difference method.
!

subroutine CalcRates(DeltaT, TsAvg, Rates)
  use gdata_utils
  implicit none

  type (GradsVar) :: TsAvg, Rates
  real :: DeltaT

  real :: f1, f2
  integer :: ix, iy, iz, it

  ! use a centered difference, uses points at t-1, t and t+1
  do ix = 1, TsAvg%Nx
    do iy = 1, TsAvg%Ny
      do iz = 1, TsAvg%Nz
        do it = 2, TsAvg%Nt-1
          f1 = (TsAvg%Vdata(it,1,1,1) + TsAvg%Vdata(it-1,1,1,1)) / 2.0
          f2 = (TsAvg%Vdata(it+1,1,1,1) + TsAvg%Vdata(it,1,1,1)) / 2.0

          Rates%Vdata(it-1,1,1,1) = (f2 - f1) / DeltaT
        end do
      end do
    end do
  end do
end subroutine
