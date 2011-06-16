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

include 'gDataTypes.h'

program main
  use GfileTypes
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  character (len=LargeString) :: Infiles, Outfiles
  character (len=LittleString) :: AvgFunc

  logical :: DoRates, DoSc

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles, OfileBases
  integer :: Nfiles, Nofiles
  integer, dimension(:), allocatable :: StmIx, StmIy

  real :: DeltaX, DeltaY, DeltaT, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
  real :: Xstart, Xinc, Ystart, Yinc
  real :: ConvertTinc
  real, dimension(:), allocatable :: ZmHeights
  real, dimension(1:1) :: DummyZcoords
  real, dimension(:), allocatable :: MinP, Xcoords, Ycoords, Zcoords

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars
  type (GradsOutDescription) :: GoutDescrip

  ! Data arrays: cloud, tempc, precipr
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of cloud, tempc, precipr in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: Cloud, TempC, Dens, PrecipR, U, V, W, CcnConc, CloudDiam, CloudConc, Press, CintLiq, TsAvg, TestSelect, Rates
  type (GradsVarLocation) :: CloudLoc, TempcLoc, DensLoc, ULoc, VLoc, WLoc, CcnLoc, CloudDiamLoc, CloudConcLoc, PreciprLoc, PressLoc, CintLiqLoc

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it

  ! Get the command line arguments
  call GetMyArgs(Infiles, Outfiles, AvgFunc, DoRates, Wthreshold, CilThresh, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')
  call String2List(Outfiles, ':', OfileBases, MaxFiles, Nofiles, 'output files')

  DoSc = .true.
  if ((AvgFunc .eq. 'wr_cloud') .or. (AvgFunc .eq. 'wr_cloud_diam')) then
    DoSc = .false.
  end if

  write (*,*) 'Time seris of average for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
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
  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''
  flush(6)

  ! Check the data description for consistency and locate the variables in the GRADS control files
  if ((AvgFunc .eq. 'sc_cloud') .or. (AvgFunc .eq. 'wr_cloud')) then
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, TempcLoc, 'tempc')
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, DensLoc, 'dn0')
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudLoc, 'cloud')
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CintLiqLoc, 'liquid_colint')
  else
    if (AvgFunc .eq. 'precipr') then
      call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PreciprLoc, 'precipr')
      call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
      call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CintLiqLoc, 'liquid_colint')
    else
      if (AvgFunc .eq. 'w_up') then
        call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, WLoc, 'w')
        call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
      else
        if ((AvgFunc .eq. 'sc_cloud_diam') .or. (AvgFunc .eq. 'wr_cloud_diam')) then
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, TempcLoc, 'tempc')
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudLoc, 'cloud')
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudDiamLoc, 'cloud_d')
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CintLiqLoc, 'liquid_colint')
        else
          if (AvgFunc .eq. 'ew_cloud') then
            call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudConcLoc, 'cloud_cm3')
            call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
          else
            if (AvgFunc .eq. 'horiz_ke') then
              call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, ULoc, 'u')
              call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, VLoc, 'v')
              call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, DensLoc, 'dn0')
              call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
            else
              if (AvgFunc .eq. 'storm_int') then
                call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, ULoc, 'u')
                call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, VLoc, 'v')
                call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
              else
                if (AvgFunc .eq. 'ccnconc') then
                  call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CcnLoc, 'ccnconcen')
                  call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
                else
                  if (AvgFunc .eq. 'test_cvs') then
                    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end if
  end if

  ! Calculate the x,y coordinates (in km) for doing selection by radius from storm center.
  ! Also for setting DeltaX and DeltaY (in m)
  allocate (Xcoords(1:Nx), Ycoords(1:Ny), Zcoords(1:Nz), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  call ConvertGridCoords(Nx, Ny, GdataDescrip(1), Xcoords, Ycoords)
  ! Using pressure for all averaging function, get the z coords from
  ! the description of the pressure variable
  do iz = 1, Nz
    Zcoords(iz) = GdataDescrip(PressLoc%Fnum)%Zcoords(iz)
  end do

  DeltaX = (Xcoords(2) - Xcoords(1)) * 1000.0
  DeltaY = (Ycoords(2) - Ycoords(1)) * 1000.0

  ! Convert string time increment to seconds
  DeltaT = ConvertTinc(GdataDescrip(1)%Tinc)

  write (*,*) 'Horizontal grid info:'
  write (*,*) '  X range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', GdataDescrip(1)%Xcoords(1), GdataDescrip(1)%Xcoords(Nx), Xcoords(1), Xcoords(Nx)
  write (*,*) '  Y range (min lat, max lat) --> (min y, max y): '
  write (*,*) '    ', GdataDescrip(1)%Ycoords(1), GdataDescrip(1)%Ycoords(Ny), Ycoords(1), Ycoords(Ny)
  write (*,*) ''
  write (*,*) 'Vertical grid info:'
  do iz = 1, Nz
    write (*,*) '  ', iz, ' --> ', Zcoords(iz)
  end do
  write (*,*) ''
  flush(6)

  ! Read in the GRADS variable data
  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', Nx
  write (*,*) '  Number of y (latitude) points:           ', Ny
  write (*,*) '  Number of z (vertical level) points:     ', Nz
  write (*,*) '  Number of t (time) points:               ', Nt
  write (*,*) '  Total number of grid variables:          ', Nvars
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''
  write (*,*) '  Grid delta x: ', DeltaX
  write (*,*) '  Grid delta y: ', DeltaY
  write (*,*) '  Time increment: ', trim(GdataDescrip(1)%Tinc), ' --> ', DeltaT
  write (*,*) ''
  flush(6)

  write (*,*) 'Locations of variables in GRADS data (file number, var number):'
  if ((AvgFunc .eq. 'sc_cloud') .or. (AvgFunc .eq. 'wr_cloud')) then
    write (*,'(a20,i3,a2,i3,a1)') 'tempc: (', TempcLoc%Fnum, ', ', TempcLoc%Vnum, ')'
    write (*,'(a20,i3,a2,i3,a1)') 'dn0: (', DensLoc%Fnum, ', ', DensLoc%Vnum, ')'
    write (*,'(a20,i3,a2,i3,a1)') 'cloud: (', CloudLoc%Fnum, ', ', CloudLoc%Vnum, ')'
    write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'
    write (*,'(a20,i3,a2,i3,a1)') 'liquid_colint: (', CintLiqLoc%Fnum, ', ', CintLiqLoc%Vnum, ')'

    ! Allocate the data arrays and read in the data from the GRADS data files
    allocate (TempC(1:Nx,1:Ny,1:Nz,1:Nt), Dens(1:Nx,1:Ny,1:Nz,1:Nt), &
              Cloud(1:Nx,1:Ny,1:Nz,1:Nt), Press(1:Nx,1:Ny,1:Nz,1:Nt), &
              CintLiq(1:Nx,1:Ny,1:1,1:Nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory allocation failed'
      stop
    end if

    ! Read in the data for the vars using the description and location information
    call ReadGradsData(GdataDescrip, 'tempc', TempcLoc, TempC, Nx, Ny, Nz, Nt)
    call ReadGradsData(GdataDescrip, 'dn0', DensLoc, Dens, Nx, Ny, Nz, Nt)
    call ReadGradsData(GdataDescrip, 'cloud', CloudLoc, Cloud, Nx, Ny, Nz, Nt)
    call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
    call ReadGradsData(GdataDescrip, 'liquid_colint', CintLiqLoc, CintLiq, Nx, Ny, 1, Nt)
  else
    if (AvgFunc .eq. 'precipr') then
      write (*,'(a20,i3,a2,i3,a1)') 'precipr: (', PreciprLoc%Fnum, ', ', PreciprLoc%Vnum, ')'
      write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'
      write (*,'(a20,i3,a2,i3,a1)') 'liquid_colint: (', CintLiqLoc%Fnum, ', ', CintLiqLoc%Vnum, ')'

      ! Allocate the data arrays and read in the data from the GRADS data files
      ! precipr is 2D variable -> Nz = 1
      allocate (PrecipR(1:Nx,1:Ny,1:1,1:Nt), Press(1:Nx,1:Ny,1:Nz,1:Nt), &
                CintLiq(1:Nx,1:Ny,1:1,1:Nt), stat=Ierror)
      if (Ierror .ne. 0) then
        write (*,*) 'ERROR: Data array memory allocation failed'
        stop
      end if

      ! Read in the data for the vars using the description and location information
      ! precipr is 2D variable -> Nz = 1
      call ReadGradsData(GdataDescrip, 'precipr', PreciprLoc, PrecipR, Nx, Ny, 1, Nt)
      call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
      call ReadGradsData(GdataDescrip, 'liquid_colint', CintLiqLoc, CintLiq, Nx, Ny, 1, Nt)
    else
      if (AvgFunc .eq. 'w_up') then
        write (*,'(a20,i3,a2,i3,a1)') 'w: (', WLoc%Fnum, ', ', WLoc%Vnum, ')'
        write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'

        ! Allocate the data arrays and read in the data from the GRADS data files
        allocate (W(1:Nx,1:Ny,1:Nz,1:Nt), Press(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
        if (Ierror .ne. 0) then
          write (*,*) 'ERROR: Data array memory allocation failed'
          stop
        end if

        ! Read in the data for the vars using the description and location information
        call ReadGradsData(GdataDescrip, 'w', WLoc, W, Nx, Ny, Nz, Nt)
        call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
      else
        if ((AvgFunc .eq. 'sc_cloud_diam') .or. (AvgFunc .eq. 'wr_cloud_diam')) then
          write (*,'(a20,i3,a2,i3,a1)') 'tempc: (', TempcLoc%Fnum, ', ', TempcLoc%Vnum, ')'
          write (*,'(a20,i3,a2,i3,a1)') 'cloud: (', CloudLoc%Fnum, ', ', CloudLoc%Vnum, ')'
          write (*,'(a20,i3,a2,i3,a1)') 'cloud_d: (', CloudDiamLoc%Fnum, ', ', CloudDiamLoc%Vnum, ')'
          write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'
          write (*,'(a20,i3,a2,i3,a1)') 'liquid_colint: (', CintLiqLoc%Fnum, ', ', CintLiqLoc%Vnum, ')'
  
          ! Allocate the data arrays and read in the data from the GRADS data files
          allocate (TempC(1:Nx,1:Ny,1:Nz,1:Nt), CloudDiam(1:Nx,1:Ny,1:Nz,1:Nt), &
                    Cloud(1:Nx,1:Ny,1:Nz,1:Nt), Press(1:Nx,1:Ny,1:Nz,1:Nt), &
                    CintLiq(1:Nx,1:Ny,1:1,1:Nt), stat=Ierror)
          if (Ierror .ne. 0) then
            write (*,*) 'ERROR: Data array memory allocation failed'
            stop
          end if

          ! Read in the data for the vars using the description and location information
          call ReadGradsData(GdataDescrip, 'tempc', TempcLoc, TempC, Nx, Ny, Nz, Nt)
          call ReadGradsData(GdataDescrip, 'cloud', CloudLoc, Cloud, Nx, Ny, Nz, Nt)
          call ReadGradsData(GdataDescrip, 'cloud_d', CloudDiamLoc, CloudDiam, Nx, Ny, Nz, Nt)
          call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
          call ReadGradsData(GdataDescrip, 'liquid_colint', CintLiqLoc, CintLiq, Nx, Ny, 1, Nt)
        else
          if (AvgFunc .eq. 'ew_cloud') then
            write (*,'(a20,i3,a2,i3,a1)') 'cloudconcen_cm3: (', CloudConcLoc%Fnum, ', ', CloudConcLoc%Vnum, ')'
            write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'
    
            ! Allocate the data arrays and read in the data from the GRADS data files
            allocate (CloudConc(1:Nx,1:Ny,1:Nz,1:Nt), Press(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
            if (Ierror .ne. 0) then
              write (*,*) 'ERROR: Data array memory allocation failed'
              stop
            end if
  
            ! Read in the data for the vars using the description and location information
            call ReadGradsData(GdataDescrip, 'cloud_cm3', CloudConcLoc, CloudConc, Nx, Ny, Nz, Nt)
            call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
          else
            if (AvgFunc .eq. 'horiz_ke') then
              write (*,'(a20,i3,a2,i3,a1)') 'u: (', ULoc%Fnum, ', ', ULoc%Vnum, ')'
              write (*,'(a20,i3,a2,i3,a1)') 'v: (', VLoc%Fnum, ', ', VLoc%Vnum, ')'
              write (*,'(a20,i3,a2,i3,a1)') 'dn0: (', DensLoc%Fnum, ', ', DensLoc%Vnum, ')'
              write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'

              ! Allocate the data arrays and read in the data from the GRADS data files
              allocate (U(1:Nx,1:Ny,1:Nz,1:Nt), V(1:Nx,1:Ny,1:Nz,1:Nt), &
                        Dens(1:Nx,1:Ny,1:Nz,1:Nt), Press(1:Nx,1:Ny,1:Nz,1:Nt), &
                        stat=Ierror)
              if (Ierror .ne. 0) then
                write (*,*) 'ERROR: Data array memory allocation failed'
                stop
              end if

              ! Read in the data for the vars using the description and location information
              call ReadGradsData(GdataDescrip, 'u', ULoc, U, Nx, Ny, Nz, Nt)
              call ReadGradsData(GdataDescrip, 'v', VLoc, V, Nx, Ny, Nz, Nt)
              call ReadGradsData(GdataDescrip, 'dn0', DensLoc, Dens, Nx, Ny, Nz, Nt)
              call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
            else
              if (AvgFunc .eq. 'storm_int') then
                write (*,'(a20,i3,a2,i3,a1)') 'u: (', ULoc%Fnum, ', ', ULoc%Vnum, ')'
                write (*,'(a20,i3,a2,i3,a1)') 'v: (', VLoc%Fnum, ', ', VLoc%Vnum, ')'
                write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'
  
                ! Allocate the data arrays and read in the data from the GRADS data files
                allocate (U(1:Nx,1:Ny,1:Nz,1:Nt), V(1:Nx,1:Ny,1:Nz,1:Nt), &
                          Press(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
                if (Ierror .ne. 0) then
                  write (*,*) 'ERROR: Data array memory allocation failed'
                  stop
                end if
  
                ! Read in the data for the vars using the description and location information
                call ReadGradsData(GdataDescrip, 'u', ULoc, U, Nx, Ny, Nz, Nt)
                call ReadGradsData(GdataDescrip, 'v', VLoc, V, Nx, Ny, Nz, Nt)
                call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
              else
                if (AvgFunc .eq. 'ccnconc') then
                  write (*,'(a20,i3,a2,i3,a1)') 'ccnconcen: (', CcnLoc%Fnum, ', ', CcnLoc%Vnum, ')'
                  write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'

                  ! Allocate the data arrays and read in the data from the GRADS data files
                  allocate (CcnConc(1:Nx,1:Ny,1:Nz,1:Nt), Press(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
                  if (Ierror .ne. 0) then
                    write (*,*) 'ERROR: Data array memory allocation failed'
                    stop
                  end if

                  ! Read in the data for the vars using the description and location information
                  call ReadGradsData(GdataDescrip, 'ccnconcen', CcnLoc, CcnConc, Nx, Ny, Nz, Nt)
                  call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
                else
                  if (AvgFunc .eq. 'test_cvs') then
                    write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'
                    ! Allocate the data arrays and read in the data from the GRADS data files
                    allocate (Press(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
                    if (Ierror .ne. 0) then
                      write (*,*) 'ERROR: Data array memory allocation failed'
                      stop
                    end if
      
                    ! Read in the data for the vars using the description and location information
                    call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end if
  end if
  write (*,*) ''
  flush(6)

  ! Allocate the output array and do the averaging
  if (AvgFunc .eq. 'test_cvs') then
    allocate (TestSelect(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Ouput data array memory allocation failed'
      stop
    end if
  else
    allocate (TsAvg(1:1,1:1,1:1,1:Nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Ouput data array memory allocation failed'
      stop
    end if
  end if

  ! Generate the storm center for all time steps
  allocate (StmIx(1:Nt), StmIy(1:Nt), MinP(1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Ouput data array memory allocation failed'
    stop
  end if
  call RecordStormCenter(Nx, Ny, Nz, Nt, Press, StmIx, StmIy, MinP)

  do it = 1, Nt
    write (*,*) 'Timestep: ', it
    write (*,'(a,i3,a,i3,a,g,a,g,a)') '  Storm Center: (', StmIx(it), ', ', StmIy(it), &
       ') --> (', Xcoords(StmIx(it)), ', ', Ycoords(StmIy(it)), ')'
    write (*,*) '  Minumum pressure: ', MinP(it)
  end do
  write (*,*) ''
  flush(6)

  ! call the averaging function

  if ((AvgFunc .eq. 'sc_cloud') .or. (AvgFunc .eq. 'wr_cloud')) then
    call DoCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, Dens, CintLiq, TsAvg)
  else
    if (AvgFunc .eq. 'precipr') then
      call DoPrecipR(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, PrecipR, CintLiq, TsAvg)
    else
      if (AvgFunc .eq. 'w_up') then
        call DoWup(Nx, Ny, Nz, Nt, DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, W, TsAvg)
      else
        if ((AvgFunc .eq. 'sc_cloud_diam') .or. (AvgFunc .eq. 'wr_cloud_diam')) then
          call DoCloudDiam(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, CloudDiam, CintLiq, TsAvg)
        else
          if (AvgFunc .eq. 'ew_cloud') then
            call DoEwCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CloudConc, TsAvg)
          else
            if (AvgFunc .eq. 'horiz_ke') then
              call DoHorizKe(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, U, V, Dens, TsAvg)
            else
              if (AvgFunc .eq. 'storm_int') then
                call DoStormInt(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, U, V, TsAvg)
              else
                if (AvgFunc .eq. 'ccnconc') then
                  call DoCcnConc(Nx, Ny, Nz, Nt, DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CcnConc, TsAvg)
                else
                  if (AvgFunc .eq. 'test_cvs') then
                    call DoTestCvs(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, TestSelect)
                  end if
                end if
              end if
            end if
          end if
        end if
      end if
    end if
  end if

  ! Convert to rates (derivative wrt time) if requested
  if (DoRates .eq. .true.) then
    write (*,*) 'Converting average data to rates'
    write (*,*) ''

    if (AvgFunc .eq. 'test_cvs') then
      write (*,*) 'WARNING: Skipping conversion due to selection of test function'
      write (*,*) ''
    else
      allocate (Rates(1:1,1:1,1:1,1:Nt-2), stat=Ierror)
      if (Ierror .ne. 0) then
        write (*,*) 'ERROR: Ouput rates data array memory allocation failed'
        stop
      end if

      ! Doing the derivative with a difference method which will reduce
      ! the number of time points by two.
      call CalcRates(1, 1, 1, Nt, DeltaT, TsAvg, Rates)
    end if
  end if

  if (AvgFunc .eq. 'test_cvs') then
    Xstart = GdataDescrip(PressLoc%Fnum)%Xcoords(1)
    Xinc = GdataDescrip(PressLoc%Fnum)%Xcoords(2) - GdataDescrip(PressLoc%Fnum)%Xcoords(1)
    Ystart = GdataDescrip(PressLoc%Fnum)%Ycoords(1)
    Yinc = GdataDescrip(PressLoc%Fnum)%Ycoords(2) - GdataDescrip(PressLoc%Fnum)%Ycoords(1)
    call BuildGoutDescrip(Nx, Ny, Nz, Nt, TestSelect, OfileBases(1), &
           GdataDescrip(PressLoc%Fnum)%UndefVal, AvgFunc, Xstart, Xinc, Ystart, Yinc, &
           GdataDescrip(PressLoc%Fnum)%Zcoords, GdataDescrip(PressLoc%Fnum)%Tstart, &
           GdataDescrip(PressLoc%Fnum)%Tinc, GoutDescrip, 'ts')
    call WriteGrads(GoutDescrip, TestSelect)
  else
    DummyZcoords(1) = 0.0
    Xstart = 0.0
    Xinc = 1.0
    Ystart = 0.0
    Yinc = 1.0
    call BuildGoutDescrip(1, 1, 1, Nt, TsAvg, OfileBases(1), GdataDescrip(1)%UndefVal, &
         AvgFunc, Xstart, Xinc, Ystart, Yinc, DummyZcoords, GdataDescrip(1)%Tstart, &
         GdataDescrip(1)%Tinc, GoutDescrip, 'ts')
    call WriteGrads(GoutDescrip, TsAvg)
    if (DoRates .eq. .true.) then
      call BuildGoutDescrip(1, 1, 1, Nt-2, Rates, OfileBases(2), GdataDescrip(1)%UndefVal, &
           AvgFunc, Xstart, Xinc, Ystart, Yinc, DummyZcoords, GdataDescrip(1)%Tstart, &
           GdataDescrip(1)%Tinc, GoutDescrip, 'dt')
      call WriteGrads(GoutDescrip, Rates)
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
    write (*,*) '            sc_cloud -> total supercooled cloud droplets'
    write (*,*) '            precipr -> total precipitation rate'
    write (*,*) '            w_up -> average w in regions of significant updrafts'
    write (*,*) '            ccnconc -> average ccn concentration'
    write (*,*) '            sc_cloud_diam -> total supercooled cloud droplet mean diameter'
    write (*,*) '            ew_cloud -> average cloud droplet concentration near eyewall region'
    write (*,*) '            wr_cloud -> total warm rain cloud droplets'
    write (*,*) '            wr_cloud_diam -> total warm rain cloud droplet mean diameter'
    write (*,*) '            horiz_ke -> total kinetic energy form horizontal winds'
    write (*,*) '            storm_int -> storm intensity metric from horizontal wind speeds'
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

  if ((AvgFunc .ne. 'sc_cloud') .and. (AvgFunc .ne. 'precipr') .and. &
      (AvgFunc .ne. 'w_up') .and. (AvgFunc .ne. 'sc_cloud_diam') .and. &
      (AvgFunc .ne. 'wr_cloud') .and. (AvgFunc .ne. 'wr_cloud_diam') .and. &
      (AvgFunc .ne. 'ew_cloud') .and. (AvgFunc .ne. 'horiz_ke') .and. &
      (AvgFunc .ne. 'storm_int') .and. (AvgFunc .ne. 'test_cvs') .and. &
      (AvgFunc .ne. 'ccnconc')) then
    write (*,*) 'ERROR: <averaging_function> must be one of:'
    write (*,*) '          sc_cloud'
    write (*,*) '          precipr'
    write (*,*) '          w_up'
    write (*,*) '          ccnconc'
    write (*,*) '          sc_cloud_diam'
    write (*,*) '          ew_cloud'
    write (*,*) '          wr_cloud'
    write (*,*) '          wr_cloud_diam'
    write (*,*) '          horiz_ke'
    write (*,*) '          storm_int'
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

subroutine DoCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, Dens, CintLiq, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  logical :: DoSc
  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: Cloud, TempC, Dens
  real, dimension(1:Nx, 1:Ny, 1:1, 1:Nt) :: CintLiq
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real, dimension(1:Nz) :: Zcoords

  real, dimension(0:Nz) :: ZmHeights
  integer :: ix, iy, iz, it, NumPoints
  logical :: InsideCylVol

  call SetZmHeights (Nz, ZmHeights)

  ! Convert the cloud mixing ratio to mass for each grid point. 
  !
  ! Mixing ratios in GRADS files are g/kg
  ! Density is in kg/m**3
  ! Heights are in m
  ! So, express the mass value in g using the formula
  !   (mix ratio) * (density) * (layer thickness) * (layer horiz area)
  ! The layer thickness for layer k is: ZmHeights(k) - ZmHeights(k-1)

  do it = 1, Nt
    ! Sum up the cloud droplet mass over the specified radial band. Only include the
    ! grid points where tempc is 0 or less (supercooled)

    TsAvg(1,1,1,it) = 0.0
    NumPoints = 0

    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          if ((InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords))  .and. (CintLiq(ix,iy,1,it) .ge. CilThresh)) then
            if (DoSc) then
              if (TempC(ix,iy,iz,it) .le. 0.0) then
                TsAvg(1,1,1,it) = TsAvg(1,1,1,it) + (Cloud(ix,iy,iz,it) * Dens(ix,iy,iz,it) * (ZmHeights(iz)-ZmHeights(iz-1)))
                NumPoints = NumPoints + 1
              end if
            else
              if (TempC(ix,iy,iz,it) .gt. 0.0) then
                TsAvg(1,1,1,it) = TsAvg(1,1,1,it) + (Cloud(ix,iy,iz,it) * Dens(ix,iy,iz,it) * (ZmHeights(iz)-ZmHeights(iz-1)))
                NumPoints = NumPoints + 1
              end if
            end if
          end if
        end do
      end do
    end do

    ! At this point TsAvg(1,1,1,it) holds g/m**2, multiply by grid cell horizontal area. Note this assumes
    ! each grid cell has the same horizontal area.

    TsAvg(1,1,1,it) = TsAvg(1,1,1,it) * DeltaX * DeltaY
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

subroutine DoWup(Nx, Ny, Nz, Nt, DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, W, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: W
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real, dimension(1:Nz) :: Zcoords

  integer ix,iy,iz,it, NumPoints
  logical :: InsideCylVol

  do it = 1, Nt
    ! Average w over regions where significant updrafts occur

    TsAvg(1,1,1,it) = 0.0
    NumPoints = 0

    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          if (InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            if (W(ix,iy,iz,it) .ge. Wthreshold) then
              TsAvg(1,1,1,it) = TsAvg(1,1,1,it) + W(ix,iy,iz,it)
              NumPoints = NumPoints + 1
            end if
          end if
        end do
      end do
    end do

    if (NumPoints .eq. 0) then
      TsAvg(1,1,1,it) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TsAvg(1,1,1,it) = TsAvg(1,1,1,it) / float(NumPoints)
      write (*,*) 'Wup: Timestep:', it, ', Number of points selected: ', NumPoints
    end if
  end do

end subroutine


!************************************************************************************
! DoCloudDiam()
!
! This routine will calculate an averaged cloud droplet diameter. Can select between
! warm rain droplets or supercooled droplets.

subroutine DoCloudDiam(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, CloudDiam, CintLiq, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  logical :: DoSc
  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: Cloud, TempC, CloudDiam
  real, dimension(1:Nx, 1:Ny, 1:1, 1:Nt) :: CintLiq
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real, dimension(1:Nz) :: Zcoords

  integer ix,iy,iz,it, NumPoints
  real SumQ, SumQD
  real MaxQ, Climit, SumD
  logical :: InsideCylVol

  do it = 1, Nt
!     ! Calculate a mass-weighted mean diameter for supercooled cloud droplets.
!     ! 
!     !    Mean diameter (TsAvg value) = Sum(cloud * cloud_d) / Sum(cloud)
!     !    where cloud and cloud_d are only included in the sum when tempc is <= 0
! 
!     SumQD = 0.0
!     SumQ = 0.0
!     NumPoints = 0
! 
!     do ix = 1, Nx
!       do iy = 1, Ny
!         do iz = 1, Nz
!           if (InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!             if (TempC(ix,iy,iz,it) .le. 0.0) then
!                SumQD = SumQD + (Cloud(ix,iy,iz,it) * CloudDiam(ix,iy,iz,it))
!                SumQ = SumQ + Cloud(ix,iy,iz,it)
!                NumPoints = NumPoints + 1
!             end if
!           end if
!         end do
!       end do
!     end do
! 
!     if (SumQ .eq. 0.0) then
!       TsAvg(1,1,1,it) = 0.0
!       write (*,*) 'WARNING: no data points selected for time step: ', it
!     else
!       TsAvg(1,1,1,it) = SumQD / SumQ
!       write (*,*) 'ScCloudDiam: Timestep:', it, ', Number of points selected: ', NumPoints
!     end if

    ! Find the max Q (mass) value and use it to filter data - select data points if
    ! the Q is within 20% of the max Q (.2 to 1.0). Then form the average of the diameters
    ! of the selected points.
!    MaxQ = 0
!    do ix = 1, Nx
!      do iy = 1, Ny
!        do iz = 1, Nz
!          if (InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!            if (Cloud(ix,iy,iz,it) .gt. MaxQ) then
!               MaxQ = Cloud(ix,iy,iz,it)
!            end if
!          end if
!        end do
!      end do
!    end do
!    Climit = 0.3*MaxQ
    
    SumQD = 0.0
    SumQ = 0.0
    SumD = 0.0
    NumPoints = 0
    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          if ((InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) .and. (CintLiq(ix,iy,1,it) .ge. CilThresh)) then
            if (DoSc) then
              if (TempC(ix,iy,iz,it) .le. 0.0) then
                 SumQD = SumQD + (Cloud(ix,iy,iz,it) * CloudDiam(ix,iy,iz,it))
                 SumQ = SumQ + Cloud(ix,iy,iz,it)
                 SumD = SumD + CloudDiam(ix,iy,iz,it)
                 NumPoints = NumPoints + 1
              end if
            else
              if (TempC(ix,iy,iz,it) .gt. 0.0) then
                 SumQD = SumQD + (Cloud(ix,iy,iz,it) * CloudDiam(ix,iy,iz,it))
                 SumQ = SumQ + Cloud(ix,iy,iz,it)
                 SumD = SumD + CloudDiam(ix,iy,iz,it)
                 NumPoints = NumPoints + 1
              end if
            end if
          end if
        end do
      end do
    end do

!    if (NumPoints .eq. 0) then
!      TsAvg(1,1,1,it) = 0.0
!      write (*,*) 'WARNING: no data points selected for time step: ', it
!    else
!      TsAvg(1,1,1,it) = SumD / float(NumPoints)
!      write (*,*) 'ScCloudDiam: Ts:', it, ', NumPoints: ', NumPoints, 'MaxQ: ', MaxQ, 'Climit: ', Climit
!    end if

     if (SumQ .eq. 0.0) then
       TsAvg(1,1,1,it) = 0.0
       write (*,*) 'WARNING: no data points selected for time step: ', it
     else
       TsAvg(1,1,1,it) = SumQD / SumQ
       write (*,*) 'ScCloudDiam: Ts:', it, ', NumPoints: ', NumPoints
     end if

  end do
end subroutine

!************************************************************************************
! DoEwCloud()
!
! This subroutine will do the average cloud droplet concentration near the eyewall
! region.

subroutine DoEwCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CloudConc, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: CloudConc
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real, dimension(1:Nz) :: Zcoords

  integer ix,iy,iz,it
  integer NumPoints
  real SumCloudConc
  logical :: InsideCylVol

  ! Calculate the average cloud droplet concentration near the eyewall region.
  ! 
  ! Call the cloud base to be between 1000m and 1200m. The z level corresponding to
  ! that is z = 4 which is 1138m. Cover from surface to the cloud base level. This
  ! results in using z = 1 to z = 4.

  do it = 1, Nt
    SumCloudConc = 0.0
    NumPoints = 0

    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          if (InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            SumCloudConc = SumCloudConc + CloudConc(ix,iy,iz,it)
            NumPoints = NumPoints + 1
          end if
        end do
      end do
    end do

    if (NumPoints .eq. 0) then
      TsAvg(1,1,1,it) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TsAvg(1,1,1,it) = SumCloudConc / float(NumPoints)
      write (*,*) 'EwCloud: Timestep:', it, ', Number of points selected: ', NumPoints
    end if
  end do
end subroutine

!************************************************************************************
! DoPrecipR()
!
! This subroutine will do the average cloud droplet concentration near the eyewall
! region.

subroutine DoPrecipR(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, PrecipR, CintLiq, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real, dimension(1:Nz) :: Zcoords
  ! PrecipR is 2D var
  real, dimension(1:Nx, 1:Ny, 1:1, 1:Nt) :: PrecipR
  real, dimension(1:Nx, 1:Ny, 1:1, 1:Nt) :: CintLiq
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg

  integer :: ix,iy,iz,it
  integer :: NumPoints
  real :: SumPrecip
  logical :: InsideCylVol

  do it = 1, Nt
    SumPrecip = 0.0
    NumPoints = 0

    do ix = 1, Nx
      do iy = 1, Ny
        if (InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, 1, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords) .and. (CintLiq(ix,iy,1,it) .ge. CilThresh)) then
          SumPrecip = SumPrecip + PrecipR(ix,iy,1,it)
          NumPoints = NumPoints + 1
        end if
      end do
    end do

    if (NumPoints .eq. 0) then
      TsAvg(1,1,1,it) = 0.0
      write (*,*) '  WARNING: no data points selected for this time step'
    else
      ! At this point TsAvg(1,1,1,it) holds mm/hr, multiply by grid cell horizontal area.
      ! Note this assumes each grid cell has the same horizontal area. What this does
      ! is convert mm/hr to kg/hr when assuming that the density of water is 1000kg/m**3.
      !   mm/hr * m**2 * 1000 kg/m**3 * 0.001 m/mm -> kg/hr
      TsAvg(1,1,1,it) = (SumPrecip / float(NumPoints)) * DeltaX * DeltaY
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

subroutine DoHorizKe(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, U, V, Dens, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: U, V, Dens
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real, dimension(1:Nz) :: Zcoords

  integer ix,iy,iz,it
  integer NumPoints
  real SumKe, CurrKe, LevThickness
  logical :: InsideCylVol

  ! KE is 1/2 * m * v^2
  !   - calculate this a every point inside the defined cylindrical volume
  !   - get m by density * volume
  !   - v^2 is based on horiz velocity so is equal to U^2 + V^2
  !
  ! Zcoords are technically the center points of the levels in the RAMS simulation. Since we don't have
  ! the level definition from the RAMS runs here, just use the difference from the i+1st z coord minus the
  ! ith z coord to approzimate the ith level thickness. This will be close enough for the measurement.

  do it = 1, Nt
    SumKe = 0.0
    NumPoints = 0

    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          if (InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            if (iz .eq. Nz) then
              ! Use the level below for this case (since no level above)
              LevThickness = Zcoords(iz) - Zcoords(iz-1)
            else
              LevThickness = Zcoords(iz+1) - Zcoords(iz)
            end if
            CurrKe = 0.5 * DeltaX * DeltaY * LevThickness * Dens(ix,iy,iz,it) * (U(ix,iy,iz,it)**2 + V(ix,iy,iz,it)**2)
            SumKe = SumKe + CurrKe
            NumPoints = NumPoints + 1
          end if
        end do
      end do
    end do

    TsAvg(1,1,1,it) = SumKe;
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

subroutine DoStormInt(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, U, V, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: U, V
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real, dimension(1:Nz) :: Zcoords

  integer ix,iy,iz,it
  integer nCat0, nCat1, nCat2, nCat3, nCat4, nCat5, NumPoints
  real Wspeed, SiMetric
  logical :: InsideCylVol

  do it = 1, Nt
    nCat0 = 0
    nCat1 = 0
    nCat2 = 0
    nCat3 = 0
    nCat4 = 0
    nCat5 = 0

    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          if (InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            ! Count up the number of grid points with wind speeds fitting each of the
            ! Saffir-Simpson categories. Then form the metric by weighting each category
            ! count.
            Wspeed = sqrt(U(ix,iy,iz,it)**2 + V(ix,iy,iz,it)**2)
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
      TsAvg(1,1,1,it) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TsAvg(1,1,1,it) = SiMetric / float(NumPoints)
      write (*,*) 'StormInt: Timestep:', it, ', Number of points selected: ', NumPoints
    end if
  end do
  
end subroutine

!************************************************************************************
! DoCcnConc()
!
! This subroutine will do the average CCN concentration.
!

subroutine DoCcnConc(Nx, Ny, Nz, Nt, DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CcnConc, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: CcnConc
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real, dimension(1:Nz) :: Zcoords

  integer ix,iy,iz,it, NumPoints
  logical :: InsideCylVol

  do it = 1, Nt
    ! Average w over regions where significant updrafts occur

    TsAvg(1,1,1,it) = 0.0
    NumPoints = 0

    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          if (InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            TsAvg(1,1,1,it) = TsAvg(1,1,1,it) + CcnConc(ix,iy,iz,it)
            NumPoints = NumPoints + 1
          end if
        end do
      end do
    end do

    if (NumPoints .eq. 0) then
      TsAvg(1,1,1,it) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TsAvg(1,1,1,it) = TsAvg(1,1,1,it) / float(NumPoints)
      write (*,*) 'Wup: Timestep:', it, ', Number of points selected: ', NumPoints
    end if
  end do

end subroutine

!************************************************************************************
! DoTestCvs()
!
! This subroutine will perform a test on the cylindrical volume selection routine.
! Just runs through the entire grid and outputs a '1' when selection occurs otherwise
! outputs a '0'. Then view the result in grads and see if selection is correct.

subroutine DoTestCvs(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, TestSelect)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: TestSelect
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real, dimension(1:Nz) :: Zcoords

  integer ix,iy,iz,it
  integer NumPoints
  logical :: InsideCylVol

  write (*,*) 'Testing cylindrical volume selection:'
  do it = 1, Nt
    NumPoints = 0
    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          if (InsideCylVol(Nx, Ny, Nz, Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
            TestSelect(ix,iy,iz,it) = 1.0
            NumPoints = NumPoints + 1
          else
            TestSelect(ix,iy,iz,it) = 0.0
          end if
        end do
      end do
    end do
    ! mark the storm center
    TestSelect(StmIx(it), StmIy(it), iz, it) = 2.0
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

subroutine CalcRates(Nx, Ny, Nz, Nt, DeltaT, TsAvg, Rates)
  implicit none

  integer Nx, Ny, Nz, Nt
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: TsAvg 
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt-2) :: Rates 
  real :: DeltaT

  real :: f1, f2
  integer :: ix, iy, iz, it

  ! use a centered difference, uses points at t-1, t and t+1
  do ix = 1, Nx
    do iy = 1, Ny
      do iz = 1, Nz
        do it = 2, Nt-1
          f1 = (TsAvg(1,1,1,it) + TsAvg(1,1,1,it-1)) / 2.0
          f2 = (TsAvg(1,1,1,it+1) + TsAvg(1,1,1,it)) / 2.0

          Rates(1,1,1,it-1) = (f2 - f1) / DeltaT
        end do
      end do
    end do
  end do
end subroutine
