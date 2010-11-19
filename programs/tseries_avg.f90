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

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: AvgFunc

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles
  integer :: Nfiles, MinLevel, MaxLevel
  integer, dimension(:), allocatable :: StmIx, StmIy

  real :: DeltaX, DeltaY, Wthreshold, MinRadius, MaxRadius
  real, dimension(:), allocatable :: ZmHeights
  real, dimension(1:1) :: DummyZcoords
  real, dimension(:), allocatable :: MinP, Xcoords, Ycoords

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars
  type (GradsOutDescription) :: GoutDescrip

  ! Data arrays: cloud, tempc, precipr
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of cloud, tempc, precipr in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: Cloud, TempC, Dens, PrecipR, W, CloudDiam, CloudConc, Press, CintLiq, TsAvg
  type (GradsVarLocation) :: CloudLoc, TempcLoc, DensLoc, WLoc, CloudDiamLoc, CloudConcLoc, PreciprLoc, PressLoc, CintLiqLoc

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBase, AvgFunc, Wthreshold, MinLevel, MaxLevel, MinRadius, MaxRadius)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Time seris of average for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  Averaging function: ', trim(AvgFunc)
  write (*,*) '  Minimum radius: ', MinRadius
  write (*,*) '  Maximum radius: ', MaxRadius
  write (*,*) '  Minimum level (for functions using level selection): ', MinLevel
  write (*,*) '  Maximum level (for functions using level selection): ', MaxLevel
  write (*,*) '  W threshold (for w_up function): ', Wthreshold
  write (*,*) ''

  ! Read the GRADS data description files and collect the information about the data
  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''
  flush(6)

  ! Check the data description for consistency and locate the variables in the GRADS control files
  if (AvgFunc .eq. 'sc_cloud') then
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, TempcLoc, 'tempc')
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, DensLoc, 'dn0')
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudLoc, 'cloud')
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CintLiqLoc, 'liquid_colint')
  else
    if (AvgFunc .eq. 'precipr') then
      call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PreciprLoc, 'precipr')
      call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
    else
      if (AvgFunc .eq. 'w_up') then
        call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, WLoc, 'w')
        call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
      else
        if (AvgFunc .eq. 'sc_cloud_diam') then
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, TempcLoc, 'tempc')
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudLoc, 'cloud')
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudDiamLoc, 'cloud_d')
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CintLiqLoc, 'liquid_colint')
        else
          if (AvgFunc .eq. 'ew_cloud') then
            call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudConcLoc, 'cloud_cm3')
            call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PressLoc, 'press')
          end if
        end if
      end if
    end if
  end if

  ! Calculate the x,y coordinates (in km) for doing selection by radius from storm center.
  ! Also for setting DeltaX and DeltaY (in m)
  allocate (Xcoords(1:Nx), Ycoords(1:Ny), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  call ConvertGridCoords(Nx, Ny, GdataDescrip(1), Xcoords, Ycoords)

  DeltaX = (Xcoords(2) - Xcoords(1)) * 1000.0
  DeltaY = (Ycoords(2) - Ycoords(1)) * 1000.0

  write (*,*) 'Horizontal grid info:'
  write (*,*) '  X range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', GdataDescrip(1)%Xcoords(1), GdataDescrip(1)%Xcoords(Nx), Xcoords(1), Xcoords(Nx)
  write (*,*) '  Y range (min lat, max lat) --> (min y, max y): '
  write (*,*) '    ', GdataDescrip(1)%Ycoords(1), GdataDescrip(1)%Ycoords(Ny), Ycoords(1), Ycoords(Ny)

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
  write (*,*) ''

  write (*,*) 'Locations of variables in GRADS data (file number, var number):'
  if (AvgFunc .eq. 'sc_cloud') then
    write (*,'(a20,i3,a2,i3,a1)') 'tempc: (', TempcLoc%Fnum, ', ', TempcLoc%Vnum, ')'
    write (*,'(a20,i3,a2,i3,a1)') 'dn0: (', DensLoc%Fnum, ', ', DensLoc%Vnum, ')'
    write (*,'(a20,i3,a2,i3,a1)') 'cloud: (', CloudLoc%Fnum, ', ', CloudLoc%Vnum, ')'
    write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'

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

      ! Allocate the data arrays and read in the data from the GRADS data files
      ! precipr is 2D variable -> Nz = 1
      allocate (PrecipR(1:Nx,1:Ny,1:1,1:Nt), Press(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
      if (Ierror .ne. 0) then
        write (*,*) 'ERROR: Data array memory allocation failed'
        stop
      end if

      ! Read in the data for the vars using the description and location information
      ! precipr is 2D variable -> Nz = 1
      call ReadGradsData(GdataDescrip, 'precipr', PreciprLoc, PrecipR, Nx, Ny, 1, Nt)
      call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
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
        if (AvgFunc .eq. 'sc_cloud_diam') then
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
          end if
        end if
      end if
    end if
  end if
  write (*,*) ''
  flush(6)

  ! Allocate the output array and do the averaging
  allocate (TsAvg(1:1,1:1,1:1,1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Ouput data array memory allocation failed'
    stop
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

  if (AvgFunc .eq. 'sc_cloud') then
    call DoScCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords, Cloud, TempC, Dens, CintLiq, TsAvg)
  else
    if (AvgFunc .eq. 'precipr') then
      call DoPrecipR(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords, PrecipR, TsAvg)
    else
      if (AvgFunc .eq. 'w_up') then
        call DoWup(Nx, Ny, Nz, Nt, DeltaX, DeltaY, Wthreshold, MinLevel, MaxLevel, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords, W, TsAvg)
      else
        if (AvgFunc .eq. 'sc_cloud_diam') then
          call DoScCloudDiam(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords, Cloud, TempC, CloudDiam, CintLiq, TsAvg)
        else
          if (AvgFunc .eq. 'ew_cloud') then
            call DoEwCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinLevel, MaxLevel, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords, CloudConc, TsAvg)
          end if
        end if
      end if
    end if
  end if

  DummyZcoords(1) = 0.0
  call BuildGoutDescrip(1, 1, 1, Nt, TsAvg, OfileBase, GdataDescrip(1)%UndefVal, AvgFunc, &
          0.0, 1.0, 0.0, 1.0, DummyZcoords, GdataDescrip(1)%Tstart, &
          GdataDescrip(1)%Tinc, GoutDescrip, 'test')

  call WriteGrads(GoutDescrip, TsAvg)

  stop
end

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   OfileBase - output GRADS file, base name for two files
!   AvgFunc - averaging function selection
!

subroutine GetMyArgs(Infiles, OfileBase, AvgFunc, Wthreshold, MinLevel, MaxLevel, MinRadius, MaxRadius)
  implicit none

  integer, parameter :: MAX_ITEMS = 5

  character (len=*) :: Infiles, OfileBase, AvgFunc
  real :: Wthreshold, MinRadius, MaxRadius
  integer :: MinLevel, MaxLevel

  integer :: iargc
  character (len=128) :: arg
  character (len=128), dimension(1:MAX_ITEMS) :: ArgList
  integer :: Nitems

  logical :: BadArgs

  if (iargc() .ne. 5) then
    write (*,*) 'ERROR: must supply exactly 5 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_file> <averaging_function>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <averaging_function>: averaging function to use on input data'
    write (*,*) '            sc_cloud -> total supercooled cloud droplets'
    write (*,*) '            precipr -> total precipitation rate'
    write (*,*) '            w_up:<min_level>:<max_level>:<w_threshold> -> average w in regions of significant updrafts'
    write (*,*) '            sc_cloud_diam -> total supercooled cloud droplet mean diameter'
    write (*,*) '            ew_cloud:<min_level>:<max_level> -> average cloud droplet concentration near eyewall region'
    write (*,*) '            Parameters:'
    write (*,*) '              <min_level>:<max_level> are for selecting the z levels. Use integers'
    write (*,*) '                (ie, the k values for the levels)'
    write (*,*) '              <w_threshold> is value to filter data, ie use data when w >= <w_threshold>'
    write (*,*) '        <min_radius> -> for selecting radial bands, this defines inner boundary'
    write (*,*) '        <max_radius> -> for selecting radial bands, this defines outer boundary'
    write (*,*) '                        NOTE: for radial band selection, you must supply the GRADS'
    write (*,*) '                        data file for pressure. The location of minimum pressure is'
    write (*,*) '                        being used to identify the storm center.'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)

  call getarg(3, arg)
  call String2List(arg, ':', ArgList, MAX_ITEMS, Nitems, 'avg function arguments')
  AvgFunc = ArgList(1)
  if (AvgFunc .eq. 'w_up') then
    read(ArgList(2), '(i)') MinLevel
    read(ArgList(3), '(i)') MaxLevel
    read(ArgList(4), '(f)') Wthreshold
  else
    if (AvgFunc .eq. 'ew_cloud') then
      read(ArgList(2), '(i)') MinLevel
      read(ArgList(3), '(i)') MaxLevel
      Wthreshold = 0.0
    else
      MinLevel = 0
      MaxLevel = 0
      Wthreshold = 0.0
    end if
  end if

  call getarg(4, arg)
  read(arg, '(f)') MinRadius

  call getarg(5, arg)
  read(arg, '(f)') MaxRadius

  BadArgs = .false.

  if ((AvgFunc .ne. 'sc_cloud') .and. (AvgFunc .ne. 'precipr') .and. &
      (AvgFunc .ne. 'w_up') .and. (AvgFunc .ne. 'sc_cloud_diam') .and. &
      (AvgFunc .ne. 'ew_cloud')) then
    write (*,*) 'ERROR: <averaging_function> must be one of:'
    write (*,*) '          sc_cloud'
    write (*,*) '          precipr'
    write (*,*) '          w_up'
    write (*,*) '          sc_cloud_diam'
    write (*,*) '          ew_cloud'
    write (*,*) ''
    BadArgs = .true.
  end if

  if ((MinRadius .lt. 0.0) .or. (MaxRadius .lt. 0.0) .or. (MaxRadius .le. MinRadius)) then
    write (*,*) 'ERROR: <min_radius> and <max_radius> must be >= 0.0, and <max_radius> must be > <min_radius>'
    write (*,*) ''
    BadArgs = .true.
  end if

  if (MaxLevel < MinLevel) then
    write (*,*) 'ERROR: <max_level> >=  <min_level>'
    write (*,*) ''
  end if

  if (BadArgs) then
    stop
  end if

  return
end subroutine

!*****************************************************************************
! InsideRadialBand()
!
! This function will compare a given grid location with the storm center
! and a given radial band from that storm center. It will return true if the
! given location falls inside the given radial band.
!

logical function InsideRadialBand(Nx, Ny, Nt, Ix, Iy, It, Rmin, Rmax, StmIx, StmIy, Xcoords, Ycoords)
  implicit none

  integer :: Nx, Ny, Nt, Ix, Iy, It
  real :: Rmin, Rmax
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords

  real :: dX, dY, Radius

  dX = Xcoords(Ix) - Xcoords(StmIx(It))
  dY = Ycoords(Iy) - Ycoords(StmIy(It))
  Radius = sqrt(dX*dX + dY*dY)

  if ((Radius .ge. Rmin) .and. (Radius .le. Rmax)) then
    InsideRadialBand = .true.
  else
    InsideRadialBand = .false.
  end if
  return
end function

!*****************************************************************************
! DoScCloud()
!
! This subroutine will perform the supercooled cloud droplet time series
! averaging.
!

subroutine DoScCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords, Cloud, TempC, Dens, CintLiq, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY, MinRadius, MaxRadius
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: Cloud, TempC, Dens
  real, dimension(1:Nx, 1:Ny, 1:1, 1:Nt) :: CintLiq
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords

  real, dimension(0:Nz) :: ZmHeights
  integer :: ix, iy, iz, it, NumPoints
  logical :: InsideRadialBand

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
        if ((InsideRadialBand(Nx, Ny, Nt, ix, iy, it, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords))  .and. (CintLiq(ix,iy,1,it) .ge. 5000.0)) then
          do iz = 1, Nz
            if (TempC(ix,iy,iz,it) .le. 0.0) then
              TsAvg(1,1,1,it) = TsAvg(1,1,1,it) + (Cloud(ix,iy,iz,it) * Dens(ix,iy,iz,it) * (ZmHeights(iz)-ZmHeights(iz-1)))
              NumPoints = NumPoints + 1
            end if
          end do
        end if
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

subroutine DoWup(Nx, Ny, Nz, Nt, DeltaX, DeltaY, Wthreshold, MinLevel, MaxLevel, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords, W, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt, MinLevel, MaxLevel
  real :: DeltaX, DeltaY, Wthreshold, MinRadius, MaxRadius
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: W
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords

  integer ix,iy,iz,it, NumPoints
  logical :: InsideRadialBand

  do it = 1, Nt
    ! Average w over regions where significant updrafts occur

    TsAvg(1,1,1,it) = 0.0
    NumPoints = 0

    do ix = 1, Nx
      do iy = 1, Ny
        if (InsideRadialBand(Nx, Ny, Nt, ix, iy, it, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords)) then
          do iz = MinLevel, MaxLevel
            if (W(ix,iy,iz,it) .ge. Wthreshold) then
              TsAvg(1,1,1,it) = TsAvg(1,1,1,it) + W(ix,iy,iz,it)
              NumPoints = NumPoints + 1
            end if
          end do
        end if
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
! DoScCloudDiam()
!
! This subroutine will do the average vertical velocity in significant supercooled
! cloud droplet regions.

subroutine DoScCloudDiam(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords, Cloud, TempC, CloudDiam, CintLiq, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY, MinRadius, MaxRadius
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: Cloud, TempC, CloudDiam
  real, dimension(1:Nx, 1:Ny, 1:1, 1:Nt) :: CintLiq
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords

  integer ix,iy,iz,it, NumPoints
  real SumQ, SumQD
  real MaxQ, Climit, SumD
  logical :: InsideRadialBand

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
!         if (InsideRadialBand(Nx, Ny, Nt, ix, iy, it, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords)) then
!           do iz = 1, Nz
!             if (TempC(ix,iy,iz,it) .le. 0.0) then
!                SumQD = SumQD + (Cloud(ix,iy,iz,it) * CloudDiam(ix,iy,iz,it))
!                SumQ = SumQ + Cloud(ix,iy,iz,it)
!                NumPoints = NumPoints + 1
!             end if
!           end do
!         end if
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
    MaxQ = 0
    do ix = 1, Nx
      do iy = 1, Ny
        if (InsideRadialBand(Nx, Ny, Nt, ix, iy, it, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords)) then
          do iz = 1, Nz
            if (Cloud(ix,iy,iz,it) .gt. MaxQ) then
               MaxQ = Cloud(ix,iy,iz,it)
            end if
          end do
        end if
      end do
    end do
    Climit = 0.3*MaxQ
    
    SumQD = 0.0
    SumQ = 0.0
    SumD = 0.0
    NumPoints = 0
    do ix = 1, Nx
      do iy = 1, Ny
!        if ((InsideRadialBand(Nx, Ny, Nt, ix, iy, it, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords)) .and. (CintLiq(ix,iy,1,it) .ge. 5000.0)) then
        if (InsideRadialBand(Nx, Ny, Nt, ix, iy, it, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords)) then
          do iz = 1, Nz
!            if ((TempC(ix,iy,iz,it) .le. 0.0) .and. (Cloud(ix,iy,iz,it) .ge. Climit)) then
            if (TempC(ix,iy,iz,it) .le. 0.0) then
               SumQD = SumQD + (Cloud(ix,iy,iz,it) * CloudDiam(ix,iy,iz,it))
               SumQ = SumQ + Cloud(ix,iy,iz,it)
               SumD = SumD + CloudDiam(ix,iy,iz,it)
               NumPoints = NumPoints + 1
            end if
          end do
        end if
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
       write (*,*) 'ScCloudDiam: Ts:', it, ', NumPoints: ', NumPoints, 'MaxQ: ', MaxQ, 'Climit: ', Climit
     end if

  end do
end subroutine

!************************************************************************************
! DoEwCloud()
!
! This subroutine will do the average cloud droplet concentration near the eyewall
! region.

subroutine DoEwCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinLevel, MaxLevel, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords, CloudConc, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt, MinLevel, MaxLevel
  real :: DeltaX, DeltaY, MinRadius, MaxRadius
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: CloudConc
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords

  integer ix,iy,iz,it
  integer NumPoints
  real SumCloudConc
  logical :: InsideRadialBand

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
        if (InsideRadialBand(Nx, Ny, Nt, ix, iy, it, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords)) then
          do iz = MinLevel, MaxLevel
             SumCloudConc = SumCloudConc + CloudConc(ix,iy,iz,it)
          end do
          NumPoints = NumPoints + 1
        end if
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

subroutine DoPrecipR(Nx, Ny, Nz, Nt, DeltaX, DeltaY, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords, PrecipR, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY, MinRadius, MaxRadius
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  ! PrecipR is 2D var
  real, dimension(1:Nx, 1:Ny, 1:1, 1:Nt) :: PrecipR
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg

  integer :: ix,iy,iz,it
  integer :: NumPoints
  real :: SumPrecip
  logical :: InsideRadialBand

  do it = 1, Nt
    SumPrecip = 0.0
    NumPoints = 0

    do ix = 1, Nx
      do iy = 1, Ny
        if (InsideRadialBand(Nx, Ny, Nt, ix, iy, it, MinRadius, MaxRadius, StmIx, StmIy, Xcoords, Ycoords)) then
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
