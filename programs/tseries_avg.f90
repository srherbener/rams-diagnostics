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
  integer :: Nfiles

  real :: DeltaX, DeltaY, ScThreshold
  real, dimension(:), allocatable :: ZmHeights
  real, dimension(1:1) :: DummyZcoords

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars
  type (GradsOutDescription) :: GoutDescrip

  ! Data arrays: cloud, tempc, precipr
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of cloud, tempc, precipr in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: Cloud, TempC, Dens, PrecipR, W, CloudDiam, CconcAz, TsAvg
  type (GradsVarLocation) :: CloudLoc, TempcLoc, DensLoc, WLoc, CloudDiamLoc, CconcAzLoc, PreciprLoc

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBase, AvgFunc, DeltaX, DeltaY, ScThreshold)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Time seris of average for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  Averaging function: ', trim(AvgFunc)
  write (*,*) '  Grid delta x: ', DeltaX
  write (*,*) '  Grid delta y: ', DeltaY
  write (*,*) '  Supercooled cloud droplet threshold (for sc_w function): ', ScThreshold
  write (*,*) ''

  ! Read the GRADS data description files and collect the information about the data
  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''

  ! Check the data description for consistency and locate the variables in the GRADS control files
  if (AvgFunc .eq. 'sc_cloud') then
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, TempcLoc, 'tempc')
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, DensLoc, 'dn0')
    call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudLoc, 'cloud')
  else
    if (AvgFunc .eq. 'precipr') then
      call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, PreciprLoc, 'precipr')
    else
      if (AvgFunc .eq. 'sc_w') then
        call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, TempcLoc, 'tempc')
        call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudLoc, 'cloud')
        call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, WLoc, 'w')
      else
        if (AvgFunc .eq. 'sc_cloud_diam') then
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, TempcLoc, 'tempc')
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudLoc, 'cloud')
          call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CloudDiamLoc, 'cloud_d')
        else
          if (AvgFunc .eq. 'ew_cloud') then
            call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, CconcAzLoc, 'cloudconcen_cm3_azavg')
          end if
        end if
      end if
    end if
  end if

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

  write (*,*) 'Locations of variables in GRADS data (file number, var number):'
  if (AvgFunc .eq. 'sc_cloud') then
    write (*,'(a20,i3,a2,i3,a1)') 'tempc: (', TempcLoc%Fnum, ', ', TempcLoc%Vnum, ')'
    write (*,'(a20,i3,a2,i3,a1)') 'dn0: (', DensLoc%Fnum, ', ', DensLoc%Vnum, ')'
    write (*,'(a20,i3,a2,i3,a1)') 'cloud: (', CloudLoc%Fnum, ', ', CloudLoc%Vnum, ')'

    ! Allocate the data arrays and read in the data from the GRADS data files
    allocate (TempC(1:Nx,1:Ny,1:Nz,1:Nt), Dens(1:Nx,1:Ny,1:Nz,1:Nt), &
              Cloud(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory allocation failed'
      stop
    end if

    ! Read in the data for the vars using the description and location information
    call ReadGradsData(GdataDescrip, 'tempc', TempcLoc, TempC, Nx, Ny, Nz, Nt)
    call ReadGradsData(GdataDescrip, 'dn0', DensLoc, Dens, Nx, Ny, Nz, Nt)
    call ReadGradsData(GdataDescrip, 'cloud', CloudLoc, Cloud, Nx, Ny, Nz, Nt)
  else
    if (AvgFunc .eq. 'precipr') then
      write (*,'(a20,i3,a2,i3,a1)') 'precipr: (', PreciprLoc%Fnum, ', ', PreciprLoc%Vnum, ')'

      ! Allocate the data arrays and read in the data from the GRADS data files
      allocate (PrecipR(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
      if (Ierror .ne. 0) then
        write (*,*) 'ERROR: Data array memory allocation failed'
        stop
      end if

      ! Read in the data for the vars using the description and location information
      call ReadGradsData(GdataDescrip, 'precipr', PreciprLoc, PrecipR, Nx, Ny, Nz, Nt)
    else
      if (AvgFunc .eq. 'sc_w') then
        write (*,'(a20,i3,a2,i3,a1)') 'tempc: (', TempcLoc%Fnum, ', ', TempcLoc%Vnum, ')'
        write (*,'(a20,i3,a2,i3,a1)') 'cloud: (', CloudLoc%Fnum, ', ', CloudLoc%Vnum, ')'
        write (*,'(a20,i3,a2,i3,a1)') 'w: (', WLoc%Fnum, ', ', WLoc%Vnum, ')'

        ! Allocate the data arrays and read in the data from the GRADS data files
        allocate (TempC(1:Nx,1:Ny,1:Nz,1:Nt), W(1:Nx,1:Ny,1:Nz,1:Nt), &
                  Cloud(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
        if (Ierror .ne. 0) then
          write (*,*) 'ERROR: Data array memory allocation failed'
          stop
        end if

        ! Read in the data for the vars using the description and location information
        call ReadGradsData(GdataDescrip, 'tempc', TempcLoc, TempC, Nx, Ny, Nz, Nt)
        call ReadGradsData(GdataDescrip, 'w', WLoc, W, Nx, Ny, Nz, Nt)
        call ReadGradsData(GdataDescrip, 'cloud', CloudLoc, Cloud, Nx, Ny, Nz, Nt)
      else
        if (AvgFunc .eq. 'sc_cloud_diam') then
          write (*,'(a20,i3,a2,i3,a1)') 'tempc: (', TempcLoc%Fnum, ', ', TempcLoc%Vnum, ')'
          write (*,'(a20,i3,a2,i3,a1)') 'cloud: (', CloudLoc%Fnum, ', ', CloudLoc%Vnum, ')'
          write (*,'(a20,i3,a2,i3,a1)') 'cloud_d: (', CloudDiamLoc%Fnum, ', ', CloudDiamLoc%Vnum, ')'
  
          ! Allocate the data arrays and read in the data from the GRADS data files
          allocate (TempC(1:Nx,1:Ny,1:Nz,1:Nt), CloudDiam(1:Nx,1:Ny,1:Nz,1:Nt), &
                    Cloud(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
          if (Ierror .ne. 0) then
            write (*,*) 'ERROR: Data array memory allocation failed'
            stop
          end if

          ! Read in the data for the vars using the description and location information
          call ReadGradsData(GdataDescrip, 'tempc', TempcLoc, TempC, Nx, Ny, Nz, Nt)
          call ReadGradsData(GdataDescrip, 'cloud', CloudLoc, Cloud, Nx, Ny, Nz, Nt)
          call ReadGradsData(GdataDescrip, 'cloud_d', CloudDiamLoc, CloudDiam, Nx, Ny, Nz, Nt)
        else
          if (AvgFunc .eq. 'ew_cloud') then
            write (*,'(a20,i3,a2,i3,a1)') 'cloudconcen_cm3: (', CconcAzLoc%Fnum, ', ', CconcAzLoc%Vnum, ')'
    
            ! Allocate the data arrays and read in the data from the GRADS data files
            allocate (CconcAz(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
            if (Ierror .ne. 0) then
              write (*,*) 'ERROR: Data array memory allocation failed'
              stop
            end if
  
            ! Read in the data for the vars using the description and location information
            call ReadGradsData(GdataDescrip, 'cloudconcen_cm3', CconcAzLoc, CconcAz, Nx, Ny, Nz, Nt)
          end if
        end if
      end if
    end if
  end if
  write (*,*) ''

  ! Allocate the output array and compute the azimuthal averaging
  allocate (TsAvg(1:1,1:1,1:1,1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Ouput data array memory allocation failed'
    stop
  end if

  ! call the averaging function

  if (AvgFunc .eq. 'sc_cloud') then
    call DoScCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, Cloud, TempC, Dens, TsAvg)
  else
    if (AvgFunc .eq. 'precipr') then
      call DoPrecipR(Nx, Ny, Nz, Nt, DeltaX, DeltaY, PrecipR, TsAvg)
    else
      if (AvgFunc .eq. 'sc_w') then
        call DoScW(Nx, Ny, Nz, Nt, DeltaX, DeltaY, ScThreshold, Cloud, TempC, W, TsAvg)
      else
        if (AvgFunc .eq. 'sc_cloud_diam') then
          call DoScCloudDiam(Nx, Ny, Nz, Nt, DeltaX, DeltaY, Cloud, TempC, CloudDiam, TsAvg)
        else
          if (AvgFunc .eq. 'ew_cloud') then
            call DoEwCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, CconcAz, TsAvg)
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

!  ! Clean up
!  deallocate (U, V, W, Press, Avar, StmIx, StmIy, MinP, AzAvg, stat=Ierror)
!  if (Ierror .ne. 0) then
!    write (*,*) 'ERROR: Data array memory de-allocation failed'
!    stop
!  end if

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

subroutine GetMyArgs(Infiles, OfileBase, AvgFunc, DeltaX, DeltaY, ScThreshold)
  implicit none

  integer, parameter :: MAX_ITEMS = 3

  character (len=*) :: Infiles, OfileBase, AvgFunc
  real :: DeltaX, DeltaY, ScThreshold

  integer :: iargc
  character (len=128) :: arg
  character (len=128), dimension(1:MAX_ITEMS) :: ArgList
  integer :: Nitems

  if (iargc() .ne. 5) then
    write (*,*) 'ERROR: must supply exactly 5 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_file> <averaging_function>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <averaging_function>: averaging function to use on input data'
    write (*,*) '            sc_cloud -> total supercooled cloud droplets'
    write (*,*) '            precipr -> total precipitation rate'
    write (*,*) '            sc_w:sc_threshold -> average vertical velocity in regions of significant'
    write (*,*) '                                 supercooled cloud droplets, sc_threshold defines limit'
    write (*,*) '                                 for selecting w - if supercooled cloud mixing ratio is'
    write (*,*) '                                 >= sc_threshold, then include w in average calculation'
    write (*,*) '            sc_cloud_diam -> total supercooled cloud droplet mean diameter'
    write (*,*) '            ew_cloud -> average cloud droplet concentration near eyewall region'
    write (*,*) '        <delta_x>: delta X of grid'
    write (*,*) '        <delta_y>: delta Y of grid'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)

  call getarg(3, arg)
  call String2List(arg, ':', ArgList, MAX_ITEMS, Nitems, 'sc_w arguments')
  AvgFunc = ArgList(1)
  if (AvgFunc .eq. 'sc_w') then
    read(ArgList(2), '(f)') ScThreshold
  else
    ScThreshold = 0.0
  end if

  call getarg(4, arg)
  read (arg, '(f)') DeltaX

  call getarg(5, arg)
  read (arg, '(f)') DeltaY


  if ((AvgFunc .ne. 'sc_cloud') .and. (AvgFunc .ne. 'precipr') .and. &
      (AvgFunc .ne. 'sc_w') .and. (AvgFunc .ne. 'sc_cloud_diam') .and. &
      (AvgFunc .ne. 'ew_cloud')) then
    write (*,*) 'ERROR: <averaging_function> must be one of:'
    write (*,*) '          sc_cloud'
    write (*,*) '          precipr'
    write (*,*) '          sc_w'
    write (*,*) '          sc_cloud_diam'
    write (*,*) '          ew_cloud'
  end if

  return
end subroutine

!*****************************************************************************
! DoScCloud()
!
! This subroutine will perform the supercooled cloud droplet time series
! averaging.
!

subroutine DoScCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, Cloud, TempC, Dens, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: Cloud, TempC, Dens
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg

  real, dimension(0:Nz) :: ZmHeights
  integer :: ix, iy, iz, it

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
    ! Sum up the cloud droplet mass over the whole domain. Only include the
    ! grid points where tempc is 0 or less (supercooled)

    TsAvg(1,1,1,it) = 0.0

    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          if (TempC(ix,iy,iz,it) .le. 0.0) then
            TsAvg(1,1,1,it) = TsAvg(1,1,1,it) + (Cloud(ix,iy,iz,it) * Dens(ix,iy,iz,it) * (ZmHeights(iz)-ZmHeights(iz-1)))
          end if
        end do
      end do
    end do

    ! At this point TsAvg(1,1,1,it) holds g/m**2, multiply by grid cell horizontal area. Note this assumes
    ! each grid cell has the same horizontal area.

    TsAvg(1,1,1,it) = TsAvg(1,1,1,it) * DeltaX * DeltaY
  end do
end subroutine

!************************************************************************************
! DoPrecipR()
!
! This subroutine will sum up the total surface precipitation.
!

subroutine DoPrecipR(Nx, Ny, Nz, Nt, DeltaX, DeltaY, PrecipR, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: PrecipR
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg

  integer :: ix,iy,iz,it

  do it = 1, Nt
    ! Sum up the precipitation rate over the whole domain.

    TsAvg(1,1,1,it) = 0.0

    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          TsAvg(1,1,1,it) = TsAvg(1,1,1,it) + PrecipR(ix,iy,iz,it)
        end do
      end do
    end do

    ! At this point TsAvg(1,1,1,it) holds mm/hr, multiply by grid cell horizontal area.
    ! Note this assumes each grid cell has the same horizontal area. What this does
    ! is convert mm/hr to kg/hr when assuming that the density of water is 1000kg/m**3.
    !   mm/hr * m**2 * 1000 kg/m**3 * 0.001 m/mm -> kg/hr

    TsAvg(1,1,1,it) = TsAvg(1,1,1,it) * DeltaX * DeltaY
  end do
end subroutine

!************************************************************************************
! DoScW()
!
! This subroutine will do the average vertical velocity in significant supercooled
! cloud droplet regions.
!

subroutine DoScW(Nx, Ny, Nz, Nt, DeltaX, DeltaY, ScThreshold, Cloud, TempC, W, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY, ScThreshold
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: Cloud, TempC, W
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg

  integer ix,iy,iz,it, NumPoints

  do it = 1, Nt
    ! Average w over regions where significant amounts supercooled cloud droplets exists

    TsAvg(1,1,1,it) = 0.0
    NumPoints = 0

    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          if (TempC(ix,iy,iz,it) .le. 0.0) then
            if (Cloud(ix,iy,iz,it) .ge. ScThreshold) then
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
    end if
  end do

end subroutine


!************************************************************************************
! DoScCloudDiam()
!
! This subroutine will do the average vertical velocity in significant supercooled
! cloud droplet regions.

subroutine DoScCloudDiam(Nx, Ny, Nz, Nt, DeltaX, DeltaY, Cloud, TempC, CloudDiam, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: Cloud, TempC, CloudDiam
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg

  integer ix,iy,iz,it
  real SumQ, SumQD

  do it = 1, Nt
    ! Calculate a mass-weighted mean diameter for supercooled cloud droplets.
    ! 
    !    Mean diameter (TsAvg value) = Sum(cloud * cloud_d) / Sum(cloud)
    !    where cloud and cloud_d are only included in the sum when tempc is <= 0

    SumQD = 0.0
    SumQ = 0.0

    do ix = 1, Nx
      do iy = 1, Ny
        do iz = 1, Nz
          if (TempC(ix,iy,iz,it) .le. 0.0) then
             SumQD = SumQD + (Cloud(ix,iy,iz,it) * CloudDiam(ix,iy,iz,it))
             SumQ = SumQ + Cloud(ix,iy,iz,it)
          end if
        end do
      end do
    end do

    if (SumQ .eq. 0.0) then
      TsAvg(1,1,1,it) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TsAvg(1,1,1,it) = SumQD / SumQ
    end if
  end do
end subroutine

!************************************************************************************
! DoEwCloud()
!
! This subroutine will do the average cloud droplet concentration near the eyewall
! region.

subroutine DoEwCloud(Nx, Ny, Nz, Nt, DeltaX, DeltaY, CconcAz, TsAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY
  real, dimension(1:Nx, 1:Ny, 1:Nz, 1:Nt) :: CconcAz
  real, dimension(1:1, 1:1, 1:1, 1:Nt) :: TsAvg

  integer ix,iy,iz,it
  integer ixStart, ixEnd, izStart, izEnd
  integer NumPoints
  real SumCloudConc

  ! Calculate the average cloud droplet concentration near the eyewall region.
  ! 
  ! The data is already an azimuthal average at different radii from the storm center.
  !
  ! Call "near the eyewall" roughly from radius 10km to radius 30km. In the data,
  ! the x axis is the radius with x = 1 being 0km, and the distance between each
  ! x value being 4.15km. Using x = 4 to x = 8 gives roughly what we want (these
  ! correspond to radius = 12.45km to radius = 29.05km).
  !
  ! Call the cloud base to be between 1000m and 1200m. The z level corresponding to
  ! that is z = 4 which is 1138m. Cover from surface to the cloud base level. This
  ! results in using z = 1 to z = 4.

  ixStart = 6
  ixEnd = 8

  izStart = 1
  izEnd = 4

  do it = 1, Nt
    SumCloudConc = 0.0
    NumPoints = 0

    do ix = ixStart, ixEnd
      do iy = 1, Ny
        do iz = izStart, izEnd
           SumCloudConc = SumCloudConc + CconcAz(ix,iy,iz,it)
           NumPoints = NumPoints + 1
        end do
      end do
    end do

    if (NumPoints .eq. 0) then
      TsAvg(1,1,1,it) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TsAvg(1,1,1,it) = SumCloudConc / float(NumPoints)
    end if
  end do
end subroutine
