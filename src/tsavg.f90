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

  real, parameter :: UndefVal = -999.0

  character (len=MediumString) :: InDir
  character (len=MediumString) :: InSuffix
  character (len=MediumString) :: OutFile
  character (len=MediumString) :: FilterFile
  character (len=LittleString) :: AvgFunc

  ! Data arrays
  ! Dims: x, y, z, t
  type (Rhdf5Var) :: U, V, AzWind, Dens, Filter, TserAvg
  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords
  type (Rhdf5Var) :: VarLon, VarLat, VarZcoords
  character (len=MediumString) :: Ufile, Vfile, AzWindFile, DensFile
  character (len=LittleString) :: Units

  integer :: ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt
  real :: DeltaX, DeltaY
  real, dimension(:), allocatable :: VarXcoords, VarYcoords

  ! Get the command line arguments
  call GetMyArgs(InDir, InSuffix, OutFile, AvgFunc, FilterFile)

  write (*,*) 'Time seris of average for RAMS data:'
  write (*,*) '  Input directory: ', trim(InDir)
  write (*,*) '  Input file suffix: ', trim(InSuffix)
  write (*,*) '  Output file:  ', trim(OutFile)
  write (*,*) '  Averaging function: ', trim(AvgFunc)
  write (*,*) '  Filter file: ', trim(FilterFile)
  write (*,*) ''
  flush(6)

  ! Check that the dimensions are consistent between the variables needed for
  ! the selected averaging function.
  !
  ! There is no associated filter with the max_azwind since a filter has already been
  ! applied by the azavg program (which created the azwind data). All other functions
  ! need the filter data.
  !
  if (AvgFunc .eq. 'max_azwind') then
    AzWindFile = trim(InDir) // '/speed_t' // trim(InSuffix)
    AzWind%vname = 'speed_t'
    call rhdf5_read_init(AzWindFile, AzWind)

    Nx = AzWind%dims(1)
    Ny = AzWind%dims(2)
    Nz = AzWind%dims(3)
    Nt = AzWind%dims(4)
  else
    ! Read in the filter and use it to check against all the other variables
    Filter%vname = 'filter'
    call rhdf5_read_init(FilterFile, Filter)
  
    Nx = Filter%dims(1)
    Ny = Filter%dims(2)
    Nz = Filter%dims(3)
    Nt = Filter%dims(4)

    if (AvgFunc .eq. 'horiz_ke') then
      DensFile = trim(InDir) // '/dn0' // trim(InSuffix)
      Dens%vname = 'dn0'
      call rhdf5_read_init(DensFile, Dens)
  
      Ufile = trim(InDir) // '/u' // trim(InSuffix)
      U%vname = 'u'
      call rhdf5_read_init(Ufile, U)
  
      Vfile = trim(InDir) // '/v' // trim(InSuffix)
      V%vname = 'v'
      call rhdf5_read_init(Vfile, V)

      if (.not. (DimsMatch(Filter, Dens) .and. DimsMatch(Filter, U) .and. DimsMatch(Filter, V))) then
        write (*,*) 'ERROR: dimensions of filter, dn0, u and v do not match'
        stop
      endif
    else if (AvgFunc .eq. 'storm_int') then
      Ufile = trim(InDir) // '/u' // trim(InSuffix)
      U%vname = 'u'
      call rhdf5_read_init(Ufile, U)
  
      Vfile = trim(InDir) // '/v' // trim(InSuffix)
      V%vname = 'v'
      call rhdf5_read_init(Vfile, V)
  
      if (.not. (DimsMatch(Filter, U) .and. DimsMatch(Filter, V))) then
        write (*,*) 'ERROR: dimensions of filter, u and v do not match'
        stop
      endif
    endif
  endif

  ! Report the dimensions
  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', Nx
  write (*,*) '  Number of y (latitude) points:           ', Ny
  write (*,*) '  Number of z (vertical level) points:     ', Nz
  write (*,*) '  Number of t (time) points:               ', Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''
  flush(6)

  ! Read in the field data, plus the t coordinates
  if (AvgFunc .eq. 'max_azwind') then
    write (*,*) 'Reading variable: speed_t'
    write (*,*) '  HDF5 file: ', trim(AzWindFile)
    write (*,*) ''
    call rhdf5_read(AzWindFile, AzWind)

    Tcoords%vname = 't_coords'
    call rhdf5_read_init(AzWindFile, Tcoords)
    call rhdf5_read(AzWindFile, Tcoords)

    Units = "m/s"
  else
    write (*,*) 'Reading variable: filter'
    write (*,*) '  HDF5 file: ', trim(FilterFile)
    write (*,*) ''
    call rhdf5_read(FilterFile, Filter)

    Tcoords%vname = 't_coords'
    call rhdf5_read_init(FilterFile, Tcoords)
    call rhdf5_read(FilterFile, Tcoords)

    if (AvgFunc .eq. 'horiz_ke') then
      write (*,*) 'Reading variable: dn0'
      write (*,*) '  HDF5 file: ', trim(DensFile)
      write (*,*) ''
      call rhdf5_read(DensFile, Dens)

      write (*,*) 'Reading variable: u'
      write (*,*) '  HDF5 file: ', trim(Ufile)
      write (*,*) ''
      call rhdf5_read(Ufile, U)

      write (*,*) 'Reading variable: v'
      write (*,*) '  HDF5 file: ', trim(Vfile)
      write (*,*) ''
      call rhdf5_read(Vfile, V)

      Units = 'Joules'
    else if (AvgFunc .eq. 'storm_int') then
      write (*,*) 'Reading variable: u'
      write (*,*) '  HDF5 file: ', trim(Ufile)
      write (*,*) ''
      call rhdf5_read(Ufile, U)

      write (*,*) 'Reading variable: v'
      write (*,*) '  HDF5 file: ', trim(Vfile)
      write (*,*) ''
      call rhdf5_read(Vfile, V)

      Units = 'int'
    endif
  endif

  write (*,*) ''
  flush(6)

  ! Allocate the output array and do the averaging
  TserAvg%vname = trim(AvgFunc)
  TserAvg%ndims = 4 
  TserAvg%dims(1) = 1
  TserAvg%dims(2) = 1
  TserAvg%dims(3) = 1
  TserAvg%dims(4) = Nt
  TserAvg%dimnames(1) = 'x' 
  TserAvg%dimnames(2) = 'y' 
  TserAvg%dimnames(3) = 'z' 
  TserAvg%dimnames(4) = 't' 
  TserAvg%units = Units
  TserAvg%descrip = 'time series averaged ' // trim(AvgFunc) 
  allocate(TserAvg%vdata(Nt))

  if (AvgFunc .eq. 'max_azwind') then
    call DoMaxAzWind(Nx, Ny, Nz, Nt, AzWind%vdata, UndefVal, TserAvg%vdata)
  else if (AvgFunc .eq. 'horiz_ke') then
    ! Need to get delta for x and y (in meters) and the z heights (in meters)
    ! for DoHorizKe to be able to calculate volume * density -> mass
    VarLon%vname = 'x_coords'
    call rhdf5_read_init(DensFile, VarLon)
    call rhdf5_read(DensFile, VarLon)

    VarLat%vname = 'y_coords'
    call rhdf5_read_init(DensFile, VarLat)
    call rhdf5_read(DensFile, VarLat)

    allocate(VarXcoords(Nx))
    allocate(VarYcoords(Ny))
    call ConvertGridCoords(Nx, Ny, Nz, Nt, VarLon%vdata, VarLat%vdata, VarXcoords, VarYcoords)

    DeltaX = VarXcoords(2) - VarXcoords(1)
    DeltaY = VarYcoords(2) - VarYcoords(1)

    VarZcoords%vname = 'z_coords'
    call rhdf5_read_init(DensFile, VarZcoords)
    call rhdf5_read(DensFile, VarZcoords)

    write (*,*) 'Horizontal grid info:'
    write (*,*) '  X range (min lon, max lon) --> (min x, max x): '
    write (*,*) '    ', VarLon%vdata(1), VarLon%vdata(Nx), VarXcoords(1), VarXcoords(Nx)
    write (*,*) '  Y range (min lat, max lat) --> (min y, max y): '
    write (*,*) '    ', VarLat%vdata(1), VarLat%vdata(Ny), VarYcoords(1), VarYcoords(Ny)
    write (*,*) ''
    write (*,*) 'Vertical grid info:'
    do iz = 1, Nz
      write (*,*) '  ', iz, ' --> ', VarZcoords%vdata(iz)
    end do
    write (*,*) ''
    flush(6)

    call DoHorizKe(Nx, Ny, Nz, Nt, Dens%vdata, U%vdata, V%vdata, Filter%vdata, DeltaX, DeltaY, VarZcoords%vdata, TserAvg%vdata)
    write (*,*) ''
    flush(6)
  else if (AvgFunc .eq. 'storm_int') then
    call DoStormInt(Nx, Ny, Nz, Nt, U%vdata, V%vdata, Filter%vdata, TserAvg%vdata)
  endif

  ! Create dummy coordinates (for GRADS sake) and write out the
  ! time series as a 4D var, (x,y,z,t), where x, y and z have
  ! dimension size of 1.
  Xcoords%vname = 'x_coords'
  Xcoords%ndims = 1
  Xcoords%dims(1) = 1
  Xcoords%dimnames(1) = 'x'
  Xcoords%units = 'degrees_east'
  Xcoords%descrip = 'longitude'
  allocate(Xcoords%vdata(1))
  Xcoords%vdata(1) = 1.0
  
  Ycoords%vname = 'y_coords'
  Ycoords%ndims = 1
  Ycoords%dims(1) = 1
  Ycoords%dimnames(1) = 'y'
  Ycoords%units = 'degrees_north'
  Ycoords%descrip = 'latitude'
  allocate(Ycoords%vdata(1))
  Ycoords%vdata(1) = 1.0
  
  Zcoords%vname = 'z_coords'
  Zcoords%ndims = 1
  Zcoords%dims(1) = 1
  Zcoords%dimnames(1) = 'z'
  Zcoords%units = 'meter'
  Zcoords%descrip = 'sigma-z'
  allocate(Zcoords%vdata(1))
  Zcoords%vdata(1) = 1.0

  ! third arg to rhdf5_write is "append" flag:
  !   0 - create new file
  !   1 - append to existing file
  write (*,*) 'Writing HDF5 output: ', trim(OutFile)
  write (*,*) ''
  call rhdf5_write(OutFile, TserAvg, 0)

  ! write out the coordinate data
  call rhdf5_write(OutFile, Xcoords, 1)
  call rhdf5_write(OutFile, Ycoords, 1)
  call rhdf5_write(OutFile, Zcoords, 1)
  call rhdf5_write(OutFile, Tcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(OutFile, Xcoords, 'x')
  call rhdf5_set_dimension(OutFile, Ycoords, 'y')
  call rhdf5_set_dimension(OutFile, Zcoords, 'z')
  call rhdf5_set_dimension(OutFile, Tcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(OutFile, TserAvg)
  
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
subroutine GetMyArgs(InDir, InSuffix, OutFile, AvgFunc, FilterFile)
  implicit none

  character (len=*) :: InDir, InSuffix, OutFile, AvgFunc, FilterFile

  integer :: iargc
  character (len=128) :: arg
  integer :: Nitems

  logical :: BadArgs

  if (iargc() .ne. 5) then
    write (*,*) 'ERROR: must supply exactly 5 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: tsavg <in_dir> <in_suffix> <out_file> <avg_function> <filter_file>'
    write (*,*) '        <in_dir>: directory where input files live'
    write (*,*) '        <in_suffix>: suffix on input file names'
    write (*,*) '        <out_file>: name of output file, HDF5 format'
    write (*,*) '        <avg_function>: averaging function to use on input data'
    write (*,*) '            horiz_ke -> total kinetic energy form horizontal winds'
    write (*,*) '            storm_int -> storm intensity metric from horizontal wind speeds'
    write (*,*) '            max_azwind -> max value of azimuthially averaged wind'
    write (*,*) '        <filter_file>: file containing the filter mask'
    stop
  end if

  call getarg(1, InDir)
  call getarg(2, InSuffix)
  call getarg(3, OutFile)
  call getarg(4, AvgFunc)
  call getarg(5, FilterFile)

  BadArgs = .false.

  if ((AvgFunc .ne. 'horiz_ke')       .and. &
      (AvgFunc .ne. 'storm_int')  .and. &
      (AvgFunc .ne. 'max_azwind')) then
    write (*,*) 'ERROR: <avg_function> must be one of:'
    write (*,*) '          horiz_ke'
    write (*,*) '          storm_int'
    write (*,*) '          max_azwind'
    write (*,*) ''
    BadArgs = .true.
  end if

  if (BadArgs) then
    stop
  end if

  return
end subroutine GetMyArgs

!************************************************************************************
! DoMaxAzWind()
!
! This subroutine will simply find the maximum wind speed in AzWind and copy that
! to TserAvg.
subroutine DoMaxAzWind(Nx, Ny, Nz, Nt, AzWind, UndefVal, TserAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(Nx,Ny,Nz,Nt) :: AzWind
  real :: UndefVal
  real, dimension(Nt) :: TserAvg

  integer :: ix,iy,iz,it

  ! dimension order is: x,y,z,t

  do it = 1, Nt
    TserAvg(it) = 0.0
    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          if ((AzWind(ix,iy,iz,it) .gt. TserAvg(it)) .and. (AzWind(ix,iy,iz,it) .ne. UndefVal)) then
            TserAvg(it) = AzWind(ix,iy,iz,it)
          endif
        end do
      end do
    end do
  end do

  return
end subroutine DoMaxAzWind

!**************************************************************************************
! DoHorizKe()
!
! This routine will calculate the total kinetic energy over the given cylindrical
! volume. Do not want average since we want the size of the storm reflected in
! this diagnostic.

subroutine DoHorizKe(Nx, Ny, Nz, Nt, Dens, U, V, Filter, DeltaX, DeltaY, Zcoords, TserAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(Nx,Ny,Nz,Nt) :: Dens, U, V, Filter
  real :: DeltaX, DeltaY
  real, dimension(Nz) :: Zcoords
  real, dimension(Nt) :: TserAvg

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

  do it = 1, Nt
    SumKe = 0.0
    NumPoints = 0

    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          if (anint(Filter(ix,iy,iz,it)) .eq. 1.0) then
            if (iz .eq. Nz) then
              ! Use the level below for this case (since no level above)
              LevThickness = Zcoords(iz) - Zcoords(iz-1)
            else
              LevThickness = Zcoords(iz+1) - Zcoords(iz)
            end if
            CurrKe = 0.5 * DeltaX * DeltaY * LevThickness * Dens(ix,iy,iz,it) * (U(ix,iy,iz,it)**2 + V(ix,iy,iz,it)**2)
            SumKe = SumKe + CurrKe
            NumPoints = NumPoints + 1
          endif
        enddo
      enddo
    enddo

    TserAvg(it) = SumKe;
    write (*,*) 'HorizKe: Timestep:', it, ', Number of points selected: ', NumPoints
  end do
  
end subroutine DoHorizKe

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

subroutine DoStormInt(Nx, Ny, Nz, Nt, U, V, Filter, TserAvg)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(Nx,Ny,Nz,Nt) :: U, V, Filter
  real, dimension(Nt) :: TserAvg

  integer ix,iy,iz,it
  integer nCat0, nCat1, nCat2, nCat3, nCat4, nCat5, NumPoints
  real Wspeed, SiMetric

  do it = 1, Nt
    nCat0 = 0
    nCat1 = 0
    nCat2 = 0
    nCat3 = 0
    nCat4 = 0
    nCat5 = 0

    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          if (anint(Filter(ix,iy,iz,it)) .eq. 1.0) then
            ! Count up the number of grid points with wind speeds fitting each of the
            ! Saffir-Simpson categories. Then form the metric by weighting each category
            ! count.
            Wspeed = sqrt(U(ix,iy,iz,it)**2 + V(ix,iy,iz,it)**2)
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
    enddo

    !Linear weighting
    NumPoints = nCat0 + nCat1 + nCat2 + nCat3 + nCat4 + nCat5
    SiMetric = float(nCat1) + (float(nCat2)*2.0) + (float(nCat3)*3.0) + (float(nCat4)*4.0) + (float(nCat5)*5.0)

    if (NumPoints .eq. 0) then
      TserAvg(it) = 0.0
      write (*,*) 'WARNING: no data points selected for time step: ', it
    else
      TserAvg(it) = SiMetric / float(NumPoints)
      write (*,*) 'StormInt: Timestep:', it, ', Number of points selected: ', NumPoints
    end if
  end do
  
end subroutine DoStormInt


!!! !*****************************************************************************
!!! ! DoCloud()
!!! !
!!! ! This subroutine will perform the cloud droplet total mass time series
!!! ! averaging. Can select between supercooled or warm rain droplets.
!!! !
!!! 
!!! subroutine DoCloud(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, Dens, CintLiq, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   logical :: DoSc
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
!!!   type (GradsVar) :: Cloud, TempC, Dens, CintLiq, TserAvg
!!!   integer, dimension(1:Cloud%Nt) :: StmIx, StmIy
!!!   real, dimension(1:Cloud%Nx) :: Xcoords
!!!   real, dimension(1:Cloud%Ny) :: Ycoords
!!!   real, dimension(1:Cloud%Nz) :: Zcoords
!!! 
!!!   real, dimension(0:Cloud%Nz) :: ZmHeights
!!!   integer :: ix, iy, iz, it, NumPoints
!!! 
!!!   call SetZmHeights (Cloud%Nz, ZmHeights)
!!! 
!!!   ! Convert the cloud mixing ratio to mass for each grid point. 
!!!   !
!!!   ! Mixing ratios in GRADS files are g/kg
!!!   ! Density is in kg/m**3
!!!   ! Heights are in m
!!!   ! So, express the mass value in g using the formula
!!!   !   (mix ratio) * (density) * (layer thickness) * (layer horiz area)
!!!   ! The layer thickness for layer k is: ZmHeights(k) - ZmHeights(k-1)
!!! 
!!!   do it = 1, Cloud%Nt
!!!     ! Sum up the cloud droplet mass over the specified radial band. Only include the
!!!     ! grid points where tempc is 0 or less (supercooled)
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, Cloud%Nz
!!!       do ix = 1, Cloud%Nx
!!!         do iy = 1, Cloud%Ny
!!!           if ((InsideCylVol(Cloud%Nx, Cloud%Ny, Cloud%Nz, Cloud%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords))  .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
!!!             if (DoSc) then
!!!               if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
!!!                 TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + (Cloud%Vdata(it,iz,ix,iy) * Dens%Vdata(it,iz,ix,iy) * (ZmHeights(iz)-ZmHeights(iz-1)))
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             else
!!!               if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
!!!                 TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + (Cloud%Vdata(it,iz,ix,iy) * Dens%Vdata(it,iz,ix,iy) * (ZmHeights(iz)-ZmHeights(iz-1)))
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     ! At this point TserAvg%Vdata(it,1,1,1) holds g/m**2, multiply by grid cell horizontal area. Note this assumes
!!!     ! each grid cell has the same horizontal area.
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) * DeltaX * DeltaY
!!!     if (NumPoints .eq. 0) then
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       write (*,*) 'ScCloud: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     endif
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoWup()
!!! !
!!! ! This subroutine will do the average vertical velocity in regions of significant
!!! ! updrafts.
!!! !
!!! 
!!! subroutine DoWup(DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, W, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: W, TserAvg
!!!   integer, dimension(1:W%Nt) :: StmIx, StmIy
!!!   real, dimension(1:W%Nx) :: Xcoords
!!!   real, dimension(1:W%Ny) :: Ycoords
!!!   real, dimension(1:W%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it, NumPoints
!!! 
!!!   do it = 1, W%Nt
!!!     ! Average w over regions where significant updrafts occur
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, W%Nz
!!!       do ix = 1, W%Nx
!!!         do iy = 1, W%Ny
!!!           if (InsideCylVol(W%Nx, W%Ny, W%Nz, W%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             if (W%Vdata(it,iz,ix,iy) .ge. Wthreshold) then
!!!               TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + W%Vdata(it,iz,ix,iy)
!!!               NumPoints = NumPoints + 1
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) / float(NumPoints)
!!!       write (*,*) 'Wup: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!   end do
!!! 
!!! end subroutine
!!! 
!!! 
!!! !************************************************************************************
!!! ! DoCloudDiam()
!!! !
!!! ! This routine will calculate an averaged cloud droplet diameter. Can select between
!!! ! warm rain droplets or supercooled droplets.
!!! 
!!! subroutine DoCloudDiam(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, CloudDiam, CintLiq, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   logical :: DoSc
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
!!!   type (GradsVar) :: Cloud, TempC, CloudDiam, CintLiq, TserAvg
!!!   integer, dimension(1:CloudDiam%Nt) :: StmIx, StmIy
!!!   real, dimension(1:CloudDiam%Nx) :: Xcoords
!!!   real, dimension(1:CloudDiam%Ny) :: Ycoords
!!!   real, dimension(1:CloudDiam%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it, NumPoints
!!!   real SumQ, SumQD
!!!   real MaxQ, Climit, SumD
!!! 
!!!   do it = 1, CloudDiam%Nt
!!!     ! Calculate a mass-weighted mean diameter for supercooled cloud droplets.
!!!     ! 
!!!     !    Mean diameter (TserAvg value) = Sum(cloud * cloud_d) / Sum(cloud)
!!!     !    where cloud and cloud_d are only included in the sum when tempc is <= 0
!!!     !
!!!     ! Find the max Q (mass) value and use it to filter data - select data points if
!!!     ! the Q is within 20% of the max Q (.2 to 1.0). Then form the average of the diameters
!!!     ! of the selected points.
!!!     SumQD = 0.0
!!!     SumQ = 0.0
!!!     SumD = 0.0
!!!     NumPoints = 0
!!!     do iz = 1, CloudDiam%Nz
!!!       do ix = 1, CloudDiam%Nx
!!!         do iy = 1, CloudDiam%Ny
!!!           if ((InsideCylVol(CloudDiam%Nx, CloudDiam%Ny, CloudDiam%Nz, CloudDiam%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
!!!             if (DoSc) then
!!!               if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
!!!                  SumQD = SumQD + (Cloud%Vdata(it,iz,ix,iy) * CloudDiam%Vdata(it,iz,ix,iy))
!!!                  SumQ = SumQ + Cloud%Vdata(it,iz,ix,iy)
!!!                  SumD = SumD + CloudDiam%Vdata(it,iz,ix,iy)
!!!                  NumPoints = NumPoints + 1
!!!               end if
!!!             else
!!!               if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
!!!                  SumQD = SumQD + (Cloud%Vdata(it,iz,ix,iy) * CloudDiam%Vdata(it,iz,ix,iy))
!!!                  SumQ = SumQ + Cloud%Vdata(it,iz,ix,iy)
!!!                  SumD = SumD + CloudDiam%Vdata(it,iz,ix,iy)
!!!                  NumPoints = NumPoints + 1
!!!               end if
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (SumQ .eq. 0.0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = SumQD / SumQ
!!!       write (*,*) 'ScCloudDiam: Ts:', it, ', NumPoints: ', NumPoints
!!!     end if
!!! 
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoCloudConc()
!!! !
!!! ! This subroutine will do the average cloud droplet concentration
!!! 
!!! subroutine DoCloudConc(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, DoSc, TempC, CloudConc, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: CloudConc, TempC, TserAvg
!!!   integer, dimension(1:CloudConc%Nt) :: StmIx, StmIy
!!!   real, dimension(1:CloudConc%Nx) :: Xcoords
!!!   real, dimension(1:CloudConc%Ny) :: Ycoords
!!!   real, dimension(1:CloudConc%Nz) :: Zcoords
!!!   logical :: DoSc
!!! 
!!!   integer ix,iy,iz,it
!!!   integer NumPoints
!!!   real SumCloudConc
!!! 
!!!   ! Calculate the average cloud droplet concentration near the eyewall region.
!!!   ! 
!!!   ! Call the cloud base to be between 1000m and 1200m. The z level corresponding to
!!!   ! that is z = 4 which is 1138m. Cover from surface to the cloud base level. This
!!!   ! results in using z = 1 to z = 4.
!!! 
!!!   do it = 1, CloudConc%Nt
!!!     SumCloudConc = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, CloudConc%Nz
!!!       do ix = 1, CloudConc%Nx
!!!         do iy = 1, CloudConc%Ny
!!!           if (InsideCylVol(CloudConc%Nx, CloudConc%Ny, CloudConc%Nz, CloudConc%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             if (DoSc) then
!!!               if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
!!!                 SumCloudConc = SumCloudConc + CloudConc%Vdata(it,iz,ix,iy)
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             else
!!!               if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
!!!                 SumCloudConc = SumCloudConc + CloudConc%Vdata(it,iz,ix,iy)
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = SumCloudConc / float(NumPoints)
!!!       write (*,*) 'CloudConc: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoPrecipR()
!!! !
!!! ! This subroutine will do the average cloud droplet concentration near the eyewall
!!! ! region.
!!! 
!!! subroutine DoPrecipR(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, PrecipR, CintLiq, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
!!!   ! PrecipR is 2D var
!!!   type (GradsVar) :: PrecipR, CintLiq, TserAvg
!!!   integer, dimension(1:PrecipR%Nt) :: StmIx, StmIy
!!!   real, dimension(1:PrecipR%Nx) :: Xcoords
!!!   real, dimension(1:PrecipR%Ny) :: Ycoords
!!!   real, dimension(1:PrecipR%Nz) :: Zcoords
!!! 
!!!   integer :: ix,iy,iz,it
!!!   integer :: NumPoints
!!!   real :: SumPrecip
!!! 
!!!   do it = 1, PrecipR%Nt
!!!     SumPrecip = 0.0
!!!     NumPoints = 0
!!! 
!!!     do ix = 1, PrecipR%Nx
!!!       do iy = 1, PrecipR%Ny
!!!         if (InsideCylVol(PrecipR%Nx, PrecipR%Ny, PrecipR%Nz, PrecipR%Nt, ix, iy, 1, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords) .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
!!!           SumPrecip = SumPrecip + PrecipR%Vdata(it,1,ix,iy)
!!!           NumPoints = NumPoints + 1
!!!         end if
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) '  WARNING: no data points selected for this time step'
!!!     else
!!!       ! At this point TserAvg%Vdata(it,1,1,1) holds mm/hr, multiply by grid cell horizontal area.
!!!       ! Note this assumes each grid cell has the same horizontal area. What this does
!!!       ! is convert mm/hr to kg/hr when assuming that the density of water is 1000kg/m**3.
!!!       !   mm/hr * m**2 * 1000 kg/m**3 * 0.001 m/mm -> kg/hr
!!!       TserAvg%Vdata(it,1,1,1) = (SumPrecip / float(NumPoints)) * DeltaX * DeltaY
!!!       write (*,*) 'Precip: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!     write (*,*) ''
!!!     flush(6)
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoCcnConc()
!!! !
!!! ! This subroutine will do the average CCN concentration.
!!! !
!!! 
!!! subroutine DoCcnConc(DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CcnConc, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: CcnConc, TserAvg
!!!   integer, dimension(1:CcnConc%Nt) :: StmIx, StmIy
!!!   real, dimension(1:CcnConc%Nx) :: Xcoords
!!!   real, dimension(1:CcnConc%Ny) :: Ycoords
!!!   real, dimension(1:CcnConc%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it, NumPoints
!!! 
!!!   do it = 1, CcnConc%Nt
!!!     ! Average w over regions where significant updrafts occur
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, CcnConc%Nz
!!!       do ix = 1, CcnConc%Nx
!!!         do iy = 1, CcnConc%Ny
!!!           if (InsideCylVol(CcnConc%Nx, CcnConc%Ny, CcnConc%Nz, CcnConc%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + CcnConc%Vdata(it,iz,ix,iy)
!!!             NumPoints = NumPoints + 1
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) / float(NumPoints)
!!!       write (*,*) 'Wup: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!   end do
!!! 
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoTestCvs()
!!! !
!!! ! This subroutine will perform a test on the cylindrical volume selection routine.
!!! ! Just runs through the entire grid and outputs a '1' when selection occurs otherwise
!!! ! outputs a '0'. Then view the result in grads and see if selection is correct.
!!! 
!!! subroutine DoTestCvs(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, TestSelect)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: TestSelect
!!!   integer, dimension(1:TestSelect%Nt) :: StmIx, StmIy
!!!   real, dimension(1:TestSelect%Nx) :: Xcoords
!!!   real, dimension(1:TestSelect%Ny) :: Ycoords
!!!   real, dimension(1:TestSelect%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it
!!!   integer NumPoints
!!! 
!!!   write (*,*) 'Testing cylindrical volume selection:'
!!!   do it = 1, TestSelect%Nt
!!!     NumPoints = 0
!!!     do iz = 1, TestSelect%Nz
!!!       do ix = 1, TestSelect%Nx
!!!         do iy = 1, TestSelect%Ny
!!!           if (InsideCylVol(TestSelect%Nx, TestSelect%Ny, TestSelect%Nz, TestSelect%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             TestSelect%Vdata(it,iz,ix,iy) = 1.0
!!!             NumPoints = NumPoints + 1
!!!           else
!!!             TestSelect%Vdata(it,iz,ix,iy) = 0.0
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!!     ! mark the storm center
!!!     TestSelect%Vdata(it,iz,StmIx(it),StmIy(it)) = 2.0
!!!     write (*,*) '  Timestep, Number of points selected: ', it, NumPoints
!!!   end do
!!! end subroutine
!!! 
!!! !******************************************************************************
!!! ! ConvertTinc()
!!! !
!!! ! This function will convert the time increment spec'd in the GRAD control
!!! ! file into a number of seconds.
!!! !
!!! 
!!! real function ConvertTinc(Tinc)
!!!   implicit none
!!! 
!!!   character (len=*) :: Tinc
!!! 
!!!   character (len=128) :: Tval, Tunits, Tfmt
!!!   integer :: i, Uloc, Tlen, Itval
!!!   
!!!   ! Walk through the Tinc string. Concatenate the numbers onto Tval and
!!!   ! the alaph characters onto Tunits. This algorithm assumes that GRADS
!!!   ! will use a format like <numeric_value><units> for Tinc where <units>
!!!   ! is a string like 'hr' or 'mn'.
!!!   Uloc = 0
!!!   Tlen = len_trim(Tinc)
!!!   do i = 1, Tlen
!!!     if ((Tinc(i:i) .ge. 'A') .and. (Tinc(i:i) .le. 'z')) then
!!!       if (Uloc .eq. 0) then
!!!         ! the first alpha character is the beginning of the spec for units
!!!         Uloc = i
!!!       end if
!!!     end if
!!!   end do
!!! 
!!!   write (Tfmt, '(a,i2.2,a,i2.2,a)') '(a', (Uloc-1), 'a' , ((Tlen-Uloc)+1), ')'
!!!   read (Tinc, Tfmt) Tval, Tunits
!!! 
!!!   read(Tval, '(i)') Itval
!!!   if (Tunits .eq. 'hr') then
!!!     ConvertTinc = float(Itval) * 3600.0
!!!   else
!!!     if (Tunits .eq. 'mn') then
!!!       ConvertTinc = float(Itval) * 60.0
!!!     else
!!!       ConvertTinc = float(Itval)
!!!     end if
!!!   end if
!!!   return
!!! end function
!!! 
!!! !******************************************************************************
!!! ! CalcRates()
!!! !
!!! ! This routine will calculate time derivatives of the input data using
!!! ! a centered difference method.
!!! !
!!! 
!!! subroutine CalcRates(DeltaT, TserAvg, Rates)
!!!   use gdata_utils
!!!   implicit none
!!! 
!!!   type (GradsVar) :: TserAvg, Rates
!!!   real :: DeltaT
!!! 
!!!   real :: f1, f2
!!!   integer :: ix, iy, iz, it
!!! 
!!!   ! use a centered difference, uses points at t-1, t and t+1
!!!   do ix = 1, TserAvg%Nx
!!!     do iy = 1, TserAvg%Ny
!!!       do iz = 1, TserAvg%Nz
!!!         do it = 2, TserAvg%Nt-1
!!!           f1 = (TserAvg%Vdata(it,1,1,1) + TserAvg%Vdata(it-1,1,1,1)) / 2.0
!!!           f2 = (TserAvg%Vdata(it+1,1,1,1) + TserAvg%Vdata(it,1,1,1)) / 2.0
!!! 
!!!           Rates%Vdata(it-1,1,1,1) = (f2 - f1) / DeltaT
!!!         end do
!!!       end do
!!!     end do
!!!   end do
!!! end subroutine

end program tsavg
