!***************************************************************
! Program to do averaging over the domain and yield a single
! data point at each time step
!
! This program will read in HDF5 data from a RAMS simulation, 
! perform an averaging function, and output the single point time
! series in HDF5 format
!

program hdata_op
  use rhdf5_utils
  use diag_utils
  implicit none

  integer, parameter :: LargeString  = 512
  integer, parameter :: MediumString = 256
  integer, parameter :: LittleString = 128
  integer, parameter :: MaxArgFields = 20

  real, parameter :: UndefVal = -999.0

  character (len=MediumString) :: InFile1
  character (len=LittleString) :: VarName1
  character (len=MediumString) :: InFile2
  character (len=LittleString) :: VarName2
  character (len=MediumString) :: OutFile
  character (len=LittleString) :: Op

  ! Data arrays
  ! Dims: x, y, z, t
  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords
  type (Rhdf5Var) :: Var1, Var2, OutVar
  character (len=LittleString) :: rh5f_facc
  integer :: rh5f_in1, rh5f_in2, rh5f_out

  integer :: i, ix, iy, iz, it, id
  integer :: Nx, Ny, Nz, Nt
  logical :: BadDims

  ! Get the command line arguments
  call GetMyArgs(InFile1, VarName1, InFile2, VarName2, OutFile, Op)

  write (*,*) 'Performing operation on HDF5 REVU data:'
  write (*,*) '  Input file 1: ', trim(InFile1)
  write (*,*) '    Variable Name: ', trim(VarName1)
  write (*,*) '  Input file 2: ', trim(InFile2)
  write (*,*) '    Variable Name: ', trim(VarName2)
  write (*,*) '  Output file:  ', trim(OutFile)
  write (*,*) '  Operator: ', trim(Op)
  write (*,*) ''
  flush(6)

  ! set up file and variable names
  Var1%vname = VarName1
  Var2%vname = VarName2

  ! Check that the dimensions are consistent between the variables needed for
  ! the selected averaging function.

  call rhdf5_read_init(InFile1, Var1)
  call rhdf5_read_init(InFile2, Var2)

  Nx = Var1%dims(1)
  Ny = Var1%dims(2)
  Nz = Var1%dims(3)
  Nt = Var1%dims(4)

  BadDims = Var2%dims(1) .ne. Nx
  BadDims = BadDims .or. (Var2%dims(2) .ne. Ny)
  BadDims = BadDims .or. (Var2%dims(3) .ne. Nz)
  BadDims = BadDims .or. (Var2%dims(4) .ne. Nt)

  if (BadDims) then
    write (*,*) 'ERROR: dimensions of variables from input files do not match'
    stop
  endif

  ! Check if the operation makes sense
  if (trim(Op) .eq. 'sub') then
    if (trim(Var1%units) .ne. trim(Var2%units)) then
      write (*,*) 'ERROR: units of variables from input files do not match'
      write (*,*) 'ERROR:   Var1: ', trim(Var1%vname), ' -> ', trim(Var1%units)
      write (*,*) 'ERROR:   Var2: ', trim(Var2%vname), ' -> ', trim(Var2%units)
      stop
    endif
  endif

  ! Set up the dimensions for reading in the input field data, one time step per read. 
  ! In other words, remove the time dimension from the input dimensions since we will 
  ! be incrementing through every time step in a loop. The time dimension is always the
  ! last dimension so what this boils down to is to decrement number of dimensions by one.
  Var1%ndims = Var1%ndims - 1
  Var2%ndims = Var2%ndims - 1
  write (*,*) 'Input variable information:'
  write (*,*) '  Number of dimensions: ', Var1%ndims
  write (*,*) '  Dimension sizes:'
  do id = 1, Var1%ndims
    write (*,*), '    ', trim(Var1%dimnames(id)), ': ', Var1%dims(id)
  enddo
  write (*,*) ''
  flush(6)

  ! Set up the dimensions for the output and allocate the output data array.
  ! Just copy the setup from Var1.
  OutVar%vname = trim(VarName1) // '_' // trim(Op) // '_' // trim(VarName2)
  if (trim(Op) .eq. 'sub') then
    OutVar%descrip = 'var1 minus var2'
  endif
  OutVar%units = Var1%units

  OutVar%ndims = Var1%ndims 
  OutVar%dims(1) = Nx
  OutVar%dims(2) = Ny
  OutVar%dims(3) = Nz
  OutVar%dimnames(1) = Var1%dimnames(1)
  OutVar%dimnames(2) = Var1%dimnames(2)
  OutVar%dimnames(3) = Var1%dimnames(3)

  allocate(OutVar%vdata(Nx*Ny*Nz))
 
  ! Report the dimensions
  write (*,*) 'Output variable information:'
  write (*,*) '  Name: ', trim(OutVar%vname)
  write (*,*) '  Units: ', trim(OutVar%units)
  write (*,*) '  Description: ', trim(OutVar%descrip)
  write (*,*) '  Number of dimensions: ', OutVar%ndims
  write (*,*) '  Dimension sizes:'
  do id = 1, OutVar%ndims
    write (*,*), '    ', trim(OutVar%dimnames(id)), ': ', OutVar%dims(id)
  enddo
  write (*,*) ''
  flush(6)

  ! Report time steps
  write (*,*) 'Number of time steps: ', Nt
  write (*,*) ''
  flush(6)

  ! Read in the input coordinates
  Xcoords%vname = 'x_coords'
  call rhdf5_read_init(InFile1, Xcoords)
  call rhdf5_read(InFile1, Xcoords)

  Ycoords%vname = 'y_coords'
  call rhdf5_read_init(InFile1, Ycoords)
  call rhdf5_read(InFile1, Ycoords)

  Zcoords%vname = 'z_coords'
  call rhdf5_read_init(InFile1, Zcoords)
  call rhdf5_read(InFile1, Zcoords)

  Tcoords%vname = 't_coords'
  call rhdf5_read_init(InFile1, Tcoords)
  call rhdf5_read(InFile1, Tcoords)

  ! Perform the operation
  rh5f_facc = 'W'
  call rhdf5_open_file(OutFile, rh5f_facc, 1, rh5f_out)

  rh5f_facc = 'R'
  call rhdf5_open_file(InFile1, rh5f_facc, 0, rh5f_in1)
  call rhdf5_open_file(InFile2, rh5f_facc, 0, rh5f_in2)

  do it = 1, Nt
    call rhdf5_read_variable(rh5f_in1, Var1%vname, Var1%ndims, it, Var1%dims, rdata=Var1%vdata)
    call rhdf5_read_variable(rh5f_in2, Var2%vname, Var2%ndims, it, Var2%dims, rdata=Var2%vdata)

    ! do the op here
    i = 0
    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          i = i + 1
          if (trim(Op) .eq. 'sub') then
            OutVar%vdata(i) = Var1%vdata(i) - Var2%vdata(i)
          endif
        enddo
      enddo
    enddo

    deallocate(Var1%vdata)
    deallocate(Var2%vdata)

    ! write out the averaged data
    call rhdf5_write_variable(rh5f_out, OutVar%vname, OutVar%ndims, it, OutVar%dims, &
       OutVar%units, OutVar%descrip, OutVar%dimnames, rdata=OutVar%vdata)
 
    ! print a message for the user on longer jobs so that it can be
    ! seen that progress is being made
    if (modulo(it,100) .eq. 0) then
      write (*,*) 'Working: Number of time steps processed so far: ', it
    endif
  enddo

  call rhdf5_close_file(rh5f_in1)
  call rhdf5_close_file(rh5f_in1)
  call rhdf5_close_file(rh5f_out)
  deallocate(OutVar%vdata)

  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''

  ! Finish off output file
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
  call rhdf5_attach_dimensions(OutFile, OutVar)
  
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
subroutine GetMyArgs(InFile1, VarName1, InFile2, VarName2, OutFile, Op)
  implicit none

  character (len=*) :: InFile1, VarName1, InFile2, VarName2, OutFile, Op

  integer :: iargc
  character (len=128) :: arg
  integer :: Nitems

  logical :: BadArgs

  if (iargc() .ne. 6) then
    write (*,*) 'ERROR: must supply exactly 6 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: hdata_op <in_file1> <var_name1> <in_file2> <var_name2> <out_file> <operator>'
    write (*,*) '        <in_file1>: input file, HDF5 format'
    write (*,*) '        <var_name1>: name of variable inside input file 1'
    write (*,*) '        <in_file2>: input file, HDF5 format'
    write (*,*) '        <var_name2>: name of variable inside input file 2'
    write (*,*) '        <out_file>: output file, HDF5 format'
    write (*,*) '        <operator>: operation to perform on data from input files'
    write (*,*) '            order of operation: <in_file1> <operator> <in_file2>'
    write (*,*) '            sub -> subtract values in <in_file2> from those in <in_file1>'
    stop
  end if

  call getarg(1, InFile1)
  call getarg(2, VarName1)
  call getarg(3, InFile2)
  call getarg(4, VarName2)
  call getarg(5, OutFile)
  call getarg(6, Op)

  BadArgs = .false.

  if (Op .ne. 'sub') then
    write (*,*) 'ERROR: <operator> must be one of:'
    write (*,*) '          sub'
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
! This subroutine will find the maximum wind speed in AzWind and copy that
! to TserAvg. The intent of this metric is to mimic the Saffir-Simpson scale
! which uses the average 10m tangential wind speed that has been sustained
! over 10 minutes.
!
! Typically, don't have enough RAMS files to do 10m time averaging so as a
! proxy use the azimuthally averaged tangetial wind speed.
!
subroutine DoMaxAzWind(Nx, Ny, Nz, AzWind, UndefVal, AzWindMax)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: AzWind
  real :: UndefVal
  real :: AzWindMax

  integer :: ix,iy,iz

  ! dimension order is: x,y,z

  ! RAMS places the first z level below the surface so use the second level
  ! to approximate the 10m winds.

  AzWindMax = 0.0
  iz = 2
  do iy = 1, Ny
    do ix = 1, Nx
      if ((AzWind(ix,iy,iz) .gt. AzWindMax) .and. (anint(AzWind(ix,iy,iz)) .ne. UndefVal)) then
        AzWindMax = AzWind(ix,iy,iz)
      endif
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

subroutine DoHorizKe(Nx, Ny, Nz, Dens, U, V, Filter, DeltaX, DeltaY, Zcoords, TserAvg)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: Dens, U, V, Filter
  real :: DeltaX, DeltaY
  real, dimension(Nz) :: Zcoords
  real :: TserAvg

  integer ix,iy,iz
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
  !
  ! The integrated kinetic engery measurement done on hurricanes is taken over the domain on the 10m
  ! level. Just use the first model level above the surface (z == 2) for this measurement which will
  ! be close enough.

  SumKe = 0.0
  NumPoints = 0

  iz = 2

  do iy = 1, Ny
    do ix = 1, Nx
      if (anint(Filter(ix,iy,iz)) .eq. 1.0) then
        if (iz .eq. Nz) then
          ! Use the level below for this case (since no level above)
          LevThickness = Zcoords(iz) - Zcoords(iz-1)
        else
          LevThickness = Zcoords(iz+1) - Zcoords(iz)
        end if
        CurrKe = 0.5 * DeltaX * DeltaY * LevThickness * Dens(ix,iy,iz) * (U(ix,iy,iz)**2 + V(ix,iy,iz)**2)
        SumKe = SumKe + CurrKe
        NumPoints = NumPoints + 1
      endif
    enddo
  enddo

  TserAvg = SumKe;
  
  return
end subroutine DoHorizKe

!**************************************************************************************
! DoHda()
!
! This routine will do horizontal domain average for all z levels.
!
subroutine DoHda(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomAvg)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(Nz) :: DomAvg

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint
  integer :: NumPoints

  do iz = 1, Nz
    DomAvg(iz) = 0.0
    NumPoints = 0

    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          DomAvg(iz) = DomAvg(iz) + Var(ix,iy,iz)
          NumPoints = NumPoints + 1
        endif
      enddo
    enddo

    if (NumPoints .eq. 0) then
      DomAvg(iz) = UndefVal
    else
      DomAvg(iz) = DomAvg(iz) / float(NumPoints)
    endif
  enddo

  return
end subroutine DoHda

!**************************************************************************************
! DoMin()
!
! This routine will return the domain minimum value. Values equal to UndefVal will
! be ignored.
!
subroutine DoMin(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomMin)
  implicit none

  real, parameter :: BigPosNum = 10.0e50

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real :: DomMin

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint

  DomMin = BigPosNum
  do iz = 1, Nz
    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          if (Var(ix,iy,iz) .lt. DomMin) then
            DomMin = Var(ix,iy,iz)
          endif
        endif
      enddo
    enddo
  enddo

  if (DomMin .eq. BigPosNum) then
    ! all entries were UndefVal
    DomMin = UndefVal
  endif

  return
end subroutine DoMin

!**************************************************************************************
! DoMax()
!
! This routine will return the domain maximum value. Values equal to UndefVal will
! be ignored.
!
subroutine DoMax(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomMax)
  implicit none

  real, parameter :: BigNegNum = -10.0e50

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real :: DomMax

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint

  DomMax = BigNegNum
  do iz = 1, Nz
    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          if (Var(ix,iy,iz) .gt. DomMax) then
            DomMax = Var(ix,iy,iz)
          endif
        endif
      enddo
    enddo
  enddo

  if (DomMax .eq. BigNegNum) then
    ! all entries were UndefVal
    DomMax = UndefVal
  endif

  return
end subroutine DoMax

!**************************************************************************************
! DoHist()
!
! This routine will do histogram binning over all of the domain.
!

subroutine DoHist(Nx, Ny, Nz, FilterNz, Nb, Var, Filter, UseFilter, UndefVal, Bins, Counts)
  implicit none

  integer :: Nx, Ny, Nz, Nb, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(Nb) :: Bins, Counts

  integer :: ib, ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint

  ! Emulate matlab histc command, ie the values in Bins are treated as the edges
  ! of the bin ranges. Do not count any values from Var that fall outside the
  ! bin range (< Bins(1) or > Bins(Nb)). The count is incremented for bin ib when:
  !
  !   Bins(ib) <= Var(ix,iy,iz) < Bins(ib+1)
  !
  ! The count in Bin(Nb), last bin, is incremented when:
  !
  !   Var(ix,iy,iz) == Bins(ib)
  !
  do ib = 1, Nb
    Counts(ib) = 0.0
  enddo

  ! Build the histogram (counts)
  do iz = 1, Nz
    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 1, Ny
      do ix = 1, Nx
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          ! Check all bins except the last.
          !
          ! Exiting out of the loop when finding the bin will help a lot when the
          ! distribution of values is biased toward smaller values. After exiting
          ! out of the loop you can either just check the last bin (which will be wasted)
          ! or put in a logical variable and check that variable saying you can skip the
          ! check of the last bin. Since you would have to check the logical variable and
          ! the last bin every time you might as well just check the last bin instead.
          do ib = 1, Nb-1
             if ((Bins(ib) .le. Var(ix,iy,iz)) .and. (Var(ix,iy,iz) .lt. Bins(ib+1))) then
                Counts(ib) = Counts(ib) + 1.0
                exit
             endif
          enddo
          ! check the last bin
          if (Bins(Nb) .eq. Var(ix,iy,iz)) then
            Counts(Nb) = Counts(Nb) + 1.0
          endif
        endif
      enddo
    enddo
  enddo

  return
end subroutine DoHist

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

subroutine DoStormInt(Nx, Ny, Nz, Speed10m, Filter, TserAvg)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: Filter
  real, dimension(Nx,Ny) :: Speed10m
  real :: TserAvg

  integer ix,iy,iz
  integer nCat0, nCat1, nCat2, nCat3, nCat4, nCat5, NumPoints
  real Wspeed, SiMetric

  nCat0 = 0
  nCat1 = 0
  nCat2 = 0
  nCat3 = 0
  nCat4 = 0
  nCat5 = 0

  iz = 2 ! use next to bottom layer in the filter

  do iy = 1, Ny
    do ix = 1, Nx
      if (anint(Filter(ix,iy,iz)) .eq. 1.0) then
        ! Count up the number of grid points with wind speeds fitting each of the
        ! Saffir-Simpson categories. Then form the metric by weighting each category
        ! count.
        Wspeed = Speed10m(ix,iy)
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

  !Linear weighting
  NumPoints = nCat0 + nCat1 + nCat2 + nCat3 + nCat4 + nCat5
  SiMetric = float(nCat1) + (float(nCat2)*2.0) + (float(nCat3)*3.0) + (float(nCat4)*4.0) + (float(nCat5)*5.0)

  if (NumPoints .eq. 0) then
    TserAvg = 0.0
  else
    TserAvg = SiMetric / float(NumPoints)
  end if
  
  return
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

end program hdata_op
