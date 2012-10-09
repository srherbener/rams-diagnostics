module diag_utils
!*********************************************************
! DATA
!*********************************************************

contains
!*********************************************************
! SUBROUTINES
!*********************************************************

!***************************************************************
! Utilities for the azimuthial averaging code.
!***************************************************************

!*************************************************************************
! AzimuthalAverage
!
! This routine will do the azimuthal averaging. It will create a 4-D array
! for its output (to keep GRADS write routine consistent/general) even though it
! really is 3-D data (radial band, z level). It will put the radial bands 
! in the x dimension and make the y dimension 1 long.
!
! This routine will step though each z value and compute the azimuthal
! average of the x-y plane at each of these points. Then it will store that
! in the output array so that you get azimuthal averaging for every original
! z level.
!

subroutine AzimuthalAverage(Nx, Ny, Nz, AvarNz, NumRbands, Avar, AzAvg, &
  RbandInc, Filter, Radius, UndefVal)

  implicit none

  integer :: Nx, Ny, Nz
  integer :: AvarNz
  integer :: NumRbands
  real, dimension(Nx,Ny,AvarNz) :: Avar
  real, dimension(NumRbands,AvarNz) :: AzAvg
  real, dimension(Nx,Ny,Nz) :: Filter
  real, dimension(Nx,Ny) :: Radius
  real :: RbandInc, UndefVal

  integer :: ix, iy, iz
  integer :: filter_z
  real, dimension(NumRbands) :: Rcounts
  integer :: ir, iRband

  ! Mask the input data (Avar) with the Filter data, ie if the corresponding
  !
  ! The storm center is taken to be the min pressure location of the
  ! x-y plane on the surface (iz .eq. 1)

  do iz = 1, AvarNz
    ! The filter data always comes in as a 3D var, and since RAMS has
    ! the first z level below the surface, the first z level in the filter
    ! tends to be all zeros (since it is typically desired to trim this
    ! level off from diagnostics). When a 2D var comes in, AvarNz will
    ! be equal to one which is the below surface level in the filter
    ! data. To work around this, when AvarNz is one (2d var), then force
    ! this code to look at the second level in the filter (first level
    ! above surface). Note that it's up to the caller to make sure
    ! that the 2nd level in the filter (first level above sfc) have
    ! good data in it.
    if (AvarNz .eq. 1) then
      filter_z = 2
    else
      filter_z = iz
    endif

    ! For the averaging
    do ir = 1, NumRbands
      Rcounts(ir) = 0.0
      AzAvg(ir,iz) = 0.0
    end do

    ! Get the averages for this x-y plane
    do iy = 1, Ny
      do ix = 1, Nx
        ! Only keep the points where Filter is equal to 1
        if ((anint(Filter(ix,iy,filter_z)) .eq. 1.0) .and. (Avar(ix,iy,iz) .ne. UndefVal)) then
           iRband = int(Radius(ix,iy) / RbandInc) + 1
           ! iRband will go from 1 to n, but n might extend beyond the last radial
           ! band due to approximations made in deriving RbandInc
           if (iRband .gt. NumRbands) then
             iRband = NumRbands
           end if

           AzAvg(iRband,iz) = AzAvg(iRband,iz) + Avar(ix,iy,iz)
           Rcounts(iRband) = Rcounts(iRband) + 1.0
        end if
      end do
    end do

    do ir = 1, NumRbands
      ! If we didn't put anything into an AzAvg slot then set it
      ! to the undefined value so the data isn't biased by trying to chose
      ! a default value.
      if (Rcounts(ir) .ne. 0.0) then
        AzAvg(ir,iz) = AzAvg(ir,iz) / Rcounts(ir)
      else
        AzAvg(ir,iz) = UndefVal
      end if
    end do
  end do

  return
end subroutine


!*****************************************************************************
! BuildCfad()
!
! This routine will sort data points out into radial bands and then do 
! the histogram binning within each radial band.
!

subroutine BuildCfad(Nx, Ny, Nz, Nt, NumRbands, NumBins, Avar, Cfad, &
   RbandInc, BinVals, Filter, Radius, UndefVal)

  implicit none

  integer :: Nx, Ny, Nz, Nt, NumRbands, NumBins
  real, dimension(Nx,Ny,Nz,Nt) :: Avar, Filter
  real, dimension(Nx,Ny,Nt) :: Radius
  real, dimension(NumRbands,NumBins,Nz,Nt) :: Cfad
  real, dimension(NumBins) :: BinVals
  real :: RbandInc, UndefVal

  integer :: ix, iy, iz, it, ib, ir
  real :: DeltaX, DeltaY, BinTotal
  integer :: iRband, iBin
  real, dimension(NumBins-1) :: BinEdges

  ! initialize the counts
  do ir = 1, NumRbands
    do ib = 1, NumBins
      do iz = 1, Nz
        do it = 1, Nt
          Cfad(ir,ib,iz,it) = 0.0;
        end do 
      end do 
    end do 
  end do 

  ! Create the points halfway in between the points given in BinVals. Use
  ! these half-way points as the boundaries between the bins.
  do ib = 1, NumBins-1
    BinEdges(ib) = (BinVals(ib) + BinVals(ib+1)) / 2.0
  enddo

  ! create the CFAD data
  write (*,*) 'Binning Data:'

  do it = 1, Nt
    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          ! Throw out undefined data points, and those screened out by Filter
          if ((anint(Filter(ix,iy,iz,it)) .eq. 1.0) .and. (Avar(ix,iy,iz,it) .ne. UndefVal)) then
             ! find the radial band for this data point
             iRband = int(Radius(ix,iy,it) / RbandInc) + 1
             ! iRband will go from 1 to n, but n might extend beyond the last radial
             ! band due to approximations made in deriving RbandInc
             if (iRband .gt. NumRbands) then
               iRband = NumRbands
             end if

             ! Find the bin for this data point. Use the values in BinEdges to decide which
             ! bin a value belongs in. Note that BinEdges(1) has the boundary between bins
             ! 1 and 2, BinEdges(2) has the boundary between bins 2 and 3, etc.
             iBin = 0
             do ib = 1, NumBins
               if (ib .lt. NumBins) then
                 if (Avar(ix,iy,iz,it) .lt. BinEdges(ib)) then
                   iBin = ib
                   exit
                 endif
               else
                 iBin = ib
               endif
             enddo

             if (iBin .eq. 0) then
               write (*,*) 'ERROR: did not find a bin for data at location:'
               write (*,*) 'ERROR:   ix, iy, iz, it: ', ix, iy, iz, it
               stop
             endif

             Cfad(iRband,iBin,iz,it) = Cfad(iRband,iBin,iz,it) + 1.0
          endif
        end do
      end do
    end do
  end do

  ! Have counts in Cfad now. Convert these to % numbers.

  do it = 1, Nt
    do iz = 1, Nz
      do ir = 1, NumRbands
        ! sum up counts in all the bins
        BinTotal = 0.0
        do ib = 1, NumBins
          BinTotal = BinTotal + Cfad(ir,ib,iz,it)
        end do
        ! convert entries to per cent values
        do ib = 1, NumBins
          ! If we selected out all data points, then the bin total will
          ! be zero. In this case, leave the zeros in the Cfad array.
          if (BinTotal .ne. 0.0) then
            Cfad(ir,ib,iz,it) = (Cfad(ir,ib,iz,it) / BinTotal) * 100.0
          endif
        end do
      end do
    end do
  end do

  write (*,*) ''
  return
end subroutine

!*******************************************************************************
! ConvertHorizVelocity()
!
! This routine will convert the horizontal velocity vectors (described in U and V)
! into tangential or radial components given the storm center location.
!

subroutine ConvertHorizVelocity(Nx, Ny, Nz, U, V, StormX, StormY, Xcoords, Ycoords, Avar, DoTangential)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: U, V, Avar
  real :: StormX, StormY
  real, dimension(Nx) :: Xcoords
  real, dimension(Ny) :: Ycoords
  logical :: DoTangential

  integer :: ix, iy, iz
  real :: StmX, StmY         ! x,y location relative to the storm center
  real :: Theta, Phi, Alpha  ! angle values, in radians
  real :: WindMag, WindX, WindY
  
  do iz = 1, Nz
    do iy = 1, Ny
      StmY = Ycoords(iy) - StormY
      do ix = 1, Nx
        StmX = Xcoords(ix) - StormX

        WindX = U(ix,iy,iz)
        WindY = V(ix,iy,iz)

        Theta = atan2(StmY, StmX) ! Angle of radius line from horizontal
        Phi = atan2(WindY, WindX) ! Angle of wind vector from horizontal
        Alpha = Phi - Theta       ! Angle of wind vector from raduis line
                                  !    radius line is the line from storm center through
                                  !    the point (StmX, StmY)

        WindMag = sqrt(WindX**2 + WindY**2)
        
        if (DoTangential) then
          ! Tangential wind
          Avar(ix,iy,iz) = WindMag * sin(Alpha)
        else
          ! Radial wind
          Avar(ix,iy,iz) = WindMag * cos(Alpha)
        end if
      end do
    end do
  end do
  
  return
end subroutine

!******************************************************************
! ConvertGridCoords()
!
! This routine will convert the longitude, latitude angle values
! in the input var to a flat plane (x and y length values)
!
subroutine ConvertGridCoords(Nx, Ny, Nz, Lon, Lat, Xcoords, Ycoords)
  implicit none

  real, parameter :: RadiusEarth = 6378.1  ! km
  real, parameter :: PI = 3.14159

  integer :: Nx, Ny, Nz
  real, dimension(Nx) :: Lon
  real, dimension(Ny) :: Lat
  real, dimension(Nx) :: Xcoords
  real, dimension(Ny) :: Ycoords

  integer :: ix, iy
  integer :: Ierror
  real :: ConvDeg2Rad;
  real :: DeltaX, DeltaY
  real :: DeltaT, DeltaP  ! Theta is longitude angle, Phi is latitude angle, radians
  real :: RadiusX
  real :: PhiN, Phi1, Phi2, Theta1, Theta2

  ! The Xcoords in the Gvar structure are longitude angles in degrees
  ! The Ycoords in the Gvar structure are latitude angles in degrees
  !
  ! To convert, figure out what DeltaX and DeltaY are according to the longitude, latitude
  ! angles. Call the lower left point of the grid (0,0) and just add in the delta values to
  ! get the remaining coordinate values. Put the x,y values in km.
  !
  ! Technically, the DeltaX values change for each unique latitude angle since DeltaX depends
  ! on the arm perpendicular to the axis of rotation of the Earth (not the center of Earth).
  ! However, since the storms are near the equator this arm doesn't change much in length.
  ! Approximate the spherical surface with a rectangle using the delta x that is an average
  ! of the deltax you find at the southernmost latidute and northernmost latitude. This will
  ! greatly simplify the code by allowing the use of just one list of x coordinates for the 
  ! entire grid.
  !
  !   DeltaX = (RadiusEarth * cos(Phi)) * DeltaTheta
  !   DeltaY = RadiusEarth * DeltaPhi
  !   angle values are in radians

  ! write a warning if the grid is located far away from the equator
  if ((abs(Lat(Ny)) .gt. 23.0) .or. (abs(Lat(1)) .gt. 23.0)) then
    write (*,*) 'Warning: extent of grid goes outside tropics, this code assumes grid is near equator'
    write (*,*) ''
  end if

  ! convert to radians
  ConvDeg2Rad = (2.0 * PI) / 360.0
  Theta1 = Lon(1) * ConvDeg2Rad
  Theta2 = Lon(2) * ConvDeg2Rad

  Phi1 = Lat(1) * ConvDeg2Rad
  Phi2 = Lat(2) * ConvDeg2Rad
  PhiN = Lat(Ny) * ConvDeg2Rad

  DeltaT = Theta2 - Theta1
  DeltaP = Phi2 - Phi1

  ! average of the arms at south lat and north lat
  RadiusX = (RadiusEarth * 0.5) * (cos(Phi1) + cos(PhiN))
  DeltaX = RadiusX * DeltaT
  DeltaY = RadiusEarth * DeltaP

  ! DeltaX and DeltaY are in units of RadiusEarth (km for now)
  do ix = 1, Nx
    Xcoords(ix) = real(ix - 1) * DeltaX
  end do

  do iy = 1, Ny
    Ycoords(iy) = real(iy - 1) * DeltaY
  end do

  return
end subroutine

!*****************************************************************************
! InsideCylVol()
!
! This function will see if a given grid location falls within a cylindrical
! volume relative to the storm center.
!

logical function InsideCylVol(Nx, Ny, Nz, Ix, Iy, Iz, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ,  StmIx, StmIy, Xcoords, Ycoords, Zcoords, Radius, Phi, Z)
  implicit none

  real, parameter :: PI = 3.141592654

  integer :: Nx, Ny, Nz, Ix, Iy, Iz
  real :: MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  integer :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real, dimension(1:Nz) :: Zcoords
  real :: Radius, Phi, Z

  real :: dX, dY

  dX = Xcoords(Ix) - Xcoords(StmIx)
  dY = Ycoords(Iy) - Ycoords(StmIy)
  Radius = sqrt(dX*dX + dY*dY)
  ! atan2 returns a value in radians between -PI and +PI
  ! convert to value between 0 and 2*PI
  Phi = atan2(dY, dX)
  if (Phi .lt. 0.0) then
    Phi = Phi + (2.0 * PI)
  end if
  Z = Zcoords(Iz)

  InsideCylVol = .false.
  if (((Radius .ge. MinR) .and. (Radius .le. MaxR)) .and. &
      ((Phi .ge. MinPhi) .and. (Phi .le. MaxPhi)) .and. &
      ((Z .ge. MinZ) .and. (Z .le. MaxZ))) then
    InsideCylVol = .true.
  end if
  return
end function

!***********************************************************
! String2List()
!
! This routine will split the input string (colon separated list)
! into the individual strings (between the colons) and load these
! into the output array.
!
! Args
!  1. input list of files, colon separated items (strings)
!  2. character used for separator
!  3. output array
!  4. maximum number of input items
!  5. actual number of input items
!  6. description of what the items are

subroutine String2List(InString, Separator, OutList, MaxItems, NumItems, ItemType)
  implicit none

  character (len=*) InString
  character Separator
  character (len=*), dimension(*) :: OutList
  character (len=*) ItemType
  integer MaxItems, NumItems

  integer i, StrStart, StrEnd
  integer BetweenSep, PrevBetweenSep

  ! walk through the string from start to end
  ! when a colon is encountered, store the string
  ! from the "start" up to the colon in the output array

  NumItems = 0
  PrevBetweenSep = 0  ! treat the beginning of the line as separator
  i = 1
  do
    if (i .gt. len_trim(InString)) exit

    ! Record if looking at the string between the separators
    ! of if looking at the separators.
    ! Remember the previous state of BetweenSep:
    !   If you just changed from separators to string, then this
    !   is the beginning of the string
    !   If you just changed from string to separators, then this
    !   is the end of the string --> then record it in the list

    if (InString(i:i) .eq. Separator) then
      BetweenSep = 0
    else
      BetweenSep = 1
    end if

    if ((PrevBetweenSep .eq. 0) .and. (BetweenSep .eq. 1)) then
      ! Start of next string
      StrStart = i
    else
      if ((PrevBetweenSep .eq. 1) .and. (BetweenSep .eq. 0)) then
        ! End of string, record it in the list
        NumItems = NumItems + 1
        StrEnd = i - 1
        if (NumItems .le. MaxItems) then
          OutList(NumItems) = InString(StrStart:StrEnd)
        else
          write (0,*) 'ERROR: Maximum number of ', ItemType, ' (', MaxItems, ') exceeded'
          stop
        end if
      end if
    end if

    i = i + 1
    PrevBetweenSep = BetweenSep
  end do

  ! may have one more string at this point
  ! treat the end of the line as a separator

  if (BetweenSep .eq. 1) then
    NumItems = NumItems + 1
    StrEnd = i - 1
    if (NumItems .le. MaxItems) then
      OutList(NumItems) = InString(StrStart:StrEnd)
    else
      write (0,*) 'ERROR: Maximum number of ', ItemType, ' (', MaxItems, ') exceeded'
      stop
    end if
  end if

  return
end subroutine

!***********************************************************************
! DimsMatch()
!
! This function will return true/false according to whether or not the
! dimensions of the two given GRADS variables match. When the variables
! are 3D then x, y, z, and t dimensions are checked. If one of the variables
! is 2D then only x, y, and t dimensions are checked.

logical function DimsMatch(Var1, Var2)
  use rhdf5_utils
  implicit none

  type (Rhdf5Var) :: Var1, Var2
  logical :: Is2D
  integer :: ix, iy, it
  integer :: jx, jy, jt

  ! Check to see if the x, y and t dimensions match
  ! For 3D variables:
  !    ndims: 4
  !    dim order: x, y, z, t
  ! For 2D variables:
  !    ndims: 3
  !    dim order: x, y, t

  if (Var1%ndims .eq. 4) then
    ix = 1
    iy = 2
    it = 4 
  else
    ix = 1
    iy = 2
    it = 3 
  endif

  if (Var2%ndims .eq. 4) then
    jx = 1
    jy = 2
    jt = 4 
  else
    jx = 1
    jy = 2
    jt = 3 
  endif

  DimsMatch = ((Var1%dims(ix) .eq. Var2%dims(jx)) .and. &
               (Var1%dims(iy) .eq. Var2%dims(jy)) .and. &
               (Var1%dims(it) .eq. Var2%dims(jt)))

  return
end function DimsMatch

!******************************************************************
! MultiDimLookup()
!
! This function will take the 1D array of real and look up a value
! as if that array were either 2D or 3D field.
!
real function MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d, Var3d)
  implicit none

  integer :: Nx, Ny, Nz, ix, iy, iz
  real, dimension(Nx,Ny), optional :: Var2d
  real, dimension(Nx,Ny,Nz), optional :: Var3d

  if (present(Var2d)) then
    MultiDimLookup = Var2d(ix,iy)
  else if (present(Var3d)) then
    MultiDimLookup = Var3d(ix,iy,iz)
  else
    print*, 'ERROR: MultiDimLookup: must use one of the optional arguments: Var2d or Var3d'
    stop 'MultiDimLookup: missing argument'
  endif

  return
end function MultiDimLookup

!******************************************************************
! MultiDimAssign()
!
! This subroutine will take the 1D array of real and look up a value
! as if that array were either 2D or 3D field.
!
subroutine MultiDimAssign(Nx, Ny, Nz, ix, iy, iz, Val, Var2d, Var3d)
  implicit none

  integer :: Nx, Ny, Nz, ix, iy, iz
  real :: Val
  real, dimension(Nx,Ny), optional :: Var2d
  real, dimension(Nx,Ny,Nz), optional :: Var3d

  if (present(Var2d)) then
    Var2d(ix,iy) = Val
  else if (present(Var3d)) then
    Var3d(ix,iy,iz) = Val
  else
    print*, 'ERROR: MultiDimAssign: must use one of the optional arguments: Var2d or Var3d'
    stop 'MultiDimAssign: missing argument'
  endif

  return
end subroutine MultiDimAssign

!******************************************************************
! CreateFilter()
!
! This subroutine will take the list of filter files and perform
! a logical intersection of the filters contained within these files.
!
subroutine CreateFilter(Nfiles, FilterFiles, Filter)
  use rhdf5_utils
  implicit none

  integer :: Nfiles
  character (len=*), dimension(Nfiles) :: FilterFiles
  type (Rhdf5Var) :: Filter, CurFilter

  integer :: i, ielem, Nelem

  ! read the first file into Filter directly, then if more files are in
  ! the list, do a logical intersection with these and what is in Filter
  Filter%vname = 'filter'
  call rhdf5_read_init(FilterFiles(1), Filter)
  call rhdf5_read(FilterFiles(1), Filter)

  ! figure out how many elements exist in Fitler and assume that
  ! all the filter files have the same number
  Nelem = 1
  do i = 1, Filter%ndims
    Nelem = Nelem * Filter%dims(i)
  enddo

  do i = 2, Nfiles
    CurFilter%vname = 'filter'
    call rhdf5_read_init(FilterFiles(i), CurFilter)
    call rhdf5_read(FilterFiles(i), CurFilter)

    ! filters are all 1's and 0's so use multiplication to do the logical intersection
    do ielem = 1, Nelem
      Filter%vdata(ielem) = Filter%vdata(ielem) * CurFilter%vdata(ielem)
    enddo

    deallocate(CurFilter%vdata)
  enddo

  return
end subroutine CreateFilter

!******************************************************************
! ConvertStormCenter()
!
! This subroutine will take the Lat/Lon values for the storm center
! and convert them to km values according to the lon/x and lat/y
! data.
!
! This routine is not very efficient so it relies on Nx, Ny
! remaining small.
!
subroutine ConvertStormCenter(Nx, Ny, Lon, Xcoords, StormLon, StormX, &
      Lat, Ycoords, StormLat, StormY)
  implicit none

  integer :: Nx, Ny
  real, dimension(Nx) :: Lon, Xcoords
  real, dimension(Ny) :: Lat, Ycoords
  real :: StormLon, StormX, StormLat, StormY

  integer :: ix, iy

  StormX = Xcoords(FindIndex(Nx, Lon, StormLon))
  StormY = Ycoords(FindIndex(Ny, Lat, StormLat))

  return
end subroutine ConvertStormCenter

!*****************************************************************
! FindIndex()
!
! This function returns the index in the Coords array where Val
! is located. This function assumes that Val exists in Coords.
!
integer function FindIndex(N, Coords, Val)
  implicit none

  real, parameter :: CloseEnough = 10.0e-4

  integer :: N
  real, dimension(N) :: Coords
  real :: Val

  integer :: i

  FindIndex = 0
  do i = 1, N
    if (abs(Coords(i) - Val) .le. CloseEnough) then
      FindIndex = i
      exit
    endif
  enddo

  return
end function FindIndex

end module diag_utils
