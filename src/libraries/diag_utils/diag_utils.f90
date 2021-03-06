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

subroutine AzimuthalAverage(Nx, Ny, Nz, FilterNz, NumRbands, NumBins, Avar, AzAvg, &
  RbandInc, Filter, Radius, HistBins, DoHist, UndefVal)

  implicit none

  integer :: Nx, Ny, Nz
  integer :: FilterNz
  integer :: NumRbands
  integer :: NumBins
  real, dimension(Nx,Ny,Nz) :: Avar
  real, dimension(NumRbands,NumBins,Nz) :: AzAvg
  real, dimension(Nx,Ny,FilterNz) :: Filter
  real, dimension(Nx,Ny) :: Radius
  real, dimension(NumBins) :: HistBins
  real :: RbandInc, UndefVal
  logical :: DoHist

  integer :: ix, iy, iz, ib
  integer :: filter_z
  real, dimension(NumRbands) :: Rcounts
  integer :: ir, iRband

  ! Mask the input data (Avar) with the Filter data, ie if the corresponding
  !
  ! The storm center is taken to be the min pressure location of the
  ! x-y plane on the surface (iz .eq. 1)

  ! zero out the output array
  do iz = 1, Nz
    do ib = 1, NumBins
      do ir = 1, NumRbands
        AzAvg(ir,ib,iz) = 0.0
      enddo
    enddo
  enddo

  do iz = 1, Nz
    ! Filter can be either 2D or 3D, use FilterNz to determine
    if (FilterNz .eq. 1) then
      filter_z = 1
    else
      filter_z = iz
    endif

    ! For the averaging
    if (.not. DoHist) then
      do ir = 1, NumRbands
        Rcounts(ir) = 0.0
      enddo
    endif

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

           if (DoHist) then
             ! Check all bins except the last.
             !    
             ! Exiting out of the loop when finding the bin will help a lot when the
             ! distribution of values is biased toward smaller values. After exiting
             ! out of the loop you can either just check the last bin (which will be wasted)
             ! or put in a logical variable and check that variable saying you can skip the
             ! check of the last bin. Since you would have to check the logical variable and
             ! the last bin every time you might as well just check the last bin instead.
             do ib = 1, NumBins-1 
                if ((HistBins(ib) .le. Avar(ix,iy,iz)) .and. (Avar(ix,iy,iz) .lt. HistBins(ib+1))) then 
                   Azavg(iRband,ib,iz) = Azavg(iRband,ib,iz) + 1.0
                   exit 
                endif
             enddo
             ! check the last bin
             if (HistBins(NumBins) .eq. Avar(ix,iy,iz)) then 
               Azavg(iRband,NumBins,iz) = Azavg(iRband,NumBins,iz) + 1.0
             endif
           else
             AzAvg(iRband,1,iz) = AzAvg(iRband,1,iz) + Avar(ix,iy,iz)
             Rcounts(iRband) = Rcounts(iRband) + 1.0
           endif
        end if
      end do
    end do

    ! If not doing histograms, then complete the average calculation.
    if (.not. DoHist) then
      do ir = 1, NumRbands
        ! If we didn't put anything into an AzAvg slot then set it
        ! to the undefined value so the data isn't biased by trying to chose
        ! a default value.
        if (Rcounts(ir) .ne. 0.0) then
          AzAvg(ir,1,iz) = AzAvg(ir,1,iz) / Rcounts(ir)
        else
          AzAvg(ir,1,iz) = UndefVal
        end if
      end do
    endif
  end do

  return
end subroutine

!*******************************************************************************
! ConvertHorizVelocity()
!
! This routine will convert the horizontal velocity vectors (described in U and V)
! into tangential or radial components given the storm center location.
!

subroutine ConvertHorizVelocity(Nx, Ny, Nz, U, V, StormX, StormY, Xcoords, Ycoords, Vt, Vr, uv2tr)
  implicit none

  integer :: Nx, Ny, Nz, uv2tr
  real, dimension(Nx,Ny,Nz) :: U, V, Vt, Vr
  real :: StormX, StormY
  real, dimension(Nx) :: Xcoords
  real, dimension(Ny) :: Ycoords

  integer :: ix, iy, iz
  real :: StmX, StmY         ! x,y location relative to the storm center
  real :: Theta, Phi, Alpha  ! angle values, in radians
  real :: WindMag, WindX, WindY

  do iz = 1, Nz
    do iy = 1, Ny
      StmY = Ycoords(iy) - StormY
      do ix = 1, Nx
        StmX = Xcoords(ix) - StormX

        if (uv2tr .eq. 1) then
          ! Convert (u,v) to (vt,vr)
          WindX = U(ix,iy,iz)
          WindY = V(ix,iy,iz)

          Theta = atan2(StmY, StmX) ! Angle of radius line from horizontal
          Phi = atan2(WindY, WindX) ! Angle of wind vector from horizontal
          Alpha = Phi - Theta       ! Angle of wind vector from raduis line
                                    !    radius line is the line from storm center through
                                    !    the point (StmX, StmY)

          WindMag = sqrt(WindX**2 + WindY**2)
        
          Vt(ix,iy,iz) = WindMag * sin(Alpha)
          Vr(ix,iy,iz) = WindMag * cos(Alpha)
        else
          ! Convert (vt,vr) to (u,v)
          WindX = Vr(ix,iy,iz)
          WindY = Vt(ix,iy,iz)

          Theta = atan2(StmY, StmX)    ! Angle of radius line from horizontal
          Alpha = atan2(WindY, WindX)  ! Angle of wind vector from raduis line
                                       !    radius line is the line from storm center through
                                       !    the point (StmX, StmY)
          Phi = Theta + Alpha          ! Angle of wind vector from horizontal

          WindMag = sqrt(WindX**2 + WindY**2)
        
          U(ix,iy,iz) = WindMag * cos(Phi)
          V(ix,iy,iz) = WindMag * sin(Phi)
        endif
      enddo
    enddo
  enddo

  
  return
end subroutine

!******************************************************************
! ConvertGridCoords()
!
! This routine will convert the longitude, latitude angle values
! in the input var to corresponding length values (km from the
! equator and prime meridian).
!
subroutine ConvertGridCoords(Nx, Ny, Lon, Lat, Xcoords, Ycoords)
  implicit none

  real, parameter :: RadiusEarth = 6378.1  ! km
  real, parameter :: PI = 3.14159

  integer :: Nx, Ny
  real, dimension(Nx) :: Lon
  real, dimension(Ny) :: Lat
  real, dimension(Nx) :: Xcoords
  real, dimension(Ny) :: Ycoords

  integer :: i
  real :: ConvDeg2Rad
  real :: AvgLat
  real :: RadiusAvgLat

  ! X corresponds to longitude
  ! Y corresponds to latitude
  !
  ! dX = [R * Cos(Lat)] * dLon
  ! dY = R * dLat
  !    Lon and Lat in radians
  !
  ! [R * Cos(Lat)] is the length from the Earth's axis to the latitude line
  ! being used to cacluate the longitude values. Assume that LON and LAT
  ! axes go through the middle of the horizontal domain. Also assume that
  ! delta-x and delta-y are uniform and equal. These assumptions are
  ! consistent with RAMS configuration. Therefore take the average of
  ! the first and last LAT entry for the latitude value used in
  ! the X calculation.
  !
  ! Define that delta lengths and angles are from the origin
  ! (intersection of equator and prime meridian) so that dX, dY
  ! dLon and dLat become equal to X, Y, Lon and Lat.

  ! convert to radians
  ConvDeg2Rad = (2.0 * PI) / 360.0

  ! average latitude for X calculation and radius corresponding to that latitude
  AvgLat = 0.5 * (Lat(1) + Lat(Ny))
  RadiusAvgLat = RadiusEarth * cos(AvgLat * ConvDeg2Rad)

  ! convert longitude to X
  do i = 1, Nx
    Xcoords(i) = RadiusAvgLat * Lon(i) * ConvDeg2Rad
  enddo

  ! convert latitude to Y
  do i = 1, Ny
    Ycoords(i) = RadiusEarth * Lat(i) * ConvDeg2Rad
  enddo

  return
end subroutine

!*****************************************************************************
! InsideCylVol()
!
! This function will see if a given grid location falls within a cylindrical
! volume relative to the storm center.
!

logical function InsideCylVol(Nx, Ny, Nz, Ix, Iy, Iz, MinR, MaxR, MinPhi, MaxPhi, &
  MinZ, MaxZ,  StmIx, StmIy, Xcoords, Ycoords, Zcoords, Radius, Phi, Z)
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
  StrStart = 0
  BetweenSep = 0
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
! TranslateField()
!
! This subroutine will subtract the given scalar from 
! all points in the given field.
!
! Handles either 2D (Nz = 1) or 3D fields.
!
subroutine TranslateField(Nx, Ny, Nz, Field, Scalar)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: Field
  real :: Scalar

  integer :: ix, iy, iz

  do iz = 1, Nz
    do iy = 1, Ny
      do ix = 1, Nx
        Field(ix,iy,iz) = Field(ix,iy,iz) - Scalar
      enddo
    enddo
  enddo

  return
end subroutine TranslateField


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

!**********************************************************************
! SetOutCoords()
!
! This routine will set the coordinate and dimension data 

subroutine SetOutCoords(Hfile, Xcoords, Ycoords, Zcoords, Tcoords)
  use rhdf5_utils
  implicit none 

  character (len=*) :: Hfile
  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords

  ! Read in longitude, latitude and height values
  Xcoords%vname = 'x_coords'
  call rhdf5_read_init(Hfile, Xcoords)
  call rhdf5_read(Hfile, Xcoords)
  
  Ycoords%vname = 'y_coords'
  call rhdf5_read_init(Hfile, Ycoords)
  call rhdf5_read(Hfile, Ycoords)

  Zcoords%vname = 'z_coords'
  call rhdf5_read_init(Hfile, Zcoords)
  call rhdf5_read(Hfile, Zcoords)

  Tcoords%vname = 't_coords'
  call rhdf5_read_init(Hfile, Tcoords)
  call rhdf5_read(Hfile, Tcoords)

  return
end subroutine SetOutCoords

!*******************************************************************
! Gsmooth2d()
!
! This routine will do Gaussian smoothing in two dimension. It assumes
! that the filter matrix is square (Npts X Npts). Because of this, it
! will split the filtering into two passes: one for horizontal the other
! for vertical.
!
subroutine Gsmooth2d(Nx, Ny, Npts, X, Sigma, Xsmooth)
  implicit none

  real, parameter :: PI = 3.14159265

  integer :: Nx, Ny, Npts
  real, dimension(Nx, Ny) :: X, Xsmooth
  real :: Sigma

  integer :: ix, iy, ig
  integer :: i_center, i_samp
  real, dimension(Npts) :: Gfilter

  real :: Gsum
  real, dimension(Nx, Ny) :: Xtemp

  ! Want to do 2D filtering, however the gaussian filter is seperable (into
  ! horizontal and vertical so it is possible to apply a horizontal filter
  ! first followed by a vertical filter. For simplicity make the 2D filter
  ! square so that the same linear filter can be applied to both horizontal
  ! and vertical directions.
  
  ! Build the linear filter
  ! For gaussian:
  !    g(x) = 1/(sqrt(2*pi) * sigma) * exp( - x**2 / (2 * sigma))
  !
  !    where x is the distance from the center of the array
  !
  ! Also, normalize the filter so that it will preserve the magnitude
  ! of the pressure.

  i_center = int(Npts / 2) + 1 ! Assumes that Npts is odd
  Gsum = 0.0
  do ig = 1, Npts
    Gfilter(ig) = 1./(sqrt(2. * PI) * Sigma) * exp(- (float(i_center-ig)**2) / (2. * (Sigma**2)))
    Gsum = Gsum + Gfilter(ig)
  enddo
  do ig = 1, Npts
    Gfilter(ig) = Gfilter(ig) / Gsum
  enddo

  ! At the edges, the filter will extend beyond the X array bounds. To handle this,
  ! just repeat the X value at the array edge.

  ! Horzontal filter - apply to X
  do ix = 1, Nx
    do iy = 1, Ny
      ! locate center of Gfilter over iy
      ! do a dot product between X (vector in y direction, across columns) and Gfilter
      Gsum = 0.0
      do ig = 1, Npts
        ! don't go past array boundaries
        i_samp = (iy + ig)  - i_center
        if (i_samp .lt. 1) then
           i_samp = 1
        else if (i_samp .gt. Ny) then
           i_samp = Ny
        endif
 
        Gsum = Gsum + (X(ix,i_samp) * Gfilter(ig))
      enddo

      Xtemp(ix,iy) = Gsum
    enddo
  enddo

  ! Vertical filter - apply to Xtemp
  do iy = 1, Ny
    do ix = 1, Nx
      ! locate center of Gfilter over ix
      ! do a dot product between X (vector in x direction, across rows) and Gfilter
      Gsum = 0.0
      do ig = 1, Npts
        ! don't go past array boundaries
        i_samp = (ix + ig)  - i_center
        if (i_samp .lt. 1) then
           i_samp = 1
        else if (i_samp .gt. Nx) then
           i_samp = Nx
        endif
 
        Gsum = Gsum + (Xtemp(i_samp,iy) * Gfilter(ig))
      enddo

      Xsmooth(ix,iy) = Gsum
    enddo
  enddo

  return
end subroutine Gsmooth2d

!*****************************************************************************
! ReadBinsFile()
!
! This subroutine will initialize a 1D array of bins. This is done by
! reading in a text file in the following format:
!
!   n
!   E1
!   E2
!   E3
!   ...
!   En
!
! First line contains number of edges (n), followed by n lines containing
! the n edges, in order, defining the bins (one edge value per line).
!
subroutine ReadBinsFile(BinsFile, Nbins, Bins)
  implicit none
  integer, parameter :: fun = 20

  integer :: Nbins
  character (len=*) :: BinsFile
  real, dimension(:), allocatable :: Bins

  integer :: i

  open(fun, file=BinsFile, status='old')

  read(fun, '(i15)') Nbins

  allocate(Bins(Nbins))
  do i = 1, Nbins
    read(fun, '(e15.7)') Bins(i)
  enddo

  close(fun)

  return
end subroutine ReadBinsFile

!*****************************************************************************
! FindBin()
!
! This function will find the bin that a given data value belongs to. Emulate
! the matlab binning where the bin values are treated like edges. Val belongs
! to a bin if:
!
!   Bins(ib) <= Val < Bins(i+1)
!
! except for the last bin where Val belongs to it if:
!
!   Bins(ib) == Val
!
integer function FindBin(Nb, Bins, Val)
  implicit none

  integer :: Nb
  real, dimension(Nb) :: Bins
  real :: Val

  integer :: ib

  ! if Val doesn't fall into any bins, FindBin will remain -1
  FindBin = -1
  do ib = 1, Nb
    if (ib .lt. Nb) then
      if ((Bins(ib) .le. Val) .and. (Val .lt. Bins(ib+1))) then
        FindBin = ib
        exit
      endif
    else
      !last bin
      if (Bins(ib) .eq. Val) then
        FindBin = ib
      endif
    endif
  enddo

  return
end function FindBin

!**********************************************************************
! FindMinPressure
!
! This routine will locate the minumum surface pressure that is 
! within the storm. There are several aids to help prevent this
! calculation from straying to false locations. Topography can
! create problems and other storm phenomena within the grid can
! also create problems. To handle these cases, the search for the
! point with minimum pressure is limited to cells that have elevation
! less than MaxElev, and fall within the box defined by MinLon, MaxLon,
! MinLat and MaxLat.
!

subroutine FindMinPressure(Nx, Ny, Press, SelectGrid, MinPressIx, MinPressIy, MinPress)
  implicit none

  integer :: Nx, Ny
  real, dimension(Nx,Ny) :: Press
  logical, dimension(Nx,Ny) :: SelectGrid
  integer :: MinPressIx, MinPressIy
  real :: MinPress

  integer :: i, j

  MinPress = 1e10 ! ridiculously large pressure
  MinPressIx = 0 
  MinPressIy = 0 

  do j = 1, Ny
    do i = 1, Nx
      if (SelectGrid(i,j) .and. (Press(i,j) .lt. MinPress)) then
        MinPress = Press(i,j)
        MinPressIx = i
        MinPressIy = j
      endif
    enddo
  enddo
  
  return
end subroutine FindMinPressure

!**********************************************************************
! CalcPolarCoords
!
! This routine will calculate polar coordinates [ radius (km), angle (radians) ]
! from the given center for each point in the horizontal domain.
!

subroutine CalcPolarCoords(Nx, Ny, Radius, Phi, X, Y, CenterIx, CenterIy)
  implicit none

  real, parameter :: PI = 3.141592654

  integer :: Nx, Ny
  real, dimension(Nx,Ny) :: Radius, Phi
  real, dimension(Nx) :: X
  real, dimension(Ny) :: Y
  integer :: CenterIx, CenterIy

  integer :: i, j
  real :: dx, dy

  do j = 1, Ny
    do i = 1, Nx
      dx = X(i) - X(CenterIx)
      dy = Y(j) - Y(CenterIy)

      Radius(i,j) = sqrt(dx*dx + dy*dy)

      ! atan2 returns a value in radians between -PI and +PI
      ! convert to value between 0 and 2*PI
      Phi(i,j) = atan2(dy, dx)
      if (Phi(i,j) .lt. 0.0) then
        Phi(i,j) = Phi(i,j) + (2.0 * PI)
      endif
      
    enddo
  enddo

  return
end subroutine CalcPolarCoords

!**********************************************************************
! FindPressureCent
!
! This routine will calculate the pressure centroid within the given
! radius.
!

subroutine FindPressureCent(Nx, Ny, Press, SelectGrid, Radius, X, Y, PressRad, PressCentIx, PressCentIy)
  implicit none

  integer :: Nx, Ny
  real, dimension(Nx,Ny) :: Press, Radius
  logical, dimension(Nx,Ny) :: SelectGrid
  real, dimension(Nx) :: X
  real, dimension(Ny) :: Y
  real :: PressRad
  integer :: PressCentIx, PressCentIy

  integer :: i, j, n
  real :: Penv, Pdiff, SumXpdiff, SumYpdiff, SumPdiff
  real :: Xcent, Ycent
  real :: DeltaX

  DeltaX = X(2) - X(1)  ! assume grid spacing is same in x and y
  
  ! Find the environmental pressure which is defined as the
  ! average pressure located at PressRad distance from the
  ! miniumum pressure location (Radius holds these distances)
  n = 0
  Penv = 0.0
  do j = 1, Ny
    do i = 1, Nx
      ! sample ~ 2*DeltaX wide band of samples for env. pressure
      ! from points in selection grid
      if (SelectGrid(i,j) .and. (abs(Radius(i,j)-PressRad) .lt. DeltaX)) then
        n = n + 1
        Penv = Penv + Press(i,j)
      endif
    enddo
  enddo

  if (n .gt. 0) then
    Penv = Penv / float(n)
  else
    print'(a)', 'ERROR: unable to calculate environmental pressure, no points selected'
    stop
  endif

  ! Calculate the centroid
  !   Pdiff = P - Penv
  !
  !    x_cent = sum(x_i * Pdiff_i) / sum(Pdiff_i)
  !    y_cent = sum(y_i * Pdiff_i) / sum(Pdiff_i)
  !
  SumPdiff = 0.0
  SumXpdiff = 0.0
  SumYpdiff = 0.0
  do j = 1, Ny
    do i = 1, Nx
      ! sample points within PressRad distance from min pressure
      ! from selection grid
      if (SelectGrid(i,j) .and. (Radius(i,j) .le. PressRad)) then
        Pdiff = Penv - Press(i,j)
        SumPdiff = SumPdiff + Pdiff
        SumXpdiff = SumXpdiff + (X(i) * Pdiff)
        SumYpdiff = SumYpdiff + (Y(j) * Pdiff)
      endif
    enddo
  enddo

  Xcent = SumXpdiff / SumPdiff
  Ycent = SumYpdiff / SumPdiff

  ! convert the centroid to the nearest grid cell x and y
  call FindCoordMatch(Nx, X, Xcent, PressCentIx)
  call FindCoordMatch(Ny, Y, Ycent, PressCentIy)

  return
end subroutine FindPressureCent

!**********************************************************************
! FindCoordMatch
!
! This routine will find which coord value that the given value
! is closest to, and return the corresponding index into the coord array.
!

subroutine FindCoordMatch(Nc, Coords, Val, Cindex)
  implicit none

  integer :: Nc, Cindex
  real, dimension(Nc) :: Coords
  real :: Val

  integer :: i
  real :: Diff, MinDiff

  Cindex = 1
  MinDiff = abs(Val - Coords(1))
  do i = 2, Nc
    Diff = abs(Val - Coords(i))
    if (Diff .lt. MinDiff) then
      Cindex = i
      MinDiff = Diff
    endif
  enddo

  return
end subroutine FindCoordMatch

!**********************************************************************
! SetSelectGrid
!
! This routine will fill in the selection grid (2D) according to the
! tests that check for each grid cells elevation and location within
! the horizontal domain.
!
subroutine SetSelectGrid(Nx, Ny, Topo, Lon, Lat, MaxElev, MinLat, MaxLat, MinLon, MaxLon, SelectGrid)
  implicit none

  integer :: Nx, Ny
  real, dimension(Nx,Ny) :: Topo
  real, dimension(Nx) :: Lon
  real, dimension(Ny) :: Lat
  real :: MaxElev, MinLat, MaxLat, MinLon, MaxLon
  logical, dimension(Nx,Ny) :: SelectGrid

  integer :: i, j

  do j = 1, Ny
    do i = 1, Nx
      ! Select if elevation less than given limit
      SelectGrid(i,j) = Topo(i,j) .lt. MaxElev

      ! Select if location is withing given limits (lat, lon)
      SelectGrid(i,j) = SelectGrid(i,j) .and. (Lon(i) .ge. MinLon)
      SelectGrid(i,j) = SelectGrid(i,j) .and. (Lon(i) .le. MaxLon)
      SelectGrid(i,j) = SelectGrid(i,j) .and. (Lat(j) .ge. MinLat)
      SelectGrid(i,j) = SelectGrid(i,j) .and. (Lat(j) .le. MaxLat)
    enddo
  enddo

  return
end subroutine SetSelectGrid

!**********************************************************************
! CalcStormSpeed
!
! This routine will calculate storm speed (keeping the x and y components
! separated) from the storm location and time vectors.
!
subroutine CalcStormSpeed(Nt, StormX, StormY, Tcoords, SpeedX, SpeedY)
  implicit none

  integer :: Nt
  real, dimension(Nt) :: StormX, StormY, Tcoords, SpeedX, SpeedY

  integer :: it
  real, dimension(Nt-1) :: IntSpeedX, IntSpeedY
  real :: dt, dx, dy

  ! Form intermediate (between time points) speeds.
  do it = 1, Nt-1
    dt = Tcoords(it+1) - Tcoords(it)           ! seconds
    dx = (StormX(it+1) - StormX(it)) * 1000.0  ! meters
    dy = (StormY(it+1) - StormY(it)) * 1000.0  ! meters
    
    IntSpeedX(it) = dx / dt;  ! m/s
    IntSpeedY(it) = dy / dt;  ! m/s
  enddo

  ! Average the intermediate speeds back onto the time points.
  !
  ! Assign the adjacent intermediate speeds at the end points
  ! to the correspoinding end point since there are no "external"
  ! points to use in the average.
  SpeedX(1)  = IntSpeedX(1)
  SpeedX(Nt) = IntSpeedX(Nt-1)

  SpeedY(1)  = IntSpeedY(1)
  SpeedY(Nt) = IntSpeedY(Nt-1)

  do it = 2, Nt-1
    SpeedX(it) = 0.5 * (IntSpeedX(it) + IntSpeedX(it-1))
    SpeedY(it) = 0.5 * (IntSpeedY(it) + IntSpeedY(it-1))
  enddo

  return
end subroutine CalcStormSpeed

!**********************************************************************
! CopyHorizData
!
! This routine will pull out the horizontal slice defined by Ilev from
! Var3d and copy this into Var2d
!
subroutine CopyHorizData(Nx, Ny, Nz, Var3d, Var2d, Ilev)
  implicit none

  integer :: Nx, Ny, Nz, Ilev
  real, dimension(Nx,Ny,Nz) :: Var3d
  real, dimension(Nx,Ny) :: Var2d

  integer :: ix, iy

  do iy = 1, Ny
    do ix = 1, Nx
      Var2d(ix,iy) = Var3d(ix,iy,Ilev)
    enddo
  enddo

  return
end subroutine CopyHorizData

end module diag_utils
