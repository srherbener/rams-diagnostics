module revu_utils

!*********************************************************
! DATA
!*********************************************************

contains
!*********************************************************
! SUBROUTINES
!*********************************************************

!***********************************************************
! GetLonLat()
!
! This routine will read in longitude and latitude values
! from the given HDF5 file. The file organization is:
!
!   Longitude values:
!     Dataset name: x_coords
!     2D data, longitude values at every horizontal grid cell
!
!   Latitude values:
!     Dataset name: y_coords
!     2D data, latitude values at every horizontal grid cell
!
! Returns two 1D arrays contating representative longitude and
! latitude values. The first row in the longitude 2D array and
! the first column in the latitude 2D array will be used for
! the output 1D arrays
!
subroutine GetLonLat(Nx, Ny, Lon, Lat, Hfile)
  use rhdf5_utils
  implicit none

  integer :: Nx, Ny
  real, dimension(Nx) :: Lon
  real, dimension(Ny) :: Lat
  character (len=*) :: Hfile

  type (Rhdf5Var) :: Lon2D, Lat2D
  integer :: i, i2d

  ! Read in the longitude and latitude values
  Lon2D%vname = 'x_coords'
  call rhdf5_read_init(Hfile, Lon2D)
  allocate(Lon2D%vdata(Nx*Ny))
  call rhdf5_read(Hfile, Lon2D)
        
  Lat2D%vname = 'y_coords'
  call rhdf5_read_init(Hfile, Lat2D)
  allocate(Lat2D%vdata(Nx*Ny))
  call rhdf5_read(Hfile, Lat2D)

  ! Pick off the first row for Lon and first column
  ! for Lat
  do i = 1, Nx
    ! First row is located in the first n entries
    Lon(i) = Lon2D%vdata(i)
  enddo

  do i = 1, Ny
    ! First column is located starting with first
    ! entry and all subsequent entries spaced Nx apart
    i2d = ((i-1) * Nx) + 1
    Lat(i) = Lat2D%vdata(i2d)
  enddo

  ! Clean up
  deallocate(Lon2D%vdata)
  deallocate(Lat2D%vdata)

  return
end subroutine GetLonLat

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

end module revu_utils
