!***************************************************************
! Program to apply filters to input data and output "mask" data
! that shows where points exist that were selected by the filtering.
!
! Args
!   1. directory containing input files
!   2. suffix to tag onto the end of input file names
!   3. output file name 
!   4. selection of averaging function
!
! Output
!   The output will be an hdf5 file containing '1's where data was
!   selected and '0's else where. The output will always be a 3D
!   field so that it can be applied downstream to either 2D or 3D fields.
!

program diag_filter
  use rhdf5_utils
  use diag_utils
  implicit none

  integer, parameter :: LargeString  = 512
  integer, parameter :: MediumString = 256
  integer, parameter :: LittleString = 128

  integer, parameter :: MaxArgFields = 20

  ! Ftype codes for the FilterDescrip type
  ! In the main loop below that is applying the filters
  ! it is advantageous to use integers (representing the filter
  ! type) versus strings. The conditional expressions in the if
  ! statements will be much faster using integers than strings.
  ! The inner statements of the loop doing the filter application
  ! can get executed a great number of times (billions) so it will
  ! make a big difference in performance using the integers.
  integer, parameter :: FTYPE_NONE      =  0
  integer, parameter :: FTYPE_CYLVOL    =  1
  integer, parameter :: FTYPE_GE        =  2
  integer, parameter :: FTYPE_GT        =  3
  integer, parameter :: FTYPE_LE        =  4
  integer, parameter :: FTYPE_LT        =  5
  integer, parameter :: FTYPE_RANGE     =  6
  integer, parameter :: FTYPE_ABS_RANGE =  7
  integer, parameter :: FTYPE_UP        =  8
  integer, parameter :: FTYPE_DN        =  9
  integer, parameter :: FTYPE_UP_DN     = 10

  integer, parameter :: MaxFilters = 10

  type FileSpec
    character (len=RHDF5_MAX_STRING) :: fname
    character (len=RHDF5_MAX_STRING) :: vname
  endtype

  type FilterDescrip
    character (len=LittleString) :: Fname
    integer :: Ftype
    type (FileSpec) :: Fvar
    real :: x1
    real :: x2
    real :: y1
    real :: y2
    real :: z1
    real :: z2
  end type FilterDescrip

  type ArgList
    character (len=LittleString) :: CombSense
    type (FileSpec) :: Output
    type (FilterDescrip), dimension(MaxFilters) :: Filters
    integer :: Nfilters
  endtype

  type (ArgList) :: Args

  integer :: i, icylvol
  integer :: imodel, OutNdims, FilterNdims
  real :: TempZ
  integer :: ix, iy, iz, it
  integer :: ih2d, ih3d
  integer :: Nx, Ny, Nz, Nt
  integer :: NhElems
  logical :: BadDims

  type (Rhdf5Var), dimension(MaxFilters) :: Vars

  character (len=RHDF5_MAX_STRING) :: FileAcc
  integer, dimension(MaxFilters) :: InFileIds
  integer :: OutFileId

  type (Rhdf5Var) :: OutFilter
  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords

  real :: DeltaX, DeltaY
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm

  type (Rhdf5Var) :: MinP, Radius, Phi, StormXindex, StormYindex
  integer :: StmIx, StmIy

  logical :: FilterVal, SelectThisPoint

  logical :: AndFilters

  logical, dimension(:,:), allocatable :: UpDnDraftMask
  logical :: DoingUpDrafts
  logical :: DoingDnDrafts
  logical :: DoingUpDnDrafts
  integer :: UDfnum

  ! Get the command line arguments
  call GetMyArgs(Args)

  ! Record which combination sense for the filters was selected
  ! GetMyArgs already checked to make sure one of 'and' or 'or' was selected for CombSense
  AndFilters = (Args%CombSense .eq. 'and')

  DoingUpDrafts = .false.
  DoingDnDrafts = .false.
  DoingUpDnDrafts = .false.

  icylvol = 0 ! if remains zero -> signifies cylvol not being used to code downstream

  write (*,*) 'Creating HDF5 data filter:'
  write (*,*) '  Output file:  ', trim(Args%Output%fname)
  write (*,*) '    Output variable name:  ', trim(Args%Output%vname)
  write (*,*) '  Data selection specs: '
  if (AndFilters) then
    write (*,*) '    Filters will be logically and-ed together'
  else
    write (*,*) '    Filters will be logically or-ed together'
  endif
  UDfnum = 0
  do i = 1, Args%Nfilters
    if (Args%Filters(i)%Ftype .eq. FTYPE_CYLVOL) then
      icylvol = i
      write (*,*) '    Cylindrical volume:'
      write (*,*) '      Storm center variable: ', trim(Args%Filters(i)%Fvar%vname)
      write (*,*) '      Storm center file: ', trim(Args%Filters(i)%Fvar%fname)
      write (*,*) '      Minimum radius: ', Args%Filters(i)%x1
      write (*,*) '      Maximum radius: ', Args%Filters(i)%x2
      write (*,*) '      Minimum angle: ', Args%Filters(i)%y1
      write (*,*) '      Maximum angle: ', Args%Filters(i)%y2
      write (*,*) '      Minimum height: ', Args%Filters(i)%z1
      write (*,*) '      Maximum height: ', Args%Filters(i)%z2
    else if (Args%Filters(i)%Ftype .eq. FTYPE_GE) then
      write (*,*) '    Greater than or equal:'
      write (*,*) '      Variable: ', trim(Args%Filters(i)%Fvar%vname)
      write (*,*) '      Variable file: ', trim(Args%Filters(i)%Fvar%fname)
      write (*,*) '      Threshold: ', Args%Filters(i)%x1
    else if (Args%Filters(i)%Ftype .eq. FTYPE_GT) then
      write (*,*) '    Greater than:'
      write (*,*) '      Variable: ', trim(Args%Filters(i)%Fvar%vname)
      write (*,*) '      Variable file: ', trim(Args%Filters(i)%Fvar%fname)
      write (*,*) '      Threshold: ', Args%Filters(i)%x1
    else if (Args%Filters(i)%Ftype .eq. FTYPE_LE) then
      write (*,*) '    Less than or equal:'
      write (*,*) '      Variable: ', trim(Args%Filters(i)%Fvar%vname)
      write (*,*) '      Variable file: ', trim(Args%Filters(i)%Fvar%fname)
      write (*,*) '      Threshold: ', Args%Filters(i)%x1
    else if (Args%Filters(i)%Ftype .eq. FTYPE_LT) then
      write (*,*) '    Less than:'
      write (*,*) '      Variable: ', trim(Args%Filters(i)%Fvar%vname)
      write (*,*) '      Variable file: ', trim(Args%Filters(i)%Fvar%fname)
      write (*,*) '      Threshold: ', Args%Filters(i)%x1
    else if (Args%Filters(i)%Ftype .eq. FTYPE_RANGE) then
      write (*,*) '    Inside range:'
      write (*,*) '      Variable: ', trim(Args%Filters(i)%Fvar%vname)
      write (*,*) '      Variable file: ', trim(Args%Filters(i)%Fvar%fname)
      write (*,*) '      Min: ', Args%Filters(i)%x1
      write (*,*) '      Max: ', Args%Filters(i)%x2
    else if (Args%Filters(i)%Ftype .eq. FTYPE_ABS_RANGE) then
      write (*,*) '    Absolute value inside range:'
      write (*,*) '      Variable: ', trim(Args%Filters(i)%Fvar%vname)
      write (*,*) '      Variable file: ', trim(Args%Filters(i)%Fvar%fname)
      write (*,*) '      Min: ', Args%Filters(i)%x1
      write (*,*) '      Max: ', Args%Filters(i)%x2
    else if (Args%Filters(i)%Ftype .eq. FTYPE_UP) then
      DoingUpDrafts = .true.
      UDfnum = i
      write (*,*) '    Updraft:'
      write (*,*) '      Variable: ', trim(Args%Filters(i)%Fvar%vname)
      write (*,*) '      Variable file: ', trim(Args%Filters(i)%Fvar%fname)
      write (*,*) '      Threshold: ', Args%Filters(i)%x1
    else if (Args%Filters(i)%Ftype .eq. FTYPE_DN) then
      DoingDnDrafts = .true.
      UDfnum = i
      write (*,*) '    Downdraft:'
      write (*,*) '      Variable: ', trim(Args%Filters(i)%Fvar%vname)
      write (*,*) '      Variable file: ', trim(Args%Filters(i)%Fvar%fname)
      write (*,*) '      Threshold: ', Args%Filters(i)%x1
    else if (Args%Filters(i)%Ftype .eq. FTYPE_UP_DN) then
      DoingUpDnDrafts = .true.
      UDfnum = i
      write (*,*) '    Updraft/Downdraft:'
      write (*,*) '      Variable: ', trim(Args%Filters(i)%Fvar%vname)
      write (*,*) '      Variable file: ', trim(Args%Filters(i)%Fvar%fname)
      write (*,*) '      Updraft Threshold: ', Args%Filters(i)%x1
      write (*,*) '      Downdraft Threshold: ', Args%Filters(i)%x2
    endif
  enddo
  write (*,*) ''
  flush(6)

  ! See if we have access to all of the required variables. Check to make sure the horizontal
  ! and time dimensions match between all vars. Compare against the model file.
  !
  ! Since we will be processing one time step at a time, the effective number of dimensions
  ! on each of the input is decreased by one. We want to eliminate the time dimension which
  ! is always the last one. This works out conveniently since all we need to do is decrement
  ! the number of dimensions by one.

  BadDims = .false.
  OutNdims = 0

  ! The number of dims for the output is determined by the maximum number of dims
  ! of all the input variables. Note that 'up', 'dn', 'up_dn' and 'cylvol' produce
  ! 2D output regarless of dimensions of the associate filter variable.
  do i = 1, Args%Nfilters
    if (i .eq. icylvol) then
      ! put radius information into Vars(i) for the dimension size check
      Vars(i)%vname = 'radius'

      ! set up the extra vars
      ! go ahead and read in the entire timeseries vars
      MinP%vname = 'min_press'
      call rhdf5_read_init(Args%Filters(i)%Fvar%fname, MinP)
      call rhdf5_read(Args%Filters(i)%Fvar%fname, MinP);

      Radius%vname = 'radius'
      call rhdf5_read_init(Args%Filters(i)%Fvar%fname, Radius)

      Phi%vname = 'phi'
      call rhdf5_read_init(Args%Filters(i)%Fvar%fname, Phi)

      StormXindex%vname = trim(Args%Filters(i)%Fvar%vname) // '_x_index'
      call rhdf5_read_init(Args%Filters(i)%Fvar%fname, StormXindex)
      call rhdf5_read(Args%Filters(i)%Fvar%fname, StormXindex);

      StormYindex%vname = trim(Args%Filters(i)%Fvar%vname) // '_y_index'
      call rhdf5_read_init(Args%Filters(i)%Fvar%fname, StormYindex)
      call rhdf5_read(Args%Filters(i)%Fvar%fname, StormYindex);
    else
      Vars(i)%vname = trim(Args%Filters(i)%Fvar%vname)
    endif

    call rhdf5_read_init(Args%Filters(i)%Fvar%fname, Vars(i))

    ! Check that the horizontal and time dimensions match with the first var. 
    if (i .gt. 1) then
      if (.not. DimsMatch(Vars(1), Vars(i))) then
        write (*,*) 'ERROR: horizontal and time dimensions of variable do not match the first filter variable: ', trim(Vars(i)%vname)
        BadDims = .true.
      endif
    endif
  enddo

  imodel = 1
  do i = 1, Args%Nfilters
    ! Prepare for reading
    Vars(i)%ndims = Vars(i)%ndims - 1

    ! Adjust the cylvol extra vars dimensions, but only for the 2D
    ! variables. The time series variables have already been read in their entirety
    ! so there is no need to condition them for reading one time step at a time.
    if (i .eq. icylvol) then
       Radius%ndims = Radius%ndims - 1
       Phi%ndims = Phi%ndims - 1
    endif

    ! Record var with maximum number of dimensions --> output dims
    ! NOTE: Vars(i) for cylvol has data corresponding to the 'radius' variable (2D)
    if ((Args%Filters(i)%Ftype .eq. FTYPE_UP) .or. (Args%Filters(i)%Ftype .eq. FTYPE_DN) .or. &
        (Args%Filters(i)%Ftype .eq. FTYPE_UP_DN) .or. (Args%Filters(i)%Ftype .eq. FTYPE_CYLVOL)) then
      FilterNdims = 2
    else
      FilterNdims = Vars(i)%ndims
    endif

    if (FilterNdims .gt. OutNdims) then
      imodel = i
      OutNdims = FilterNdims

      ! record dimension sizes
      Nx = Vars(i)%dims(1)
      Ny = Vars(i)%dims(2)
      if (Vars(i)%ndims .eq. 2) then
        Nz = 1
        Nt = Vars(i)%dims(3)
      else
        Nz = Vars(i)%dims(3)
        Nt = Vars(i)%dims(4)
      endif
    endif
  enddo

  if (BadDims) then
    stop
  endif

  ! Set dimensions of output
  !
  ! It is possible that a filter type that forces 2D output is using a 3D variable
  ! ('up' uses variable 'w' for example). In this case change the Zcoords to a
  ! size of 1 using the value in Zcoords(2) (first level above surface) as the
  ! coordinate value.
  call SetOutCoords(Args%Filters(imodel)%Fvar%fname, Xcoords, Ycoords, Zcoords, Tcoords)
  if ((OutNdims .eq. 2) .and. (Vars(imodel)%ndims .eq. 3)) then
    TempZ = Zcoords%vdata(2)

    deallocate(Zcoords%vdata)
    allocate(Zcoords%vdata(1))

    Zcoords%dims(1) = 1
    Zcoords%vdata(1) = TempZ
    Nz = 1
  endif

  ! Set the output dimensions and coordinates to those of the selected input var
  
  ! Convert lat (x coords) and lon (y coords) to distances in km
  allocate(XcoordsKm(Nx))
  allocate(YcoordsKm(Ny))
  call ConvertGridCoords(Nx, Ny, Xcoords%vdata, Ycoords%vdata, XcoordsKm, YcoordsKm)

  DeltaX = (XcoordsKm(2) - XcoordsKm(1)) * 1000.0
  DeltaY = (YcoordsKm(2) - YcoordsKm(1)) * 1000.0

  write (*,*) 'Horizontal grid info:'
  write (*,*) '  X range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', Xcoords%vdata(1), Xcoords%vdata(Nx), XcoordsKm(1), XcoordsKm(Nx)
  write (*,*) '  Y range (min lat, max lat) --> (min y, max y): '
  write (*,*) '    ', Ycoords%vdata(1), Ycoords%vdata(Ny), YcoordsKm(1), YcoordsKm(Ny)
  write (*,*) ''
  write (*,*) '  Delta x of domain:     ', DeltaX
  write (*,*) '  Delta y of domain:     ', DeltaY
  write (*,*) ''
  write (*,*) 'Vertical grid info:'
  do iz = 1, Nz
    write (*,*) '  ', iz, ' --> ', Zcoords%vdata(iz)
  end do
  write (*,*) ''
  flush(6)

  ! Stats
  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', Nx
  write (*,*) '  Number of y (latitude) points:           ', Ny
  write (*,*) '  Number of z (vertical level) points:     ', Nz
  write (*,*) '  Number of t (time) points:               ', Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''
  write (*,*) '  Grid delta x: ', DeltaX
  write (*,*) '  Grid delta y: ', DeltaY
  write (*,*) ''
  flush(6)

  ! Prepare the output filter, only set up for the data field since
  ! the time dimension will be built as we go through the loop below.
  OutFilter%vname = 'filter'
  OutFilter%ndims = 3
  OutFilter%dims(1) = Nx
  OutFilter%dims(2) = Ny
  OutFilter%dims(3) = Nz
  OutFilter%dimnames(1) = 'x'
  OutFilter%dimnames(2) = 'y'
  OutFilter%dimnames(3) = 'z'
  OutFilter%units = 'logical'
  OutFilter%descrip = 'data selection locations'
  allocate(OutFilter%vdata(Nx*Ny*Nz))

  ! Check up and down draft filtering
  if (DoingUpDrafts .and. DoingDnDrafts) then
    write (*,*) 'ERROR: cannot use both updraft and downdraft filtering'
    stop
  else if (DoingUpDrafts .and. DoingUpDnDrafts) then
    write (*,*) 'ERROR: cannot use both updraft and updraft/downdraft filtering'
    stop
  else if (DoingDnDrafts .and. DoingUpDnDrafts) then
    write (*,*) 'ERROR: cannot use both downdraft and updraft/downdraft filtering'
    stop
  endif

  if (DoingUpDrafts .or. DoingDnDrafts .or. DoingUpDnDrafts) then
    allocate(UpDnDraftMask(Nx,Ny))
  endif

  ! Open the input files and the output file
  FileAcc = 'R'
  do i = 1, Args%Nfilters
    call rhdf5_open_file(Args%Filters(i)%Fvar%fname, FileAcc, 0, InFileIds(i))
    write (*,*) 'Reading HDF5 file: ', trim(Args%Filters(i)%Fvar%fname)
  enddo
  write (*,*) ''

  FileAcc = 'W'
  call rhdf5_open_file(Args%Output%fname, FileAcc, 1, OutFileId)
  write (*,*) 'Writing HDF5 file: ', trim(Args%Output%fname)
  write (*,*) ''

  ! Do the filtering one time step at a time.
  !
  ! The necessary input files have been opened, data buffers allocated,
  ! plus the time dimensions have been "removed" from the descriptions
  ! of the input variables.
  !

  NhElems = Nx * Ny
  do it = 1, Nt
    ! read the input vars
    do i = 1, Args%Nfilters
      if (i .eq. icylvol) then
        call rhdf5_read_variable(InFileIds(i), Radius%vname, Radius%ndims, it, Radius%dims, rdata=Radius%vdata)
        call rhdf5_read_variable(InFileIds(i), Phi%vname, Phi%ndims, it, Phi%dims, rdata=Phi%vdata)

        StmIx = nint(StormXindex%vdata(it))
        StmIy = nint(StormYindex%vdata(it))
      else
        call rhdf5_read_variable(InFileIds(i), Vars(i)%vname, Vars(i)%ndims, it, Vars(i)%dims, rdata=Vars(i)%vdata)
      endif
    enddo

    ! If doing up- or down-draft filtering, build the mask
    ! Note that last argument is a code telling subrouting which mask to build:
    !    1 - updraft
    !    2 - downdraft
    !    3 - up and down drafts
    if (DoingUpDrafts) then
      call BuildUpDnMask(Nx, Ny, Vars(UDfnum)%dims(3), Vars(UDfnum)%vdata, UpDnDraftMask, Args%Filters(UDfnum)%x1, Args%Filters(UDfnum)%x2, 1)
    elseif (DoingDnDrafts) then
      call BuildUpDnMask(Nx, Ny, Vars(UDfnum)%dims(3), Vars(UDfnum)%vdata, UpDnDraftMask, Args%Filters(UDfnum)%x1, Args%Filters(UDfnum)%x2, 2)
    elseif (DoingUpDnDrafts) then
      call BuildUpDnMask(Nx, Ny, Vars(UDfnum)%dims(3), Vars(UDfnum)%vdata, UpDnDraftMask, Args%Filters(UDfnum)%x1, Args%Filters(UDfnum)%x2, 3)
    endif

    ! do the filter selection
    !
    ! For execution speed, want to walk thorugh Vars in the same order as
    ! contiguous memory (which will minimize cache misses). Since the array
    ! is column major, this means the looping from outer to inner should be
    ! done as: z, y, x.
    ! 
    ! Calling subroutines (MultiDimLookup, MultiDimAssign) for every point in domain seems wasteful
    ! in terms of execution speed, so keep track of the proper linear index for 3D and 2D variables.
    ! Use these index variables to reference linear arrays in the variables.
    !
    ! WARNING: this algorithm depends on the looping order going from outside to inside to be: z, y, x

    ih3d = 0
    do iz = 1, Nz
      ih2d = 0
      do iy = 1, Ny
        do ix = 1, Nx
          ih3d = ih3d + 1
          ih2d = ih2d + 1

          ! Either and-ing the filter values or or-ing them. The logical var AndFilters is
          ! true if and-ing and false if or-ing. Want SelectThisPoint to be initialized
          ! to true if and-ing and false if or-ing (ie, to the value in AndFilters).
          SelectThisPoint = AndFilters

          do i = 1, Args%Nfilters
            ! Apply the current filter
            if (Args%Filters(i)%Ftype .eq. FTYPE_CYLVOL) then
              ! select if inside the cylindrical volume
              FilterVal =                 (Radius%vdata(ih2d) .ge. Args%Filters(i)%x1)  ! x1 -> min radius
              FilterVal = FilterVal .and. (Radius%vdata(ih2d) .le. Args%Filters(i)%x2)  ! x2 -> max radius

              FilterVal = FilterVal .and. (Phi%vdata(ih2d) .ge. Args%Filters(i)%y1)     ! y1 -> min phi
              FilterVal = FilterVal .and. (Phi%vdata(ih2d) .le. Args%Filters(i)%y2)     ! y2 -> max phi

              FilterVal = FilterVal .and. (Zcoords%vdata(iz) .ge. Args%Filters(i)%z1)   ! z1 -> min height
              FilterVal = FilterVal .and. (Zcoords%vdata(iz) .le. Args%Filters(i)%z2)   ! z2 -> max height

            else if (Args%Filters(i)%Ftype .eq. FTYPE_GT) then
              ! select if var > threshold
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (Vars(i)%vdata(ih2d) .gt. Args%Filters(i)%x1)
              else
                FilterVal = (Vars(i)%vdata(ih3d) .gt. Args%Filters(i)%x1)
              endif

            else if (Args%Filters(i)%Ftype .eq. FTYPE_GE) then
              ! select if var >= threshold
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (Vars(i)%vdata(ih2d) .ge. Args%Filters(i)%x1)
              else
                FilterVal = (Vars(i)%vdata(ih3d) .ge. Args%Filters(i)%x1)
              endif

            else if (Args%Filters(i)%Ftype .eq. FTYPE_LT) then
              ! select if var < threshold
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (Vars(i)%vdata(ih2d) .lt. Args%Filters(i)%x1)
              else
                FilterVal = (Vars(i)%vdata(ih3d) .lt. Args%Filters(i)%x1)
              endif

            else if (Args%Filters(i)%Ftype .eq. FTYPE_LE) then
              ! select if var <= threshold
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (Vars(i)%vdata(ih2d) .le. Args%Filters(i)%x1)
              else
                FilterVal = (Vars(i)%vdata(ih3d) .le. Args%Filters(i)%x1)
              endif

            else if (Args%Filters(i)%Ftype .eq. FTYPE_RANGE) then
              ! select if min <= var <= max
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (Vars(i)%vdata(ih2d) .ge. Args%Filters(i)%x1) .and. (Vars(i)%vdata(ih2d) .le. Args%Filters(i)%x2)
              else
                FilterVal = (Vars(i)%vdata(ih3d) .ge. Args%Filters(i)%x1) .and. (Vars(i)%vdata(ih3d) .le. Args%Filters(i)%x2)
              endif

            else if (Args%Filters(i)%Ftype .eq. FTYPE_ABS_RANGE) then
              ! select if min <= abs(var) <= max
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (abs(Vars(i)%vdata(ih2d)) .ge. Args%Filters(i)%x1) .and. (abs(Vars(i)%vdata(ih2d)) .le. Args%Filters(i)%x2)
              else
                FilterVal = (abs(Vars(i)%vdata(ih3d)) .ge. Args%Filters(i)%x1) .and. (abs(Vars(i)%vdata(ih3d)) .le. Args%Filters(i)%x2)
              endif

            else if ((Args%Filters(i)%Ftype .eq. FTYPE_UP) .or. (Args%Filters(i)%Ftype .eq. FTYPE_DN) .or. (Args%Filters(i)%Ftype .eq. FTYPE_UP_DN)) then
              ! select according to up down draft mask
              FilterVal = UpDnDraftMask(ix,iy)
            endif

            ! Combine the current filter value with the overall selection according to CombSense
            if (AndFilters) then
              SelectThisPoint = SelectThisPoint .and. FilterVal
            else
              ! or-ing the filters
              SelectThisPoint = SelectThisPoint .or. FilterVal
            endif
          enddo

          if (SelectThisPoint) then
            OutFilter%vdata(ih3d) = 1.0
          else
            OutFilter%vdata(ih3d) = 0.0
          endif
        enddo
      enddo
    enddo

    ! Write the filter data, and storm info if doing cylvol selection, to the output file
    call rhdf5_write_variable(OutFileId, OutFilter%vname, OutFilter%ndims, it, OutFilter%dims, &
      OutFilter%units, OutFilter%descrip, OutFilter%dimnames, rdata=OutFilter%vdata)

    ! cleanup
    do i = 1, Args%Nfilters
      if (i .eq. icylvol) then
        deallocate(Radius%vdata)
        deallocate(Phi%vdata)
      else
        deallocate(Vars(i)%vdata)
      endif
    enddo

    ! Write out status to screen every 50 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,50) .eq. 0) then
      write (*,*) 'Working: Timestep: ', it

      if (icylvol .gt. 0) then
        write (*,'(a,i3,a,i3,a,g15.2,a,g15.2,a)') '    Storm Center: (', StmIx, ', ', StmIy, &
         ') --> (', XcoordsKm(StmIx), ', ', YcoordsKm(StmIy), ')'
        write (*,*) '   Minumum pressure: ', MinP%vdata(it)
      endif

      write (*,*) ''
    endif
  enddo

  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''

  ! write out the coordinate data
  call rhdf5_write(Args%Output%fname, Xcoords, 1)
  call rhdf5_write(Args%Output%fname, Ycoords, 1)
  call rhdf5_write(Args%Output%fname, Zcoords, 1)
  call rhdf5_write(Args%Output%fname, Tcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(Args%Output%fname, Xcoords, 'x')
  call rhdf5_set_dimension(Args%Output%fname, Ycoords, 'y')
  call rhdf5_set_dimension(Args%Output%fname, Zcoords, 'z')
  call rhdf5_set_dimension(Args%Output%fname, Tcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(Args%Output%fname, OutFilter)

  ! cleanup
  call rhdf5_close_file(OutFileId)
  do i = 1, Args%Nfilters
    call rhdf5_close_file(InFileIds(i))
  enddo

  stop

contains
!**********************************************************************
! Subroutines go here. Want these contained within the main program
! so that the interfaces to these subroutine are 'explicit' (ie, like
! ANSI C external declarations). This is needed to get the passing
! of allocatable arrays working properly.
!**********************************************************************

!**********************************************************************
! GetMyArgs()
!

 subroutine GetMyArgs(Args)
  use getoptions

  implicit none

  type (ArgList) :: Args

  character :: optval
  logical :: BadArgs
  character (len=MediumString), dimension(MaxArgFields) :: ArgList
  integer :: Nfields
  integer :: ifilt

  ! Default values
  Args%CombSense = 'and'

  ! Initialization
  Args%Output%fname = 'none'
  Args%Output%vname = 'none'

  do ifilt = 1, MaxFilters
    Args%Filters(i)%Fvar%fname = 'none'
    Args%Filters(i)%Fvar%vname = 'none'
  enddo 

  BadArgs = .false.

  ifilt = 0

  ! loop through all of the command line tokens
  ! optarg is a character string variable that the getoptions module supplies
  !   optarg gets set to the argument value for an option that uses an argument
  ! getopt returns single character:
  !    '>' finished
  !    '!' unknown option
  !    '.' command line argument (no option)
  do
    optval = getopt('c:')

    select case (optval)
      case ('>') ! finished
        exit

      case ('!') ! unrecognized argument
        write(*,*) 'ERROR: unknow option: ', trim(optarg)
        write(*,*) ''
        BadArgs = .true.

      case ('c')
        Args%CombSense = optarg
        if ((Args%CombSense .ne. 'and') .and. (Args%CombSense .ne. 'or')) then
          write (*,*) 'ERROR: <comb_sense> must be one of: "and", "or"'
          BadArgs = .true.
        endif

      case ('.')  ! file and filter specs
        call String2List(optarg, ':', ArgList, MaxArgFields, Nfields, 'file/filter spec')
        select case (ArgList(1))
          case ('out')
            Args%Output%fname = trim(ArgList(2))
            Args%Output%vname = trim(ArgList(3))

          case ('filter')
            ifilt = ifilt + 1
            select case (ArgList(2))
              case ('cylvol')
                if (Nfields .lt. 10) then
                  write (*,*) 'ERROR: not enough arguments to fully specify the cylvol filter'
                  BadArgs = .true.
                else
                  ! have enough args
                  Args%Filters(ifilt)%Fname      = ArgList(2)
                  Args%Filters(ifilt)%Ftype      = FTYPE_CYLVOL
                  Args%Filters(ifilt)%Fvar%vname = ArgList(3)
                  Args%Filters(ifilt)%Fvar%fname = ArgList(4)
        
                  read(ArgList(5), '(f15.7)') Args%Filters(ifilt)%x1
                  read(ArgList(6), '(f15.7)') Args%Filters(ifilt)%x2
                  read(ArgList(7), '(f15.7)') Args%Filters(ifilt)%y1
                  read(ArgList(8), '(f15.7)') Args%Filters(ifilt)%y2
                  read(ArgList(9), '(f15.7)') Args%Filters(ifilt)%z1
                  read(ArgList(10), '(f15.7)') Args%Filters(ifilt)%z2
                endif
        
                if ((Args%Filters(ifilt)%x1 .lt. 0.0) .or. (Args%Filters(ifilt)%x2 .lt. 0.0) .or. (Args%Filters(ifilt)%x2 .le. Args%Filters(ifilt)%x1)) then
                  write (*,*) 'ERROR: <x1> and <x2> must be >= 0.0, and <x2> must be > <x1>'
                  write (*,*) ''
                  BadArgs = .true.
                endif
                if ((Args%Filters(ifilt)%y1 .lt. 0.0) .or. (Args%Filters(ifilt)%y2 .lt. 0.0) .or. (Args%Filters(ifilt)%y2 .le. Args%Filters(ifilt)%y1)) then
                  write (*,*) 'ERROR: <y1> and <y2> must be >= 0.0, and <y2> must be > <y1>'
                  write (*,*) ''
                  BadArgs = .true.
                endif
                if ((Args%Filters(ifilt)%z1 .lt. 0.0) .or. (Args%Filters(ifilt)%z2 .lt. 0.0) .or. (Args%Filters(ifilt)%z2 .le. Args%Filters(ifilt)%z1)) then
                  write (*,*) 'ERROR: <z1> and <z2> must be >= 0.0, and <z2> must be > <z1>'
                  write (*,*) ''
                  BadArgs = .true.
                endif

              case ('ge')
                if (Nfields .lt. 5) then
                  write (*,*) 'ERROR: not enough arguments to fully specify the ge filter'
                  BadArgs = .true.
                else
                  ! have enough args
                  Args%Filters(ifilt)%Fname      = ArgList(2)
                  Args%Filters(ifilt)%Ftype      = FTYPE_GE
                  Args%Filters(ifilt)%Fvar%vname = ArgList(3)
                  Args%Filters(ifilt)%Fvar%fname = ArgList(4)
        
                  read(ArgList(5),  '(f15.7)') Args%Filters(ifilt)%x1
                  Args%Filters(ifilt)%x2 = 0
                  Args%Filters(ifilt)%y1 = 0
                  Args%Filters(ifilt)%y2 = 0
                  Args%Filters(ifilt)%z1 = 0
                  Args%Filters(ifilt)%z2 = 0
                endif

              case ('gt')
                if (Nfields .lt. 5) then
                  write (*,*) 'ERROR: not enough arguments to fully specify the gt filter'
                  BadArgs = .true.
                else
                  ! have enough args
                  Args%Filters(ifilt)%Fname      = ArgList(2)
                  Args%Filters(ifilt)%Ftype      = FTYPE_GT
                  Args%Filters(ifilt)%Fvar%vname = ArgList(3)
                  Args%Filters(ifilt)%Fvar%fname = ArgList(4)
        
                  read(ArgList(5),  '(f15.7)') Args%Filters(ifilt)%x1
                  Args%Filters(ifilt)%x2 = 0
                  Args%Filters(ifilt)%y1 = 0
                  Args%Filters(ifilt)%y2 = 0
                  Args%Filters(ifilt)%z1 = 0
                  Args%Filters(ifilt)%z2 = 0
                endif

              case ('le')
                if (Nfields .lt. 5) then
                  write (*,*) 'ERROR: not enough arguments to fully specify the le filter'
                  BadArgs = .true.
                else
                  ! have enough args
                  Args%Filters(ifilt)%Fname      = ArgList(2)
                  Args%Filters(ifilt)%Ftype      = FTYPE_LE
                  Args%Filters(ifilt)%Fvar%vname = ArgList(3)
                  Args%Filters(ifilt)%Fvar%fname = ArgList(4)
        
                  read(ArgList(5),  '(f15.7)') Args%Filters(ifilt)%x1
                  Args%Filters(ifilt)%x2 = 0
                  Args%Filters(ifilt)%y1 = 0
                  Args%Filters(ifilt)%y2 = 0
                  Args%Filters(ifilt)%z1 = 0
                  Args%Filters(ifilt)%z2 = 0
                endif

              case ('lt')
                if (Nfields .lt. 5) then
                  write (*,*) 'ERROR: not enough arguments to fully specify the lt filter'
                  BadArgs = .true.
                else
                  ! have enough args
                  Args%Filters(ifilt)%Fname      = ArgList(2)
                  Args%Filters(ifilt)%Ftype      = FTYPE_LT
                  Args%Filters(ifilt)%Fvar%vname = ArgList(3)
                  Args%Filters(ifilt)%Fvar%fname = ArgList(4)
        
                  read(ArgList(5),  '(f15.7)') Args%Filters(ifilt)%x1
                  Args%Filters(ifilt)%x2 = 0
                  Args%Filters(ifilt)%y1 = 0
                  Args%Filters(ifilt)%y2 = 0
                  Args%Filters(ifilt)%z1 = 0
                  Args%Filters(ifilt)%z2 = 0
                endif

              case ('range')
                if (Nfields .lt. 6) then
                  write (*,*) 'ERROR: not enough arguments to fully specify the range filter'
                  BadArgs = .true.
                else
                  ! have enough args
                  Args%Filters(ifilt)%Fname      = ArgList(2)
                  Args%Filters(ifilt)%Ftype      = FTYPE_RANGE
                  Args%Filters(ifilt)%Fvar%vname = ArgList(3)
                  Args%Filters(ifilt)%Fvar%fname = ArgList(4)
        
                  read(ArgList(5),  '(f15.7)') Args%Filters(ifilt)%x1
                  read(ArgList(6),  '(f15.7)') Args%Filters(ifilt)%x2
                  Args%Filters(ifilt)%y1 = 0
                  Args%Filters(ifilt)%y2 = 0
                  Args%Filters(ifilt)%z1 = 0
                  Args%Filters(ifilt)%z2 = 0
                endif

              case ('abs_range')
                if (Nfields .lt. 6) then
                  write (*,*) 'ERROR: not enough arguments to fully specify the abs_range filter'
                  BadArgs = .true.
                else
                  ! have enough args
                  Args%Filters(ifilt)%Fname      = ArgList(2)
                  Args%Filters(ifilt)%Ftype      = FTYPE_ABS_RANGE
                  Args%Filters(ifilt)%Fvar%vname = ArgList(3)
                  Args%Filters(ifilt)%Fvar%fname = ArgList(4)
        
                  read(ArgList(5),  '(f15.7)') Args%Filters(ifilt)%x1
                  read(ArgList(6),  '(f15.7)') Args%Filters(ifilt)%x2
                  Args%Filters(ifilt)%y1 = 0
                  Args%Filters(ifilt)%y2 = 0
                  Args%Filters(ifilt)%z1 = 0
                  Args%Filters(ifilt)%z2 = 0
                endif

              case ('up')
                if (Nfields .lt. 5) then
                  write (*,*) 'ERROR: not enough arguments to fully specify the updraft filter'
                  BadArgs = .true.
                else
                  ! have enough args
                  Args%Filters(ifilt)%Fname      = ArgList(2)
                  Args%Filters(ifilt)%Ftype      = FTYPE_UP
                  Args%Filters(ifilt)%Fvar%vname = ArgList(3)
                  Args%Filters(ifilt)%Fvar%fname = ArgList(4)
        
                  read(ArgList(5), '(f15.7)') Args%Filters(ifilt)%x1
                  Args%Filters(ifilt)%x2 = 0
                  Args%Filters(ifilt)%y1 = 0
                  Args%Filters(ifilt)%y2 = 0
                  Args%Filters(ifilt)%z1 = 0
                  Args%Filters(ifilt)%z2 = 0
                endif

              case ('dn')
                if (Nfields .lt. 5) then
                  write (*,*) 'ERROR: not enough arguments to fully specify the downdraft filter'
                  BadArgs = .true.
                else
                  ! have enough args
                  Args%Filters(ifilt)%Fname      = ArgList(2)
                  Args%Filters(ifilt)%Ftype      = FTYPE_DN
                  Args%Filters(ifilt)%Fvar%vname = ArgList(3)
                  Args%Filters(ifilt)%Fvar%fname = ArgList(4)
        
                  read(ArgList(5), '(f15.7)') Args%Filters(ifilt)%x1
                  Args%Filters(ifilt)%x2 = 0
                  Args%Filters(ifilt)%y1 = 0
                  Args%Filters(ifilt)%y2 = 0
                  Args%Filters(ifilt)%z1 = 0
                  Args%Filters(ifilt)%z2 = 0
                endif

              case ('up_dn')
                if (Nfields .lt. 6) then
                  write (*,*) 'ERROR: not enough arguments to fully specify the updraft/downdraft filter'
                  BadArgs = .true.
                else
                  ! have enough args
                  Args%Filters(ifilt)%Fname      = ArgList(2)
                  Args%Filters(ifilt)%Ftype      = FTYPE_UP_DN
                  Args%Filters(ifilt)%Fvar%vname = ArgList(3)
                  Args%Filters(ifilt)%Fvar%fname = ArgList(4)
        
                  read(ArgList(5), '(f15.7)') Args%Filters(ifilt)%x1
                  read(ArgList(6), '(f15.7)') Args%Filters(ifilt)%x2
                  Args%Filters(ifilt)%y1 = 0
                  Args%Filters(ifilt)%y2 = 0
                  Args%Filters(ifilt)%z1 = 0
                  Args%Filters(ifilt)%z2 = 0
                endif

              case default
                write (*,*) 'ERROR: <ftype>, ', trim(ArgList(2)), ', must be one of:'
                write (*,*) '          cylvol'
                write (*,*) '          ge'
                write (*,*) '          gt'
                write (*,*) '          le'
                write (*,*) '          lt'
                write (*,*) '          range'
                write (*,*) '          abs_range'
                write (*,*) '          up'
                write (*,*) '          dn'
                write (*,*) '          up_dn'
                write (*,*) ''
                BadArgs = .true.
            endselect

          case default
            write(*,*) 'ERROR: unknown spec type: ', trim(ArgList(1))
            write(*,*) 'ERROR:   must use one of "out" or "filter"'
            write(*,*) ''
            BadArgs = .true.
        endselect

    endselect
  enddo

  Args%Nfilters = ifilt

  if (Args%Nfilters .eq. 0) then
    write (*,*) 'ERROR: must specify at least one <filter_spec>'
    write (*,*) ''
    BadArgs = .true.
  endif

  if (BadArgs) then
    write (*,*) 'USAGE: diag_filter [-c <cobm_sense>]  <out_file> <filter_spec> [<filter_spec>...]'
    write (*,*) '        -c: combine filters using <comb_sense> logical operator'
    write (*,*) '            <comb_sense> must be one of "and", "or"'
    write (*,*) '              and: multiple filter specs are and-ed together (default)' 
    write (*,*) '              or: multiple filter specs are or-ed together' 
    write (*,*) '        <filter_spec>: filter:<ftype>:<vname>:<vfile>:<x1>:<x2>:<y1>:<y2>:<z1>:<z2>'
    write (*,*) ''
    write (*,*) '            <ftype>:'
    write (*,*) '                gt: select point if value >  <x1>'
    write (*,*) '                ge: select point if value >= <x1>'
    write (*,*) '                lt: select point if value <  <x1>'
    write (*,*) '                le: select point if value <= <x1>'
    write (*,*) '                range: select point if <x1> <= value <= <x2>'
    write (*,*) '                abs_range: select point if <x1> <= abs(value) <= <x2>'
    write (*,*) '                up: select entire column if any value in column >=  <x1>'
    write (*,*) '                dn: select entire column if any value in column <=  <x1>'
    write (*,*) '                up_dn: select entire column if any value in column >= <x1> or <=  <x2>'
    write (*,*) '                cylvol: select point if it is located within a cylindrical volume'
    write (*,*) '                   <x1> -> minimum radius (in km)'
    write (*,*) '                   <x2> -> maximum radius (in km)'
    write (*,*) '                   <y1> -> minimum angle (phi, in radians)'
    write (*,*) '                   <y2> -> maximum angle (phi, in radians)'
    write (*,*) '                   <z1> -> minimum height (in m)'
    write (*,*) '                   <z2> -> maximum height (in m)'
    write (*,*) ''
    write (*,*) '            <vname>: name of variable inside HDF5 file'
    write (*,*) '            <vfile>: complete path to file containing variable' 
    write (*,*) ''
    write (*,*) '          Note that more than one <filter_spec> can be specified. When this is done the data selection'
    write (*,*) '          becomes the intersection or union of all the filter specs depending on the <comb_sense> spec.'
    write (*,*) ''
    write (*,*) '        <out_file>: out:<file_name>:<variable_name>'
    write (*,*) '          <file_name> is the complete path to the file'
    write (*,*) '          <variable_name> is the name of the variable inside the file'
    write (*,*) ''

    stop
  end if

  return
end subroutine GetMyArgs

!**********************************************************************
! BuildUpDnMask()
!
! This routine will construct a mask (horizontal 2d field) that shows
! where updfraft and/or downdraft speed surpasses a threshold.
!
! UpDnCode:
!   1 - check for updrafts only
!   2 - check for downdrafts only
!   3 - check for updrafts and downdrafts
!
!   Note: Need to keep the list of legal values in sync with caller.
!
! For updrafts only (UpDnCode == 1) and downdrafts only (UpDnCode == 2), the
! associated threshold is passed in Thresh1. For up and down drafs
! (UpDnCode == 3) the updraft threshold is in Thresh1 and the downdraft
! threshold is in Thresh2.

subroutine BuildUpDnMask(Nx, Ny, Nz, Wvar, Wmask, Thresh1, Thresh2, UpDnCode)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: Wvar
  logical, dimension(Nx,Ny) :: Wmask
  real :: Thresh1, Thresh2
  integer :: UpDnCode

  integer :: ix, iy, iz

  ! For execution speed, want to walk thorugh Wvar in the same order as
  ! contiguous memory (which will minimize cache misses). Since the array
  ! is column major, this means the looping from outer to inner should be
  ! done as: z, y, x.
  !
  ! Want to look at columns and set Wmask(ix,iy) to .true. if any w value
  ! in the column surpasses the corresponding threshold. Putting z on the
  ! outer loop makes this a bit awkward but it is still do-able.

  do iy = 1, Ny
    do ix = 1, Nx
      Wmask(ix,iy) = .false.
    enddo
  enddo

  do iz = 1, Nz
    do iy = 1, Ny
      do ix = 1, Nx
        ! if already found a qualifying w value, skip checking this column
        if (Wmask(ix,iy) .eqv. .false.) then
          if ((UpDnCode .eq. 1) .and. (Wvar(ix,iy,iz) .ge. Thresh1)) then
            Wmask(ix,iy) = .true.
          else if ((UpDnCode .eq. 2) .and. (Wvar(ix,iy,iz) .le. Thresh1)) then
            Wmask(ix,iy) = .true.
          else if ((UpDnCode .eq. 3) .and. ((Wvar(ix,iy,iz) .ge. Thresh1) .or. (Wvar(ix,iy,iz) .le. Thresh2))) then
            Wmask(ix,iy) = .true.
          endif
        endif
      enddo
    enddo
  enddo

end subroutine BuildUpDnMask

end program diag_filter
