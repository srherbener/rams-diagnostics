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

  integer, parameter :: MaxFilters = 10

  type FilterDescrip
    character (len=LittleString) :: Ftype
    character (len=LittleString) :: Vname
    character (len=LittleString) :: Vfprefix
    real :: x1
    real :: x2
    real :: y1
    real :: y2
    real :: z1
    real :: z2
  end type FilterDescrip

  character (len=LargeString) :: InDir, InSuffix, OutFile
  type (FilterDescrip), dimension(MaxFilters) :: Filters
  character (len=LargeString) :: Mvname, Mfprefix, CombSense

  integer :: Nfilters

  integer :: i, ipress
  integer :: ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt
  logical :: BadDims
  logical :: DoingCylVol

  type (Rhdf5Var), dimension(MaxFilters) :: Vars
  character (len=RHDF5_MAX_STRING), dimension(MaxFilters) :: InFiles
  character (len=RHDF5_MAX_STRING) :: Mfile
  character (len=RHDF5_MAX_STRING) :: FsFilterFile
  type (Rhdf5Var) :: Mvar
  type (Rhdf5Var) :: FsFilterVar

  character (len=RHDF5_MAX_STRING) :: FileAcc
  integer, dimension(MaxFilters) :: InFileIds
  integer :: OutFileId
  integer :: rh5f_fsf__fid

  type (Rhdf5Var) :: OutFilter
  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords

  real :: DeltaX, DeltaY
  real :: Xstart, Xinc, Ystart, Yinc
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm

  type (Rhdf5Var) :: Radius, MinP, StormX, StormY, MaxRadius, StormXidx, StormYidx
  real :: Rval, Pval, Zval, MaxRval
  integer :: StmIx, StmIy

  logical :: FilterVal, SelectThisPoint
  logical :: UseFsFilter
  type (FilterDescrip) :: FsFilter

  logical, dimension(:,:), allocatable :: UpDnDraftMask
  logical :: DoingUpDrafts
  logical :: DoingDnDrafts
  logical :: DoingUpDnDrafts
  integer :: UDfnum

  ! Get the command line arguments
  call GetMyArgs(LargeString, MaxFilters, InDir, InSuffix, OutFile, Mvname, Mfprefix, FsFilter, CombSense, Filters, Nfilters)

  DoingCylVol = .false.
  DoingUpDrafts = .false.
  DoingDnDrafts = .false.
  DoingUpDnDrafts = .false.
  write (*,*) 'Creating HDF5 data filter:'
  write (*,*) '  Input directory: ', trim(InDir)
  write (*,*) '  Input file name suffix: ', trim(InSuffix)
  write (*,*) '  Output file name:  ', trim(OutFile)
  write (*,*) '  Model variable name:  ', trim(Mvname)
  write (*,*) '  Model file prefix: ', trim(Mfprefix)
  write (*,*) '  Find storm filter:'
  write (*,*) '    Filter type:  ', trim(FsFilter%Ftype)
  write (*,*) '    Variable name:  ', trim(FsFilter%Vname)
  write (*,*) '    File prefix: ', trim(FsFilter%Vfprefix)
  write (*,*) '    Treshold: ', FsFilter%x1
  write (*,*) '  Data selection specs: '
  write (*,*) '    Combination sense: ', trim(CombSense)
  do i = 1, Nfilters
    if (Filters(i)%Ftype .eq. 'cylvol') then
      DoingCylVol = .true.
      ipress = i
      write (*,*) '    Cylindrical volume:'
      write (*,*) '      Pressure variable: ', trim(Filters(i)%Vname)
      write (*,*) '      Pressure file prefix: ', trim(Filters(i)%Vfprefix)
      write (*,*) '      Minimum radius: ', Filters(i)%x1
      write (*,*) '      Maximum radius: ', Filters(i)%x2
      write (*,*) '      Minimum angle: ', Filters(i)%y1
      write (*,*) '      Maximum angle: ', Filters(i)%y2
      write (*,*) '      Minimum height: ', Filters(i)%z1
      write (*,*) '      Maximum height: ', Filters(i)%z2

      MaxRval = Filters(i)%x2
    else if (Filters(i)%Ftype .eq. 'ge') then
      write (*,*) '    Greater than or equal:'
      write (*,*) '      Variable: ', trim(Filters(i)%Vname)
      write (*,*) '      Variable file prefix: ', trim(Filters(i)%Vfprefix)
      write (*,*) '      Threshold: ', Filters(i)%x1
    else if (Filters(i)%Ftype .eq. 'gt') then
      write (*,*) '    Greater than:'
      write (*,*) '      Variable: ', trim(Filters(i)%Vname)
      write (*,*) '      Variable file prefix: ', trim(Filters(i)%Vfprefix)
      write (*,*) '      Threshold: ', Filters(i)%x1
    else if (Filters(i)%Ftype .eq. 'le') then
      write (*,*) '    Less than or equal:'
      write (*,*) '      Variable: ', trim(Filters(i)%Vname)
      write (*,*) '      Variable file prefix: ', trim(Filters(i)%Vfprefix)
      write (*,*) '      Threshold: ', Filters(i)%x1
    else if (Filters(i)%Ftype .eq. 'lt') then
      write (*,*) '    Less than:'
      write (*,*) '      Variable: ', trim(Filters(i)%Vname)
      write (*,*) '      Variable file prefix: ', trim(Filters(i)%Vfprefix)
      write (*,*) '      Threshold: ', Filters(i)%x1
    else if (Filters(i)%Ftype .eq. 'range') then
      write (*,*) '    Inside range:'
      write (*,*) '      Variable: ', trim(Filters(i)%Vname)
      write (*,*) '      Variable file prefix: ', trim(Filters(i)%Vfprefix)
      write (*,*) '      Min: ', Filters(i)%x1
      write (*,*) '      Max: ', Filters(i)%x2
    else if (Filters(i)%Ftype .eq. 'abs_range') then
      write (*,*) '    Absolute value inside range:'
      write (*,*) '      Variable: ', trim(Filters(i)%Vname)
      write (*,*) '      Variable file prefix: ', trim(Filters(i)%Vfprefix)
      write (*,*) '      Min: ', Filters(i)%x1
      write (*,*) '      Max: ', Filters(i)%x2
    else if (Filters(i)%Ftype .eq. 'up') then
      DoingUpDrafts = .true.
      UDfnum = i
      write (*,*) '    Updraft:'
      write (*,*) '      Variable: ', trim(Filters(i)%Vname)
      write (*,*) '      Variable file prefix: ', trim(Filters(i)%Vfprefix)
      write (*,*) '      Threshold: ', Filters(i)%x1
    else if (Filters(i)%Ftype .eq. 'dn') then
      DoingDnDrafts = .true.
      UDfnum = i
      write (*,*) '    Downdraft:'
      write (*,*) '      Variable: ', trim(Filters(i)%Vname)
      write (*,*) '      Variable file prefix: ', trim(Filters(i)%Vfprefix)
      write (*,*) '      Threshold: ', Filters(i)%x1
    else if (Filters(i)%Ftype .eq. 'up_dn') then
      DoingUpDnDrafts = .true.
      UDfnum = i
      write (*,*) '    Updraft/Downdraft:'
      write (*,*) '      Variable: ', trim(Filters(i)%Vname)
      write (*,*) '      Variable file prefix: ', trim(Filters(i)%Vfprefix)
      write (*,*) '      Updraft Threshold: ', Filters(i)%x1
      write (*,*) '      Downdraft Threshold: ', Filters(i)%x2
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

  ! Model var is expected to be 3d
  Mfile = trim(InDir) // '/' // trim(Mfprefix) // trim(InSuffix)
  Mvar%vname = trim(Mvname)
  call rhdf5_read_init(Mfile, Mvar)

  if (Mvar%ndims .ne. 4) then
    write (*,*) 'ERROR: model variable must be a 3d variable: ', trim(Mvname)
    stop
  else
    Nx = Mvar%dims(1)
    Ny = Mvar%dims(2)
    Nz = Mvar%dims(3)
    Nt = Mvar%dims(4)
  endif

  do i = 1, Nfilters
    InFiles(i) = trim(InDir) // '/' // trim(Filters(i)%Vfprefix) // trim(InSuffix)
    Vars(i)%vname = trim(Filters(i)%Vname)

    call rhdf5_read_init(InFiles(i), Vars(i))

    ! Check that the horizontal and time dimensions match with the model var. 
    if (.not. DimsMatch(Mvar, Vars(i))) then
      write (*,*) 'ERROR: horizontal and time dimensions of variable do not match the model variable: ', trim(Vars(i)%vname)
      BadDims = .true.
    else
      ! Prepare for reading
      Vars(i)%ndims = Vars(i)%ndims - 1
    endif
  enddo

  ! Check to see if we are using the find storm "hints". It is only used by
  ! the cylvol filter in the call to FindStormCenter with the purpose of keeping the
  ! the search for the storm center over the ocean.
  UseFsFilter = .false.
  if (DoingCylVol) then
    if (trim(FsFilter%Ftype) .ne. 'none') then
      UseFsFilter = .true.
      FsFilterFile = trim(InDir) // '/' // trim(FsFilter%Vfprefix) // trim(InSuffix)
      FsFilterVar%vname = trim(FsFilter%Vname)
      call rhdf5_read_init(FsFilterFile, FsFilterVar)

      ! check dimensions and prepare for reading
      if ((trim(FsFilterVar%Vname) .eq. 'sfclat') .or. (trim(FsFilterVar%Vname) .eq. 'sfclon')) then
        ! just check horizontal dims (not time)
        if ((FsFilterVar%dims(1) .ne. Mvar%dims(1)) .or. (FsFilterVar%dims(2) .ne. Mvar%dims(2))) then
          write (*,*) 'ERROR: horizontal and time dimensions of find storm variable do not match the model variable: ', trim(FsFilterVar%vname)
          BadDims = .true.
        endif
      else
        if (.not. DimsMatch(Mvar, FsFilterVar)) then
          write (*,*) 'ERROR: horizontal and time dimensions of find storm variable do not match the model variable: ', trim(FsFilterVar%vname)
          BadDims = .true.
        endif
      endif
    endif
  endif

  if (BadDims) then
    stop
  endif

  ! Set the output dimensions and coordinates to those of the selected input var
  call SetOutCoords(Mfile, Xcoords, Ycoords, Zcoords, Tcoords)
  
  ! Convert lat (x coords) and lon (y coords) to distances in km
  allocate(XcoordsKm(Nx))
  allocate(YcoordsKm(Ny))
  call ConvertGridCoords(Nx, Ny, Nz, Xcoords%vdata, Ycoords%vdata, XcoordsKm, YcoordsKm)

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

  if (DoingCylVol) then
    ! prepare for writing the radius values into the output file.
    Radius%vname = 'radius'
    Radius%ndims = 2
    Radius%dims(1) = Nx
    Radius%dims(2) = Ny
    Radius%dimnames(1) = 'x'
    Radius%dimnames(2) = 'y'
    Radius%units = 'km'
    Radius%descrip = 'radius from storm center'
    allocate(Radius%vdata(Nx*Ny))

    ! ndims == 0 means this var will be a time series, f(t)
    MinP%vname = 'min_press'
    MinP%ndims = 0
    MinP%units = 'mb'
    MinP%descrip = 'minimum SLP of storm'
    allocate(MinP%vdata(1))

    StormX%vname = 'min_press_xloc'
    StormX%ndims = 0
    StormX%units = 'deg lon'
    StormX%descrip = 'longitude location of minimum SLP of storm'
    allocate(StormX%vdata(1))
    
    StormY%vname = 'min_press_yloc'
    StormY%ndims = 0
    StormY%units = 'deg lat'
    StormY%descrip = 'latitude location of minimum SLP of storm'
    allocate(StormY%vdata(1))

    StormXidx%vname = 'min_press_x_index'
    StormXidx%ndims = 0
    StormXidx%units = 'lon index'
    StormXidx%descrip = 'longitude index of minimum SLP of storm'
    allocate(StormXidx%vdata(1))
    
    StormYidx%vname = 'min_press_y_index'
    StormYidx%ndims = 0
    StormYidx%units = 'lat index'
    StormYidx%descrip = 'latitude index of minimum SLP of storm'
    allocate(StormYidx%vdata(1))

    MaxRadius%vname = 'max_radius'
    MaxRadius%ndims = 0 
    MaxRadius%units = 'km'
    MaxRadius%descrip = 'maximum radius across domain and time'
    allocate(MaxRadius%vdata(1))
    MaxRadius%vdata(1) = MaxRval ! use the maximum radius spec from the command line
  endif

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
  do i = 1, Nfilters
    call rhdf5_open_file(InFiles(i), FileAcc, 0, InFileIds(i))
    write (*,*) 'Reading HDF5 file: ', trim(InFiles(i))
  enddo
  if (UseFsFilter) then
    call rhdf5_open_file(FsFilterFile, FileAcc, 0, rh5f_fsf__fid)
    write (*,*) 'Reading HDF5 file: ', trim(FsFilterFile)
  endif
  write (*,*) ''

  FileAcc = 'W'
  call rhdf5_open_file(OutFile, FileAcc, 1, OutFileId)
  write (*,*) 'Writing HDF5 file: ', trim(OutFile)
  write (*,*) ''

  ! Do the filtering one time step at a time.
  !
  ! The necessary input files have been opened, data buffers allocated,
  ! plus the time dimensions have been "removed" from the descriptions
  ! of the input variables.
  !

  do it = 1, Nt
    ! read the input vars
    do i = 1, Nfilters
      call rhdf5_read_variable(InFileIds(i), Vars(i)%vname, Vars(i)%ndims, it, Vars(i)%dims, rdata=Vars(i)%vdata)
    enddo

    ! In the case that we are doing cylvol, the above loop will have read in the pressure
    ! data, and ipress will point to the index of the pressure data. Find the storm center
    ! to prepare for the InsideCylVol call.
    if (DoingCylVol) then
      ! zero out radius
      do iy = 1, Ny
        do ix = 1, Nx
          call MultiDimAssign(Nx, Ny, Nz, ix, iy, iz, 0.0, Var2d=Radius%vdata)
        enddo
      enddo

      ! Find the storm center
      if (UseFsFilter) then
        if ((trim(FsFilterVar%Vname) .eq. 'sfclat') .or. (trim(FsFilterVar%Vname) .eq. 'sfclon')) then
          call rhdf5_read_variable(rh5f_fsf__fid, FsFilterVar%vname, FsFilterVar%ndims, 0, FsFilterVar%dims, rdata=FsFilterVar%vdata)
        else
          call rhdf5_read_variable(rh5f_fsf__fid, FsFilterVar%vname, FsFilterVar%ndims, it, FsFilterVar%dims, rdata=FsFilterVar%vdata)
        endif
      endif
      call FindStormCenter(Nx, Ny, Vars(ipress)%vdata, UseFsFilter, FsFilter, FsFilterVar%vdata, StmIx, StmIy, MinP%vdata(1))

      ! Record the lat,lon of the storm center
      StormX%vdata(1) = Xcoords%vdata(StmIx)
      StormY%vdata(1) = Ycoords%vdata(StmIy)
      StormXidx%vdata(1) = float(StmIx)
      StormYidx%vdata(1) = float(StmIy)

      ! clean up
      if (UseFsFilter) then
        deallocate(FsFilterVar%vdata)
      endif
    endif

    ! If doing up- or down-draft filtering, build the mask
    if (DoingUpDrafts .or. DoingDnDrafts .or. DoingUpDnDrafts) then
      do iy = 1, Ny
        do ix = 1, Nx
          UpDnDraftMask(ix,iy) = .false.

          do iz = 1, Nz
            if (DoingUpDrafts) then
              if (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(UDfnum)%vdata) .ge. Filters(UDfnum)%x1) then
                UpDnDraftMask(ix,iy) = .true.
                exit
              endif
            else if (DoingDnDrafts) then
              ! doing down drafts
              if (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(UDfnum)%vdata) .le. Filters(UDfnum)%x1) then
                UpDnDraftMask(ix,iy) = .true.
                exit
              endif
            else
              ! doing up/down drafts
              if ((MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(UDfnum)%vdata) .ge. Filters(UDfnum)%x1) .or. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(UDfnum)%vdata) .le. Filters(UDfnum)%x2)) then
                UpDnDraftMask(ix,iy) = .true.
                exit
              endif
            endif
          enddo
        enddo
      enddo
    endif

    ! do the selection
    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          if (CombSense .eq. 'and') then
            SelectThisPoint = .true.
          else if (CombSense .eq. 'or') then
            SelectThisPoint = .false.
          endif

          do i = 1, Nfilters
            ! select if inside the cylindrical volume
            if (Filters(i)%Ftype .eq. 'cylvol') then
              FilterVal = InsideCylVol(Nx, Ny, Nz, ix, iy, iz, Filters(i)%x1, Filters(i)%x2, Filters(i)%y1, Filters(i)%y2, &
                    Filters(i)%z1, Filters(i)%z2, StmIx, StmIy, XcoordsKm, YcoordsKm, Zcoords%vdata, Rval, Pval, Zval) 
            endif

            ! select if var > threshold
            if (Filters(i)%Ftype .eq. 'gt') then
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .gt. Filters(i)%x1)
              else
                FilterVal = (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .gt. Filters(i)%x1)
              endif
            endif

            ! select if var >= threshold
            if (Filters(i)%Ftype .eq. 'ge') then
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .ge. Filters(i)%x1)
              else
                FilterVal = (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .ge. Filters(i)%x1)
              endif
            endif

            ! select if var < threshold
            if (Filters(i)%Ftype .eq. 'lt') then
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .lt. Filters(i)%x1)
              else
                FilterVal = (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .lt. Filters(i)%x1)
              endif
            endif

            ! select if var <= threshold
            if (Filters(i)%Ftype .eq. 'le') then
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .le. Filters(i)%x1)
              else
                FilterVal = (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .le. Filters(i)%x1)
              endif
            endif

            ! select if min <= var <= max
            if (Filters(i)%Ftype .eq. 'range') then
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .ge. Filters(i)%x1) .and. &
                            (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .le. Filters(i)%x2)
              else
                FilterVal = (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .ge. Filters(i)%x1) .and. &
                            (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .le. Filters(i)%x2)
              endif
            endif

            ! select if min <= abs(var) <= max
            if (Filters(i)%Ftype .eq. 'abs_range') then
              if (Vars(i)%ndims .eq. 2) then
                FilterVal = (abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata)) .ge. Filters(i)%x1) .and. &
                            (abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata)) .le. Filters(i)%x2)
              else
                FilterVal = (abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata)) .ge. Filters(i)%x1) .and. &
                            (abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata)) .le. Filters(i)%x2)
              endif
            endif

            ! select according to up down draft mask
            if ((Filters(i)%Ftype .eq. 'up') .or. (Filters(i)%Ftype .eq. 'dn') .or. (Filters(i)%Ftype .eq. 'up_dn')) then
              FilterVal = UpDnDraftMask(ix,iy)
            endif

            ! combine the current filter value with the overall selection according to CombSense
            if (CombSense .eq. 'and') then
              SelectThisPoint = SelectThisPoint .and. FilterVal
            else if (CombSense .eq. 'or') then
              SelectThisPoint = SelectThisPoint .or. FilterVal
            end if
          enddo

          if (DoingCylVol) then
            ! save the radius value
            ! note that this re-writes the same horizontal radius values for each z level
            call MultiDimAssign(Nx, Ny, Nz, ix, iy, iz, Rval, Var2d=Radius%vdata)
          endif

          if (SelectThisPoint) then
            call MultiDimAssign(Nx, Ny, Nz, ix, iy, iz, 1.0, Var3d=OutFilter%vdata)
          else
            call MultiDimAssign(Nx, Ny, Nz, ix, iy, iz, 0.0, Var3d=OutFilter%vdata)
          endif
        enddo
      enddo
    enddo

    ! Write the filter data, and storm info if doing cylvol selection, to the output file
    call rhdf5_write_variable(OutFileId, OutFilter%vname, OutFilter%ndims, it, OutFilter%dims, &
      OutFilter%units, OutFilter%descrip, OutFilter%dimnames, rdata=OutFilter%vdata)
    if (DoingCylVol) then
      call rhdf5_write_variable(OutFileId, Radius%vname, Radius%ndims, it, Radius%dims, &
        Radius%units, Radius%descrip, Radius%dimnames, rdata=Radius%vdata)
      call rhdf5_write_variable(OutFileId, MinP%vname, MinP%ndims, it, MinP%dims, &
        MinP%units, MinP%descrip, MinP%dimnames, rdata=MinP%vdata)
      call rhdf5_write_variable(OutFileId, StormX%vname, StormX%ndims, it, StormX%dims, &
        StormX%units, StormX%descrip, StormX%dimnames, rdata=StormX%vdata)
      call rhdf5_write_variable(OutFileId, StormY%vname, StormY%ndims, it, StormY%dims, &
        StormY%units, StormY%descrip, StormY%dimnames, rdata=StormY%vdata)
      call rhdf5_write_variable(OutFileId, StormXidx%vname, StormXidx%ndims, it, StormXidx%dims, &
        StormXidx%units, StormXidx%descrip, StormXidx%dimnames, rdata=StormXidx%vdata)
      call rhdf5_write_variable(OutFileId, StormYidx%vname, StormYidx%ndims, it, StormYidx%dims, &
        StormYidx%units, StormYidx%descrip, StormYidx%dimnames, rdata=StormYidx%vdata)
      ! This repeats the max radius value for each time step, however this keeps the size of
      ! the data array consistent with the size of t_coords
      call rhdf5_write_variable(OutFileId, MaxRadius%vname, MaxRadius%ndims, it, MaxRadius%dims, &
        MaxRadius%units, MaxRadius%descrip, MaxRadius%dimnames, rdata=MaxRadius%vdata)
    endif

    ! cleanup
    do i = 1, Nfilters
      deallocate(Vars(i)%vdata)
    enddo

    ! Write out status to screen every 100 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,100) .eq. 0) then
      write (*,*) 'Working: Timestep: ', it

      if (DoingCylVol) then
        write (*,'(a,i3,a,i3,a,g,a,g,a)') '    Storm Center: (', StmIx, ', ', StmIy, &
         ') --> (', XcoordsKm(StmIx), ', ', YcoordsKm(StmIy), ')'
        write (*,*) '   Minumum pressure: ', MinP%vdata(1)
      endif

      write (*,*) ''
    endif
  enddo

  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''

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
  call rhdf5_attach_dimensions(OutFile, OutFilter)

  if (DoingCylVol) then
    ! Attach dims to the auxillary data
    call rhdf5_attach_dimensions(OutFile, MinP)
    call rhdf5_attach_dimensions(OutFile, StormX)
    call rhdf5_attach_dimensions(OutFile, StormY)
    call rhdf5_attach_dimensions(OutFile, StormXidx)
    call rhdf5_attach_dimensions(OutFile, StormYidx)
    call rhdf5_attach_dimensions(OutFile, Radius)
    call rhdf5_attach_dimensions(OutFile, MaxRadius)
  endif

  ! cleanup
  call rhdf5_close_file(OutFileId)
  do i = 1, Nfilters
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

 subroutine GetMyArgs(Nstr, Nfilt, InDir, InSuffix, OutFile, Mvname, Mfprefix, FsFilter, CombSense, Filters, NumFilters)
  implicit none

  integer, parameter :: MaxFields = 15

  character (len=Nstr) :: InDir, InSuffix, OutFile, Mvname, Mfprefix, CombSense
  integer :: Nstr, Nfilt, NumFilters
  type (FilterDescrip), dimension(Nfilt) :: Filters
  type (FilterDescrip) :: FsFilter

  integer :: iargc, i, j, Nargs, Nfld
  character (len=Nstr) :: Arg
  character (len=Nstr), dimension(MaxFields) :: Fields

  logical :: BadArgs

  BadArgs = .false.
  InDir = ''
  InSuffix = ''
  OutFile = ''
  Mvname = ''
  Mfprefix = ''
  NumFilters = 0

  ! walk through the list of arguments
  !   see the USAGE string at the end of this routine
  !
  !   first arg --> InDir
  !   second arg --> InSuffix
  !   third arg --> OutFile
  !
  !   remaining arguments are the filter specs

  i = 1
  Nargs = iargc()
  do while (i .le. Nargs)
    call getarg(i, Arg)

    if (i .eq. 1) then
      InDir = Arg
      i = i + 1
    else if (i .eq. 2) then
      InSuffix = Arg
      i = i + 1
    else if (i .eq. 3) then
      OutFile = Arg
      i = i + 1
    else if (i .eq. 4) then
      call String2List(Arg, ':', Fields, MaxFields, Nfld, 'model spec')
      i = i + 1

      Mvname   = Fields(1)
      Mfprefix = Fields(2)
    else if (i .eq. 5) then
      call String2List(Arg, ':', Fields, MaxFields, Nfld, 'find storm')
      i = i + 1

      FsFilter%Ftype    = Fields(1)
      FsFilter%Vname    = Fields(2)
      FsFilter%Vfprefix = Fields(3)
      read(Fields(4), '(f)') FsFilter%x1
      ! for now just support 'le', 'ge', 'lt', and 'gt'
      FsFilter%x2 = 0
      FsFilter%y1 = 0
      FsFilter%y2 = 0
      FsFilter%z1 = 0
      FsFilter%z2 = 0
    else if (i .eq. 6) then
      CombSense = Arg

      if ((CombSense .ne. 'and') .and. (CombSense .ne. 'or')) then
        write (*,*) 'ERROR: <comb_sense> must be one of: "and", "or"'
        BadArgs = .true.
      endif
      i = i + 1
    else
      call String2List(Arg, ':', Fields, MaxFields, Nfld, 'filter spec')
      i = i + 1

      !************
      !* CYLVOL
      !************
      if (Fields(1) .eq. 'cylvol') then
        NumFilters = NumFilters + 1

        if (Nfld .lt. 9) then
          write (*,*) 'ERROR: not enough arguments to fully specify the cylvol filter'
          BadArgs = .true.
        else
          ! have enough args
          Filters(NumFilters)%Ftype    = Fields(1)
          Filters(NumFilters)%Vname    = Fields(2)
          Filters(NumFilters)%Vfprefix = Fields(3)

          read(Fields(4), '(f)') Filters(NumFilters)%x1
          read(Fields(5), '(f)') Filters(NumFilters)%x2
          read(Fields(6), '(f)') Filters(NumFilters)%y1
          read(Fields(7), '(f)') Filters(NumFilters)%y2
          read(Fields(8), '(f)') Filters(NumFilters)%z1
          read(Fields(9), '(f)') Filters(NumFilters)%z2
        endif

        if ((Filters(NumFilters)%x1 .lt. 0.0) .or. (Filters(NumFilters)%x2 .lt. 0.0) .or. (Filters(NumFilters)%x2 .le. Filters(NumFilters)%x1)) then
          write (*,*) 'ERROR: <x1> and <x2> must be >= 0.0, and <x2> must be > <x1>'
          write (*,*) ''
          BadArgs = .true.
        end if
        if ((Filters(NumFilters)%y1 .lt. 0.0) .or. (Filters(NumFilters)%y2 .lt. 0.0) .or. (Filters(NumFilters)%y2 .le. Filters(NumFilters)%y1)) then
          write (*,*) 'ERROR: <y1> and <y2> must be >= 0.0, and <y2> must be > <y1>'
          write (*,*) ''
          BadArgs = .true.
        end if
        if ((Filters(NumFilters)%z1 .lt. 0.0) .or. (Filters(NumFilters)%z2 .lt. 0.0) .or. (Filters(NumFilters)%z2 .le. Filters(NumFilters)%z1)) then
          write (*,*) 'ERROR: <z1> and <z2> must be >= 0.0, and <z2> must be > <z1>'
          write (*,*) ''
          BadArgs = .true.
        end if
      !************
      !* GE
      !************
      else if (Fields(1) .eq. 'ge') then
        NumFilters = NumFilters + 1

        if (Nfld .lt. 4) then
          write (*,*) 'ERROR: not enough arguments to fully specify the ge filter'
          BadArgs = .true.
        else
          ! have enough args
          Filters(NumFilters)%Ftype    = Fields(1)
          Filters(NumFilters)%Vname    = Fields(2)
          Filters(NumFilters)%Vfprefix = Fields(3)

          read(Fields(4),  '(f)') Filters(NumFilters)%x1
          Filters(NumFilters)%x2 = 0
          Filters(NumFilters)%y1 = 0
          Filters(NumFilters)%y2 = 0
          Filters(NumFilters)%z1 = 0
          Filters(NumFilters)%z2 = 0
        endif
      !************
      !* GT
      !************
      else if (Fields(1) .eq. 'gt') then
        NumFilters = NumFilters + 1

        if (Nfld .lt. 4) then
          write (*,*) 'ERROR: not enough arguments to fully specify the gt filter'
          BadArgs = .true.
        else
          ! have enough args
          Filters(NumFilters)%Ftype    = Fields(1)
          Filters(NumFilters)%Vname    = Fields(2)
          Filters(NumFilters)%Vfprefix = Fields(3)

          read(Fields(4),  '(f)') Filters(NumFilters)%x1
          Filters(NumFilters)%x2 = 0
          Filters(NumFilters)%y1 = 0
          Filters(NumFilters)%y2 = 0
          Filters(NumFilters)%z1 = 0
          Filters(NumFilters)%z2 = 0
        endif
      !************
      !* LE
      !************
      else if (Fields(1) .eq. 'le') then
        NumFilters = NumFilters + 1

        if (Nfld .lt. 4) then
          write (*,*) 'ERROR: not enough arguments to fully specify the le filter'
          BadArgs = .true.
        else
          ! have enough args
          Filters(NumFilters)%Ftype    = Fields(1)
          Filters(NumFilters)%Vname    = Fields(2)
          Filters(NumFilters)%Vfprefix = Fields(3)

          read(Fields(4),  '(f)') Filters(NumFilters)%x1
          Filters(NumFilters)%x2 = 0
          Filters(NumFilters)%y1 = 0
          Filters(NumFilters)%y2 = 0
          Filters(NumFilters)%z1 = 0
          Filters(NumFilters)%z2 = 0
        endif
      !************
      !* LT
      !************
      else if (Fields(1) .eq. 'lt') then
        NumFilters = NumFilters + 1

        if (Nfld .lt. 4) then
          write (*,*) 'ERROR: not enough arguments to fully specify the lt filter'
          BadArgs = .true.
        else
          ! have enough args
          Filters(NumFilters)%Ftype    = Fields(1)
          Filters(NumFilters)%Vname    = Fields(2)
          Filters(NumFilters)%Vfprefix = Fields(3)

          read(Fields(4),  '(f)') Filters(NumFilters)%x1
          Filters(NumFilters)%x2 = 0
          Filters(NumFilters)%y1 = 0
          Filters(NumFilters)%y2 = 0
          Filters(NumFilters)%z1 = 0
          Filters(NumFilters)%z2 = 0
        endif
      !************
      !* RANGE
      !************
      else if (Fields(1) .eq. 'range') then
        NumFilters = NumFilters + 1

        if (Nfld .lt. 5) then
          write (*,*) 'ERROR: not enough arguments to fully specify the range filter'
          BadArgs = .true.
        else
          ! have enough args
          Filters(NumFilters)%Ftype    = Fields(1)
          Filters(NumFilters)%Vname    = Fields(2)
          Filters(NumFilters)%Vfprefix = Fields(3)

          read(Fields(4),  '(f)') Filters(NumFilters)%x1
          read(Fields(5),  '(f)') Filters(NumFilters)%x2
          Filters(NumFilters)%y1 = 0
          Filters(NumFilters)%y2 = 0
          Filters(NumFilters)%z1 = 0
          Filters(NumFilters)%z2 = 0
        endif
      !************
      !* ABS_RANGE
      !************
      else if (Fields(1) .eq. 'abs_range') then
        NumFilters = NumFilters + 1

        if (Nfld .lt. 5) then
          write (*,*) 'ERROR: not enough arguments to fully specify the abs_range filter'
          BadArgs = .true.
        else
          ! have enough args
          Filters(NumFilters)%Ftype    = Fields(1)
          Filters(NumFilters)%Vname    = Fields(2)
          Filters(NumFilters)%Vfprefix = Fields(3)

          read(Fields(4),  '(f)') Filters(NumFilters)%x1
          read(Fields(5),  '(f)') Filters(NumFilters)%x2
          Filters(NumFilters)%y1 = 0
          Filters(NumFilters)%y2 = 0
          Filters(NumFilters)%z1 = 0
          Filters(NumFilters)%z2 = 0
        endif
      else if (Fields(1) .eq. 'up') then
        NumFilters = NumFilters + 1

        if (Nfld .lt. 4) then
          write (*,*) 'ERROR: not enough arguments to fully specify the updraft filter'
          BadArgs = .true.
        else
          ! have enough args
          Filters(NumFilters)%Ftype    = Fields(1)
          Filters(NumFilters)%Vname    = Fields(2)
          Filters(NumFilters)%Vfprefix = Fields(3)

          read(Fields(4), '(f)') Filters(NumFilters)%x1
          Filters(NumFilters)%x2 = 0
          Filters(NumFilters)%y1 = 0
          Filters(NumFilters)%y2 = 0
          Filters(NumFilters)%z1 = 0
          Filters(NumFilters)%z2 = 0
        endif
      else if (Fields(1) .eq. 'dn') then
        NumFilters = NumFilters + 1

        if (Nfld .lt. 4) then
          write (*,*) 'ERROR: not enough arguments to fully specify the downdraft filter'
          BadArgs = .true.
        else
          ! have enough args
          Filters(NumFilters)%Ftype    = Fields(1)
          Filters(NumFilters)%Vname    = Fields(2)
          Filters(NumFilters)%Vfprefix = Fields(3)

          read(Fields(4), '(f)') Filters(NumFilters)%x1
          Filters(NumFilters)%x2 = 0
          Filters(NumFilters)%y1 = 0
          Filters(NumFilters)%y2 = 0
          Filters(NumFilters)%z1 = 0
          Filters(NumFilters)%z2 = 0
        endif
      else if (Fields(1) .eq. 'up_dn') then
        NumFilters = NumFilters + 1

        if (Nfld .lt. 5) then
          write (*,*) 'ERROR: not enough arguments to fully specify the updraft/downdraft filter'
          BadArgs = .true.
        else
          ! have enough args
          Filters(NumFilters)%Ftype    = Fields(1)
          Filters(NumFilters)%Vname    = Fields(2)
          Filters(NumFilters)%Vfprefix = Fields(3)

          read(Fields(4), '(f)') Filters(NumFilters)%x1
          read(Fields(5), '(f)') Filters(NumFilters)%x2
          Filters(NumFilters)%y1 = 0
          Filters(NumFilters)%y2 = 0
          Filters(NumFilters)%z1 = 0
          Filters(NumFilters)%z2 = 0
        endif
      else
        write (*,*) 'ERROR: <ftype>, ', trim(Fields(1)), ', must be one of:'
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
      endif
    endif
  enddo
  
  if (NumFilters .eq. 0) then
    write (*,*) 'ERROR: must specify at least one <filter_spec>'
    write (*,*) ''
    BadArgs = .true.
  endif

  if (BadArgs) then
    write (*,*) 'USAGE: diag_filter <in_dir> <in_suffix> <out_file> <var3d_model> <2d_find_storm_filter> <comb_sense> <filter_spec> [<filter_spec>...]'
    write (*,*) '        <in_dir>: directory containing input files'
    write (*,*) '        <in_suffix>: suffix to tag onto the end of input file names'
    write (*,*) '           Note: input file names become: <in_dir>/<var_name><in_suffix>'
    write (*,*) '        <out_file>: output hdf5 file name'
    write (*,*) '        <var3d_model>: <vname>:<vfprefix>'
    write (*,*) '            <vname>: hdf5 file name of variable that will serve as a model (same dimensions)'
    write (*,*) '                     for the filter output'
    write (*,*) '            <vfprefix>: file name prefix (goes with <in_suffix> to form the whole file name)' 
    write (*,*) '        <2d_find_storm_filter>: <ftype>:<vname>:<vfprefix>:<val>'
    write (*,*) '            <ftype>:'
    write (*,*) '                gt: select point if value >  <val>'
    write (*,*) '                ge: select point if value >= <val>'
    write (*,*) '                lt: select point if value <  <val>'
    write (*,*) '                le: select point if value <= <val>'
    write (*,*) '            <vname>: hdf5 variable'
    write (*,*) '            <vfprefix>: file name prefix (goes with <in_suffix> to form the whole file name)' 
    write (*,*) '            <val>'
    write (*,*) '        <comb_sense>: one of "and", "or"'
    write (*,*) '            and: multiple filter specs are and-ed together' 
    write (*,*) '            or: multiple filter specs are or-ed together' 
    write (*,*) '        <filter_spec>: <ftype>:<vname>:<vfprefix>:<x1>:<x2>:<y1>:<y2>:<z1>:<z2>'
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
    write (*,*) '            <vfprefix>: file name prefix (goes with <in_suffix> to form the whole file name)' 
    write (*,*) ''
    write (*,*) '        Note that more than one <filter_spec> can be specified. When this is done the data selection'
    write (*,*) '        becomes the intersection or union of all the filter specs depending on the <comb_sense> spec.'
    write (*,*) ''

    stop
  end if

  return
end subroutine GetMyArgs

!**********************************************************************
! FindStormCenter
!
! This routine will locate the storm center using the simple hueristic
! of the center being where the minimum surface pressure exists.
!
! Argument it holds the time step that you want to analyze. (iStmCtr, jStmCtr)
! hold the grid position of the minumum pressure value on the first vertical
! level (iz = 1).
!

subroutine FindStormCenter(Nx, Ny, Press, UseFsFilter, FsFilter, Fs, StmCtrX, StmCtrY, MinP)
  implicit none

  integer :: Nx, Ny
  real, dimension(Nx,Ny) :: Press, Fs
  type(FilterDescrip) :: FsFilter
  integer :: StmCtrX, StmCtrY
  real :: MinP
  logical :: UseFsFilter

  integer :: ix, iy
  logical :: SelectPoint

  MinP = 1e10 ! ridiculously large pressure
  StmCtrX = 0 
  StmCtrY = 0 

  do ix = 1, Nx
    do iy = 1, Ny
      ! Figure out if we want to consider this point
      SelectPoint = .true.
      if (UseFsFilter) then
        SelectPoint = .false.
        if (FsFilter%Ftype .eq. 'gt') then
          if (Fs(ix,iy) .gt. FsFilter%x1) then
            SelectPoint = .true.
          endif
        elseif (FsFilter%Ftype .eq. 'ge') then
          if (Fs(ix,iy) .ge. FsFilter%x1) then
            SelectPoint = .true.
          endif
        elseif (FsFilter%Ftype .eq. 'lt') then
          if (Fs(ix,iy) .lt. FsFilter%x1) then
            SelectPoint = .true.
          endif
        elseif (FsFilter%Ftype .eq. 'le') then
          if (Fs(ix,iy) .le. FsFilter%x1) then
            SelectPoint = .true.
          endif
        endif
      endif

      if (SelectPoint) then 
        if (Press(ix,iy) .lt. MinP) then
          MinP = Press(ix,iy)
          StmCtrX = ix
          StmCtrY = iy
        endif
      endif
    enddo
  enddo
  
  return
end subroutine FindStormCenter

end program diag_filter
