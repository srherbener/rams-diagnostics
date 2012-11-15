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
  character (len=LargeString) :: Mvname, Mfprefix

  integer :: Nfilters

  integer :: i, ipress
  integer :: ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt
  logical :: BadDims
  logical :: DoingCylVol

  type (Rhdf5Var), dimension(MaxFilters) :: Vars
  character (len=RHDF5_MAX_STRING), dimension(MaxFilters) :: InFiles
  character (len=RHDF5_MAX_STRING) :: Mfile
  type (Rhdf5Var) :: Mvar

  character (len=RHDF5_MAX_STRING) :: FileAcc
  integer, dimension(MaxFilters) :: InFileIds
  integer :: OutFileId

  type (Rhdf5Var) :: OutFilter
  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords, Dcoords

  real :: DeltaX, DeltaY
  real :: Xstart, Xinc, Ystart, Yinc
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm

  type (Rhdf5Var) :: Radius, MinP, StormX, StormY, MaxRadius
  real :: Rval, Pval, Zval
  integer :: StmIx, StmIy

  logical :: SelectThisPoint

  ! Get the command line arguments
  call GetMyArgs(LargeString, MaxFilters, InDir, InSuffix, OutFile, Mvname, Mfprefix, Filters, Nfilters)

  DoingCylVol = .false.
  write (*,*) 'Creating HDF5 data filter:'
  write (*,*) '  Input directory: ', trim(InDir)
  write (*,*) '  Input file name suffix: ', trim(InSuffix)
  write (*,*) '  Output file name:  ', trim(OutFile)
  write (*,*) '  Model variable name:  ', trim(Mvname)
  write (*,*) '  Model file prefix: ', trim(Mfprefix)
  write (*,*) '  Data selection specs: '
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
      if (Vars(i)%ndims .eq. 2) then
        allocate(Vars(i)%vdata(Nx*Ny))
      else
        allocate(Vars(i)%vdata(Nx*Ny*Nz))
      endif
    endif
  enddo

  if (BadDims) then
    stop
  endif

  ! Set the output dimensions and coordinates to those of the selected input var
  call SetOutCoords(Mfile, Xcoords, Ycoords, Zcoords, Tcoords)
  
  ! In order for GRADS to be able to read in the output HDF5 file, every variable
  ! needs to be four dimensional: (x,y,z,t) regardless if that is the true nature
  ! of the variable. To accommodate this, fill in missing dimensions with the following
  ! dummy coordinate so that the variariable meets the (x,y,z,t) structure requirement
  ! of GRADS. For example, say a variable is really (x,y,t), then add a z dimension
  ! with a size of one to it so it will be compatible with GRADS, (x,y,z,t).
  Dcoords%vname = 'd_coords'
  Dcoords%ndims = 1
  Dcoords%dims(1) = 1
  Dcoords%dimnames(1) = 'd'
  Dcoords%units = 'dummy'
  Dcoords%descrip = 'dummy coordinates'
  allocate(Dcoords%vdata(1))
  Dcoords%vdata = 0.0

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
    Radius%ndims = 3
    Radius%dims(1) = Nx
    Radius%dims(2) = Ny
    Radius%dims(3) = 1
    Radius%dimnames(1) = 'x'
    Radius%dimnames(2) = 'y'
    Radius%dimnames(3) = 'd'
    Radius%units = 'km'
    Radius%descrip = 'radius from storm center'
    allocate(Radius%vdata(Nx*Ny))

    ! Generate the storm center for all time steps
    MinP%vname = 'min_press'
    MinP%ndims = 3
    MinP%dims(1) = 1
    MinP%dims(2) = 1
    MinP%dims(3) = 1
    MinP%dimnames(1) = 'd'
    MinP%dimnames(2) = 'd'
    MinP%dimnames(3) = 'd'
    MinP%units = 'mb'
    MinP%descrip = 'minimum SLP of storm'
    allocate(MinP%vdata(1))

    StormX%vname = 'min_press_xloc'
    StormX%ndims = 3
    StormX%dims(1) = 1
    StormX%dims(2) = 1
    StormX%dims(3) = 1
    StormX%dimnames(1) = 'd'
    StormX%dimnames(2) = 'd'
    StormX%dimnames(3) = 'd'
    StormX%units = 'deg lon'
    StormX%descrip = 'longitude location of minimum SLP of storm'
    allocate(StormX%vdata(1))
    
    StormY%vname = 'min_press_yloc'
    StormY%ndims = 3
    StormY%dims(1) = 1
    StormY%dims(2) = 1
    StormY%dims(3) = 1
    StormY%dimnames(1) = 'd'
    StormY%dimnames(2) = 'd'
    StormY%dimnames(3) = 'd'
    StormY%units = 'deg lat'
    StormY%descrip = 'latitude location of minimum SLP of storm'
    allocate(StormY%vdata(1))

    MaxRadius%vname = 'max_radius'
    MaxRadius%ndims = 4
    MaxRadius%dims(1) = 1
    MaxRadius%dims(2) = 1
    MaxRadius%dims(3) = 1
    MaxRadius%dims(4) = 1
    MaxRadius%dimnames(1) = 'd'
    MaxRadius%dimnames(2) = 'd'
    MaxRadius%dimnames(3) = 'd'
    MaxRadius%dimnames(4) = 'd'
    MaxRadius%units = 'km'
    MaxRadius%descrip = 'maximum radius across domain and time'
    allocate(MaxRadius%vdata(1))
    MaxRadius%vdata(1) = -1.0 ! not expecting negative radii
  endif

  ! Open the input files and the output file
  FileAcc = 'R'
  do i = 1, Nfilters
    call rhdf5_open_file(InFiles(i), FileAcc, 0, InFileIds(i))
    write (*,*) 'Reading HDF5 file: ', trim(InFiles(i))
  enddo
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
      call FindStormCenter(Nx, Ny, Vars(ipress)%vdata, StmIx, StmIy, MinP%vdata(1))

      ! Record the lat,lon of the storm center
      StormX%vdata(1) = Xcoords%vdata(StmIx)
      StormY%vdata(1) = Ycoords%vdata(StmIy)
    endif

    ! do the selection
    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          SelectThisPoint = .true.
          do i = 1, Nfilters
            ! select if inside the cylindrical volume
            if (Filters(i)%Ftype .eq. 'cylvol') then
              SelectThisPoint = SelectThisPoint .and. InsideCylVol(Nx, Ny, Nz, ix, iy, iz, &
                    Filters(i)%x1, Filters(i)%x2, Filters(i)%y1, Filters(i)%y2, &
                    Filters(i)%z1, Filters(i)%z2, StmIx, StmIy, &
                    XcoordsKm, YcoordsKm, Zcoords%vdata, Rval, Pval, Zval) 
            endif

            ! select if var > threshold
            if (Filters(i)%Ftype .eq. 'gt') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .gt. Filters(i)%x1)
              else
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .gt. Filters(i)%x1)
              endif
            endif

            ! select if var >= threshold
            if (Filters(i)%Ftype .eq. 'ge') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .ge. Filters(i)%x1)
              else
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .ge. Filters(i)%x1)
              endif
            endif

            ! select if var < threshold
            if (Filters(i)%Ftype .eq. 'lt') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .lt. Filters(i)%x1)
              else
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .lt. Filters(i)%x1)
              endif
            endif

            ! select if var <= threshold
            if (Filters(i)%Ftype .eq. 'le') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .le. Filters(i)%x1)
              else
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .le. Filters(i)%x1)
              endif
            endif

            ! select if min <= var <= max
            if (Filters(i)%Ftype .eq. 'range') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  ((MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .ge. Filters(i)%x1) .and. &
                   (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata) .le. Filters(i)%x2))
              else
                SelectThisPoint = SelectThisPoint .and. &
                  ((MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .ge. Filters(i)%x1) .and. &
                   (MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata) .le. Filters(i)%x2))
              endif
            endif

            ! select if min <= abs(var) <= max
            if (Filters(i)%Ftype .eq. 'abs_range') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  ((abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata)) .ge. Filters(i)%x1) .and. &
                   (abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var2d=Vars(i)%vdata)) .le. Filters(i)%x2))
              else
                SelectThisPoint = SelectThisPoint .and. &
                  ((abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata)) .ge. Filters(i)%x1) .and. &
                   (abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Vars(i)%vdata)) .le. Filters(i)%x2))
              endif
            endif
          enddo

          if (SelectThisPoint) then
            call MultiDimAssign(Nx, Ny, Nz, ix, iy, iz, 1.0, Var3d=OutFilter%vdata)
            if (DoingCylVol) then
              ! save the radius value
              ! note that this re-writes the same horizontal radius values for each z level
              call MultiDimAssign(Nx, Ny, Nz, ix, iy, iz, Rval, Var2d=Radius%vdata)
              if (Rval .gt. MaxRadius%vdata(1)) then
                MaxRadius%vdata(1) = Rval
              endif
            endif
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
    endif

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

  ! Write out coordinates and attach dimensions

  ! write out the coordinate data
  call rhdf5_write(OutFile, Xcoords, 1)
  call rhdf5_write(OutFile, Ycoords, 1)
  call rhdf5_write(OutFile, Zcoords, 1)
  call rhdf5_write(OutFile, Tcoords, 1)
  call rhdf5_write(OutFile, Dcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(OutFile, Xcoords, 'x')
  call rhdf5_set_dimension(OutFile, Ycoords, 'y')
  call rhdf5_set_dimension(OutFile, Zcoords, 'z')
  call rhdf5_set_dimension(OutFile, Tcoords, 't')
  call rhdf5_set_dimension(OutFile, Dcoords, 'd')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(OutFile, OutFilter)

  if (DoingCylVol) then
    ! Write out the maximum radius value
    call rhdf5_write(OutFile, MaxRadius, 1)

    ! Attach dims to the auxillary data
    call rhdf5_attach_dimensions(OutFile, MinP)
    call rhdf5_attach_dimensions(OutFile, StormX)
    call rhdf5_attach_dimensions(OutFile, StormY)
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

 subroutine GetMyArgs(Nstr, Nfilt, InDir, InSuffix, OutFile, Mvname, Mfprefix, Filters, NumFilters)
  implicit none

  integer, parameter :: MaxFields = 15

  character (len=Nstr) :: InDir, InSuffix, OutFile, Mvname, Mfprefix
  integer :: Nstr, Nfilt, NumFilters
  type (FilterDescrip), dimension(Nfilt) :: Filters

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
    else
      call String2List(Arg, ':', Fields, MaxFields, Nfld, 'filter spec')
      i = i + 1

      !************
      !* CYLVOL
      !************
      if (Fields(1) .eq. 'cylvol') then
        NumFilters = NumFilters + 1

        if (Nfld .lt. 9) then
          write (*,*) 'ERROR: not enough arguments to fully specify the clyvol filter'
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
      else
        write (*,*) 'ERROR: <ftype>, ', trim(Fields(1)), ', must be one of:'
        write (*,*) '          cylvol'
        write (*,*) '          ge'
        write (*,*) '          gt'
        write (*,*) '          le'
        write (*,*) '          lt'
        write (*,*) '          range'
        write (*,*) '          abs_range'
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
    write (*,*) 'USAGE: diag_filter <in_dir> <in_suffix> <out_file> <var3d_model> <filter_spec> [<filter_spec>...]'
    write (*,*) '        <in_dir>: directory containing input files'
    write (*,*) '        <in_suffix>: suffix to tag onto the end of input file names'
    write (*,*) '           Note: input file names become: <in_dir>/<var_name><in_suffix>'
    write (*,*) '        <out_file>: output hdf5 file name'
    write (*,*) '        <var3d_model>: <vname>:<vfprefix>'
    write (*,*) '            <vname>: hdf5 file name of variable that will serve as a model (same dimensions)'
    write (*,*) '                     for the filter output'
    write (*,*) '            <vfprefix>: file name prefix (goes with <in_suffix> to form the whole file name)' 
    write (*,*) '        <filter_spec>: <ftype>:<vname>:<vfprefix>:<x1>:<x2>:<y1>:<y2>:<z1>:<z2>'
    write (*,*) ''
    write (*,*) '            <ftype>:'
    write (*,*) '                gt: select point if value >  <x1>'
    write (*,*) '                ge: select point if value >= <x1>'
    write (*,*) '                lt: select point if value <  <x1>'
    write (*,*) '                le: select point if value <= <x1>'
    write (*,*) '                range: select point if <x1> <= value <= <x2>'
    write (*,*) '                abs_range: select point if <x1> <= abs(value) <= <x2>'
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
    write (*,*) '        becomes the intersection of all the filter specs. Eg., specifying both cylvol and tempc less that or equal '
    write (*,*) '        to zero will select only the data that is both inside the cylindrical volume and in regions <= 0 deg C'
    write (*,*) ''

    stop
  end if

  return
end subroutine GetMyArgs

!**********************************************************************
! SetOutCoords()
!
! This routine will set the coordinate and dimension data 

subroutine SetOutCoords(Hfile, Xcoords, Ycoords, Zcoords, Tcoords)
  use rhdf5_utils
  use diag_utils
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

subroutine FindStormCenter(Nx, Ny, Press, StmCtrX, StmCtrY, MinP)
  implicit none

  integer :: Nx, Ny
  real, dimension(Nx,Ny) :: Press
  integer :: StmCtrX, StmCtrY
  real :: MinP

  integer :: ix, iy

  MinP = 1e10 ! ridiculously large pressure
  StmCtrX = 0 
  StmCtrY = 0 

  do ix = 1, Nx
    do iy = 1, Ny
      if (Press(ix,iy) .lt. MinP) then
        MinP = Press(ix,iy)
        StmCtrX = ix
        StmCtrY = iy
      end if
    end do
  end do
  
  return
end subroutine FindStormCenter

end program diag_filter
