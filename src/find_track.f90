!***************************************************************
! Program to find the storm track.
!
! Args
!   1. directory containing input files
!   2. suffix to tag onto the end of input file names
!   3. output file name 
!   4. selection of averaging function
!
! Output
!   1. minimum pressure
!   2. longitude of minimum pressure location
!   3. latitude of minimum pressure location
!   4. radius of every horizontal point (distance from storm center)
!

program find_track
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

  type PressureSpec
    character (len=LittleString) :: Vname
    character (len=LittleString) :: Vfprefix
    integer :: Npts
    real :: Sigma
  end type PressureSpec

  character (len=LargeString) :: InDir, InSuffix, OutFile
  type (FilterDescrip), dimension(MaxFilters) :: Filters
  type (PressureSpec) :: Pspecs

  integer :: Nfilters

  integer :: i
  integer :: ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt
  logical :: BadDims

  type (Rhdf5Var), dimension(MaxFilters) :: Vars
  character (len=RHDF5_MAX_STRING), dimension(MaxFilters) :: InFiles
  type (Rhdf5Var) :: Pvar
  character (len=RHDF5_MAX_STRING) :: Pfile
  integer :: PfileId

  character (len=RHDF5_MAX_STRING) :: FileAcc
  integer, dimension(MaxFilters) :: InFileIds
  integer :: OutFileId

  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm

  real :: DeltaX, DeltaY
  real :: Rval

  type (Rhdf5Var) :: Radius, MinP, StormX, StormY
  integer :: StmIx, StmIy

  logical, dimension(:,:), allocatable :: DataSelect
  logical :: SelectThisPoint

  ! Get the command line arguments
  call GetMyArgs(LargeString, MaxFilters, InDir, InSuffix, OutFile, Pspecs, Filters, Nfilters)

  write (*,*) 'Finding storm track:'
  write (*,*) '  Input directory: ', trim(InDir)
  write (*,*) '  Input file name suffix: ', trim(InSuffix)
  write (*,*) '  Output file name:  ', trim(OutFile)
  write (*,*) '  Pressure specs:'
  write (*,*) '    Variable:  ', trim(Pspecs%Vname)
  write (*,*) '    Variable file prefix:  ', trim(Pspecs%Vfprefix)
  write (*,*) '    Gaussian smoothing:'
  write (*,*) '      Number of points:  ', Pspecs%Npts
  write (*,*) '      Sigma:  ', Pspecs%Sigma
  do i = 1, Nfilters
    if (i .eq. 1) then
      write (*,*) '  Filter specs: '
    endif

    if (Filters(i)%Ftype .eq. 'ge') then
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
  ! and time dimensions match between all vars. Compare against the pressure file.
  !
  ! Since we will be processing one time step at a time, the effective number of dimensions
  ! on each of the input is decreased by one. We want to eliminate the time dimension which
  ! is always the last one. This works out conveniently since all we need to do is decrement
  ! the number of dimensions by one.

  BadDims = .false.

  ! Always have pressure var. Use x,y,t to set output dims (output will be 2D).
  Pfile = trim(InDir) // '/' // trim(Pspecs%Vfprefix) // trim(InSuffix)
  Pvar%vname = trim(Pspecs%Vname)
  call rhdf5_read_init(Pfile, Pvar)

  Nx = Pvar%dims(1)
  Ny = Pvar%dims(2)
  Nz = 1
  if (Pvar%ndims .eq. 3) then
    ! 2D field
    Nt = Pvar%dims(3)
  else
    ! 3D field
    Nt = Pvar%dims(4)
  endif

  ! Prepare for reading
  Pvar%ndims = Pvar%ndims - 1
  if (Pvar%ndims .eq. 2) then
    allocate(Pvar%vdata(Pvar%dims(1)*Pvar%dims(2)))
  else
    allocate(Pvar%vdata(Pvar%dims(1)*Pvar%dims(2)*Pvar%dims(3)))
  endif

  do i = 1, Nfilters
    InFiles(i) = trim(InDir) // '/' // trim(Filters(i)%Vfprefix) // trim(InSuffix)
    Vars(i)%vname = trim(Filters(i)%Vname)

    call rhdf5_read_init(InFiles(i), Vars(i))


    ! check dimensions and prepare for reading
    if ((trim(Vars(i)%Vname) .eq. 'sfclat') .or. (trim(Vars(i)%Vname) .eq. 'sfclon')) then 
      ! just check horizontal dims (not time)
      if ((Vars(i)%dims(1) .ne. Pvar%dims(1)) .or. (Vars(i)%dims(2) .ne. Pvar%dims(2))) then 
        write (*,*) 'ERROR: horizontal dimensions of filter variable do not match the pressure variable: ', trim(Vars(i)%vname)
        BadDims = .true.
      else 
        ! Prepare for reading
        allocate(Vars(i)%vdata(Nx*Ny))
      endif
    else 
      if (.not. DimsMatch(Pvar, Vars(i))) then
        write (*,*) 'ERROR: horizontal and time dimensions of variable do not match the pressure variable: ', trim(Vars(i)%vname)
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
    endif
  enddo

  if (BadDims) then
    stop
  endif

  ! Set the output dimensions and coordinates to those of the selected input var
  call SetOutCoords(Pfile, Xcoords, Ycoords, Zcoords, Tcoords)
  
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

  ! Prepare output variables

  Radius%vname = 'radius'
  Radius%ndims = 2
  Radius%dims(1) = Nx
  Radius%dims(2) = Ny
  Radius%dimnames(1) = 'x'
  Radius%dimnames(2) = 'y'
  Radius%units = 'km'
  Radius%descrip = 'radius from storm center'
  allocate(Radius%vdata(Nx*Ny))

  ! Generate the storm center for all time steps
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

  ! Open the input files and the output file
  FileAcc = 'R'
  call rhdf5_open_file(Pfile, FileAcc, 0, PfileId)
  write (*,*) 'Reading HDF5 file: ', trim(Pfile)
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

  allocate(DataSelect(Nx,Ny))
  do it = 1, Nt

    ! Do the selection. If no filters, fill up the selection array with all .true.
    if (Nfilters .gt. 0) then

      do i = 1, Nfilters
        if ((trim(Vars(i)%Vname) .eq. 'sfclat') .or. (trim(Vars(i)%Vname) .eq. 'sfclon')) then 
          call rhdf5_read_variable(InFileIds(i), Vars(i)%vname, Vars(i)%ndims, 0, Vars(i)%dims, rdata=Vars(i)%vdata)
        else
          call rhdf5_read_variable(InFileIds(i), Vars(i)%vname, Vars(i)%ndims, it, Vars(i)%dims, rdata=Vars(i)%vdata)
        endif
      enddo

      ! Use iz = 1 for 2D data, iz = 2 (nearest level to surface above ground) for 3D data
      do iy = 1, Ny
        do ix = 1, Nx
          SelectThisPoint = .true.
          do i = 1, Nfilters
            ! select if var > threshold
            if (Filters(i)%Ftype .eq. 'gt') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, 1, Var2d=Vars(i)%vdata) .gt. Filters(i)%x1)
              else
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, 2, Var3d=Vars(i)%vdata) .gt. Filters(i)%x1)
              endif
            endif

            ! select if var >= threshold
            if (Filters(i)%Ftype .eq. 'ge') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, 1, Var2d=Vars(i)%vdata) .ge. Filters(i)%x1)
              else
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, 2, Var3d=Vars(i)%vdata) .ge. Filters(i)%x1)
              endif
            endif

            ! select if var < threshold
            if (Filters(i)%Ftype .eq. 'lt') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, 1, Var2d=Vars(i)%vdata) .lt. Filters(i)%x1)
              else
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, 2, Var3d=Vars(i)%vdata) .lt. Filters(i)%x1)
              endif
            endif

            ! select if var <= threshold
            if (Filters(i)%Ftype .eq. 'le') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, 1, Var2d=Vars(i)%vdata) .le. Filters(i)%x1)
              else
                SelectThisPoint = SelectThisPoint .and. &
                  (MultiDimLookup(Nx, Ny, Nz, ix, iy, 2, Var3d=Vars(i)%vdata) .le. Filters(i)%x1)
              endif
            endif

            ! select if min <= var <= max
            if (Filters(i)%Ftype .eq. 'range') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  ((MultiDimLookup(Nx, Ny, Nz, ix, iy, 1, Var2d=Vars(i)%vdata) .ge. Filters(i)%x1) .and. &
                   (MultiDimLookup(Nx, Ny, Nz, ix, iy, 1, Var2d=Vars(i)%vdata) .le. Filters(i)%x2))
              else
                SelectThisPoint = SelectThisPoint .and. &
                  ((MultiDimLookup(Nx, Ny, Nz, ix, iy, 2, Var3d=Vars(i)%vdata) .ge. Filters(i)%x1) .and. &
                   (MultiDimLookup(Nx, Ny, Nz, ix, iy, 2, Var3d=Vars(i)%vdata) .le. Filters(i)%x2))
              endif
            endif

            ! select if min <= abs(var) <= max
            if (Filters(i)%Ftype .eq. 'abs_range') then
              if (Vars(i)%ndims .eq. 2) then
                SelectThisPoint = SelectThisPoint .and. &
                  ((abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, 1, Var2d=Vars(i)%vdata)) .ge. Filters(i)%x1) .and. &
                   (abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, 1, Var2d=Vars(i)%vdata)) .le. Filters(i)%x2))
              else
                SelectThisPoint = SelectThisPoint .and. &
                  ((abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, 2, Var3d=Vars(i)%vdata)) .ge. Filters(i)%x1) .and. &
                   (abs(MultiDimLookup(Nx, Ny, Nz, ix, iy, 2, Var3d=Vars(i)%vdata)) .le. Filters(i)%x2))
              endif
            endif
          enddo

          DataSelect(ix,iy) = SelectThisPoint

        enddo
      enddo
    else
      do iy = i, Ny
        do ix = i, Nx
          DataSelect(ix,iy) = .true.
        enddo
      enddo
    endif

    call rhdf5_read_variable(PfileId, Pvar%vname, Pvar%ndims, it, Pvar%dims, rdata=Pvar%vdata)
    call FindStormCenter(Nx, Ny, Pvar%vdata, DataSelect, StmIx, StmIy, MinP%vdata(1))

    ! Record the lat,lon of the storm center
    StormX%vdata(1) = Xcoords%vdata(StmIx)
    StormY%vdata(1) = Ycoords%vdata(StmIy)

    ! calculate the radius
    do iy = 1, Ny
      DeltaY = YcoordsKm(iy) - YcoordsKm(StmIy)
      do ix = 1, Nx
        DeltaX = XcoordsKm(ix) - XcoordsKm(StmIx)
        Rval = sqrt(DeltaX**2 + DeltaY**2)
        call MultiDimAssign(Nx, Ny, Nz, ix, iy, 1, Rval, Var2d=Radius%vdata)
      enddo
    enddo

    ! Write the filter data, and storm info if doing cylvol selection, to the output file
    call rhdf5_write_variable(OutFileId, Radius%vname, Radius%ndims, it, Radius%dims, &
      Radius%units, Radius%descrip, Radius%dimnames, rdata=Radius%vdata)
    call rhdf5_write_variable(OutFileId, MinP%vname, MinP%ndims, it, MinP%dims, &
      MinP%units, MinP%descrip, MinP%dimnames, rdata=MinP%vdata)
    call rhdf5_write_variable(OutFileId, StormX%vname, StormX%ndims, it, StormX%dims, &
      StormX%units, StormX%descrip, StormX%dimnames, rdata=StormX%vdata)
    call rhdf5_write_variable(OutFileId, StormY%vname, StormY%ndims, it, StormY%dims, &
      StormY%units, StormY%descrip, StormY%dimnames, rdata=StormY%vdata)

    ! Write out status to screen every 100 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,100) .eq. 0) then
      write (*,*) 'Working: Timestep: ', it

      write (*,'(a,i3,a,i3,a,g,a,g,a)') '    Storm Center: (', StmIx, ', ', StmIy, &
       ') --> (', XcoordsKm(StmIx), ', ', YcoordsKm(StmIy), ')'
      write (*,*) '   Minumum pressure: ', MinP%vdata(1)

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

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(OutFile, Xcoords, 'x')
  call rhdf5_set_dimension(OutFile, Ycoords, 'y')
  call rhdf5_set_dimension(OutFile, Zcoords, 'z')
  call rhdf5_set_dimension(OutFile, Tcoords, 't')

  ! attach the dimension specs to the output variable

  ! Attach dims to the auxillary data
  call rhdf5_attach_dimensions(OutFile, MinP)
  call rhdf5_attach_dimensions(OutFile, StormX)
  call rhdf5_attach_dimensions(OutFile, StormY)
  call rhdf5_attach_dimensions(OutFile, Radius)

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

 subroutine GetMyArgs(Nstr, Nfilt, InDir, InSuffix, OutFile, Pspecs, Filters, NumFilters)
  implicit none

  integer, parameter :: MaxFields = 15

  character (len=Nstr) :: InDir, InSuffix, OutFile, Pvname, Pfprefix
  integer :: Nstr, Nfilt, NumFilters
  type (FilterDescrip), dimension(Nfilt) :: Filters
  type (PressureSpec) :: Pspecs

  integer :: iargc, i, j, Nargs, Nfld
  character (len=Nstr) :: Arg
  character (len=Nstr), dimension(MaxFields) :: Fields

  logical :: BadArgs

  BadArgs = .false.
  InDir = ''
  InSuffix = ''
  OutFile = ''
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
      call String2List(Arg, ':', Fields, MaxFields, Nfld, 'pressure spec')
      i = i + 1

      Pspecs%Vname = Fields(1)
      Pspecs%Vfprefix = Fields(2)
      read(Fields(3), '(i)') Pspecs%Npts
      read(Fields(4), '(f)') Pspecs%Sigma
    else
      call String2List(Arg, ':', Fields, MaxFields, Nfld, 'filter spec')
      i = i + 1

      !************
      !* GE
      !************
      if (Fields(1) .eq. 'ge') then
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
      else if (Fields(1) .eq. 'none') then
        ! skip over this, just here to allow zero filters to be spec'd
      else
        write (*,*) 'ERROR: <ftype>, ', trim(Fields(1)), ', must be one of:'
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
  
  if (BadArgs) then
    write (*,*) 'USAGE: find_track <in_dir> <in_suffix> <out_file> <filter_spec> [<filter_spec>...]'
    write (*,*) '        <in_dir>: directory containing input files'
    write (*,*) '        <in_suffix>: suffix to tag onto the end of input file names'
    write (*,*) '           Note: input file names become: <in_dir>/<var_name><in_suffix>'
    write (*,*) '        <out_file>: output hdf5 file name'
    write (*,*) '        <filter_spec>: <ftype>:<vname>:<vfprefix>:<x1>:<x2>:<y1>:<y2>:<z1>:<z2>'
    write (*,*) ''
    write (*,*) '            <ftype>:'
    write (*,*) '                gt: select point if value >  <x1>'
    write (*,*) '                ge: select point if value >= <x1>'
    write (*,*) '                lt: select point if value <  <x1>'
    write (*,*) '                le: select point if value <= <x1>'
    write (*,*) '                range: select point if <x1> <= value <= <x2>'
    write (*,*) '                abs_range: select point if <x1> <= abs(value) <= <x2>'
    write (*,*) ''
    write (*,*) '            <vname>: name of variable inside HDF5 file'
    write (*,*) '            <vfprefix>: file name prefix (goes with <in_suffix> to form the whole file name)' 
    write (*,*) ''
    write (*,*) '        Note that more than one <filter_spec> can be specified. When this is done the data selection'
    write (*,*) '        becomes the intersection of all the filter specs.'
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

subroutine FindStormCenter(Nx, Ny, Press, DataSelect, StmCtrX, StmCtrY, MinP)
  implicit none

  integer :: Nx, Ny
  real, dimension(Nx,Ny) :: Press
  logical, dimension(Nx,Ny) :: DataSelect
  integer :: StmCtrX, StmCtrY
  real :: MinP
  logical :: UseFsFilter

  integer :: ix, iy

  MinP = 1e10 ! ridiculously large pressure
  StmCtrX = 0 
  StmCtrY = 0 

  do ix = 1, Nx
    do iy = 1, Ny
      ! DataSelect holds true for points that we want to consider
      ! for minimum pressure
      if (DataSelect(ix,iy)) then 
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

end program find_track
