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

program sig_proc
  use rhdf5_utils
  use diag_utils
  implicit none

  integer, parameter :: LargeString  = 512
  integer, parameter :: MediumString = 256
  integer, parameter :: LittleString = 128

  integer, parameter :: MaxPspecs = 10

  type PspecDescrip
    character (len=LittleString) :: Ptype
    integer :: i1
    integer :: i2
    real :: s1
    real :: s2
  end type PspecDescrip

  character (len=LargeString) :: InFile, InVarName, OutFile
  type (PspecDescrip), dimension(MaxPspecs) :: Pspecs

  integer :: Npspecs

  integer :: i
  integer :: ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt

  type (Rhdf5Var) :: InVar

  character (len=RHDF5_MAX_STRING) :: FileAcc
  integer :: InFileId
  integer :: OutFileId

  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm

  real :: DeltaX, DeltaY

  type (Rhdf5Var) :: OutVar

  real, dimension(:,:), allocatable :: InField, OutField

  ! Get the command line arguments
  call GetMyArgs(LargeString, MaxPspecs, InFile, InVarName, OutFile, Pspecs, Npspecs)

  ! For now, only supporting doing one processing step per run
  if (Npspecs .ne. 1) then
    write (*,*) 'ERROR: For now, must specify exactly one processing step per run'
    stop
  endif

  write (*,*) 'Signal Processing:'
  write (*,*) '  Input file: ', trim(InFile)
  write (*,*) '  Input variable name: ', trim(InVarName)
  write (*,*) '  Output file:  ', trim(OutFile)
  do i = 1, Npspecs
    if (i .eq. 1) then
      write (*,*) '  Processing specs: '
    endif

    if (Pspecs(i)%Ptype .eq. 'gauss2d') then
      write (*,*) '    Gaussian 2D smoothing:'
      write (*,*) '      Number of points: ', Pspecs(i)%i1
      write (*,*) '      Sigma: ', Pspecs(i)%s1
    endif
  enddo
  write (*,*) ''
  flush(6)

  ! Since we will be processing one time step at a time, the effective number of dimensions
  ! on each of the input is decreased by one. We want to eliminate the time dimension which
  ! is always the last one. This works out conveniently since all we need to do is decrement
  ! the number of dimensions by one.

  InVar%vname = trim(InVarName)
  call rhdf5_read_init(InFile, InVar)

  Nx = InVar%dims(1)
  Ny = InVar%dims(2)
  if (InVar%ndims .eq. 3) then
    ! 2D field
    Nz = 1
    Nt = InVar%dims(3)
  else
    ! 3D field
    Nz = InVar%dims(3)
    Nt = InVar%dims(4)
  endif

  ! Prepare for reading
  InVar%ndims = InVar%ndims - 1
  if (InVar%ndims .eq. 2) then
    allocate(InVar%vdata(Nx*Ny))
  else
    allocate(InVar%vdata(Nx*Ny*Nz))
  endif

  ! Set the output dimensions and coordinates to those of the selected input var
  call SetOutCoords(InFile, Xcoords, Ycoords, Zcoords, Tcoords)
  
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

  ! Prepare output variable

  OutVar%vname = InVar%vname
  OutVar%units = InVar%units
  OutVar%descrip = InVar%descrip
  OutVar%ndims = InVar%ndims
  OutVar%dims(1) = Nx
  OutVar%dims(2) = Ny
  OutVar%dimnames(1) = 'x'
  OutVar%dimnames(2) = 'y'
  if (OutVar%ndims .gt. 2) then
    OutVar%dims(3) = Nz
    OutVar%dimnames(3) = 'z'
  endif
  allocate(OutVar%vdata(Nx*Ny*Nz)) ! Nz will be 1 when have a 2D field

  ! allocate scratch arrays
  allocate(InField(Nx,Ny))
  allocate(OutField(Nx,Ny))

  ! Open the input files and the output file
  FileAcc = 'R'
  call rhdf5_open_file(InFile, FileAcc, 0, InFileId)
  write (*,*) 'Reading HDF5 file: ', trim(InFile)
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

    call rhdf5_read_variable(InFileId, InVar%vname, InVar%ndims, it, InVar%dims, rdata=InVar%vdata)

    do iz = 1, Nz
      ! Copy horizontal slice to scratch var (InData), do the processing, copy result
      ! into horizontal slice of output.

      ! copy from input var
      do iy = 1, Ny
        do ix = 1, Nx
          InField(ix,iy) = MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=InVar%vdata)
        enddo
      enddo

      ! For now, only support doing one pspec
      if (Pspecs(1)%Ptype .eq. 'gauss2d') then
        call Gsmooth2d(Nx, Ny, Pspecs(1)%i1, InField, Pspecs(1)%s1, OutField)
      endif

      ! copy to output var
      do iy = 1, Ny
        do ix = 1, Nx
          call MultiDimAssign(Nx, Ny, Nz, ix, iy, iz, OutField(ix,iy), Var3d=OutVar%vdata)
        enddo
      enddo
    enddo

    ! Write the filter data, and storm info if doing cylvol selection, to the output file
    call rhdf5_write_variable(OutFileId, OutVar%vname, OutVar%ndims, it, OutVar%dims, &
      OutVar%units, OutVar%descrip, OutVar%dimnames, rdata=OutVar%vdata)

    ! Write out status to screen every 100 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,100) .eq. 0) then
      write (*,*) 'Working: Timestep: ', it
    endif
  enddo
  write (*,*) ''

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
  call rhdf5_attach_dimensions(OutFile, OutVar)

  ! cleanup
  call rhdf5_close_file(OutFileId)
  call rhdf5_close_file(InFileId)

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

 subroutine GetMyArgs(MaxStr, MaxPspecs, InFile, InVarName, OutFile, Pspecs, NumPspecs)
  implicit none

  integer, parameter :: MaxFields = 15

  character (len=MaxStr) :: InFile, InVarName, OutFile, Pvname, Pfprefix
  integer :: MaxStr, MaxPspecs, NumPspecs
  type (PspecDescrip), dimension(MaxPspecs) :: Pspecs

  integer :: iargc, i, j, Nargs, Nfld
  character (len=MaxStr) :: Arg
  character (len=MaxStr), dimension(MaxFields) :: Fields

  logical :: BadArgs

  BadArgs = .false.
  InFile = ''
  InVarName = ''
  OutFile = ''
  NumPspecs = 0


  ! parse args
  i = 1
  Nargs = iargc()
  do while (i .le. Nargs)
    call getarg(i, Arg)

    if (i .eq. 1) then
      InFile = Arg
      i = i + 1
    else if (i .eq. 2) then
      InVarName = Arg
      i = i + 1
    else if (i .eq. 3) then
      OutFile = Arg
      i = i + 1
    else
      call String2List(Arg, ':', Fields, MaxFields, Nfld, 'sproc spec')
      i = i + 1

      !****************************
      !* Gaussian smoothing, 2D
      !****************************
      if (Fields(1) .eq. 'gauss2d') then
        NumPspecs = NumPspecs + 1

        if (Nfld .lt. 3) then
          write (*,*) 'ERROR: not enough arguments to fully specify the 2D Gaussian smoothing function'
          BadArgs = .true.
        else
          ! have enough args
          Pspecs(NumPspecs)%Ptype    = Fields(1)
          read(Fields(2),  '(i)') Pspecs(NumPspecs)%i1
          read(Fields(3),  '(f)') Pspecs(NumPspecs)%s1
        endif
      else
        write (*,*) 'ERROR: <ptype>, ', trim(Fields(1)), ', must be one of:'
        write (*,*) '          gauss2d'
        write (*,*) ''
        BadArgs = .true.
      endif
    endif
  enddo

  if (NumPspecs .lt. 1) then
    write(*,*) 'ERROR: Must give at least one processing spec'
    BadArgs = .true.
  endif
  
  if (BadArgs) then
    write (*,*) 'USAGE: sig_proc <in_file> <in_var> <out_file> <proc_spec> [<proc_spec>...]'
    write (*,*) '        <in_file>: input HDF5 file'
    write (*,*) '        <in_var>: variable (dataset) name inside the input HDF5 file'
    write (*,*) '        <out_file>: output HDF5 file name'
    write (*,*) '        <proc_spec>:<ptype>:<v1>:<v2>:<v3>:<v4>>'
    write (*,*) ''
    write (*,*) '            <ptype>:'
    write (*,*) '                gauss2d: 2D Gaussian smoothing'
    write (*,*) '                  <v1> is number of points (integer)'
    write (*,*) '                  <v2> is sigma (std dev) of filter (real)'
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
! Gaussian smoothing will be applied to the pressure field in order to help
! prevent mistakenly using topological features as the storm center.
!

subroutine FindStormCenter(Nx, Ny, Press, DataSelect, Npts, Sigma, StmCtrX, StmCtrY, MinP)
  implicit none

  integer :: Nx, Ny, Npts
  real, dimension(Nx,Ny) :: Press
  logical, dimension(Nx,Ny) :: DataSelect
  integer :: StmCtrX, StmCtrY
  real :: MinP, Sigma
  logical :: UseFsFilter

  integer :: ix, iy
  real, dimension(Nx,Ny) :: Psmooth

  ! apply 2D Gaussian smoothing to the pressure field
  call Gsmooth2d(Nx, Ny, Npts, Press, Sigma, Psmooth)

  MinP = 1e10 ! ridiculously large pressure
  StmCtrX = 0 
  StmCtrY = 0 

  do ix = 1, Nx
    do iy = 1, Ny
      ! DataSelect holds true for points that we want to consider
      ! for minimum pressure
      if (DataSelect(ix,iy)) then 
        if (Psmooth(ix,iy) .lt. MinP) then
          MinP = Psmooth(ix,iy)
          StmCtrX = ix
          StmCtrY = iy
        endif
      endif
    enddo
  enddo
  
  return
end subroutine FindStormCenter

end program sig_proc
