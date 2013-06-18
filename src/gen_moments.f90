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

program gen_moments
  use rhdf5_utils
  use diag_utils
  implicit none

  integer, parameter :: LargeString  = 512
  integer, parameter :: MediumString = 256
  integer, parameter :: LittleString = 128

  character (len=LargeString) :: InFile, InVarName, OutFile
  integer :: TsStart, TsEnd

  integer :: ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt
  integer :: Ntsteps

  type (Rhdf5Var) :: InVar

  character (len=RHDF5_MAX_STRING) :: FileAcc

  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm

  real :: DeltaX, DeltaY

  integer :: InFileId, OutFileId

  type (Rhdf5Var) :: VarM1, VarM2, VarM3, VarM4 ! first four moments
  type (Rhdf5Var) :: Npts
  real :: Vdiff, Vterm

  ! Get the command line arguments
  call GetMyArgs(LargeString, InFile, InVarName, OutFile, TsStart, TsEnd)

  write (*,*) 'Generating moments:'
  write (*,*) '  Input file: ', trim(InFile)
  write (*,*) '  Input variable name: ', trim(InVarName)
  write (*,*) '  Output file:  ', trim(OutFile)
  write (*,*) '  Beginning time step: ', TsStart
  write (*,*) '  Ending time step: ', TsEnd
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

  ! Prepare output variables
  Npts%vname = 'num_points'
  Npts%units = 'count'
  Npts%descrip = 'number of points selected'
  Npts%ndims = 1
  Npts%dims(1) = Nz
  Npts%dimnames(1) = 'z'
  allocate(Npts%vdata(Nz))

  VarM1%vname = trim(InVar%vname) // '_M1'
  VarM1%units = trim(InVar%units)
  VarM1%descrip = 'first moment ' // trim(InVar%descrip)
  VarM1%ndims = 1
  VarM1%dims(1) = Nz
  VarM1%dimnames(1) = 'z'
  allocate(VarM1%vdata(Nz))

  VarM2%vname = trim(InVar%vname) // '_M2'
  VarM2%units = '(' // trim(InVar%units) // ')^2'
  VarM2%descrip = 'second moment ' // trim(InVar%descrip)
  VarM2%ndims = 1
  VarM2%dims(1) = Nz
  VarM2%dimnames(1) = 'z'
  allocate(VarM2%vdata(Nz))

  VarM3%vname = trim(InVar%vname) // '_M3'
  VarM3%units = '(' // trim(InVar%units) // ')^3'
  VarM3%descrip = 'third moment ' // trim(InVar%descrip)
  VarM3%ndims = 1
  VarM3%dims(1) = Nz
  VarM3%dimnames(1) = 'z'
  allocate(VarM3%vdata(Nz))

  VarM4%vname = trim(InVar%vname) // '_M4'
  VarM4%units = '(' // trim(InVar%units) // ')^4'
  VarM4%descrip = 'fourth moment ' // trim(InVar%descrip)
  VarM4%ndims = 1
  VarM4%dims(1) = Nz
  VarM4%dimnames(1) = 'z'
  allocate(VarM4%vdata(Nz))

  ! Open the input files and the output file
  FileAcc = 'R'
  call rhdf5_open_file(InFile, FileAcc, 0, InFileId)
  write (*,*) 'Reading HDF5 file: ', trim(InFile)
  write (*,*) ''

  FileAcc = 'W'
  call rhdf5_open_file(OutFile, FileAcc, 1, OutFileId)
  write (*,*) 'Writing HDF5 file: ', trim(OutFile)
  write (*,*) ''

  ! The necessary input files have been opened, data buffers allocated,
  ! plus the time dimensions have been "removed" from the descriptions
  ! of the input variables.
  !
  ! Need to make two passes through the data:
  !   Pass 1 -> compute the mean (1st moment)
  !   Pass 2 -> compute the higher moments

  ! Pass 1 -> mean (first moment)
  do iz = 1, Nz
    VarM1%vdata(iz) = 0.0
    Npts%vdata(iz) = 0.0
  enddo
  write (*,*) 'Pass 1 - calculating mean'
  do it = TsStart, TsEnd
    call rhdf5_read_variable(InFileId, InVar%vname, InVar%ndims, it, InVar%dims, rdata=InVar%vdata)

    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          VarM1%vdata(iz) = VarM1%vdata(iz) + MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=InVar%vdata)
          Npts%vdata(iz) = Npts%vdata(iz) + 1.0
        enddo
      enddo
    enddo
  enddo

  do iz = 1, Nz 
    VarM1%vdata(iz) = VarM1%vdata(iz) / Npts%vdata(iz)
  enddo
  write (*,*) '  Done!'
  write (*,*) ''

  ! Pass 2 --> higher moments
  Ntsteps = 0
  write (*,*) 'Pass 2 - higher moments'
  do it = TsStart, TsEnd
    call rhdf5_read_variable(InFileId, InVar%vname, InVar%ndims, it, InVar%dims, rdata=InVar%vdata)

    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          ! Form (X-Xmean)
          Vdiff = MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=InVar%vdata) - VarM1%vdata(iz)

          ! second moment: sum up (X-Xmean)^2
          Vterm = Vdiff * Vdiff
          VarM2%vdata(iz) = VarM2%vdata(iz) + Vterm

          ! third moment: sum up (X-Xmean)^3
          Vterm = Vterm * Vdiff
          VarM3%vdata(iz) = VarM3%vdata(iz) + Vterm

          ! fourth moment: sum up (X-Xmean)^4
          Vterm = Vterm * Vdiff
          VarM4%vdata(iz) = VarM4%vdata(iz) + Vterm
        enddo
      enddo
    enddo

    ! Write out status to screen every 100 timesteps so that the user can see that a long
    ! running job is progressing okay.
    Ntsteps = Ntsteps + 1
    if (modulo(Ntsteps,100) .eq. 0) then
      write (*,*) 'Working: Timestep: ', it
    endif
  enddo

  do iz = 1, Nz 
    VarM2%vdata(iz) = VarM2%vdata(iz) / Npts%vdata(iz)
    VarM3%vdata(iz) = VarM3%vdata(iz) / Npts%vdata(iz)
    VarM4%vdata(iz) = VarM4%vdata(iz) / Npts%vdata(iz)
  enddo
  write (*,*) '  Done!'
  write (*,*) ''

! Write out the moment data
  call rhdf5_write_variable(OutFileId, Npts%vname, Npts%ndims, 0, Npts%dims, &
     Npts%units, Npts%descrip, Npts%dimnames, rdata=Npts%vdata)
  call rhdf5_write_variable(OutFileId, VarM1%vname, VarM1%ndims, 0, VarM1%dims, &
     VarM1%units, VarM1%descrip, VarM1%dimnames, rdata=VarM1%vdata)
  call rhdf5_write_variable(OutFileId, VarM2%vname, VarM2%ndims, 0, VarM2%dims, &
     VarM2%units, VarM2%descrip, VarM2%dimnames, rdata=VarM2%vdata)
  call rhdf5_write_variable(OutFileId, VarM3%vname, VarM3%ndims, 0, VarM3%dims, &
     VarM3%units, VarM3%descrip, VarM3%dimnames, rdata=VarM3%vdata)
  call rhdf5_write_variable(OutFileId, VarM4%vname, VarM4%ndims, 0, VarM4%dims, &
     VarM4%units, VarM4%descrip, VarM4%dimnames, rdata=VarM4%vdata)

  write (*,*) ''

  write (*,*) 'Finished: Total number of time steps processed: ', Ntsteps
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
  call rhdf5_attach_dimensions(OutFile, Npts)
  call rhdf5_attach_dimensions(OutFile, VarM1)
  call rhdf5_attach_dimensions(OutFile, VarM2)
  call rhdf5_attach_dimensions(OutFile, VarM3)
  call rhdf5_attach_dimensions(OutFile, VarM4)

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

 subroutine GetMyArgs(MaxStr, InFile, InVarName, OutFile, TsStart, TsEnd)
  implicit none

  integer :: MaxStr, TsStart, TsEnd
  character (len=MaxStr) :: InFile, InVarName, OutFile

  integer :: iargc, i, j, Nargs
  character (len=MaxStr) :: Arg

  logical :: BadArgs

  BadArgs = .false.

  ! parse args
  i = 1
  Nargs = iargc()

  if (Nargs .eq. 5) then

    call getarg(1, Arg)
    InFile = Arg

    call getarg(2, Arg)
    InVarName = Arg

    call getarg(3, Arg)
    OutFile = Arg

    call getarg(4, Arg)
    read(Arg, '(i)') TsStart

    call getarg(5, Arg)
    read(Arg, '(i)') TsEnd

    if (TsEnd .lt. TsStart) then
      write (*,*) 'ERROR: <ts_end> must be greater than or equal to <ts_start>'
      BadArgs = .true.
    endif
  else
    write (*,*) 'ERROR: Must supply exactly 5 arguments'
    BadArgs = .true.
  endif
  
  if (BadArgs) then
    write (*,*) 'USAGE: gen_moments <in_file> <in_var> <out_file> <ts_start> <ts_end>'
    write (*,*) '        <in_file>: input HDF5 file'
    write (*,*) '        <in_var>: variable (dataset) name inside the input HDF5 file'
    write (*,*) '        <out_file>: output HDF5 file name'
    write (*,*) '        <ts_start>: beginning time step number'
    write (*,*) '        <ts_end>: end time step number'
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

end program gen_moments
