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
  integer, parameter :: MaxVars = 5

  type VarSpec
    character (len=LittleString) :: Vname
    character (len=LittleString) :: Vfprefix
  end type VarSpec

  character (len=LargeString) :: InDir, InSuffix, OutFile
  integer :: TsStart, TsEnd
  type (VarSpec), dimension(MaxVars) :: VarList
  integer :: Nvars

  logical :: BadDims, ZdimMatches

  integer :: iv, ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt
  integer :: Ntsteps

  character (len=LargeString), dimension(MaxVars) :: InFiles, InVarNames
  type (Rhdf5Var), dimension(MaxVars) :: InVars
  integer, dimension(MaxVars) :: InFileIds

  character (len=RHDF5_MAX_STRING) :: FileAcc

  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm

  real :: DeltaX, DeltaY

  integer :: OutFileId
  type (Rhdf5Var) :: OutVar, Npts
  character (len=RHDF5_MAX_STRING) :: OutVarName, OutVarUnits, OutVarDescrip

  real, dimension(:,:), allocatable :: Means
  real :: Vdiff, Vterm

  ! Get the command line arguments
  call GetMyArgs(LargeString, MaxVars, InDir, InSuffix, OutFile, TsStart, TsEnd, VarList, Nvars)
  do iv = 1, Nvars
    InFiles(iv) = trim(InDir) // '/' // trim(VarList(iv)%Vfprefix) // trim(InSuffix)
    InVarNames(iv) = trim(VarList(iv)%Vname)
  enddo

  write (*,*) 'Generating moments:'
  write (*,*) '  Input variables:'
  do iv = 1, Nvars
    write (*,*) '    File: ', trim(InFiles(iv))
    write (*,*) '    Variable: ', trim(InVarNames(iv))
  enddo
  write (*,*) '  Output file:  ', trim(OutFile)
  write (*,*) '  Beginning time step: ', TsStart
  write (*,*) '  Ending time step: ', TsEnd
  write (*,*) ''
  flush(6)

  ! Since we will be processing one time step at a time, the effective number of dimensions
  ! on each of the input is decreased by one. We want to eliminate the time dimension which
  ! is always the last one. This works out conveniently since all we need to do is decrement
  ! the number of dimensions by one.

  BadDims = .false.
  do iv = 1, Nvars
    InVars(iv)%vname = trim(InVarNames(iv))
    call rhdf5_read_init(InFiles(iv), InVars(iv))

    if (iv .eq. 1) then
      ! If first variable, record the dimensions
      Nx = InVars(iv)%dims(1)
      Ny = InVars(iv)%dims(2)
      if (InVars(iv)%ndims .eq. 3) then
        ! 2D field
        Nz = 1
        Nt = InVars(iv)%dims(3)
      else
        ! 3D field
        Nz = InVars(iv)%dims(3)
        Nt = InVars(iv)%dims(4)
      endif
    else
      ! If beyond the first variable, check to make sure its dimensions match
      ! The DimsMatch() function only checks horzontal and time dimensions. We also want to
      ! force the z dimensions to match as well.
      if (Nz .eq. 1) then
        ! 2D vars
        ZdimMatches = InVars(iv)%ndims .eq. 3
      else
        ! 3D vars
        ZdimMatches = Nz .eq. InVars(iv)%dims(3)
      endif

      if (DimsMatch(InVars(1), InVars(iv)) .and. ZdimMatches) then
      else
        write (*,*) 'ERROR: dimensions of variable do not match the first variable: ', trim(InVars(iv)%vname)
        BadDims = .true.
      endif
    endif
  enddo

  if (BadDims) then
    stop
  endif

  ! Prepare for reading
  ! Create names for ouput variable
  do iv = 1, Nvars
    InVars(iv)%ndims = InVars(iv)%ndims - 1
    if (InVars(iv)%ndims .eq. 2) then
      allocate(InVars(iv)%vdata(Nx*Ny))
    else
      allocate(InVars(iv)%vdata(Nx*Ny*Nz))
    endif
 
    if (iv .eq. 1) then
      OutVarName = trim(InVars(iv)%vname)
      OutVarUnits = '(' // trim(InVars(iv)%units) // ')'
      OutVarDescrip = 'moment: (' // trim(InVars(iv)%vname) // ')'
    else
      OutVarName = trim(OutVarName) // '-' // trim(InVars(iv)%vname)
      OutVarUnits = trim(OutVarUnits) // '(' // trim(InVars(iv)%units) // ')'
      OutVarDescrip = trim(OutVarDescrip) // '(' // trim(InVars(iv)%vname) // ')'
    endif
  enddo

  ! Set the output dimensions and coordinates to those of the selected input var
  call SetOutCoords(InFiles(1), Xcoords, Ycoords, Zcoords, Tcoords)

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

  OutVar%vname = trim(OutVarName)
  OutVar%units = trim(OutVarUnits)
  OutVar%descrip = trim(OutVarDescrip)
  OutVar%ndims = 1
  OutVar%dims(1) = Nz
  OutVar%dimnames(1) = 'z'
  allocate(OutVar%vdata(Nz))

  ! Open the input files and the output file
  FileAcc = 'R'
  do iv = 1, Nvars
    call rhdf5_open_file(InFiles(iv), FileAcc, 0, InFileIds(iv))
    write (*,*) 'Reading HDF5 file: ', trim(InFiles(iv))
  enddo
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
  !   Pass 2 -> compute the higher moments, if asked for

  ! Pass 1 -> mean (first moment)
  !   record the means of each variable separately, however the number of points will
  !   be the same for each variable so only count those up during the first variable
  !   Count the number of points in a loop (as opposed to multiplying Nx * Ny * NumberTimeSteps)
  !   so that filtering can be added later on.
  allocate(Means(Nvars, Nz))
  do iz = 1, Nz
    Npts%vdata(iz) = 0.0
    do iv = 1, Nvars
      Means(iv,iz) = 0.0
    enddo
  enddo
  write (*,*) 'Pass 1 - calculating mean'
  Ntsteps = 0
  do it = TsStart, TsEnd
    do iv = 1, Nvars
      call rhdf5_read_variable(InFileIds(iv), InVars(iv)%vname, InVars(iv)%ndims, it, InVars(iv)%dims, rdata=InVars(iv)%vdata)

      ! calculte mean for each level
      do iz = 1, Nz
        do iy = 1, Ny
          do ix = 1, Nx
            Means(iv,iz) = Means(iv,iz) + MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=InVars(iv)%vdata)
            ! Only count Npts during first variable
            if (iv .eq. 1) then
              Npts%vdata(iz) = Npts%vdata(iz) + 1.0
            endif
          enddo
        enddo
      enddo

      ! free up the memory that rhdf5_read_variable allocated
      deallocate(InVars(iv)%vdata)
    enddo

    ! Write out status to screen every 100 timesteps so that the user can see that a long
    ! running job is progressing okay.
    Ntsteps = Ntsteps + 1
    if (modulo(Ntsteps,100) .eq. 0) then
      write (*,*) 'Working: Timestep: ', it
    endif
  enddo

  do iv = 1, Nvars
    do iz = 1, Nz 
      Means(iv,iz) = Means(iv,iz) / Npts%vdata(iz)
    enddo
  enddo
  write (*,*) '  Done!'
  write (*,*) ''

  ! Pass 2 --> higher moments
  !   If there is only one variable, copy the mean values into OutVar
  !   Otherwise go through and compute the higher moments
  if (Nvars .eq. 1) then
    do iz = 1, Nz
      OutVar%vdata(iz) = Means(1,iz)
    enddo 
  else
    do iz = 1, Nz
      OutVar%vdata(iz) = 0.0
    enddo 

    Ntsteps = 0
    write (*,*) 'Pass 2 - higher moments'
    do it = TsStart, TsEnd
      ! Read in all variables
      do iv = 1, Nvars
        call rhdf5_read_variable(InFileIds(iv), InVars(iv)%vname, InVars(iv)%ndims, it, InVars(iv)%dims, rdata=InVars(iv)%vdata)
      enddo
  
      do iz = 1, Nz
        do iy = 1, Ny
          do ix = 1, Nx

            ! Form (V1 - V1bar)*(V2 - V2bar)*...
            Vterm = 1.0
            do iv = 1, Nvars
              Vdiff = MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=InVars(iv)%vdata) - Means(iv,iz)
              Vterm = Vterm * Vdiff
            enddo
  
            OutVar%vdata(iz) = OutVar%vdata(iz) + Vterm
          enddo
        enddo
      enddo

      ! free up the memory that rhdf5_read_variable allocated
      do iv = 1, Nvars
        deallocate(InVars(iv)%vdata)
      enddo
  
      ! Write out status to screen every 100 timesteps so that the user can see that a long
      ! running job is progressing okay.
      Ntsteps = Ntsteps + 1
      if (modulo(Ntsteps,100) .eq. 0) then
        write (*,*) 'Working: Timestep: ', it
      endif
    enddo
  
    do iz = 1, Nz 
      OutVar%vdata(iz) = OutVar%vdata(iz) / Npts%vdata(iz)
    enddo
  endif
  write (*,*) '  Done!'
  write (*,*) ''

  ! Write out the moment data
  call rhdf5_write_variable(OutFileId, Npts%vname, Npts%ndims, 0, Npts%dims, &
     Npts%units, Npts%descrip, Npts%dimnames, rdata=Npts%vdata)
  call rhdf5_write_variable(OutFileId, OutVar%vname, OutVar%ndims, 0, OutVar%dims, &
     OutVar%units, OutVar%descrip, OutVar%dimnames, rdata=OutVar%vdata)

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
  call rhdf5_attach_dimensions(OutFile, OutVar)

  ! cleanup
  call rhdf5_close_file(OutFileId)
  do iv = 1, Nvars
    call rhdf5_close_file(InFileIds(iv))
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

 subroutine GetMyArgs(MaxStr, MaxVars, InDir, InSuffix, OutFile, TsStart, TsEnd, VarList, NumVars)
  implicit none

  integer, parameter :: MaxFields = 15

  integer :: MaxStr, MaxVars, TsStart, TsEnd, NumVars
  character (len=MaxStr) :: InDir, InSuffix, OutFile
  type (VarSpec), dimension(MaxVars) :: VarList

  integer :: iargc, i, Nargs, Nfld
  character (len=MaxStr) :: Arg
  character (len=MaxStr), dimension(MaxFields) :: Fields

  logical :: BadArgs

  BadArgs = .false.
  InDir = ''
  InSuffix = ''
  OutFile = ''
  TsStart = 0
  TsEnd = 0
  NumVars = 0

  ! parse args
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
      read(Arg, '(i)') TsStart
      i = i + 1
    else if (i .eq. 5) then
      read(Arg, '(i)') TsEnd
      i = i + 1
    else
      call String2List(Arg, ':', Fields, MaxFields, Nfld, 'var spec')
      i = i + 1

      NumVars = NumVars + 1
      if (NumVars .gt. MaxVars) then
        write (*,'(a,i1,a)') 'ERROR: Exceeded maximum number of times <var_spec> can be used (', MaxVars, ')'
        BadArgs = .true.
      else
        VarList(NumVars)%Vname = Fields(1)
        VarList(NumVars)%Vfprefix = Fields(2)
      endif
    endif
  enddo

  if (TsEnd .lt. TsStart) then
    write (*,*) 'ERROR: <ts_end> must be greater than or equal to <ts_start>'
    BadArgs = .true.
  endif
  
  if (NumVars .lt. 1) then
    write (*,*) 'ERROR: must specify at least one <var_spec>'
    BadArgs = .true.
  endif
  
  if (BadArgs) then
    write (*,*) 'USAGE: gen_moments <in_dir> <in_suffix> <out_file> <ts_start> <ts_end> <var_spec> [<var_spec>...]'
    write (*,*) '        <in_dir>: directory containing input HDF5 files'
    write (*,*) '        <in_suffix>: suffix attached to all input HDF5 files'
    write (*,*) '        <out_file>: output HDF5 file name'
    write (*,*) '        <ts_start>: beginning time step number'
    write (*,*) '        <ts_end>: end time step number'
    write (*,*) '        <var_spec> = <var_name>:<file_prefix>'
    write (*,*) '          <var_name>: name of variable inside input HDF5 file'
    write (*,*) '          <file_prefix>: leading name of input HDF5 file'
    write (*,*) '              input file name gets built by joining: <file_prefix><in_suffix>'
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
