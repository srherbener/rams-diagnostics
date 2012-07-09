!***************************************************************
! Program to take a list of files and join the HDF5 data within
! them together.
!
! Args
!   1. input HDF5 file names (control files, colon separated list)
!   2. output HDF5 file name
!

program join_hdata
  use rhdf5_utils
  use diag_utils
  implicit none

  integer, parameter :: LargeString  = 512
  integer, parameter :: MediumString = 256
  integer, parameter :: LittleString = 128
  integer, parameter :: MaxFiles     = 10

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OutFile
  character (len=LittleString) :: VarName

  character (len=MediumString), dimension(MaxFiles) :: InFileList
  integer :: Nfiles

  integer, dimension(MaxFiles) :: InSizes
  integer :: Nx, Ny, Nz, Nt
  integer :: InNx, InNy, InNz, InNt

  ! Data arrays
  ! OutVar dims: x, y, z, t
  ! Invars, array where each element has a data array with dims: x, y, z, t
  type (Rhdf5Var), dimension(MaxFiles) :: InVars
  type (Rhdf5Var), dimension(RHDF5_MAX_DIMS) :: Coords
  type (Rhdf5Var) :: InTcoords
  type (Rhdf5Var) :: OutVar

  integer :: i, iin, iout, itin, itout
  integer :: Ndims
  character (len=LittleString), dimension(RHDF5_MAX_DIMS) :: Cnames
  character (len=LittleString), dimension(RHDF5_MAX_DIMS) :: GradsCnames
  
  real :: Xstart, Xinc, Ystart, Yinc
  logical :: BadData

  ! Get the command line arguments
  call GetMyArgs(Infiles, OutFile, VarName)
  call String2List(Infiles, ':', InFileList, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Joining HDF5 data:'
  write (*,*) '  HDF5 input files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(InFileList(i))
  end do
  write (*,*) '  Variable: ', trim(VarName)
  write (*,*) ''

  ! For this function we need to locate the variable in every input file since
  ! we want to join these together. Use arrayed versions of the GradsControlFiles type
  ! and GradVar type. Walk though each individual input file given and collect the variable
  ! information from just that particular file. Store each result in the elements of the
  ! arrayed versions of GradsControlFiles and GradsVar types.
  do i = 1, Nfiles
    InVars(i)%vname = VarName
    call rhdf5_read_init(InFileList(i), InVars(i))

    ! Do some dimension checking
    ! Assume that the variable dimension orders are:
    !   1D: x, t
    !   2D: x, y, t
    !   2D: x, y, t
    !   3D: x, y, z, t
    !   Note that t is last and it is the dimension that we are joining 
    !   the pieces together on.
    !
    ! Grab the dimensions of the current input file. Record these if the first
    ! file, and if not the first file compare the dimensions to the ones just recorded.
    ! Record the total data size (number of elements) for each of the input files.
    ! Sum up the number of time steps in Nt however so that the recored dimensions can
    ! be used to describe the output data size.
    !
    ! Then use the recorded dimensions for figuring the size of the output variable.

    if (InVars(i)%ndims .eq. 2) then
      ! 1D data plus time
      InNx = InVars(i)%dims(1)
      InNy = 1 ! use one so that the total number of elements gets calculated correctly
      InNz = 1
      InNt = InVars(i)%dims(2)
    elseif (InVars(i)%ndims .eq. 3) then
      ! 2D data plus time
      InNx = InVars(i)%dims(1)
      InNy = InVars(i)%dims(2)
      InNz = 1 ! use one so that the total number of elements gets calculated correctly
      InNt = InVars(i)%dims(3)
    else
      ! 3D data plus time
      InNx = InVars(i)%dims(1)
      InNy = InVars(i)%dims(2)
      InNz = InVars(i)%dims(3)
      InNt = InVars(i)%dims(4)
    endif

    InSizes(i) = InNx * InNy * InNz * InNt

    if (i .eq. 1) then
      ! if first variable, then record the dimensions for comparing
      ! and initialze the coordinate datasets
      Nx = InNx
      Ny = InNy
      Nz = InNz
      Nt = InNt
    else
      ! if not first variable
      !   compare all but the time dimensions (since these can change)
      !     unused dimensions (Nz for 2D data for example) will always be set to 1 and
      !     will compare okay
      !   accumulate the number of time steps
      !
      BadData = .false.
      if (Nx .ne. InNx) then
        write (*,*) 'ERROR: number of x points in GRADS control files do not match'
        BadData = .true.
      endif
      if (Ny .ne. InNy) then
        write (*,*) 'ERROR: number of y points in GRADS control files do not match'
        BadData = .true.
      endif
      if (Nz .ne. InNz) then
        write (*,*) 'ERROR: number of z points in GRADS control files do not match'
        BadData = .true.
      endif
      if (BadData) then
        write (*,*) 'ERROR:  Control file: ', trim(InFileList(i))
        stop
      endif

      ! Accumulate the number of time steps so we know how many time steps will exist
      ! in the output data
      Nt = Nt + InNt
    endif
  enddo

  Ndims = InVars(1)%ndims

  write (*,*) 'Variable dimensions after joining files:'
  write (*,*) '  Ndims: ', Ndims
  write (*,*) '  Nx: ', Nx
  write (*,*) '  Ny: ', Ny
  write (*,*) '  Nz: ', Nz
  write (*,*) '  Nt: ', Nt
  write (*,*) ''
  flush(6)

  ! Generate the names of the dimension coordinate variables
  do i = 1, Ndims
    Cnames(i) = trim(InVars(1)%dimnames(i)) // '_coords'
    if (i .eq. 1) then
      GradsCnames(i) = 'x'
    else if (i .eq. 2) then
      GradsCnames(i) = 'y'
    else if (i .eq. 3) then
      GradsCnames(i) = 'z'
    else
      GradsCnames(i) = 't'
    endif
  enddo 
  ! the last dimension is always 't'
  GradsCnames(Ndims) = 't'

  ! Read in all but the last dimension from the first input file
  !   Last dimension is 't' which will be handled below
  do i = 1, Ndims-1
    Coords(i)%vname = Cnames(i)
    call rhdf5_read_init(InFileList(1), Coords(i))
    allocate (Coords(i)%vdata(Coords(i)%dims(1)))
    call rhdf5_read(InFileList(1), Coords(i))
  enddo

  ! Initialize the output variable, go through the list of input variables and read
  ! in their data and append it to the output variable
  ! We have verified that the x, y, z dimensions of the input data in all HDF5 input
  ! files match so it's okay use InVars(1) as a representative.
  !
  OutVar%vname = VarName
  OutVar%ndims = InVars(1)%ndims
  OutVar%dims = InVars(1)%dims
  OutVar%dimnames = InVars(1)%dimnames
  ! The last dimension size will be wrong since we copied it from InVars(1). It needs
  ! to be replaced with Nt.
  OutVar%dims(Invars(1)%ndims) = Nt
  OutVar%units = Invars(1)%units
  OutVar%descrip = Invars(1)%descrip
  allocate(OutVar%vdata(Nx*Ny*Nz*Nt))

  ! Set up the 't' coordinate variable based on contents of first input file
  ! except make the size equal to Nt
  Coords(Ndims)%vname = Cnames(Ndims)
  call rhdf5_read_init(InFileList(1), Coords(Ndims))
  Coords(Ndims)%dims(1) = Nt
  allocate (Coords(Ndims)%vdata(Nt))

  iout = 0
  itout = 0
  do i = 1, Nfiles
    ! append the variable data
    allocate(InVars(i)%vdata(InSizes(i)))
    call rhdf5_read(InFileList(i), InVars(i))
    
    do iin = 1, InSizes(i)
      iout = iout + 1
      OutVar%vdata(iout) = InVars(i)%vdata(iin)
    end do

    deallocate(InVars(i)%vdata)

    ! append the 't' coordinate variable data
    InTcoords%vname = Cnames(Ndims)
    call rhdf5_read_init(InFileList(i), InTcoords)
    allocate(InTcoords%vdata(InTcoords%dims(1)))
    call rhdf5_read(InFileList(i), InTcoords)

    do itin = 1, InTcoords%dims(1)
      itout = itout + 1
      Coords(Ndims)%vdata(itout) = InTcoords%vdata(itin)
    enddo

    deallocate(InTcoords%vdata)
    
  end do

  !Write out the joined data
  call rhdf5_write(OutFile, OutVar, 0)

  ! Write out the coordinate data and mark as such
  ! variable
  do i = 1, Ndims
    call rhdf5_write(OutFile, Coords(i), 1)
    call rhdf5_set_dimension(OutFile, Coords(i), GradsCnames(i))
  enddo

  ! attach the coordinates to the output variable
  call rhdf5_attach_dimensions(OutFile, OutVar)

  deallocate(OutVar%vdata)
!  deallocate(InTcoords)
  do i = 1, Ndims
    deallocate(Coords(i)%vdata)
  enddo

  stop

contains

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   OutFile - output GRADS file, base name for two files
!

subroutine GetMyArgs(Infiles, OutFile, VarName)
  implicit none

  integer, parameter :: MAX_ITEMS = 5

  character (len=*) :: Infiles, OutFile, VarName

  integer :: iargc
  character (len=128) :: arg
  character (len=128), dimension(MAX_ITEMS) :: ArgList
  integer :: Nitems

  logical :: BadArgs

  if (iargc() .ne. 3) then
    write (*,*) 'ERROR: must supply exactly 3 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_file> <variable_name>'
    write (*,*) '        <in_data_files>: HDF5 format, control file, colon separated list'
    write (*,*) '        <out_data_file>: HDF5 format'
    write (*,*) '        <variable_name>: HDF5 variable name, this program will look for this variable in all of the <in_data_files>'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OutFile)
  call getarg(3, VarName)

  return
end subroutine

end program join_hdata
