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
  integer :: Nfiles, Nt
  integer, dimension(MaxFiles) :: Ntsteps

  integer :: rh5_file, rh5_out_file
  character (len=RHDF5_MAX_STRING) :: rh5_file_acc

  ! Data arrays
  ! Invars, array where each element has a data array with dims: x, y, z, t
  type (Rhdf5Var), dimension(RHDF5_MAX_DIMS) :: Coords
  type (Rhdf5Var) :: InTcoords

  integer :: i, id, it, itin, iout, itout
  integer :: ndims
  integer, dimension(RHDF5_MAX_DIMS) :: dims, check_dims
  integer :: InSize
  real, dimension(:), allocatable :: VarData
  character (len=RHDF5_MAX_STRING) :: units, descrip
  character (len=RHDF5_MAX_STRING), dimension(RHDF5_MAX_DIMS) :: dimnames
  character (len=LittleString), dimension(RHDF5_MAX_DIMS) :: Cnames
  
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
    ! Read the file meta-data
    rh5_file_acc = 'R'
    call rhdf5_open_file(InFileList(i), rh5_file_acc, 0, rh5_file)
    call rhdf5_read_variable_init(rh5_file, VarName, ndims, 0, dims, units, descrip, dimnames)
    call rhdf5_close_file(rh5_file)

    ! record the number of times steps for each input file
    Ntsteps(i) = dims(ndims)

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

    if (i .eq. 1) then
      ! if first variable, then record the dimensions for comparing
      ! and figure out the size of the spatial field data
      InSize = 1
      do id = 1, ndims-1
        InSize = InSize * dims(id)
        check_dims(id) = dims(id)
      enddo
      Nt = dims(ndims)
    else
      ! if not first variable
      !   compare all but the time dimensions (since these can change)
      BadData = .false.
      do id = 1, ndims - 1
        if (dims(id) .ne. check_dims(id)) then
          write (*,*) 'ERROR: dimensions in the input files do not match:'
          write (*,*) 'ERROR:   input file with mis-match: ', trim(InFileList(i))
          write (*,*) 'ERROR:   dimension number with mis-match: ', id
          write (*,*) 'ERROR:   this file dimension size: ', dims(id)
          write (*,*) 'ERROR:   first file dimension size: ', check_dims(id)
          BadData = .true.
        endif
      enddo

      if (BadData) then
        write (*,*) 'ERROR:  Input file: ', trim(InFileList(i))
        stop
      endif

      ! calculate total number of time steps
      Nt = Nt + dims(ndims)
    endif
  enddo

  write (*,*) 'Variable information:'
  write (*,*) '  Variable name: ', trim(VarName)
  write (*,*) '  Variable units: ', trim(units)
  write (*,*) '  Variable description: ', trim(descrip)
  write (*,*) '  Number of dimensions: ', ndims-1
  write (*,*) '  Dimension sizes:'
  do id = 1, ndims-1
    write (*,*) '    ', trim(dimnames(id)), ': ', dims(id)
  enddo
  write (*,*) 'Number of time steps per input file:'
  do i = 1, Nfiles
    write (*,*) '    File number: ', i, ': ', Ntsteps(i)
  enddo
  write (*,*) ''
  write (*,*) 'Total number of time steps: ', Nt
  write (*,*) ''
  flush(6)

  ! Generate the names of the dimension coordinate variables
  do i = 1, ndims
    Cnames(i) = trim(dimnames(i)) // '_coords'
  enddo 

  ! Read in all but the last dimension from the first input file
  !   Last dimension is 't' which will be handled below
  do i = 1, ndims-1
    Coords(i)%vname = Cnames(i)
    call rhdf5_read_init(InFileList(1), Coords(i))
    call rhdf5_read(InFileList(1), Coords(i))
  enddo

  ! Set up the 't' coordinate variable based on contents of first input file
  ! except make the size equal to Nt
  Coords(ndims)%vname = Cnames(ndims)
  call rhdf5_read_init(InFileList(1), Coords(ndims))
  Coords(ndims)%dims(1) = Nt
  allocate (Coords(ndims)%vdata(Nt))

  ! Copy the input data fields and t coordinates to the output
  iout = 0
  itout = 0
  rh5_file_acc = 'W'
  call rhdf5_open_file(OutFile, rh5_file_acc, 1, rh5_out_file)

  do i = 1, Nfiles
    ! append the variable data
    rh5_file_acc = 'R'
    call rhdf5_open_file(InFileList(i), rh5_file_acc, 0, rh5_file)
    do it = 1, Ntsteps(i)
      call rhdf5_read_variable(rh5_file, VarName, ndims-1, it, dims, rdata=VarData)

      ! copy to the next time step in the output file
      iout = iout + 1
      call rhdf5_write_variable(rh5_out_file, VarName, ndims-1, iout, dims, &
        units, descrip, dimnames, rdata=VarData)

      if (modulo(iout,100) .eq. 0) then
        write (*,*) 'Working: Number of time steps transferred so far: ', iout
      endif

      deallocate(VarData)
    enddo
    call rhdf5_close_file(rh5_file);

    ! append the 't' coordinate variable data
    InTcoords%vname = Cnames(ndims)
    call rhdf5_read_init(InFileList(i), InTcoords)
    call rhdf5_read(InFileList(i), InTcoords)

    do itin = 1, InTcoords%dims(1)
      itout = itout + 1
      Coords(ndims)%vdata(itout) = InTcoords%vdata(itin)
    enddo

    deallocate(InTcoords%vdata)
  end do

  call rhdf5_close_file(rh5_out_file);
  write (*,*) 'Finished: Total number of time steps transferred: ', iout
  write (*,*) ''

  ! Write out the coordinate data and mark as such
  ! variable
  do i = 1, ndims
    call rhdf5_write(OutFile, Coords(i), 1)
    call rhdf5_set_dimension(OutFile, Coords(i), dimnames(i))
  enddo

  ! attach the coordinates to the output variable
  !call rhdf5_attach_dimensions(OutFile, OutVar)
  !deallocate(OutVar%vdata)
  rh5_file_acc = 'RW' 
  call rhdf5_open_file(OutFile, rh5_file_acc, 1, rh5_out_file)
  call rhdf5_attach_dims_to_var(rh5_out_file, VarName)
  call rhdf5_close_file(rh5_out_file)

  do i = 1, ndims
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
