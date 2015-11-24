!***************************************************************
! Program to do averaging over the domain and yield a single
! data point at each time step
!
! This program will read in HDF5 data from a RAMS simulation, 
! perform an averaging function, and output the single point time
! series in HDF5 format
!

program hdata_op
  use rhdf5_utils
  use diag_utils
  implicit none

  ! Use integer representatives for the different operations so that
  ! the if statements in the main loop below can use interger comparisons
  ! instead of string comparisons. This will help speed up that loop.
  integer, parameter :: OP_ADD  = 1
  integer, parameter :: OP_SUB  = 2
  integer, parameter :: OP_MULT = 3
  integer, parameter :: OP_DIV  = 4
  integer, parameter :: OP_AND  = 5
  integer, parameter :: OP_OR   = 6
  integer, parameter :: OP_INV  = 7
  integer, parameter :: OP_ABS  = 8

  integer, parameter :: LargeString  = 512
  integer, parameter :: MediumString = 256
  integer, parameter :: LittleString = 128
  integer, parameter :: MaxArgFields = 20

  real, parameter :: UndefVal = -999.0

  character (len=MediumString) :: InFile1
  character (len=LittleString) :: VarName1
  character (len=MediumString) :: InFile2
  character (len=LittleString) :: VarName2
  character (len=MediumString) :: OutFile
  character (len=LittleString) :: OutVarName
  character (len=LittleString) :: OpName
  integer :: Op

  ! Data arrays
  ! Dims: x, y, z, t
  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords
  type (Rhdf5Var) :: Var1, Var2, OutVar
  character (len=RHDF5_MAX_STRING) :: rh5f_facc
  integer :: rh5f_in1, rh5f_in2, rh5f_out

  integer :: i, ix, iy, iz, it, id
  integer :: Nx, Ny, Nz, Nt, Nelems
  logical :: BadDims

  ! Get the command line arguments
  call GetMyArgs(InFile1, VarName1, InFile2, VarName2, OutFile, OutVarName, OpName, Op)

  write (*,*) 'Performing operation on HDF5 REVU data:'
  write (*,*) '  Input file 1: ', trim(InFile1)
  write (*,*) '    Variable Name: ', trim(VarName1)
  write (*,*) '  Input file 2: ', trim(InFile2)
  write (*,*) '    Variable Name: ', trim(VarName2)
  write (*,*) '  Output file:  ', trim(OutFile)
  write (*,*) '    Variable Name: ', trim(OutVarName)
  write (*,*) '  Operator: ', trim(OpName)
  write (*,*) ''
  flush(6)

  ! set up file and variable names
  Var1%vname = VarName1
  Var2%vname = VarName2

  ! Check that the dimensions are consistent between the variables needed for
  ! the selected averaging function.

  call rhdf5_read_init(InFile1, Var1)
  call rhdf5_read_init(InFile2, Var2)

  Nx = Var1%dims(1)
  Ny = Var1%dims(2)
  if (Var1%ndims .eq. 3) then
    Nz = 1
    Nt = Var1%dims(3)
  else
    Nz = Var1%dims(3)
    Nt = Var1%dims(4)
  endif

  Nelems = Nx * Ny * Nz ! number of elements in spatial field

  BadDims = Var2%dims(1) .ne. Nx
  BadDims = BadDims .or. (Var2%dims(2) .ne. Ny)
  if (Var1%ndims .eq. 3) then
    BadDims = BadDims .or. (Var2%dims(3) .ne. Nt)
  else
    BadDims = BadDims .or. (Var2%dims(3) .ne. Nz)
    BadDims = BadDims .or. (Var2%dims(4) .ne. Nt)
  endif

  if (BadDims) then
    write (*,*) 'ERROR: dimensions of variables from input files do not match'
    stop
  endif

  ! Check if the operation makes sense
  if ((Op .eq. OP_ADD) .or. (Op .eq. OP_SUB)) then
    if (trim(Var1%units) .ne. trim(Var2%units)) then
      write (*,*) 'ERROR: units of variables from input files do not match'
      write (*,*) 'ERROR:   Var1: ', trim(Var1%vname), ' -> ', trim(Var1%units)
      write (*,*) 'ERROR:   Var2: ', trim(Var2%vname), ' -> ', trim(Var2%units)
      stop
    endif
  endif

  ! Set up the dimensions for reading in the input field data, one time step per read. 
  ! In other words, remove the time dimension from the input dimensions since we will 
  ! be incrementing through every time step in a loop. The time dimension is always the
  ! last dimension so what this boils down to is to decrement number of dimensions by one.
  Var1%ndims = Var1%ndims - 1
  Var2%ndims = Var2%ndims - 1
  write (*,*) 'Input variable information:'
  write (*,*) '  Number of dimensions: ', Var1%ndims
  write (*,*) '  Dimension sizes:'
  do id = 1, Var1%ndims
    write (*,*), '    ', trim(Var1%dimnames(id)), ': ', Var1%dims(id)
  enddo
  write (*,*) ''
  flush(6)

  ! Set up the dimensions for the output and allocate the output data array.
  OutVar%vname = trim(OutVarName)
  OutVar%descrip = Var1%descrip
  OutVar%units = Var1%units
  OutVar%ndims = Var1%ndims 
  OutVar%dims(1) = Nx
  OutVar%dims(2) = Ny
  OutVar%dims(3) = Nz
  OutVar%dimnames(1) = Var1%dimnames(1)
  OutVar%dimnames(2) = Var1%dimnames(2)
  OutVar%dimnames(3) = Var1%dimnames(3)

  allocate(OutVar%vdata(Nx*Ny*Nz))
 
  ! Report the dimensions
  write (*,*) 'Output variable information:'
  write (*,*) '  Name: ', trim(OutVar%vname)
  write (*,*) '  Units: ', trim(OutVar%units)
  write (*,*) '  Description: ', trim(OutVar%descrip)
  write (*,*) '  Number of dimensions: ', OutVar%ndims
  write (*,*) '  Dimension sizes:'
  do id = 1, OutVar%ndims
    write (*,*), '    ', trim(OutVar%dimnames(id)), ': ', OutVar%dims(id)
  enddo
  write (*,*) ''
  flush(6)

  ! Report time steps
  write (*,*) 'Number of time steps: ', Nt
  write (*,*) ''
  flush(6)

  ! Read in the input coordinates
  Xcoords%vname = 'x_coords'
  call rhdf5_read_init(InFile1, Xcoords)
  call rhdf5_read(InFile1, Xcoords)

  Ycoords%vname = 'y_coords'
  call rhdf5_read_init(InFile1, Ycoords)
  call rhdf5_read(InFile1, Ycoords)

  Zcoords%vname = 'z_coords'
  call rhdf5_read_init(InFile1, Zcoords)
  call rhdf5_read(InFile1, Zcoords)

  Tcoords%vname = 't_coords'
  call rhdf5_read_init(InFile1, Tcoords)
  call rhdf5_read(InFile1, Tcoords)

  ! Perform the operation
  rh5f_facc = 'W'
  call rhdf5_open_file(OutFile, rh5f_facc, 1, rh5f_out)

  rh5f_facc = 'R'
  call rhdf5_open_file(InFile1, rh5f_facc, 0, rh5f_in1)
  call rhdf5_open_file(InFile2, rh5f_facc, 0, rh5f_in2)

  do it = 1, Nt
    call rhdf5_read_variable(rh5f_in1, Var1%vname, Var1%ndims, it, Var1%dims, rdata=Var1%vdata)
    call rhdf5_read_variable(rh5f_in2, Var2%vname, Var2%ndims, it, Var2%dims, rdata=Var2%vdata)

    ! do the op here
    do i = 1, Nelems
      if (Op .eq. OP_ADD) then
        OutVar%vdata(i) = Var1%vdata(i) + Var2%vdata(i)
      else if (Op .eq. OP_SUB) then
        OutVar%vdata(i) = Var1%vdata(i) - Var2%vdata(i)
      else if (Op .eq. OP_MULT) then
        OutVar%vdata(i) = Var1%vdata(i) * Var2%vdata(i)
      else if (Op .eq. OP_DIV) then
        if (Var2%vdata(i) .eq. 0.0) then
          OutVar%vdata(i) = UndefVal
        else
          OutVar%vdata(i) = Var1%vdata(i) / Var2%vdata(i)
        endif
      else if (Op .eq. OP_AND) then
        if ((anint(Var1%vdata(i)) .eq. 1.0) .and. (anint(Var2%vdata(i)) .eq. 1.0)) then
          OutVar%vdata(i) = 1.0
        else
          OutVar%vdata(i) = 0.0
        endif
      else if (Op .eq. OP_OR) then
        if ((anint(Var1%vdata(i)) .eq. 0.0) .and. (anint(Var2%vdata(i)) .eq. 0.0)) then
          OutVar%vdata(i) = 0.0
        else
          OutVar%vdata(i) = 1.0
        endif
      else if (Op .eq. OP_INV) then
        if ((anint(Var1%vdata(i)) .eq. 0.0)) then
          OutVar%vdata(i) = 1.0
        else
          OutVar%vdata(i) = 0.0
        endif
      else if (Op .eq. OP_ABS) then
        OutVar%vdata(i) = abs(Var1%vdata(i))
      endif
    enddo

    deallocate(Var1%vdata)
    deallocate(Var2%vdata)

    ! write out the result
    call rhdf5_write_variable(rh5f_out, OutVar%vname, OutVar%ndims, it, OutVar%dims, &
       OutVar%units, OutVar%descrip, OutVar%dimnames, rdata=OutVar%vdata)
 
    ! print a message for the user on longer jobs so that it can be
    ! seen that progress is being made
    if (modulo(it,100) .eq. 0) then
      write (*,*) 'Working: Number of time steps processed so far: ', it
    endif
  enddo

  call rhdf5_close_file(rh5f_in1)
  call rhdf5_close_file(rh5f_in2)
  call rhdf5_close_file(rh5f_out)
  deallocate(OutVar%vdata)

  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''

  ! Finish off output file
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
  
  stop

contains
!**********************************************************************
! SUBROUTINES
!**********************************************************************

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the command line arguments
!
subroutine GetMyArgs(InFile1, VarName1, InFile2, VarName2, OutFile, OutVarName, OpName, Op)
  implicit none

  character (len=*) :: InFile1, VarName1, InFile2, VarName2, OutFile, OutVarName, OpName
  integer :: Op

  integer :: iargc
  character (len=128) :: arg
  integer :: Nitems

  logical :: BadArgs

  if (iargc() .ne. 7) then
    write (*,*) 'ERROR: must supply exactly 7 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: hdata_op <in_file1> <var_name1> <in_file2> <var_name2> <out_file> <out_var_name> <operator>'
    write (*,*) '        <in_file1>: input file, HDF5 format'
    write (*,*) '        <var_name1>: name of variable inside input file 1'
    write (*,*) '        <in_file2>: input file, HDF5 format'
    write (*,*) '        <var_name2>: name of variable inside input file 2'
    write (*,*) '        <out_file>: output file, HDF5 format'
    write (*,*) '        <out_var_name>: name of variable in output file'
    write (*,*) '        <operator>: operation to perform on data from input files'
    write (*,*) '            order of operation: <in_file1> <operator> <in_file2>'
    write (*,*) '            add  -> add values in <in_file1> to those in <in_file2>'
    write (*,*) '            sub  -> subtract values in <in_file2> from those in <in_file1>'
    write (*,*) '            mult -> multiply values in <in_file1> by those in <in_file2>'
    write (*,*) '            div  -> divide values in <in_file1> by those in <in_file2>'
    write (*,*) '            and  -> logical and values in <in_file1> with those in <in_file2>'
    write (*,*) '            or   -> logical or values in <in_file1> with those in <in_file2>'
    write (*,*) '            inv  -> logical inversion of values in <in_file1>'
    write (*,*) '            abs  -> absolute value of values in <in_file1>'
    stop
  end if

  call getarg(1, InFile1)
  call getarg(2, VarName1)
  call getarg(3, InFile2)
  call getarg(4, VarName2)
  call getarg(5, OutFile)
  call getarg(6, OutVarName)
  call getarg(7, OpName)

  BadArgs = .false.

  if (OpName .eq. 'add') then
    Op = OP_ADD
  else if (OpName .eq. 'sub') then
    Op = OP_SUB
  else if (OpName .eq. 'mult') then
    Op = OP_MULT
  else if (OpName .eq. 'div') then
    Op = OP_DIV
  else if (OpName .eq. 'and') then
    Op = OP_AND
  else if (OpName .eq. 'or') then
    Op = OP_OR
  else if (OpName .eq. 'inv') then
    Op = OP_INV
  else if (OpName .eq. 'abs') then
    Op = OP_ABS
  else
    write (*,*) 'ERROR: <operator> must be one of:'
    write (*,*) '          add'
    write (*,*) '          sub'
    write (*,*) '          mult'
    write (*,*) '          div'
    write (*,*) '          and'
    write (*,*) '          or'
    write (*,*) '          inv'
    write (*,*) '          abs'
    write (*,*) ''
    BadArgs = .true.
  end if

  if (BadArgs) then
    stop
  end if

  return
end subroutine GetMyArgs

end program hdata_op
