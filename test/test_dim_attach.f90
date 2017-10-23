! test hdf5 storage order

program test_dim_attach
  use rhdf5_utils
  implicit none

  character (len=RHDF5_MAX_STRING) :: H5File, H5Var
  character (len=RHDF5_MAX_STRING) :: rh5_file_acc
  character (len=RHDF5_MAX_STRING), dimension(RHDF5_MAX_DIMS) :: dimnames
  character (len=RHDF5_MAX_STRING) :: units, descrip
  integer :: rh5_file

  call GetArgs(H5File, H5Var)
  print*, 'Updating HDF5 file: ', trim(H5File)
  print*, '  Attaching dimensions to variable: ', trim(H5Var)

  rh5_file_acc = 'RW'
  call rhdf5_open_file(H5File, rh5_file_acc, 0, rh5_file)

  call rhdf5_set_var_to_dim(rh5_file, 'x_coords', 'x')
  call rhdf5_set_var_to_dim(rh5_file, 'y_coords', 'y')
  call rhdf5_set_var_to_dim(rh5_file, 'z_coords', 'z')
  call rhdf5_set_var_to_dim(rh5_file, 't_coords', 't')

  call rhdf5_attach_dims_to_var(rh5_file, trim(H5Var))

  call rhdf5_close_file(rh5_file)

  stop

contains
!**************************************************
subroutine GetArgs(H5File, H5Var)
  implicit none

  character (len=RHDF5_MAX_STRING) :: H5File, H5Var

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 2) then
     print*, 'ERROR: must supply exactly 1 argument'
     print*, ''
     print*, 'USAGE: test_dim_attach <hdf5_file> <hdf5_variable>'
     print*, ''
   stop
  endif

  call getarg(1, H5File)
  call getarg(2, H5Var)

  return
end subroutine GetArgs

end program test_dim_attach
