! test hdf5 storage order

program testh5
  use rhdf5_utils
  implicit none

  integer, parameter :: Nx = 3
  integer, parameter :: Ny = 2
  integer, parameter :: Nz = 2
  integer, parameter :: BigNx = 3000
  !integer, parameter :: BigNx = 6
  integer, parameter :: BigNy = 200
  integer, parameter :: BigNz = 65
  !integer, parameter :: BigNitems = 300
  !integer, parameter :: BigNitems = 150000
  integer, parameter :: BigNitems = 10
  !integer, parameter :: BigNitems = 5000

  integer, dimension(Nx,Ny,Nz) :: A, B, C
  integer, dimension(Nx,Ny) :: Array2D
  real, dimension(BigNx,BigNy,BigNz) :: BigVar
 
  integer :: i, ix, iy, iz
  integer, dimension(3) :: dims
  integer :: rh5_file
  character (len=HDF5_MAX_STRING) :: rh5_file_name
  character (len=HDF5_MAX_STRING) :: rh5_file_acc
  character (len=HDF5_MAX_STRING), dimension(HDF5_MAX_DIMS) :: dimnames
  integer :: TestNum
  real :: tstart, tend, ttot
  character (len=128) :: Vname

  call GetArgs(TestNum)

print*,'DEBUG: sizeof(BigNx): ', sizeof(BigNx)
print*,'DEBUG: sizeof(tstart): ', sizeof(tstart)
print*,'DEBUG: sizeof(Vname): ', sizeof(Vname)
print*,'DEBUG: sizeof(BigVar): ', sizeof(BigVar)

  if (TestNum .eq. 1) then
    dims(1) = Nx
    dims(2) = Ny
    dims(3) = Nz
  
    i = 1
    do iz = 1,Nz
      do iy = 1,Ny
        do ix = 1,Nx
          A(ix,iy,iz) = i
          B(ix,iy,iz) = i + 12
          C(ix,iy,iz) = i + 24
          i = i + 1
        enddo
      enddo
    enddo
  
    i = 1
    do iy = 1,Ny
      do ix = 1,Nx
        Array2D(ix,iy) = i
        i = i + 1
      enddo
    enddo
  
    print*,'A before hdf5 write (linear storage): '
    call DumpLinearStorage(A,Nx*Ny*Nz)
    print*,'B before hdf5 write (linear storage): '
    call DumpLinearStorage(B,Nx*Ny*Nz)
    print*,'C before hdf5 write (linear storage): '
    call DumpLinearStorage(C,Nx*Ny*Nz)
    print*,'Array2D before hdf5 write (linear storage): '
    call DumpLinearStorage(Array2D,Nx*Ny)
  
    ! write into an HDF5 file
    rh5_file_name = 'test_write.h5'
    rh5_file_acc = 'W'
    call rhdf5_open_file(rh5_file_name, rh5_file_acc, 1, rh5_file)
  
    ! Write out the A, 5th argument is the time step. Setting to 0 tells rhdf5_write_variable
    ! that these variables do not have a time dimension, just (x,y,z).
    ! subroutine rhdf5_write_variable(rh5_file, vloc, vname, ndims, itstep, &
    !    dims, units, descrip, rdata, idata, sdata, ssize)
    dimnames(1) = 'x'
    dimnames(2) = 'y'
    dimnames(3) = 'z'
    call rhdf5_write_variable(rh5_file, 'test_4D', 3, 1, dims, 'test', 'test', dimnames, idata=A)
    call rhdf5_write_variable(rh5_file, 'test_4D', 3, 2, dims, 'test', 'test', dimnames, idata=B)
    call rhdf5_write_variable(rh5_file, 'test_4D', 3, 3, dims, 'test', 'test', dimnames, idata=C)
  
    call rhdf5_write_variable(rh5_file, 'test_3D', 3, 0, dims, 'test', 'test', dimnames, idata=A)
  
    call rhdf5_write_variable(rh5_file, 'test_2D', 2, 0, dims, 'test', 'test', dimnames, idata=Array2D)
  
    call rhdf5_close_file(rh5_file)
  elseif (TestNum .eq. 2) then
    ! load up the array and write it out a bunch of times to get a big file created

    i = 1
    do iz = 1,BigNz
      do iy = 1,BigNy
        do ix = 1,BigNx
          BigVar(ix,iy,iz) = i
          i = i + 1
        enddo
      enddo
    enddo
    dims(1) = BigNx
    dims(2) = BigNy
    dims(3) = BigNz
    dimnames(1) = 'x'
    dimnames(2) = 'y'
    dimnames(3) = 'z'

    ! write into an HDF5 file
    rh5_file_name = 'test_big.h5'
    rh5_file_acc = 'W'
    call rhdf5_open_file(rh5_file_name, rh5_file_acc, 1, rh5_file)
  
    print*, 'Writing file:'
    call cpu_time(tstart)
    do i = 1, BigNitems
      write(Vname, '(a6,i0)') 'BigVar', i
      call rhdf5_write_variable(rh5_file, Vname, 3, 0, dims, 'test', 'test', dimnames, rdata=BigVar)
    enddo
    call cpu_time(tend)
    ttot = tend - tstart

    print*, ''
    print*, 'Total elapsed time (s): ', ttot
    print*, 'Average write time (s): ', ttot/real(BigNitems);

    call rhdf5_close_file(rh5_file)

  endif

end program testh5

!**************************************************
subroutine GetArgs(TestNum)
  implicit none

  integer :: TestNum

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 1) then
     print*, 'ERROR: must supply exactly 1 argument'
     print*, ''
     print*, 'USAGE: testh5 <test_num>'
     print*, '       <test_num>: 1 --> run short test to check file format'
     print*, '                   2 --> run big file test to check performance'
     print*, ''
   stop
  endif

  call getarg(1, arg)
  read(arg,'(i)') TestNum
  if ((TestNum .lt. 1) .and. (TestNum .gt. 2)) then
    print*, 'ERROR: <test_num> must be integer between 1 and 2'
    stop 'GetArgs'
  endif

  return
end subroutine GetArgs

!********************************************
subroutine DumpLinearStorage(A,N)
  implicit none

  integer :: N
  integer, dimension(N) :: A

  integer i

  do i = 1,N
    write(*,'(a3,i5,a4,i5)') 'i: ', i, ' -> ', A(i)
  enddo

  return
end subroutine DumpLinearStorage
