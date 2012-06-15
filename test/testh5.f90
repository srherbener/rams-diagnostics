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
  integer, parameter :: ssize = 15

  integer, dimension(Nx,Ny,Nz) :: A, B, C
  integer, dimension(Nx,Ny) :: Array2D
  real, dimension(BigNx,BigNy,BigNz) :: BigVar
  integer, dimension(Nx,Ny,Nz) :: Iarray
  real, dimension(Nx,Ny,Nz) :: Rarray
  character (len=ssize), dimension(Nx,Ny,Nz) :: Sarray
  character (len=ssize) :: Carray
 
  integer :: i, ix, iy, iz
  integer, dimension(3) :: dims
  integer :: ndims
  integer :: rh5_file
  character (len=RHDF5_MAX_STRING) :: rh5_file_name
  character (len=RHDF5_MAX_STRING) :: rh5_file_acc
  character (len=RHDF5_MAX_STRING), dimension(RHDF5_MAX_DIMS) :: dimnames
  character (len=RHDF5_MAX_STRING) :: units, descrip
  integer :: TestNum
  real :: tstart, tend, ttot
  character (len=128) :: Vname

  call GetArgs(TestNum)

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

  elseif (TestNum .eq. 3) then
    ! try reading
    rh5_file_name = 'test_read.h5'
    rh5_file_acc = 'R'
    call rhdf5_open_file(rh5_file_name, rh5_file_acc, 0, rh5_file)

    print*, 'Reading IntData:'
    call rhdf5_read_variable(rh5_file,'IntData', ndims, dims, units, descrip, dimnames, idata=Iarray)
    call DumpVarInfo(ndims, dims, units, descrip, dimnames, RHDF5_MAX_STRING)

    print*, 'Reading RealData:'
    call rhdf5_read_variable(rh5_file,'RealData', ndims, dims, units, descrip, dimnames, rdata=Rarray)
    call DumpVarInfo(ndims, dims, units, descrip, dimnames, RHDF5_MAX_STRING)

    print*, 'Reading StringData:'
    call rhdf5_read_variable(rh5_file,'StringData', ndims, dims, units, descrip, dimnames, sdata=Sarray, ssize=ssize)
    call DumpVarInfo(ndims, dims, units, descrip, dimnames, RHDF5_MAX_STRING)
    print*, '  ssize: ', ssize

    print*, 'Reading CharData:'
    call rhdf5_read_variable(rh5_file,'CharData', ndims, dims, units, descrip, dimnames, cdata=Carray)
    call DumpVarInfo(ndims, dims, units, descrip, dimnames, RHDF5_MAX_STRING)

    print*, '3D Data Arrays:'
    print*, ''
    write(*,'(a5,a5,a5,a10,a10,a20)') 'ix', 'iy', 'iz', 'Iarray', 'Rarray', 'Sarray'
    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          write(*,'(i5,i5,i5,i10,f10.2,a20)') ix, iy, iz, &
            Iarray(ix,iy,iz), Rarray(ix,iy,iz), trim(Sarray(ix,iy,iz))
        enddo
      enddo
    enddo
    print*, ''

    print*, 'Char Data Arrays:'
    write(*,'(a2,a20)') '  ', Carray
    print*, ''

    call rhdf5_close_file(rh5_file)

  elseif (TestNum .eq. 4) then
    i = 1
    do iz = 1,Nz
      do iy = 1,Ny
        do ix = 1,Nx
          Iarray(ix,iy,iz) = i
          Rarray(ix,iy,iz) = float(i)
          write(Sarray(ix,iy,iz), '(a8,i0.3)') 'Sarray: ', i
          i = i + 1
        enddo
      enddo
    enddo
    write(Carray, '(a8,i0.3)') 'Carray: ', i

    ! write into an HDF5 file
    rh5_file_name = 'test_wr.h5'
    rh5_file_acc = 'W'
    call rhdf5_open_file(rh5_file_name, rh5_file_acc, 1, rh5_file)

    dims(1) = Nx
    dims(2) = Ny
    dims(3) = Nz

    dimnames(1) = 'x'
    dimnames(2) = 'y'
    dimnames(3) = 'z'
  
    ! write out one of each supported type: integer, real, string, character
    call rhdf5_write_variable(rh5_file, 'IntData', 3, 0, dims, 'test_ui', 'test_di', dimnames, idata=Iarray)
    call rhdf5_write_variable(rh5_file, 'RealData', 3, 0, dims, 'test_ur', 'test_dr', dimnames, rdata=Rarray)
    call rhdf5_write_variable(rh5_file, 'StringData', 3, 0, dims, 'test_us', 'test_ds', dimnames, sdata=Sarray, ssize=ssize)

    dims(1) = ssize
    dims(2) = 0
    dims(3) = 0
    call rhdf5_write_variable(rh5_file, 'CharData', 1, 0, dims, 'test_uc', 'test_dc', dimnames, cdata=Carray)

    call rhdf5_close_file(rh5_file)

    ! now try reading
    ! blank out the arrays
    do iz = 1,Nz
      do iy = 1,Ny
        do ix = 1,Nx
          Iarray(ix,iy,iz) = -1
          Rarray(ix,iy,iz) = -1.0
          write(Sarray(ix,iy,iz), '(a5)') 'EMPTY'
          i = i + 1
        enddo
      enddo
    enddo
    write(Carray, '(a5)') 'EMPTY'

    ndims = -1
    dims(1) = -1
    dims(2) = -1
    dims(3) = -1
    dimnames(1) = 'Empty'
    dimnames(2) = 'Empty'
    dimnames(3) = 'Empty'
    units = 'Empty'
    descrip = 'Empty'

    rh5_file_name = 'test_wr.h5'
    rh5_file_acc = 'R'
    call rhdf5_open_file(rh5_file_name, rh5_file_acc, 0, rh5_file)

    print*, 'Reading IntData:'
    call rhdf5_read_variable(rh5_file,'IntData', ndims, dims, units, descrip, dimnames, idata=Iarray)
    call DumpVarInfo(ndims, dims, units, descrip, dimnames, RHDF5_MAX_STRING)

    print*, 'Reading RealData:'
    call rhdf5_read_variable(rh5_file,'RealData', ndims, dims, units, descrip, dimnames, rdata=Rarray)
    call DumpVarInfo(ndims, dims, units, descrip, dimnames, RHDF5_MAX_STRING)

    print*, 'Reading StringData:'
    call rhdf5_read_variable(rh5_file,'StringData', ndims, dims, units, descrip, dimnames, sdata=Sarray, ssize=ssize)
    call DumpVarInfo(ndims, dims, units, descrip, dimnames, RHDF5_MAX_STRING)
    print*, '  ssize: ', ssize

    print*, 'Reading CharData:'
    call rhdf5_read_variable(rh5_file,'CharData', ndims, dims, units, descrip, dimnames, cdata=Carray)
    call DumpVarInfo(ndims, dims, units, descrip, dimnames, RHDF5_MAX_STRING)

    print*, '3D Data Arrays:'
    print*, ''
    write(*,'(a5,a5,a5,a10,a10,a20)') 'ix', 'iy', 'iz', 'Iarray', 'Rarray', 'Sarray'
    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          write(*,'(i5,i5,i5,i10,f10.2,a20)') ix, iy, iz, &
            Iarray(ix,iy,iz), Rarray(ix,iy,iz), trim(Sarray(ix,iy,iz))
        enddo
      enddo
    enddo
    print*, ''

    print*, 'Char Data Arrays:'
    write(*,'(a2,a20)') '  ', Carray
    print*, ''

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
     print*, '       <test_num>: 1 --> run short write test to check file format'
     print*, '                   2 --> run big file write test to check performance'
     print*, '                   3 --> read test'
     print*, '                   4 --> write then read test'
     print*, ''
   stop
  endif

  call getarg(1, arg)
  read(arg,'(i)') TestNum
  if ((TestNum .lt. 1) .and. (TestNum .gt. 4)) then
    print*, 'ERROR: <test_num> must be integer between 1 and 4'
    stop 'GetArgs'
  endif

  return
end subroutine GetArgs

!********************************************
subroutine DumpVarInfo(ndims, dims, units, descrip, dimnames, ssize)
  implicit none

  integer :: ndims
  integer :: ssize
  integer, dimension(ndims) :: dims
  character (len=ssize) :: units
  character (len=ssize) :: descrip
  character (len=ssize), dimension(ndims) :: dimnames

  integer :: i

  print*, 'Variable information:'
  print*, '  ndims: ', ndims
  print*, '  dims: '
  do i = 1,ndims
    print*, '    i, dims(i): ', i, dims(i)
  enddo
  print*, '  units: ', trim(units)
  print*, '  descrip: ', trim(descrip)
  print*, '  dimnames: '
  do i = 1,ndims
    print*, '    i, dimnames(i): ', i, trim(dimnames(i))
  enddo
  print*, ''

  return
end subroutine DumpVarInfo

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
