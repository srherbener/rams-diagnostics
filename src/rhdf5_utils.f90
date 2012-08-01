!
! Copyright (C) 1991-2004  ; All Rights Reserved ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
!############################# Change Log ##################################
! 2.2.0
!###########################################################################

Module rhdf5_utils

! max limits for arrays, strings, etc, keep these in sync with like named
! defines in rhdf5_f2c.c
integer, parameter :: RHDF5_MAX_STRING = 128
integer, parameter :: RHDF5_MAX_DIMS   =  10

! integer coding for HDF5 types, keep these in sync with like named defines
! in rhdf5_f2c.c
integer, parameter :: RHDF5_TYPE_STRING  = 0
integer, parameter :: RHDF5_TYPE_INTEGER = 1
integer, parameter :: RHDF5_TYPE_FLOAT   = 2
integer, parameter :: RHDF5_TYPE_CHAR    = 3

! structure to hold variable information
type Rhdf5Var
  character (len=RHDF5_MAX_STRING) :: vname
  integer :: ndims
  integer, dimension(RHDF5_MAX_DIMS) :: dims
  character (len=RHDF5_MAX_STRING) :: units
  character (len=RHDF5_MAX_STRING) :: descrip
  character (len=RHDF5_MAX_STRING), dimension(RHDF5_MAX_DIMS) :: dimnames
  real, dimension(:), allocatable :: vdata
end type Rhdf5Var

Contains

!********************************************************************************
! Routines for HDF5 REVU IO
!********************************************************************************

!********************************************************************************
! High level routines
!********************************************************************************

!**************************************************************
! rhdf5_read_init()
!
! This routine will locate the given variable in the given HDF5 file
! and read in information about the variable to prepare for a
! subsequent read.
!
! Splitting up the read into two pieces allows the caller to
! quickly get size information about variables without having
! to read in all the data. One purpose for this is get the
! dimensions so that the data buffer can be allocated.
! 
subroutine rhdf5_read_init(fname, rvar)
  implicit none

  character (len=RHDF5_MAX_STRING) :: fname
  type (Rhdf5Var) :: rvar

  integer :: fileid
  character (len=RHDF5_MAX_STRING):: facc

  facc = 'R'
  call rhdf5_open_file(fname, facc, 0, fileid)

  call rhdf5_read_variable_init(fileid, rvar%vname, rvar%ndims, rvar%dims, rvar%units, &
    rvar%descrip, rvar%dimnames)

  call rhdf5_close_file(fileid)
  return
end subroutine rhdf5_read_init

!**************************************************************
! rhdf5_read()
!
! This routine will read the given variable out of the given HDF5 file.
!
! The allocation for the memory buffer (vdata) is done here, the caller
! is responsible for deallocating this memory when done with the
! variable.
! 
subroutine rhdf5_read(fname, rvar)
  implicit none

  character (len=RHDF5_MAX_STRING) :: fname
  type (Rhdf5Var) :: rvar

  integer :: fileid
  character (len=RHDF5_MAX_STRING):: facc
  integer :: i, nelem

  ! Calculate the number of elements needed for the variable data
  ! based on the dimensions of that variable.
  nelem = 1
  do i = 1, rvar%ndims
    nelem = nelem * rvar%dims(i)
  enddo

  allocate(rvar%vdata(nelem))

  facc = 'R'
  call rhdf5_open_file(fname, facc, 0, fileid)

  call rhdf5_read_variable(fileid, rvar%vname, rvar%ndims, rvar%dims, rdata=rvar%vdata)

  call rhdf5_close_file(fileid)
  return
end subroutine rhdf5_read

!**************************************************************
! rhdf5_write()
!
! This routine will write the given variable out of the given HDF5 file.
! The caller must set all of the elements in the Rhdf5Var structure before
! calling this routine.
subroutine rhdf5_write(fname, rvar, append)
  implicit none

  character (len=RHDF5_MAX_STRING) :: fname
  type (Rhdf5Var) :: rvar
  integer :: append

  integer :: fileid
  character (len=RHDF5_MAX_STRING):: facc

  if (append .eq. 1) then
    facc = 'RW'
  else
    facc = 'W'
  endif
  call rhdf5_open_file(fname, facc, 1, fileid)

  ! Fourth arg is 'itstep' (time step) which is used for writing into an extendable
  ! dimension. We already have the complete data (no need to extend) so set this to zero.
  call rhdf5_write_variable(fileid, rvar%vname, rvar%ndims, 0, rvar%dims, rvar%units, &
    rvar%descrip, rvar%dimnames, rdata=rvar%vdata)

  call rhdf5_close_file(fileid)
  return
end subroutine rhdf5_write

!**************************************************************
! rhdf5_set_dimension()
!
! This routine will set up a given variable as a "dimension"
! for GRADS. In the HDF5 lingo, this routine will set the
! dimension scale on the given variable.
!
! The value in dname will be set in a new attribute on
! the given variable called "axis" of which GRADS uses
! to figure out which dimension this variable represents.
! (Analagous to the [XYZT]DEF statements in the GRADS
! descriptor file.) Because of this, dname must be one
! of 'x', 'y', 'z', and 't'.
!
subroutine rhdf5_set_dimension(fname, rvar, dname)
  implicit none

  character (len=*) :: fname
  character (len=*) :: dname
  type (Rhdf5Var) :: rvar

  integer :: fileid
  character (len=RHDF5_MAX_STRING):: facc

  facc = 'RW'
  call rhdf5_open_file(fname, facc, 1, fileid)

  call rhdf5_set_var_to_dim(fileid, rvar%vname, dname)

  call rhdf5_close_file(fileid)
  return
end subroutine rhdf5_set_dimension

!**************************************************************
! rhdf5_attach_dimensions()
!
! This routine will attach the dimension specs to the given
! variable. The code assumes that for each dimension in rvar,
! there exists a dataset in the file that contains the coordinate
! values for that dimension. The name of the coordinate variable
! is given by the suffix "_coords" appended to the dimension name. 
!
! Example:
!  If the variable is U and it has dimensions (x,y,z,t), ie.
!  dimnames contains the strings: "x", "y", "z" and "t", then
!  the followind datasets will be expected in the file for
!  attachement:
!    x_coords
!    y_coords
!    z_coords
!    t_coords
!
subroutine rhdf5_attach_dimensions(fname, rvar)
  implicit none

  character (len=*) :: fname
  type (Rhdf5Var) :: rvar

  integer :: fileid
  character (len=RHDF5_MAX_STRING):: facc

  facc = 'RW'
  call rhdf5_open_file(fname, facc, 1, fileid)

  call rhdf5_attach_dims_to_var(fileid, rvar%vname)

  call rhdf5_close_file(fileid)
  return
end subroutine rhdf5_attach_dimensions

!********************************************************************
! rhdf5_write_attribute()
!
! This routine will create and populate an attribute for an hdf5 file group.
!
subroutine rhdf5_write_attribute(fname, vname, aname, cval, ival, rval)
  implicit none

  character (len=*) :: fname
  character (len=*) :: vname
  character (len=*) :: aname
  character (len=*), optional :: cval
  integer, optional :: ival
  real, optional :: rval

  integer :: fileid
  integer :: dsetid
  character (len=RHDF5_MAX_STRING):: facc
  integer :: hdferr

  facc = 'RW'
  call rhdf5_open_file(fname, facc, 1, fileid)

  call rh5d_open(fileid, trim(vname)//char(0), dsetid, hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_write_attribute: cannot open dataset for variable: ',trim(vname)
    stop 'hdf5_write_attribute: variable open error'
  endif

  ! Figure out the type
  !    type codes for rh5a_write_anyscalar are:
  !      1 - string
  !      2 - integer
  !      3 - real
  if (present(cval)) then
    call rh5a_write_anyscalar(dsetid, trim(aname)//char(0), trim(cval)//char(0), RHDF5_TYPE_STRING, hdferr)
  elseif (present(ival)) then
    call rh5a_write_anyscalar(dsetid, trim(aname)//char(0), ival, RHDF5_TYPE_INTEGER, hdferr)
  elseif (present(rval)) then
    call rh5a_write_anyscalar(dsetid, trim(aname)//char(0), rval, RHDF5_TYPE_FLOAT, hdferr)
  else
    print*,'ERROR: hdf5_write_attribute: must use one of the arguments cval, ival or rval'
    stop 'hdf5_write_attribute: bad argument assignment'
  endif

  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_write_attribute: cannot write attribute: ',trim(aname)
    stop 'hdf5_write_attribute: bad attribute write'
  endif

  call rh5d_close(dsetid, hdferr)
  call rhdf5_close_file(fileid)

  return
end subroutine rhdf5_write_attribute

!********************************************************************************
! Low level routines
!********************************************************************************

!*******************************************************************
! rhdf5_open_file()
!
! This routine will open an HDF5 file specified by the fname, facc, fdelete
! arguments, and pass back the file id number in the fid argument.
!
!  facc: 'R'  --> read
!        'RW' --> read/write
!        'W'  --> write
!
!  fdelete: 0 --> do not delete if file exists
!           1 --> delete if opening in write mode
!
! If read mode
!    file exists --> open (rh5f_open)
!    file does not exist --> error
!
! If read/write mode
!    file exists --> open (rh5f_open)
!    file does not exist --> error
!
! If write mode
!    file exists:
!      fdelete flag is 1 --> truncate file (rh5f_create)
!      fdelete fla0 is 0 --> error
!    file does not exists --> create (rh5f_create)
!
! For the 'f_flgs' argument to rh5f_open and rh5f_create routines
! (f_flgs is 2nd argument to both routines)
!   If calling rh5f_open:
!      read only  --> set f_flgs to 1
!      read/write --> set f_flgs to 2
!   If calling rh5f_create
!      truncate the file if file exists --> set f_flgs to 1
!      fail if file exists              --> set f_flgs to 2
!
subroutine rhdf5_open_file(fname,facc,fdelete,fid)
  implicit none

  character (len=RHDF5_MAX_STRING) :: fname
  character (len=RHDF5_MAX_STRING) :: facc
  integer :: fdelete
  integer :: fid

  character (len=RHDF5_MAX_STRING) :: command
  integer :: hdferr
  logical :: exists

  ! Check for existence of HDF5 file.
  inquire(file=trim(fname),exist=exists)

  ! Create a new file or open an existing file.
  if (trim(facc) .eq. 'R') then
    ! read only mode
    if (.not.exists) then
      print*, 'ERROR: rhdf5_open_file: HDF5 file does not exist: ', trim(fname)
      stop 'rhdf5_open_file: File does not exist'
    else
      call rh5f_open(trim(fname)//char(0), 1, fid, hdferr)
      if (hdferr < 0) then
        print*,'ERROR: rhdf5_open_file: Cannot open HDF5 file in read mode: ', trim(fname)
        stop 'rhdf5_open_file: Cannot open file'
      endif
    endif
  elseif (trim(facc) .eq. 'RW') then
    ! read/write mode
    if (.not.exists) then
      print*, 'ERROR: rhdf5_open_file: HDF5 file does not exist: ', trim(fname)
      stop 'rhdf5_open_file: File does not exist'
    else
      call rh5f_open(trim(fname)//char(0), 2, fid, hdferr)
      if (hdferr < 0) then
        print*, 'ERROR: rhdf5_open_file: Cannot open HDF5 file in read/write mode: ', trim(fname)
        stop 'rhdf5_open_file: Cannot open file'
      endif
    endif
  elseif (trim(facc) .eq. 'W') then
    if (.not.exists) then
      call rh5f_create(trim(fname)//char(0), 2, fid, hdferr)
    else
      if(fdelete .eq. 0) then
         print*, 'ERROR: rhdf5_open_file: File exists when delete flag is set to zero: ', trim(fname)
         stop 'rhdf5_open_file: File exists'
      else
         command = 'rm -f '//trim(fname)//char(0)
         call system(trim(command))
         call rh5f_create(trim(fname)//char(0), 1, fid, hdferr)
      endif
    endif
    if(hdferr < 0) then
      print*, 'ERROR: rhdf5_open_file: Cannot create HDF5 file in write mode:',hdferr
      stop 'rhdf5_open_file: Cannot create file'
    endif
  endif

  return
end subroutine rhdf5_open_file

!*****************************************************************
! rhdf5_close_file()
!
! This routine will close the hdf5 file.
subroutine rhdf5_close_file(fid)
  implicit none

  integer :: fid

  integer :: hdferr

  call rh5f_close(fid, hdferr)

  return
end subroutine rhdf5_close_file

!********************************************************************
! rhdf5_write_variable()
!
! This routine will create and populate the dataset and attributes
! for an hdf5 file variable.
!
! We are going through the C interface since the HDF5 Fortran interface does not work with the
! pgf90 compiler. This interface will store the array data in row major fashion (first dimension
! changes the slowest when deteriming the linear storage order for the array data) whereas Fortran
! stores array data in column major fashion (first dimension changes the fastest).
!
! The data coming into this routine is ordered (x,y,z) and is in column major fashion.
!
! This routine will handle two styles of variable storage.
!   1. Store the variable once as is. This is intended for grid information such as lat/lon/levels values.
!      It is not necessary to store the same thing over for each time step, just once when setting up
!      the grid for the first time. Setting the itstep argument to zero denotes this type of variable.
!      This variable will be stored with fixed dimensions in the HDF5 file.
!
!   2. Store the variable each time step. This is intended for RAMS variable (2D or 3D field) storage
!      where for each time step a version of the field is stored. The storage will be optimized toward
!      reading out the variable in a time step by time step fashion. This routine will allow for an
!      unlimited number of time steps as it will create the HDF5 dataset with the time dimension
!      having an unlimited max dimension. Along with this HDF5 chunking will be used so that we don't have
!      to preallocate contiguous file space for all time steps before we start writing in data. Making
!      the chunk size match the size of the data coming into this routine (x,y,z) will allow for writing
!      each time step's data into one chunk (which is a contiguous file space). Because of this scheme, we
!      want time to be the slowest changing dimension which makes it so that for each time step only one
!      chunk needs to be read out of the file to get the (x,y,z) data. Note that if time was not the slowest
!      changing dimension, then multiple chunks would have to read out of the file in order to fill in one
!      time step's (x,y,z) data.
!
! Since we want to make t the slowest changing dimension (first dimension due to the
! C interface with the row major array storage) just tagging on the (x,y,z) data that
! came into this routine would make for an awkward storage in the HDF5 file, where you
! are taking memory containing (x,y,z) in column major order and writing that out to
! the HDF5 file as if it were row major order.
!
! MATLAB and IDL use a clever trick to deal with this issue. Leave the linear storage alone,
! and simply tell the C HDF5 write routine that the dimensions for (x,y,z) go in the reverse order.
! This automatically changes the perspective on the linear storage to row major and the data in the
! HDF5 file will be consistent. That is, the value at a given x,y,z location will remain at that
! location in the file making it possible for a C routine, using row major order, to read in the
! data and have it match what was written out by this routine.
!
! For example if the data coming in has (x,y,z) with dimension sizes (4,3,2) then send the data as
! is to the C hdf5 write routine, but tell that routine that the order of the data is now (z,y,x) 
! with dimension sizes (2,3,4).
!
! If attaching time as a variable always place it first in the variable order. It is the slowest
! changing variable so placing it first is consistent with row major ordering.
!
! Note that sdata is an array of strings where each string gets cast into a C style string
! terminated with a null byte. cdata on the other hand gets cast into a character array with
! a fixed length (RHDF5_MAX_STRING).
!
subroutine rhdf5_write_variable(id, vname, ndims, itstep, dims, units, descrip, dimnames, &
  rdata, idata, sdata, ssize, cdata)
  implicit none

  integer :: id
  character (len=*) :: vname
  integer :: ndims
  integer :: itstep
  integer, dimension(ndims) :: dims
  character (len=*) :: units
  character (len=*) :: descrip
  character (len=RHDF5_MAX_STRING), dimension(ndims) :: dimnames
  real, dimension(*), optional :: rdata
  integer, dimension(*), optional :: idata
  character (len=*), dimension(*), optional :: sdata
  character (len=*), optional :: cdata
  integer, optional :: ssize

  integer :: hdferr
  integer :: dtype
  integer :: dsetid
  integer :: deflvl
  integer :: dset_ndims
  integer, dimension(ndims+1) :: dset_dims 
  integer :: i
  integer, dimension(ndims+1) :: ext_dims
  integer, dimension(ndims+1) :: chunk_dims
  character (len=RHDF5_MAX_STRING) :: dnstring
  character (len=RHDF5_MAX_STRING) :: arrayorg
  integer :: dsize
  integer :: idummy   ! to get the size of an integer
  real :: rdummy      ! to get the size of a real
  character :: cdummy ! to get the size of a character

  ! Check for valid argument values
  if (itstep .lt. 0) then
    print*,'ERROR: hdf5_write_variable: itstep argument must be >= zero'
    stop 'hdf5_write_variable: bad variable write'
  endif

  if (present(sdata)) then
    dtype = RHDF5_TYPE_STRING
    if (.not.present(ssize)) then
      print*,'ERROR: hdf5_write_variable: must use "ssize" argument when using "sdata" argument'
      stop 'hdf5_write_variable: bad variable write'
    endif
    dsize = ssize
  elseif (present(idata)) then
    dtype = RHDF5_TYPE_INTEGER
    dsize = sizeof(idummy)
  elseif (present(rdata)) then
    dtype = RHDF5_TYPE_FLOAT
    dsize = sizeof(rdummy)
  elseif (present(cdata)) then
    dtype = RHDF5_TYPE_CHAR
    dsize = sizeof(cdummy)
  else
    print*,'ERROR: hdf5_write_variable: must use one of the "rdata", "idata", "sdata" arguments'
    stop 'hdf5_write_variable: bad variable write'
  endif

  ! If ndims is zero then the data coming in is time data. This is okay as long as itstep is not zero
  ! as well (which would mean that there is no data to write). Print a warning and return if ndims
  ! and itstep are both zero.
  if ((ndims .eq. 0) .and. (itstep .eq. 0)) then
    print*,'WARNING: hdf5_write_variable: Both ndims and itstep are zero, skipping write for variable: ', trim(vname)
    return
  endif

  ! Array organization is always row major for now
  arrayorg = 'row major'

  ! Figure out what the "real" dimensions for the data will be. Takes into account that the x,y,z
  ! dimensions need to be reversed for row major ordering
  call rhdf5_build_dims_for_write(itstep, ndims, dims, dimnames, &
    dset_ndims, dset_dims, ext_dims, chunk_dims, dnstring)
  call rhdf5_adjust_chunk_sizes(dset_ndims, chunk_dims, dsize)

  ! write out the data
  !deflvl = 6
  deflvl = 1
  if (present(sdata)) then
    ! Use the original ndims and dims values for the call to rhd5_add_null_bytes(). Since
    ! this routine treats the sdata array as a 1D array with a total number of elements
    ! given by the product of the sizes in dims, it doesn't matter what order the dimension
    ! sizes are listed in dims. We also want to use dims and ndims since they describe what
    ! is in sdata.
    !
    call rhdf5_add_null_bytes(sdata, ssize, ndims, dims)
    call rh5d_setup_and_write(id, trim(vname)//char(0), dtype, ssize, dset_ndims, dset_dims, ext_dims, &
      chunk_dims, deflvl, dsetid, sdata, hdferr)
  elseif (present(idata)) then
    call rh5d_setup_and_write(id, trim(vname)//char(0), dtype, 0, dset_ndims, dset_dims, ext_dims, &
      chunk_dims, deflvl, dsetid, idata, hdferr)
  elseif (present(rdata)) then
    call rh5d_setup_and_write(id, trim(vname)//char(0), dtype, 0, dset_ndims, dset_dims, ext_dims, &
      chunk_dims, deflvl, dsetid, rdata, hdferr)
  elseif (present(cdata)) then
    call rh5d_setup_and_write(id, trim(vname)//char(0), dtype, 0, dset_ndims, dset_dims, ext_dims, &
      chunk_dims, deflvl, dsetid, cdata, hdferr)
  endif
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_write_variable: cannot write data for variable: ',trim(vname)
    stop 'hdf5_write_variable: bad variable write'
  endif

  ! write out the attributes
  call rh5a_write_anyscalar(dsetid, 'units'//char(0),       trim(units)//char(0), &
    RHDF5_TYPE_STRING, hdferr)
  call rh5a_write_anyscalar(dsetid, 'long_name'//char(0), trim(descrip)//char(0), &
    RHDF5_TYPE_STRING, hdferr)
  call rh5a_write_anyscalar(dsetid, 'ArrayOrg'//char(0),    trim(arrayorg)//char(0), &
    RHDF5_TYPE_STRING, hdferr)
  call rh5a_write_anyscalar(dsetid, 'DimNames'//char(0),    trim(dnstring)//char(0), &
    RHDF5_TYPE_STRING, hdferr)
  if (trim(vname) .eq. 'time') then
    ! helper attribute for ncdump
    call rh5a_write_anyscalar(dsetid, 'C_format'//char(0), '%c'//char(0), RHDF5_TYPE_STRING, hdferr)
  endif

  ! close the dataset
  call rh5d_close(dsetid, hdferr)

  return
end subroutine rhdf5_write_variable

!*****************************************************************************
! rhdf5_build_dims_for_write()
!
! This routine will figure out the dimensions we want for the upcoming write.
!
! If itstep is > 0, then the time dimension will be tagged onto the front of
! dims creating one extra dimension to the data.
!
! ext_dims gets set to show which dimensions, if any, are extendable (only
! extending time dimension for now).
!
! chunk_dims shows how to organize chunking in the HDF5 file.
!
subroutine rhdf5_build_dims_for_write(itstep, ndims, dims, dimnames, &
  dset_ndims, dset_dims, ext_dims, chunk_dims, dnstring)

  implicit none

  integer :: itstep, ndims
  integer, dimension(ndims) :: dims
  character (len=RHDF5_MAX_STRING), dimension(ndims) :: dimnames
  integer :: dset_ndims
  integer, dimension(ndims+1) :: dset_dims, ext_dims, chunk_dims
  character (len=RHDF5_MAX_STRING) :: dnstring

  integer :: i, irev
  integer :: dummy_int
  real :: dummy_real
  character :: dummy_char

  ! The array ext_dims holds a 1 (extendable) or 0 (not extendable) to mark which,
  ! if any, dimensions of the variable are to be extendable. For now we just want
  ! the time dimension, if it exists, to be extendable.
  !
  ! The array chunk_dims holds a description of how big to make the chunks in the HDF file. We want
  ! one chunk per time step. Always make the chunk size match the size of the field (set the individual
  ! chunk sizes to the size of the dimensions of the field). If the variable has a time dimension
  ! set the corresponding chunk_size for the time dimension (first dimension) to one which keep the
  ! overall chunk size matching the size of the field.
  !
  ! Always place time as the first dimension. Then reverse the order of the field dimensions which
  ! will automatically make the data appear to be row major ordered.
  !
  if (itstep .eq. 0) then
    dset_ndims = ndims
    do i = 1, ndims
      ! reverse the ordering
      irev = (ndims - i) + 1
      dset_dims(irev) = dims(i)
      ext_dims(irev) = 0
      chunk_dims(irev) = dims(i)
      if (i .eq. 1) then
        dnstring = trim(dimnames(irev))
      else
        dnstring = trim(dnstring) // ' ' // trim(dimnames(irev))
      endif
    enddo
  elseif (itstep .gt. 0) then
    dset_ndims = ndims + 1
    dset_dims(1) = itstep
    ext_dims(1) = 1
    chunk_dims(1) = 1
    dnstring = 't'
    do i = 1, ndims
      ! reverse the ordering, time is in the first slot
      ! so shift this up by one
      irev = (ndims - i) + 2
      dset_dims(irev) = dims(i)
      ext_dims(irev) = 0
      chunk_dims(irev) = dims(i)
      ! note that the dimnames array is not to be shifted
      dnstring = trim(dnstring) // ' ' // trim(dimnames(irev-1))
    enddo
  endif

  return
end subroutine rhdf5_build_dims_for_write

!********************************************************************
! rhdf5_adjust_chunk_sizes()
!
! This routine will attempt to keep the chunk size at optimal
! values for the dataset.
!
! First of all, the HDF5 interface uses a default chunk cache (hash
! table) with 521 entries where each entry is 1MB in size. For small
! datasets (data buffer for each write under 512KB) there is no
! need to change the chunk sizes since they default to the data
! dimension sizes making the entire data buffer one chunk.
!
! For large datasets (data buffer for each write > 512KB) diffferent
! sources have different recommendations for the optimal chunk size.
! The HDF5 documentation recommendation is: < 512KB. They 
! all agree that the thing that really matters is the access order.
!
! Since we are just streaming out the data time step by time step
! with no subset selection, we need to make the chunks fit in a
! contiguous line through the memory. We also would like to make
! all the chunks add up to the total data buffer size without
! any extra padding so that we keep the file size at a minimum.
! Should be able to accomplish this by dividing up each dimension
! starting with the first (slowest changing) dimension into smaller
! pieces (that don't leave any remainder) until the chunk size is
! at the right size. For exmaple if you have 4-byte data organized as
! (20, 50, 1000), 4 million bytes, then allow the chunk size
! (2, 50, 1000), 400 thousand bytes, but don't allow (3, 50, 1000),
! 600 thousand bytes (closer to 512KB) since 3 does not divide into
! 20 evenly resulting in unused space in the dataset.
! 
subroutine rhdf5_adjust_chunk_sizes(ndims, chunk_dims, dsize)
  implicit none
  integer, parameter :: KB_512 = 524288

  integer :: ndims, dsize
  integer, dimension(ndims) :: chunk_dims

  integer :: i, chunk_size, new_dim, shrink
  integer :: rhdf5_find_even_divisor

  chunk_size = dsize
  do i = 1, ndims
    chunk_size = chunk_size * chunk_dims(i)
  enddo

  ! only adjust if the default chunk size is greater than 512KB
  if (chunk_size .gt. KB_512) then
    shrink = ceiling(float(chunk_size) / float(KB_512))

    ! shrink is the factor that we want to reduce the chunk size by,
    ! but we want the actual reduction to be an even divisor of the
    ! given chunk dimensions so that the chunks fit over the data
    ! exactly with no wasted space.
    !
    ! Walk through the chunk dimensions. Keep setting the current
    ! chunk dimension to 1 (and dividing the shrink factor by that
    ! dimension) while the shrink factor remains larger than that
    ! dimension. As soon as the dimension is greater than the shrink
    ! factor, divide up the dimension with the nearest even divisor
    ! less than the current dimension divided by the shrink factor.
    !
    ! The intent of this algorithm is to yield exact fitting tiles
    ! (chunks) over the data that have their size as close to 512KB
    ! as possible without going over 512KB.
    do i = 1, ndims
      if (shrink .gt. chunk_dims(i)) then
        shrink = ceiling(float(shrink) / float(chunk_dims(i)))
        chunk_dims(i) = 1
      elseif (shrink .gt. 1) then
        new_dim = floor(float(chunk_dims(i))/float(shrink))
        chunk_dims(i) = rhdf5_find_even_divisor(chunk_dims(i), new_dim)
        shrink = 1
      endif
    enddo
  endif

  return
end subroutine rhdf5_adjust_chunk_sizes

!********************************************************************
! rhdf5_add_null_bytes()
!
! This routine will walk through the sdata array adding null bytes
! to the ends of all the strings.
subroutine rhdf5_add_null_bytes(sdata,ssize,ndims,dims)
  implicit none

  integer :: ndims, ssize
  integer, dimension(*) :: dims
  character (len=ssize), dimension(*) :: sdata

  integer :: i, ntot

  ! No matter what the dimensions of sdata really are, the strings are all lined up in contiguous
  ! memory, each one spaced apart by ssize bytes. Therefore we can treat sdata as a 1D array
  ! of strings with length ssize, and just step through as many strings as the actual dimensions
  ! dicatate.
  ntot = 1
  do i = 1, ndims
    ntot = ntot * dims(i)
  enddo

  do i = 1, ntot
    sdata(i) = trim(sdata(i)) // char(0)
  enddo

  return
end subroutine rhdf5_add_null_bytes

!********************************************************************
! rhdf5_find_even_divisor()
!
! This function will find the nearest divisor of num which is
! <= div that divides num evenly.
!
integer function rhdf5_find_even_divisor(num, div)
  implicit none

  integer :: num, div

  integer :: i

  rhdf5_find_even_divisor = 1
  do i = div, 1, -1
    if (mod(num, i) .eq. 0) then
      rhdf5_find_even_divisor = i
      exit
    endif
  enddo

  return
end function rhdf5_find_even_divisor

!********************************************************************
! rhdf5_read_variable_init()
!
! This routine read the dimension and attribute info for an hdf5 file variable.
!
! This routine needs to be called before hdf5_read_variable().
!
subroutine rhdf5_read_variable_init(id, vname, ndims, dims, units, descrip, dimnames)
  implicit none

  integer :: id
  character (len=*) :: vname
  integer :: ndims
  integer, dimension(RHDF5_MAX_DIMS) :: dims
  character (len=RHDF5_MAX_STRING) :: units
  character (len=RHDF5_MAX_STRING) :: descrip
  character (len=RHDF5_MAX_STRING), dimension(RHDF5_MAX_DIMS) :: dimnames

  integer :: hdferr
  integer :: dsetid
  character (len=RHDF5_MAX_STRING) :: dnstring
  character (len=RHDF5_MAX_STRING) :: arrayorg
  character (len=RHDF5_MAX_STRING) :: stemp

  ! open dataset
  call rh5d_open(id, trim(vname)//char(0), dsetid, hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_read_variable_init: cannot open dataset for variable: ',trim(vname)
    stop 'hdf5_read_variable: bad variable read'
  endif

  ! dimensions
  call rh5d_read_get_dims(dsetid, ndims, dims, hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_read_variable_init: cannot obtain dimension information for variable: ',trim(vname)
    stop 'hdf5_read_variable: bad variable read'
  endif

  ! read the attributes
  call rh5a_read_anyscalar(dsetid, 'units'//char(0),       stemp, RHDF5_TYPE_STRING, hdferr)
  call rhdf5_c2f_string(stemp, units, RHDF5_MAX_STRING)
  call rh5a_read_anyscalar(dsetid, 'long_name'//char(0), stemp, RHDF5_TYPE_STRING, hdferr)
  call rhdf5_c2f_string(stemp, descrip, RHDF5_MAX_STRING)
  call rh5a_read_anyscalar(dsetid, 'ArrayOrg'//char(0),    stemp, RHDF5_TYPE_STRING, hdferr)
  call rhdf5_c2f_string(stemp, arrayorg, RHDF5_MAX_STRING)
  call rh5a_read_anyscalar(dsetid, 'DimNames'//char(0),    stemp, RHDF5_TYPE_STRING, hdferr)
  call rhdf5_c2f_string(stemp, dnstring, RHDF5_MAX_STRING)
  call rhdf5_read_dim_name_string(ndims, dimnames, dnstring)

  ! If arrayorg is "row major" then reverse the dimensions
  if (trim(arrayorg) .eq. 'row major') then
    call rhdf5_reverse_dims(ndims, dims, dimnames)
  endif

  ! close the dataset
  call rh5d_close(dsetid, hdferr)

  return
end subroutine rhdf5_read_variable_init

!********************************************************************
! rhdf5_read_variable()
!
! This routine read the dataset for an hdf5 file variable.
!
! The caller is responsible for allocating the memory for the data buffer.
!
! Note that sdata is an array of strings where each string gets cast into a C style string
! terminated with a null byte. cdata on the other hand gets cast into a character array with
! a fixed length (RHDF5_MAX_STRING).
!
subroutine rhdf5_read_variable(id, vname, ndims, dims, rdata, idata, sdata, ssize, cdata)
  implicit none

  integer :: id
  character (len=*) :: vname
  integer :: ndims
  integer, dimension(ndims) :: dims
  real, dimension(:), optional, allocatable :: rdata
  integer, dimension(:), optional, allocatable :: idata
  character (len=*), dimension(:), optional, allocatable :: sdata
  integer, optional :: ssize
  character, dimension(:), optional, allocatable :: cdata

  integer :: hdferr
  integer :: dsetid
  character (len=RHDF5_MAX_STRING) :: dnstring
  character (len=RHDF5_MAX_STRING) :: arrayorg
  character (len=RHDF5_MAX_STRING) :: stemp
  integer :: dsize

  ! open dataset
  call rh5d_open(id, trim(vname)//char(0), dsetid, hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_read_variable: cannot open dataset for variable: ',trim(vname)
    stop 'hdf5_read_variable: bad variable read'
  endif

  ! read the data
  if (present(sdata)) then
    call rh5d_read(dsetid, RHDF5_TYPE_STRING, ssize, sdata, hdferr)
    call rhdf5_convert_c2f_strings(sdata,ssize,ndims,dims)
  elseif (present(idata)) then
    call rh5d_read(dsetid, RHDF5_TYPE_INTEGER, dsize, idata, hdferr)
  elseif (present(rdata)) then
    call rh5d_read(dsetid, RHDF5_TYPE_FLOAT, dsize, rdata, hdferr)
  elseif (present(cdata)) then
    call rh5d_read(dsetid, RHDF5_TYPE_CHAR, dsize, cdata, hdferr)
  else
    print*,'ERROR: hdf5_read_variable: must use one of the "rdata", "idata", "sdata" arguments'
    stop 'hdf5_read_variable: bad variable read'
  endif
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_read_variable: cannot read data for variable: ',trim(vname)
    stop 'hdf5_read_variable: bad variable read'
  endif

  ! close the dataset
  call rh5d_close(dsetid, hdferr)

  return
end subroutine rhdf5_read_variable

!********************************************************************
! rhdf5_convert_c2f_strings()
!
! This routine will walk through the sdata and convert the C strings
! (null terminated) to FORTRAN strings.
subroutine rhdf5_convert_c2f_strings(sdata,ssize,ndims,dims)
  implicit none

  integer :: ndims, ssize
  integer, dimension(*) :: dims
  character (len=ssize), dimension(*) :: sdata

  integer :: i, ntot
  character (len=ssize) :: stemp

  ! No matter what the dimensions of sdata really are, the strings are all lined up in contiguous
  ! memory, each one spaced apart by ssize bytes. Therefore we can treat sdata as a 1D array
  ! of strings with length ssize, and just step through as many strings as the actual dimensions
  ! dicatate.
  ntot = 1
  do i = 1, ndims
    ntot = ntot * dims(i)
  enddo

  do i = 1, ntot
    call rhdf5_c2f_string(sdata(i), stemp, ssize)
    sdata(i) = stemp
  enddo

  return
end subroutine rhdf5_convert_c2f_strings

!***********************************************************************
! rhdf5_c2f_string()
!
! This routine will convert the C style string (null terminated) in
! cstring and convert it to a FORTRAN style string in fstring.

subroutine rhdf5_c2f_string(cstring, fstring, ssize)
  implicit none

  integer :: ssize
  character (len=ssize) :: cstring
  character (len=ssize) :: fstring

  integer :: i
  integer :: null_pos

  do i = 1, ssize
    if (ichar(cstring(i:i)) .eq. 0) then
      null_pos = i
      exit
    endif
  enddo

  if (null_pos .gt. 1) then
    fstring = cstring(1:null_pos-1)
  else
    fstring = ''
  endif

  return
end subroutine rhdf5_c2f_string


!***********************************************************************
! rhdf5_c2f_string()
!
! This routine will reverse the dimensions and dimension names. This is
! intended to be used for reading HDF5 files with row major storage.
!
subroutine rhdf5_reverse_dims(ndims, dims, dimnames)
  implicit none

  integer :: ndims
  integer, dimension(ndims) :: dims
  character (len=RHDF5_MAX_STRING), dimension(ndims) :: dimnames

  integer :: i
  integer :: irev
  integer, dimension(ndims) :: temp_dims
  character (len=RHDF5_MAX_STRING), dimension(ndims) :: temp_dimnames

  ! copy dims into temp_dims in reverse order
  ! then copy temp_dims back into dims
  do i = 1, ndims
    irev = (ndims - i) + 1
    temp_dims(irev) = dims(i)
    temp_dimnames(irev) = trim(dimnames(i))
  enddo

  do i = 1, ndims
    dims(i) = temp_dims(i)
    dimnames(i) = trim(temp_dimnames(i))
  enddo

  return
end subroutine rhdf5_reverse_dims

!***********************************************************************
! rhdf5_read_dim_name_string()
!
! This routine will read the dimension name string (dnstring) and load
! up the entries in the dimnames array. The different names in dnstring
! are space separated.
!
subroutine rhdf5_read_dim_name_string(ndims, dimnames, dnstring)
  implicit none

  integer :: ndims
  character (len=RHDF5_MAX_STRING), dimension(ndims) :: dimnames
  character (len=RHDF5_MAX_STRING) :: dnstring

  integer :: i
  integer :: j

  ! initialize dimnames to null strings
  do j = 1, ndims
    dimnames(j) = ''
  enddo

  j = 1 ! index into dimnames
  do i = 1, len_trim(dnstring)
    if (dnstring(i:i) .ne. ' ') then
       ! keep appending dnstring(i) to dimnames(j)
       dimnames(j) = trim(dimnames(j)) // dnstring(i:i)
    else
       ! hit a space, bump j to next entry
       j = j + 1
    endif
  enddo

  return
end subroutine rhdf5_read_dim_name_string

!***********************************************************************
! rhdf5_set_var_to_dim()
!
! This routine will set a "dimension scale" on the given variable so that
! variable can be used as dimensions coordinate values in the hdf5 file.
!
subroutine rhdf5_set_var_to_dim(fileid, vname, dname)
  implicit none

  integer :: fileid
  character (len=*) :: vname
  character (len=*) :: dname

  integer :: hdferr
  integer :: dsetid

  call rh5d_open(fileid, trim(vname)//char(0), dsetid, hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_set_var_to_dim: cannot open dataset for variable: ',trim(vname)
    stop 'hdf5_set_var_to_dim: bad variable read'
  endif

  ! set the scale on the variable
  ! add a new attribute "axis" to the variable for GRADS
  call rh5ds_set_scale(dsetid, trim(dname)//char(0), hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_set_var_to_dim: cannot set coordinate variable to a dimension'
    stop 'hdf5_set_var_to_dim: bad variable set dimension'
  endif

  call rh5a_write_anyscalar(dsetid, 'axis'//char(0), trim(dname)//char(0), RHDF5_TYPE_STRING, hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_set_var_to_dim: cannot add "axis" attribute to variable: ', trim(vname)
    stop 'hdf5_set_var_to_dim: bad variable set dimension'
  endif

  ! close the dataset
  call rh5d_close(dsetid, hdferr)

  return
end subroutine rhdf5_set_var_to_dim

!***********************************************************************
! rhdf5_attach_dims_to_var()
!
! This routine will attach dimensions specs to the given variable. This
! is done via the HDF5 dimension scaling feature, and the purpose of
! doing the attachment is so that GRADS can read in the resulting
! HDF5 file directly without the use of a descriptor file.
!
subroutine rhdf5_attach_dims_to_var(fileid, vname)
  implicit none

  integer :: fileid
  character (len=*) :: vname

  integer :: i
  integer :: ndims
  integer, dimension(RHDF5_MAX_DIMS) :: dims
  character (len=RHDF5_MAX_STRING) :: units
  character (len=RHDF5_MAX_STRING) :: descrip
  character (len=RHDF5_MAX_STRING), dimension(RHDF5_MAX_DIMS) :: dimnames

  integer :: hdferr
  integer :: dsetid
  integer :: dsclid
  character (len=RHDF5_MAX_STRING) :: dnstring
  character (len=RHDF5_MAX_STRING) :: stemp
  character (len=RHDF5_MAX_STRING) :: coord_name

  ! open dataset
  call rh5d_open(fileid, trim(vname)//char(0), dsetid, hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_attach_dims_to_var: cannot open dataset for variable: ',trim(vname)
    stop 'hdf5_attach_dims_to_var: bad variable read'
  endif

  ! dimensions
  call rh5d_read_get_dims(dsetid, ndims, dims, hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_attach_dims_to_var: cannot obtain dimension information for variable: ',trim(vname)
    stop 'hdf5_attach_dims_to_var: bad variable read'
  endif

  ! read the DimNames attribute
  call rh5a_read_anyscalar(dsetid, 'DimNames'//char(0),    stemp, RHDF5_TYPE_STRING, hdferr)
  call rhdf5_c2f_string(stemp, dnstring, RHDF5_MAX_STRING)
  call rhdf5_read_dim_name_string(ndims, dimnames, dnstring)

  ! There is no need to reverse the dimensions since we are going immediately back
  ! into the C interface to do the dimension scale attach.
  !
  ! Assume that the coordinate information is stored in datasets with names
  ! that reflect the dimension names as shown below:
  !     dim name            dataset name with coordinates
  !       x                     x_coords
  !       y                     y_coords
  !       z                     z_coords
  !       t                     t_coords
  !      <n>                    <n>_coords  (in general)

  do i = 1, ndims
    coord_name = trim(dimnames(i)) // '_coords'

    call rh5d_open(fileid, trim(coord_name)//char(0), dsclid, hdferr)
    if (hdferr .ne. 0) then
      print*,'ERROR: hdf5_attach_dims_to_var: cannot open dataset for coordinates: ',trim(coord_name)
      stop 'hdf5_attach_dims_to_var: bad variable read'
    endif

    ! do the attach, C indices start with zero so send (i-1) to this routine
    call rh5ds_attach_scale(dsetid, dsclid, i-1, hdferr)
    if (hdferr .ne. 0) then
      print*,'ERROR: hdf5_attach_dims_to_var: cannot attach coordinate data to variable:'
      print*,'ERROR:  Variable dataset: ',trim(vname)
      print*,'ERROR:  Coordinate dataset: ',trim(coord_name)
      stop 'hdf5_attach_dims_to_var: bad variable attach'
    endif

    call rh5d_close(dsclid, hdferr)
  enddo


  ! close the dataset
  call rh5d_close(dsetid, hdferr)

  return
end subroutine rhdf5_attach_dims_to_var

end module
