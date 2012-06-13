/*
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
*/

// Using these C wrapper routines to be able to utilize the C API for HDF5. Unfortunately,
// the fortran .mod files that came with our installation of HDF5 do not work with the
// pgf90 compiler.

#include "rhdf5_snames.h"
#include <string.h>
#include "hdf5.h"

////////////////////////////////////////////////////////////////////////////
// Routines for HDF5 REVU IO. 
////////////////////////////////////////////////////////////////////////////

///////////// FILE ROUTINES //////////////////////////
void rh5f_open(char *fname, int *f_flgs, int *fileid, int *hdferr)
{
unsigned flags;

if(*f_flgs == 1) flags = H5F_ACC_RDONLY;
if(*f_flgs == 2) flags = H5F_ACC_RDWR;

*fileid = H5Fopen(fname, flags, H5P_DEFAULT);

if (*fileid < 0)
  {
  *hdferr = *fileid;
  }
else
  {
  *hdferr = 0;
  }

return;
}

void rh5f_create(char *fname, int *f_flgs, int *fileid, int *hdferr)
{
unsigned flags;

if(*f_flgs == 1) flags = H5F_ACC_TRUNC;
if(*f_flgs == 2) flags = H5F_ACC_EXCL ;

*fileid = H5Fcreate(fname, flags, H5P_DEFAULT, H5P_DEFAULT);

if (*fileid < 0)
  {
  *hdferr = *fileid;
  }
else
  {
  *hdferr = 0;
  }

return;
}

void rh5f_close(int *fileid, int *hdferr)
{

*hdferr = H5Fclose(*fileid);

return;
}

////////////////// GROUP ROUTINES ////////////////

void rh5g_create(int *id, char *name, int *gid, int *hdferr)
{

#ifdef H5_USE_16_API
// 1.6 API
// 3rd arg is size hint
*gid = H5Gcreate(*id, name, 0);
#else
// 1.8 API
// 3rd arg is link creation property list
// 4th arg is group creation property list
// 5th art is group access property list (must be H5P_DEFAULT as of 4/19/12)
*gid = H5Gcreate(*id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

if (*gid < 0)
  {
  *hdferr = *gid;
  }
else
  {
  *hdferr = 0;
  }
 
return;
}

void rh5g_open(int *id, char *name, int *gid, int *hdferr)
{

#ifdef H5_USE_16_API
// 1.6 API
*gid = H5Gopen(*id, name);
#else
// 1.8 API
// 3rd arg is group access property list
*gid = H5Gopen(*id, name, H5P_DEFAULT);
#endif

if (*gid < 0)
  {
  *hdferr = *gid;
  }
else
  {
  *hdferr = 0;
  }
 
return;
}

void rh5g_close(int *grpid, int *hdferr)
{

*hdferr = H5Gclose(*grpid);

return;
}

////////////// LINK (GROUP, DATASET, etc) ROUTINES ///////////////

void rh5l_exists(int *id, char *name, int *exists)
{

// Third arg is "link access property list id"
*exists = H5Lexists(*id, name, H5P_DEFAULT);

return;
}

////////////// DATASET ROUTINES ////////////////

void rh5d_open(int *id, char *name, int *dsetid, int *hdferr)
{

// Open the existing dataset
#ifdef H5_USE_16_API
// 1.6 API
*dsetid = H5Dopen(*id, name);
#else
// 1.8 API
*dsetid = H5Dopen(*id, name, H5P_DEFAULT);
#endif

if (*dsetid < 0)
  {
  *hdferr = *dsetid;
  }
else
  {
  *hdferr = 0;
  }

return;
}

//**********************************************************************
// rh5d_setup_and_write()
//
// This routine will either create or extend a dataset and then
// perform the write into that dataset.
//
void rh5d_setup_and_write(int *id, char *name, int *dtype, int *ssize,
   int *ndims, int *dims, int *ext_dims, int *chunk_dims,
   int *deflvl, int *dsetid, void *data, int *hdferr)
{
int i;
int cur_ndims;
hsize_t cur_dims[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t max_dims[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t new_dims[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t mem_dims[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t chunk_sizes[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
hsize_t hs_count[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t hs_start[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hid_t mtypid;
hid_t mspcid;
hid_t fspcid;
hid_t dcplid;

*hdferr = 0;

// set the memory type (translate integer code in dtype to HDF5 type)
mtypid = rf2c_set_mem_type(*dtype, *ssize);
if (mtypid < 0)
  {
  printf("ERROR: rh5d_setup_append: unrecognized dtype: %d\n", *dtype);
  *hdferr = -1;
  return;
  }

// Check to see if the dataset exists aleady. If so extend it and do
// the write. If not, then create it and do the write.
if (H5Lexists(*id, name, H5P_DEFAULT))
  {
  // dataset exists --> extend and write
  //
  // Open the existing dataset
#ifdef H5_USE_16_API
  // 1.6 API
  *dsetid = H5Dopen(*id, name);
#else
  // 1.8 API
  *dsetid = H5Dopen(*id, name, H5P_DEFAULT);
#endif

  // Get the currently existing dimension information from the dataset
  // Close the dataspce after grabbing the current dimension information
  // since we need to reopen the dataspace after extending the dataset.
  fspcid = H5Dget_space(*dsetid);
  if (fspcid < 0)
    {
    printf("ERROR: rh5d_setup_append: Cannot get file space id from dataset\n");
    *hdferr = -1;
    return;
    }
  cur_ndims = H5Sget_simple_extent_dims(fspcid, cur_dims, max_dims);
  H5Sclose(fspcid);
  
  if (cur_ndims != *ndims)
    {
    printf("ERROR: rh5d_setup_append: number of dimensions in dataset (%d) does not match specified number of dimensions (ndims): %d\n",cur_ndims, *ndims);
    *hdferr = -1;
    return;
    }
  
  // Get the new dimension information from the input arguments
  for (i = 0; i<*ndims; i++)
    {
    new_dims[i] = dims[i];
    hs_count[i] = cur_dims[i];  // hs_* vars are for the hyperslab selection later on
    hs_start[i] = 0;
    if (ext_dims[i] == 1)
      {
      mem_dims[i] = 1; // this is the dimension we are extended and assume for now
                       // that we are extending by just one element
  
      if (new_dims[i] < cur_dims[i])
        {
        printf("ERROR: rh5d_setup_append: Shrinking an extendible dimension is not supported\n");
        printf("ERROR:   Current dimension number: %d, Current dimension size: %d\n", i+1, (int) cur_dims[i]);
        printf("ERROR:   New dimension number:     %d, New dimension size:     %d\n", i+1, (int) new_dims[i]);
        *hdferr = -1;
        return;
        }
      }
    else
      {
      mem_dims[i] = dims[i];
      }
    }
  
  // Create a memory dataspace to control the selection from the data buffer (input arg "data")
  // We are assuming that the entries in the data buffer are organized so that the dimension
  // that we are extending is a size of one so if we set up the memory space with that in mind
  // the selection will happen automatically without the need to apply any selection mechanism
  // on the memory space.
  mspcid = H5Screate_simple(cur_ndims, mem_dims, NULL);
  
  // Extend the dataset and re-open the dataspace
  H5Dset_extent(*dsetid, new_dims);
  fspcid = H5Dget_space(*dsetid);
  
  if ((mspcid < 0) || (fspcid < 0))
    {
    *hdferr = -1;
    return;
    }
  
  // Create a hyperslab to select just the region where the new data will be stored
  // into the file. That is, we are appending the new data just after the end of the
  // existing data so select out the existing data. To do this first select the entire
  // dataspace (extended) and subtract out the dataspace prior to the extension. The
  // select operator H5S_SELECT_NOTB is used to do the subtraction.
  H5Sselect_all(fspcid);
  *hdferr = H5Sselect_hyperslab(fspcid, H5S_SELECT_NOTB, hs_start, NULL, hs_count, NULL);
  
  if (*hdferr != 0)
    {
    return;
    }

  // write and free up resources
  *hdferr = H5Dwrite(*dsetid, mtypid, mspcid, fspcid, H5P_DEFAULT, data);

  if (*dtype == 1)
    {
    H5Tclose(mtypid);
    }
  H5Sclose(mspcid);
  H5Sclose(fspcid);
  }
else
  {
  // dataset does not exist --> create and write
  //
  // H5Screate_simple wants all of the elements of the dims and maxdims arrays (2nd, 3rd args)
  // to be set to a non-zero value, even the elements that are not being used (since they
  // are beyond what ndims specifies). Copy the used contents of dims (the first ndims entries)
  // into an internal array which is initialized to all 1's so that the caller doesn't need
  // to bother with this detail.
  //
  // H5Pset_chunk wants elements not used in the chunk size array (2nd arg) to be set to
  // zero and the used elements set to non-zero values.
  for (i=0 ; i<*ndims; i++)
    {
    cur_dims[i] = dims[i];
    if (ext_dims[i] == 0)
      {
      max_dims[i] = dims[i];
      }
    else
      {
      max_dims[i] = H5S_UNLIMITED;
      }
    chunk_sizes[i] = chunk_dims[i];
    }
  
  // create the file dataspace
  mspcid = H5Screate_simple(*ndims, cur_dims, max_dims);
  if (mspcid < 0)
    {
    printf("ERROR: rh5d_create: Cannot create memory space for dataset\n");
    *hdferr = -1;
    return;
    }
  
  // create the dataset
  // dataset creation properties
  //
  
  dcplid = H5Pcreate(H5P_DATASET_CREATE);
  if (dcplid < 0)
    {
    printf("ERROR: rh5d_create: Cannot make create properties for dataset\n");
    *hdferr = -1;
    return;
    }
  
  // Set up for data compression
  //   enable chunking of data
  //   enable shuffling of data
  //   enable deflation (compression) of data
  //
  H5Pset_chunk(dcplid, *ndims, chunk_sizes);
  H5Pset_shuffle(dcplid);
  H5Pset_deflate(dcplid, *deflvl);
  
  // dataset creation properties have been successfully built
  
  // dataset
#ifdef H5_USE_16_API
  // 1.6 API
  // 3rd arg is type of the data elements - always using float for now
  *dsetid = H5Dcreate(*id, name, mtypid, mspcid, dcplid);
#else
  // 1.8 API
  // 3rd arg is type of the data elements - always using float for now
  // 5th arg is the link creation property list (eg, automatically build intermediate groups)
  // 7th arg is the data access property list
  *dsetid = H5Dcreate(*id, name, mtypid, mspcid, H5P_DEFAULT, dcplid, H5P_DEFAULT);
#endif
  if (*dsetid < 0)
    {
    printf("ERROR: rh5d_create: Cannot make create dataset\n");
    *hdferr = -1;
    return;
    }
  
  // write and free up resources
  *hdferr = H5Dwrite(*dsetid, mtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  if (*dtype == 1)
    {
    H5Tclose(mtypid);
    }
  H5Pclose(dcplid);
  H5Sclose(mspcid);
  }

return;
}

//**********************************************************************
void rh5d_write(int *dsetid, int *mtypid, int *mspcid, int *fspcid, int *xfplid, void *data, int *hdferr)
  {
  *hdferr = H5Dwrite(*dsetid, *mtypid, *mspcid, *fspcid, *xfplid, data);

  return;
  }

//**********************************************************************
void rh5d_close(int *dsetid, int *hdferr)
{
*hdferr = H5Dclose(*dsetid);

return;
}

//////////////////// DATASPACE ROUTINES ////////////////////

//**********************************************************
void rh5s_close(int *dspcid, int *hdferr)
{
*hdferr = H5Sclose(*dspcid);

return;
}

//////////////////// PROPERTY ROUTINES ////////////////////

//**********************************************************
void rh5p_close(int *propid, int *hdferr)
{
*hdferr = H5Pclose(*propid);

return;
}

//////////////////// DATATYPE ROUTINES ////////////////////

//**********************************************************
void rh5t_close(int *dtypid, int *hdferr)
{
*hdferr = H5Tclose(*dtypid);

return;
}

//////////////////// ATTRIBUTE ROUTINES ////////////////////

//**********************************************************
// rh5a_write_anyscalar()
//
// interface to external FORTRAN caller
void rh5a_write_anyscalar(int *id, char *name, void *value, int *vtype, int *hdferr)
{
// value is a void pointer so the caller can use any type (char, int, float for now). vtype
// is used to denote what type value is. The encoding for vtype is:
//   1 -> C style string
//   2 -> integer
//   3 -> float

hid_t a_id, amem_id, atype;
htri_t attr_exists;
int ssize;

*hdferr = 0;

// Set up memory type
if (*vtype == 1)
  {
  ssize = strlen(value) + 1;
  }
else
  {
  ssize = 0;
  }

atype = rf2c_set_mem_type(*vtype, ssize);

if (atype < 0)
  {
  printf("ERROR: rh5a_write_anyscalar: unrecognized vtype: %d\n", *vtype);
  *hdferr = -1;
  return;
  }

if (*hdferr == 0)
  {
  // got a recognized type so keep going

  attr_exists = H5Aexists(*id,name);

  if (attr_exists)
    {
    // attribute exists already --> update mode

    a_id = H5Aopen(*id, name, H5P_DEFAULT);
    }
  else
    {
    // attribute does not exist --> create mode

    // scalar dataspace to hold the string for the file write
    amem_id = H5Screate(H5S_SCALAR);
  
    // create the attribute and save the result for returning the status in hdferr
    // set hdferr here since we will be closing the a_id pointer at the end
#ifdef H5_USE_16_API
    // 1.6 API
    // 5th arg is attribute creation property list (must be H5P_DEFAULT as of 4/19/12)
    a_id = H5Acreate(*id, name, atype, amem_id, H5P_DEFAULT);
#else
    // 1.8 API
    // 5th arg is attribute creation property list (must be H5P_DEFAULT as of 4/19/12)
    // 6th arg is attribute access property list (must be H5P_DEFAULT as of 4/19/12)
    a_id = H5Acreate(*id, name, atype, amem_id, H5P_DEFAULT, H5P_DEFAULT);
#endif
    }

  if (a_id < 0)
    {
    *hdferr = a_id;
    }

  // If we have an error don't attempt the write, but do continue on to
  // the close statments.

  if (*hdferr == 0)
    {
    // The third arg to H5Awrite is a pointer to the data. value is that already,
    // but this does put the onus on the caller to make sure that value is really
    // pointing to the type that matches what's in vtype
    H5Awrite(a_id, atype, value);
    }

  // clean up the attribute related pointers
  if (*vtype == 1)
    {
    // Only call the close when using the string type (which is a complex structure)
    H5Tclose(atype);
    }
  if (attr_exists == 0)
    {
    // if created the attribute, then close the memory space id
    H5Sclose(amem_id);
    }
  H5Aclose(a_id);
  }

return;
}


//////////// INTERNAL ROUTINES //////////////////

//*******************************************************************
// rf2c_set_mem_type()
//
// This routine translates the integer coded memory type to the
// corresponding HDF5 type.
//
// Note that ssize input argument is only used for string types.

hid_t rf2c_set_mem_type(int dtype, int ssize)
{
hid_t mtype;

switch (dtype)
  { 
  case 1:
    // define type for the string
    // C style string, null byte terminated, size taken from input argument ssize
    mtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(mtype, ssize);
    H5Tset_strpad(mtype, H5T_STR_NULLTERM);
    break;

  case 2:
    // integer
    mtype = H5T_NATIVE_INT;
    break;

  case 3:
    // float
    mtype = H5T_NATIVE_FLOAT;
    break;

  case 4:
    // char
    mtype = H5T_NATIVE_UCHAR;
    break;

  default:
    // unrecognized type, send back -1 to inform caller
    mtype = -1;
    break;
  }

return(mtype);
}
