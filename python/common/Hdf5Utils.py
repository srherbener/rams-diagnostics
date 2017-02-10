##############################################################
# HDF5 file utilities
#
import h5py
import numpy as np

# classes

# Base class Dset which will be used by DimCoards and DsetCoards
class BaseDset:
    '''Class for creating a dataset in HDF5 file'''
    def __init__(self, name, ndims, dims):
        # name can be a path
        # ndims is number of dimensions
        # dims is list containing sizes of dimensions
        self.name = name
        self.ndims = ndims
        self.dims = dims

    # Build will:
    #    Build the dataset in the file
    #    Return handle to dataset
    def Build(self, Fid):
        return Fid.create_dataset(self.name, self.dims)

    # Write will:
    #   Look up dataset in the file
    #   Load variable data (Var) into dataset
    def Write(self, Fid, Var):
        Dset = Fid[self.name]
        Dset[...] = Var
         

class DimCoards(BaseDset):
    '''Class for creating dimension datasets in HDF5 file using COARDS convention'''
    def __init__(self, name, ndims, dims, kind, tstring='0000-00-00 00:00:00 00:00'):
        BaseDset.__init__(self, name, ndims, dims)
        self.kind = kind
        self.tstring = tstring
        self.units = { 'x': 'degrees_east',
                       'y': 'degrees_north',
                       'z': 'meters',
                       't': 'seconds since {0:s}'.format(self.tstring) }
        self.longnames = { 'x': 'longitude',
                           'y': 'latitude',
                           'z': 'sigma-z',
                           't': 'simulation time' }


    ######################################################
    # Create dimension
    #
    # In order to be able to read into GRADS, 'x', 'y', 'z',
    # and 't' dimensions need to exist. These may not all
    # be used by the datasets in the file that are going to
    # have dimensions attached to them. create_scale
    # will be called using the dataset object from one or
    # more of these datasets. Call create_scale from the dataset
    # just created for the coordinate variable so that
    # create_scale is called at least once for each coordinate
    # dataset.
    #
    def Create(self, Fid, Var):
        # Build dataset and write Var data into it
        Dset = BaseDset.Build(self, Fid)
        BaseDset.Write(self, Fid, Var)

        # Turn into dimension scale, and attach COARDS attributes
        Dset.dims.create_scale(Dset, self.kind)
        Dset.attrs.create('axis', np.string_(self.kind))
        Dset.attrs.create('long_name', np.string_(self.longnames[self.kind]))
        Dset.attrs.create('units', np.string_(self.units[self.kind]))


class DsetCoards(BaseDset):
    '''Class for creating datasets in HDF5 file using COARDS convention'''
    def __init__(self, name, ndims, dims):
        BaseDset.__init__(self, name, ndims, dims)

    ######################################################
    # Create dataset with dimensions
    #
    # Aruguments (*args) consist of elements which are
    # DimCoards objects.
    #
    def Create(self, Fid, *args):
        # Build dataset
        Dset = BaseDset.Build(self, Fid)

        # attach dimensions to dataset (create dimensions as needed)
        for i, arg in enumerate(args):
            DimDset = Fid[arg.name]
            Dset.dims.create_scale(DimDset, arg.kind)
            Dset.dims[i].attach_scale(DimDset)


