##############################################################
# HDF5 file utilities
#
import h5py
import numpy as np

# classes

# Base class Dset which will be used by DimCoards and DsetCoards
class BaseDset:
    '''Class for creating a dataset in HDF5 file'''
    def __init__(self, name, ndims, dims, chunks=()):
        # name can be a path
        # ndims is number of dimensions
        # dims is list containing sizes of dimensions
        self.name = name
        self.ndims = ndims
        self.dims = dims
        self.chunks = chunks

    def Create(self, Fid):
        if (self.chunks):
            Dset = Fid.create_dataset(self.name, self.dims,
              chunks=self.chunks, compression="gzip", compression_opts=6, shuffle=True)
        else:
            Dset = Fid.create_dataset(self.name, self.dims)

        return Dset


class DimCoards(BaseDset):
    '''Class for creating dimension datasets in HDF5 file using COARDS convention'''
    def __init__(self, name, ndims, dims, kind, tstring='0000-00-00 00:00:00 00:00', chunks=()):
        BaseDset.__init__(self, name, ndims, dims, chunks=chunks)
        self.kind = kind
        self.tstring = tstring
        self.units = { 'x': 'degrees_east',
                       'y': 'degrees_north',
                       'z': 'meters',
                       'p': 'millibars',
                       't': 'seconds since {0:s}'.format(self.tstring) }
        self.longnames = { 'x': 'longitude',
                           'y': 'latitude',
                           'z': 'sigma-z',
                           'p': 'pressure',
                           't': 'simulation time' }


    def Build(self, Fid, Var):
        # build the dataset
        Dset = BaseDset.Create(self, Fid)
        Dset[...] = Var
 
        # Turn into dimension scale, and attach COARDS attributes
        Dset.dims.create_scale(Dset, self.kind)
        Dset.attrs.create('axis', np.string_(self.kind))
        Dset.attrs.create('long_name', np.string_(self.longnames[self.kind]))
        Dset.attrs.create('units', np.string_(self.units[self.kind]))

        return Dset


class DsetCoards(BaseDset):
    '''Class for creating datasets in HDF5 file using COARDS convention'''
    def __init__(self, name, ndims, dims, chunks=()):
        BaseDset.__init__(self, name, ndims, dims, chunks=chunks)

    def AttachDims(self, Fid, *args):
        Dset = Fid[self.name]
        for i, arg in enumerate(args):
            DimDset = Fid[arg.name]
            Dset.dims.create_scale(DimDset, arg.kind)
            Dset.dims[i].attach_scale(DimDset)

    def Build(self, Fid, Var, *args):
        # build the dataset
        Dset = BaseDset.Create(self, Fid)
        Dset[...] = Var

        # Attach the dimensions, given by args
        self.AttachDims(Fid, Dset, *args)

