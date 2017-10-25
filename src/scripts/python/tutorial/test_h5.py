#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))

import h5py
import numpy as np
import Hdf5Utils as h5u
import copy

TimeBaseString = '2006-08-22 06:00:00 00:00'

X = [ 1, 2, 3 ]
Y = [ 1, 2, 3, 4 ] 
Z = [ 1, 2, 3, 4, 5 ] 
T = np.arange(10) * 3600

Nx = len(X)
Ny = len(Y)
Nz = len(Z)
Nt = len(T)

Ndims = 3
Dims = [ Nx, Ny, Nz ]
RevDims = [ Nz, Ny, Nx ]

Var1 = np.ones(Dims)
Var2 = np.ones(RevDims)
Var3 = np.ones([ Nt ])

count = 0
for k in range(Nz):
    for j in range(Ny):
        for i in range(Nx):
            count += 1
            Var1[i,j,k] = count
            Var2[k,j,i] = count

for i in range(Nt):
    Var3[i] = Nt - i

# set up dimension vars
Xdim = h5u.DimCoards('/x_coords', 1, [ Nx ], 'x')
Ydim = h5u.DimCoards('/y_coords', 1, [ Ny ], 'y')
Zdim = h5u.DimCoards('/z_coords', 1, [ Nz ], 'z')
Tdim = h5u.DimCoards('/t_coords', 1, [ Nt ], 't', tstring=TimeBaseString)

Var1Dset = h5u.DsetCoards('/Var1', Ndims, Dims)
Var2Dset = h5u.DsetCoards('/Var2', Ndims, RevDims)
Var3Dset = h5u.DsetCoards('/Var3', 1, [ Nt ])

Fname = 'test.h5'

Fid = h5py.File(Fname, mode='w')
Xdim.Create(Fid, X)
Ydim.Create(Fid, Y)
Zdim.Create(Fid, Z)
Tdim.Create(Fid, T)

Var1Dset.Create(Fid, Var1, Xdim, Ydim, Zdim)
Var2Dset.Create(Fid, Var2, Zdim, Ydim, Xdim)
Var3Dset.Create(Fid, Var3, Tdim)

Fid.close()

# try reading
Fid = h5py.File(Fname, mode='r')

Dset = Fid['/Var1']
print('Var1:\n', Dset[...])

Dset = Fid['/Var2']
print('Var2:\n', Dset[...])

Dset = Fid['/Var3']
print('Var3:\n', Dset[...])

Fid.close()
