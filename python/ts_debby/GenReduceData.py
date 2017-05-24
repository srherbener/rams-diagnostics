#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
import numpy as np
import ConfigTsd as conf
import Hdf5Utils as h5u
from skimage.measure import block_reduce

Tstring = conf.SetTimeString()

SimList = [
    'TSD_SAL_DUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_NONSAL_NODUST'
    ]
Nsims = len(SimList)

FileList = [
    # 2D vars
#    [ 'HDF5/<SIM>/HDF5/sea_press-<SIM>-AS-2006-08-20-120000-g3.h5', '/sea_press', 'HDF5/<SIM>/HDF5/sea_press_lite-<SIM>-AS-2006-08-20-120000-g3.h5', '/sea_press'       ],

    # Vertical coords = sigma-z
    [ 'HDF5/<SIM>/HDF5/u-<SIM>-AS-2006-08-20-120000-g3.h5',       '/u',       'HDF5/<SIM>/HDF5/u_lite-<SIM>-AS-2006-08-20-120000-g3.h5',       '/u'       ],
    [ 'HDF5/<SIM>/HDF5/v-<SIM>-AS-2006-08-20-120000-g3.h5',       '/v',       'HDF5/<SIM>/HDF5/v_lite-<SIM>-AS-2006-08-20-120000-g3.h5',       '/v'       ],
    [ 'HDF5/<SIM>/HDF5/tempc-<SIM>-AS-2006-08-20-120000-g3.h5',   '/tempc',   'HDF5/<SIM>/HDF5/tempc_lite-<SIM>-AS-2006-08-20-120000-g3.h5',   '/tempc'   ],
    [ 'HDF5/<SIM>/HDF5/theta-<SIM>-AS-2006-08-20-120000-g3.h5',   '/theta',   'HDF5/<SIM>/HDF5/theta_lite-<SIM>-AS-2006-08-20-120000-g3.h5',   '/theta'   ],
    [ 'HDF5/<SIM>/HDF5/theta_e-<SIM>-AS-2006-08-20-120000-g3.h5', '/theta_e', 'HDF5/<SIM>/HDF5/theta_e_lite-<SIM>-AS-2006-08-20-120000-g3.h5', '/theta_e' ],
    [ 'HDF5/<SIM>/HDF5/vapor-<SIM>-AS-2006-08-20-120000-g3.h5',   '/vapor',   'HDF5/<SIM>/HDF5/vapor_lite-<SIM>-AS-2006-08-20-120000-g3.h5',   '/vapor'   ],

    # Vertical coords = pressure
    [ 'HDF5/<SIM>/HDF5/u-<SIM>-AP-2006-08-20-120000-g3.h5',       '/u',       'HDF5/<SIM>/HDF5/u_lite-<SIM>-AP-2006-08-20-120000-g3.h5',       '/u'       ],
    [ 'HDF5/<SIM>/HDF5/v-<SIM>-AP-2006-08-20-120000-g3.h5',       '/v',       'HDF5/<SIM>/HDF5/v_lite-<SIM>-AP-2006-08-20-120000-g3.h5',       '/v'       ],
    [ 'HDF5/<SIM>/HDF5/tempc-<SIM>-AP-2006-08-20-120000-g3.h5',   '/tempc',   'HDF5/<SIM>/HDF5/tempc_lite-<SIM>-AP-2006-08-20-120000-g3.h5',   '/tempc'   ],
    [ 'HDF5/<SIM>/HDF5/theta-<SIM>-AP-2006-08-20-120000-g3.h5',   '/theta',   'HDF5/<SIM>/HDF5/theta_lite-<SIM>-AP-2006-08-20-120000-g3.h5',   '/theta'   ],
    [ 'HDF5/<SIM>/HDF5/theta_e-<SIM>-AP-2006-08-20-120000-g3.h5', '/theta_e', 'HDF5/<SIM>/HDF5/theta_e_lite-<SIM>-AP-2006-08-20-120000-g3.h5', '/theta_e' ],
    [ 'HDF5/<SIM>/HDF5/vapor-<SIM>-AP-2006-08-20-120000-g3.h5',   '/vapor',   'HDF5/<SIM>/HDF5/vapor_lite-<SIM>-AP-2006-08-20-120000-g3.h5',   '/vapor'   ],
    ]
Nfiles = len(FileList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

RfactX = 15
RfactY = 15

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Reducing data for simulation: {0:s}".format(Sim))
    print("  Reduction factors (x,y): {0:.2f}, {1:.2f}".format(RfactX, RfactY))
    print("")

    for ifile in range(Nfiles):
        Ifname = FileList[ifile][0].replace("<SIM>", Sim)
        Ivname = FileList[ifile][1]
        Ofname = FileList[ifile][2].replace("<SIM>", Sim)
        Ovname = FileList[ifile][3]

        # Read in the coordinates first and get the sizes of their corresponding
        # dimensions. Determine the largest evenly divisible, by x or y factors,
        # sizes for selecting the field that block_reduce will be run on.
        Ifile = h5py.File(Ifname, mode='r')

        print("  Reading {0:s} ({1:s})".format(Ifname, Xname))
        X = Ifile[Xname][...]
        Nx = len(X)

        print("  Reading {0:s} ({1:s})".format(Ifname, Yname))
        Y = Ifile[Yname][...]
        Ny = len(Y)

        print("  Reading {0:s} ({1:s})".format(Ifname, Zname))
        Z = Ifile[Zname][...]
        Nz = len(Z)

        print("  Reading {0:s} ({1:s})".format(Ifname, Tname))
        T = Ifile[Tname][...]
        Nt = len(T)
        print("")

        # Compute max size that is evenly divisible by corresponding factor.
        # Split the difference between the beginning and end of the actual size
        # and the max evenly divisible size for the selection. Note that // is
        # integer divide.
        Xsize = (Nx // RfactX) * RfactX
        Ysize = (Ny // RfactY) * RfactY

        # Form start and finish indices for x and y
        X1 = (Nx - Xsize) // 2
        X2 = X1 + Xsize
        Y1 = (Ny - Ysize) // 2
        Y2 = Y1 + Ysize

        # Reduce the x and y dimensions. Use block_reduce which requires at least a 2d
        # input array. Use meshgrid to form x and y into 2d grids and run block_reduce
        # on those.
        X = X[X1:X2]
        Y = Y[Y1:Y2]

        XG, YG = np.meshgrid(X, Y)

        Xreduce = block_reduce(XG, block_size=(1, RfactX), func=np.mean)
        Yreduce = block_reduce(YG, block_size=(RfactY, 1), func=np.mean)

        Xreduce = np.squeeze(Xreduce[0,:])
        Yreduce = np.squeeze(Yreduce[:,0])

        NxReduce = len(Xreduce)
        NyReduce = len(Yreduce)

        # Write out reduced x and y dims, plus the others.
        Ofile = h5py.File(Ofname, mode='w')

        Xdim = h5u.DimCoards(Xname, 1, [ NxReduce ], 'x')
        Xdim.Build(Ofile, Xreduce)
        Ydim = h5u.DimCoards(Yname, 1, [ NyReduce ], 'y')
        Ydim.Build(Ofile, Yreduce)
        Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
        Zdim.Build(Ofile, Z)
        Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
        Tdim.Build(Ofile, T)


        # Do one time step at a time to help reduce memory requirements. Use block_reduce to
        # to form the coarse grid by averaging the neigboring fine grid cells.
        print("  Reading {0:s} ({1:s})".format(Ifname, Ivname))

        # Create the output dataset
        OutDset = h5u.DsetCoards(Ivname, 4, [ Nt, Nz, NyReduce, NxReduce ], chunks=( 1, Nz, NyReduce, NxReduce ))
        OutDset.Create(Ofile)

        for it in range(Nt):
            # Input field is (t,z,y,x). Select along the x, y dimensions so that
            # their resultant sizes are evenly divisible by teh x, y reduction factors.
            Var = Ifile[Ivname][it,:,Y1:Y2,X1:X2]

            VarReduce = block_reduce(Var, block_size=(1, RfactY, RfactX), func=np.mean)

            Ofile[Ovname][it,:,:,:] = VarReduce

            if ((it % 10) == 0):
                print("    Processing time step: {0:d}".format(it))

        Ifile.close()
        print("    Finished: number of time steps processed: {0:d}".format(it))
        print("")


        print("  Writing {0:s} ({1:s})".format(Ofname, Ovname))
        OutDset.AttachDims(Ofile, Tdim, Zdim, Ydim, Xdim)
        print("")

        Ofile.close()

