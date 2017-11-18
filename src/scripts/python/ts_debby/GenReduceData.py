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
#    [ 'SIMS/<SIM>/HDF5/topo-<SIM>-AS-2006-08-20-120000-g3.h5', '/topo', 'SIMS/<SIM>/HDF5/topo_lite-<SIM>-AS-2006-08-20-120000-g3.h5', '/topo'       ],
#    [ 'SIMS/<SIM>/HDF5/sea_press-<SIM>-AS-2006-08-20-120000-g3.h5', '/sea_press', 'SIMS/<SIM>/HDF5/sea_press_lite-<SIM>-AS-2006-08-20-120000-g3.h5', '/sea_press'       ],
#    [ 'SIMS/<SIM>/HDF5/sst-<SIM>-AS-2006-08-20-120000-g3.h5', '/sst', 'SIMS/<SIM>/HDF5/sst_lite-<SIM>-AS-2006-08-20-120000-g3.h5', '/sst'       ],
#    [ 'SIMS/<SIM>/HDF5/pcprate-<SIM>-AS-2006-08-20-120000-g3.h5', '/pcprate', 'SIMS/<SIM>/HDF5/pcprate_lite-<SIM>-AS-2006-08-20-120000-g3.h5', '/pcprate'       ],

    # Vertical coords = sigma-z
#    [ 'SIMS/<SIM>/HDF5/u-<SIM>-AS-2006-08-20-120000-g3.h5',       '/u',       'SIMS/<SIM>/HDF5/u_lite-<SIM>-AS-2006-08-20-120000-g3.h5',       '/u'       ],
#    [ 'SIMS/<SIM>/HDF5/v-<SIM>-AS-2006-08-20-120000-g3.h5',       '/v',       'SIMS/<SIM>/HDF5/v_lite-<SIM>-AS-2006-08-20-120000-g3.h5',       '/v'       ],
#    [ 'SIMS/<SIM>/HDF5/tempc-<SIM>-AS-2006-08-20-120000-g3.h5',   '/tempc',   'SIMS/<SIM>/HDF5/tempc_lite-<SIM>-AS-2006-08-20-120000-g3.h5',   '/tempc'   ],
#    [ 'SIMS/<SIM>/HDF5/theta-<SIM>-AS-2006-08-20-120000-g3.h5',   '/theta',   'SIMS/<SIM>/HDF5/theta_lite-<SIM>-AS-2006-08-20-120000-g3.h5',   '/theta'   ],
#    [ 'SIMS/<SIM>/HDF5/theta_e-<SIM>-AS-2006-08-20-120000-g3.h5', '/theta_e', 'SIMS/<SIM>/HDF5/theta_e_lite-<SIM>-AS-2006-08-20-120000-g3.h5', '/theta_e' ],
#    [ 'SIMS/<SIM>/HDF5/vapor-<SIM>-AS-2006-08-20-120000-g3.h5',   '/vapor',   'SIMS/<SIM>/HDF5/vapor_lite-<SIM>-AS-2006-08-20-120000-g3.h5',   '/vapor'   ],

#    [ 'SIMS/<SIM>/HDF5/press-<SIM>-AS-2006-08-20-120000-g3.h5',   '/press',   'SIMS/<SIM>/HDF5/press_lite-<SIM>-AS-2006-08-20-120000-g3.h5',   '/press'   ],
#    [ 'SIMS/<SIM>/HDF5/relhum-<SIM>-AS-2006-08-20-120000-g3.h5',   '/relhum',   'SIMS/<SIM>/HDF5/relhum_lite-<SIM>-AS-2006-08-20-120000-g3.h5',   '/relhum'   ],

    # Vertical coords = pressure
#    [ 'SIMS/<SIM>/HDF5/u-<SIM>-AP-2006-08-20-120000-g3.h5',       '/u',       'SIMS/<SIM>/HDF5/u_lite-<SIM>-AP-2006-08-20-120000-g3.h5',       '/u'       ],
#    [ 'SIMS/<SIM>/HDF5/v-<SIM>-AP-2006-08-20-120000-g3.h5',       '/v',       'SIMS/<SIM>/HDF5/v_lite-<SIM>-AP-2006-08-20-120000-g3.h5',       '/v'       ],
#    [ 'SIMS/<SIM>/HDF5/tempc-<SIM>-AP-2006-08-20-120000-g3.h5',   '/tempc',   'SIMS/<SIM>/HDF5/tempc_lite-<SIM>-AP-2006-08-20-120000-g3.h5',   '/tempc'   ],
#    [ 'SIMS/<SIM>/HDF5/theta-<SIM>-AP-2006-08-20-120000-g3.h5',   '/theta',   'SIMS/<SIM>/HDF5/theta_lite-<SIM>-AP-2006-08-20-120000-g3.h5',   '/theta'   ],
#    [ 'SIMS/<SIM>/HDF5/theta_e-<SIM>-AP-2006-08-20-120000-g3.h5', '/theta_e', 'SIMS/<SIM>/HDF5/theta_e_lite-<SIM>-AP-2006-08-20-120000-g3.h5', '/theta_e' ],
#    [ 'SIMS/<SIM>/HDF5/vapor-<SIM>-AP-2006-08-20-120000-g3.h5',   '/vapor',   'SIMS/<SIM>/HDF5/vapor_lite-<SIM>-AP-2006-08-20-120000-g3.h5',   '/vapor'   ],
#    [ 'SIMS/<SIM>/HDF5/relhum-<SIM>-AP-2006-08-20-120000-g3.h5',  '/relhum',  'SIMS/<SIM>/HDF5/relhum_lite-<SIM>-AP-2006-08-20-120000-g3.h5',  '/relhum'  ],
#    [ 'SIMS/<SIM>/HDF5/vapliqt-<SIM>-AP-2006-08-20-120000-g3.h5', '/vapliqt', 'SIMS/<SIM>/HDF5/vapliqt_lite-<SIM>-AP-2006-08-20-120000-g3.h5', '/vapliqt' ],
#    [ 'SIMS/<SIM>/HDF5/vapicet-<SIM>-AP-2006-08-20-120000-g3.h5', '/vapicet', 'SIMS/<SIM>/HDF5/vapicet_lite-<SIM>-AP-2006-08-20-120000-g3.h5', '/vapicet' ],

#    [ 'SIMS/<SIM>/HDF5/ice2raint-<SIM>-AP-2006-08-20-120000-g3.h5', '/ice2raint', 'SIMS/<SIM>/HDF5/ice2raint_lite-<SIM>-AP-2006-08-20-120000-g3.h5', '/ice2raint' ],
#    [ 'SIMS/<SIM>/HDF5/melticet-<SIM>-AP-2006-08-20-120000-g3.h5',  '/melticet',  'SIMS/<SIM>/HDF5/melticet_lite-<SIM>-AP-2006-08-20-120000-g3.h5',  '/melticet'  ],
#    [ 'SIMS/<SIM>/HDF5/rain2icet-<SIM>-AP-2006-08-20-120000-g3.h5', '/rain2icet', 'SIMS/<SIM>/HDF5/rain2icet_lite-<SIM>-AP-2006-08-20-120000-g3.h5', '/rain2icet' ],
#    [ 'SIMS/<SIM>/HDF5/rimecldt-<SIM>-AP-2006-08-20-120000-g3.h5',  '/rimecldt',  'SIMS/<SIM>/HDF5/rimecldt_lite-<SIM>-AP-2006-08-20-120000-g3.h5',  '/rimecldt'  ],
#    [ 'SIMS/<SIM>/HDF5/cld2raint-<SIM>-AP-2006-08-20-120000-g3.h5', '/cld2raint', 'SIMS/<SIM>/HDF5/cld2raint_lite-<SIM>-AP-2006-08-20-120000-g3.h5', '/cld2raint' ],

    [ 'SIMS/<SIM>/HDF5/cloud-<SIM>-AP-2006-08-20-120000-g3.h5',      '/cloud',           'SIMS/<SIM>/HDF5/cloud_lite-<SIM>-AP-2006-08-20-120000-g3.h5',      '/cloud'      ],
    [ 'SIMS/<SIM>/HDF5/cloud_diam-<SIM>-AP-2006-08-20-120000-g3.h5', '/cloud_diam',      'SIMS/<SIM>/HDF5/cloud_diam_lite-<SIM>-AP-2006-08-20-120000-g3.h5', '/cloud_diam' ],
    [ 'SIMS/<SIM>/HDF5/cloud_num-<SIM>-AP-2006-08-20-120000-g3.h5',  '/cloud_concen_kg', 'SIMS/<SIM>/HDF5/cloud_num_lite-<SIM>-AP-2006-08-20-120000-g3.h5',       '/cloud_num'  ],

    [ 'SIMS/<SIM>/HDF5/rain-<SIM>-AP-2006-08-20-120000-g3.h5',      '/rain',           'SIMS/<SIM>/HDF5/rain_lite-<SIM>-AP-2006-08-20-120000-g3.h5',      '/rain'      ],
    [ 'SIMS/<SIM>/HDF5/rain_diam-<SIM>-AP-2006-08-20-120000-g3.h5', '/rain_diam',      'SIMS/<SIM>/HDF5/rain_diam_lite-<SIM>-AP-2006-08-20-120000-g3.h5', '/rain_diam' ],
    [ 'SIMS/<SIM>/HDF5/rain_num-<SIM>-AP-2006-08-20-120000-g3.h5',  '/rain_concen_kg', 'SIMS/<SIM>/HDF5/rain_num_lite-<SIM>-AP-2006-08-20-120000-g3.h5',       '/rain_num'  ],

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

        # Get the dimensions of the input dataset. 2D spatial variables will be (t,y,x)
        # and 3D spatial variables will be (t,z,y,x).
        InDset = Ifile[Ivname]
        InDims = InDset.shape

        Ndims = len(InDims)
        if ((Ndims != 3) and (Ndims != 4)):
            print("WARNING: Only know how to do 2D and 3D spatial fields, skipping this variable")
            Ifile.close()
            break

        if (Ndims == 3):
            (Nt, Ny, Nx) = InDims
            Nz = 1
        elif (Ndims == 4):
            (Nt, Nz, Ny, Nx) = InDims

        print("  Reading {0:s} ({1:s})".format(Ifname, Xname))
        X = Ifile[Xname][...]

        print("  Reading {0:s} ({1:s})".format(Ifname, Yname))
        Y = Ifile[Yname][...]

        print("  Reading {0:s} ({1:s})".format(Ifname, Zname))
        Z = Ifile[Zname][...]

        print("  Reading {0:s} ({1:s})".format(Ifname, Tname))
        T = Ifile[Tname][...]
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

        Xdim = h5u.DimCoards(Xname, 1, Xreduce.shape, 'x')
        Xdim.Build(Ofile, Xreduce)
        Ydim = h5u.DimCoards(Yname, 1, Yreduce.shape, 'y')
        Ydim.Build(Ofile, Yreduce)
        Zdim = h5u.DimCoards(Zname, 1, Z.shape, 'z')
        Zdim.Build(Ofile, Z)
        Tdim = h5u.DimCoards(Tname, 1, T.shape, 't', tstring=Tstring)
        Tdim.Build(Ofile, T)

        # Create the output dataset
        if (Ndims == 3):
            OutDims = (Nt, NyReduce, NxReduce)
            OutChunks = (1, NyReduce, NxReduce)
        elif (Ndims == 4):
            OutDims = (Nt, Nz, NyReduce, NxReduce)
            OutChunks = (1, Nz, NyReduce, NxReduce)

        OutDset = h5u.DsetCoards(Ovname, Ndims, OutDims, chunks=OutChunks)
        OutDset.Create(Ofile)


        # Do one time step at a time to help reduce memory requirements. Use block_reduce to
        # to form the coarse grid by averaging the neigboring fine grid cells.
        print("  Reading {0:s} ({1:s})".format(Ifname, Ivname))

        for it in range(Nt):
            # Input field is either (t,y,x) or (t,z,y,x). Select along the x, y dimensions so that
            # their resultant sizes are evenly divisible by the x, y reduction factors.
            if (Ndims == 3):
                Var = InDset[it,Y1:Y2,X1:X2]
                VarReduce = block_reduce(Var, block_size=(RfactY, RfactX), func=np.mean)
                Ofile[Ovname][it,:,:] = VarReduce
            elif (Ndims == 4):
                Var = InDset[it,:,Y1:Y2,X1:X2]
                VarReduce = block_reduce(Var, block_size=(1, RfactY, RfactX), func=np.mean)
                Ofile[Ovname][it,:,:,:] = VarReduce

            if ((it % 10) == 0):
                print("    Processing time step: {0:d}".format(it))

        Ifile.close()
        print("    Finished: number of time steps processed: {0:d}".format(it+1))
        print("")


        print("  Writing {0:s} ({1:s})".format(Ofname, Ovname))
        if (Ndims == 3):
            OutDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
        elif (Ndims == 4):
            OutDset.AttachDims(Ofile, Tdim, Zdim, Ydim, Xdim)
        print("")

        Ofile.close()

