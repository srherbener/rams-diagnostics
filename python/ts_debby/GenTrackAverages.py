#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
import numpy as np
import ConfigTsd as conf
import Hdf5Utils as h5u


Tstring = conf.SetTimeString()
RhoAir = conf.SetRhoAir()

SimList = [
    'TSD_SAL_DUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_NONSAL_NODUST'
    ]
Nsims = len(SimList)

# For now, need to run either the set with full horizontal data, or the set with reduced
# horizontal data. Can't run both at the same time.
TrackList = [
    # PTRACK
    [   "DIAGS/ptrack_avgs_<SIM>.h5",
        [
#        # full horizontal data
#        [ 'XsectionData/ptrack_theta_<SIM>.h5', '/theta',  'pre_sal', '/ps_theta', 'z' ],
#        [ 'XsectionData/ptrack_theta_<SIM>.h5', '/theta',  'sal',     '/s_theta',  'z' ],
#
#        [ 'XsectionData/ptrack_theta_e_<SIM>.h5', '/theta_e',  'pre_sal', '/ps_theta_e', 'z' ],
#        [ 'XsectionData/ptrack_theta_e_<SIM>.h5', '/theta_e',  'sal',     '/s_theta_e',  'z' ],
#
#        [ 'XsectionData/ptrack_tempc_<SIM>.h5',   '/tempc',  'pre_sal', '/ps_tempc',   'z' ],
#        [ 'XsectionData/ptrack_tempc_<SIM>.h5',   '/tempc',  'sal',     '/s_tempc',    'z' ],
#        [ 'XsectionData/ptrack_tempc_p_<SIM>.h5', '/tempc',  'pre_sal', '/ps_tempc_p', 'p' ],
#        [ 'XsectionData/ptrack_tempc_p_<SIM>.h5', '/tempc',  'sal',     '/s_tempc_p',  'p' ],

        # reduced (lite) horizontal data
#        [ 'XsectionData/ptrack_tempc_nv_lite_<SIM>.h5',   '/tempc',  'pre_sal', '/ps_tempc_nv',   'z' ],
#        [ 'XsectionData/ptrack_tempc_nv_lite_<SIM>.h5',   '/tempc',  'sal',     '/s_tempc_nv',    'z' ],
        [ 'XsectionData/ptrack_tempc_nv_lite_p_<SIM>.h5', '/tempc',  'pre_sal', '/ps_tempc_nv_p', 'p' ],
        [ 'XsectionData/ptrack_tempc_nv_lite_p_<SIM>.h5', '/tempc',  'sal',     '/s_tempc_nv_p',  'p' ],

#        [ 'XsectionData/ptrack_vapor_nv_lite_<SIM>.h5',   '/vapor',  'pre_sal', '/ps_vapor_nv',   'z' ],
#        [ 'XsectionData/ptrack_vapor_nv_lite_<SIM>.h5',   '/vapor',  'sal',     '/s_vapor_nv',    'z' ],
        [ 'XsectionData/ptrack_vapor_nv_lite_p_<SIM>.h5', '/vapor',  'pre_sal', '/ps_vapor_nv_p', 'p' ],
        [ 'XsectionData/ptrack_vapor_nv_lite_p_<SIM>.h5', '/vapor',  'sal',     '/s_vapor_nv_p',  'p' ],
        ]
    ],

    # STRACK
    [   "DIAGS/strack_avgs_<SIM>.h5",
        [
        [ 'XsectionData/strack_tempc_nv_lite_<SIM>.h5',   '/tempc',  'pre_sal', '/ps_tempc_nv',   'z' ],
        [ 'XsectionData/strack_tempc_nv_lite_<SIM>.h5',   '/tempc',  'sal',     '/s_tempc_nv',    'z' ],
        [ 'XsectionData/strack_tempc_nv_lite_p_<SIM>.h5', '/tempc',  'pre_sal', '/ps_tempc_nv_p', 'p' ],
        [ 'XsectionData/strack_tempc_nv_lite_p_<SIM>.h5', '/tempc',  'sal',     '/s_tempc_nv_p',  'p' ],

        [ 'XsectionData/strack_theta_nv_lite_<SIM>.h5',   '/theta',  'pre_sal', '/ps_theta_nv',   'z' ],
        [ 'XsectionData/strack_theta_nv_lite_<SIM>.h5',   '/theta',  'sal',     '/s_theta_nv',    'z' ],
        [ 'XsectionData/strack_theta_nv_lite_p_<SIM>.h5', '/theta',  'pre_sal', '/ps_theta_nv_p', 'p' ],
        [ 'XsectionData/strack_theta_nv_lite_p_<SIM>.h5', '/theta',  'sal',     '/s_theta_nv_p',  'p' ],

        [ 'XsectionData/strack_theta_e_nv_lite_<SIM>.h5',   '/theta_e',  'pre_sal', '/ps_theta_e_nv',   'z' ],
        [ 'XsectionData/strack_theta_e_nv_lite_<SIM>.h5',   '/theta_e',  'sal',     '/s_theta_e_nv',    'z' ],
        [ 'XsectionData/strack_theta_e_nv_lite_p_<SIM>.h5', '/theta_e',  'pre_sal', '/ps_theta_e_nv_p', 'p' ],
        [ 'XsectionData/strack_theta_e_nv_lite_p_<SIM>.h5', '/theta_e',  'sal',     '/s_theta_e_nv_p',  'p' ],

        [ 'XsectionData/strack_vapor_nv_lite_<SIM>.h5',   '/vapor',  'pre_sal', '/ps_vapor_nv',   'z' ],
        [ 'XsectionData/strack_vapor_nv_lite_<SIM>.h5',   '/vapor',  'sal',     '/s_vapor_nv',    'z' ],
        [ 'XsectionData/strack_vapor_nv_lite_p_<SIM>.h5', '/vapor',  'pre_sal', '/ps_vapor_nv_p', 'p' ],
        [ 'XsectionData/strack_vapor_nv_lite_p_<SIM>.h5', '/vapor',  'sal',     '/s_vapor_nv_p',  'p' ],

        ]
    ],
    ]
Nsets = len(TrackList)

PsapStart = 10  # PSAP, sim time in hours
PsapEnd   = 30

SapStart = 40  # SAP, sim time in hours
SapEnd   = 60

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'  # vertical coords: height
Pname = '/p_coords'  # vertical coords: pressure
Tname = '/t_coords'

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Creating track averages for simulation: {0:s}".format(Sim))
    print('')

    for iset in range(Nsets):
        OfileTemplate = TrackList[iset][0]
        VarList = TrackList[iset][1]

        Nvars = len(VarList)
        
        NoX = True
        NoY = True
        NoZ = True
        NoP = True
        NoT = True

        OutFname = OfileTemplate.replace("<SIM>", Sim)
        Ofile = h5py.File(OutFname, mode='w')

        for ivar in range(Nvars):
            InFname    = VarList[ivar][0].replace("<SIM>", Sim)
            InVname    = VarList[ivar][1]
            AvgPeriod  = VarList[ivar][2]
            OutVname   = VarList[ivar][3]
            VcoordType = VarList[ivar][4]
    
            # Read input data. If this is the first set, read in the coordinates
            # and build the dimensions.
            print("  Reading {0:s} ({1:s})".format(InFname, InVname))
            print("    Averaging time period: {0:s}".format(AvgPeriod))
            print("    Vertical coordinate type: {0:s}".format(VcoordType))
            Ifile = h5py.File(InFname, mode='r')
            Var = np.squeeze(Ifile[InVname][...]) # Y-dim is size 1, reduce Var to (t,z,x)
    
            if (NoX):
                print("    Reading {0:s} ({1:s})".format(InFname, Xname))
                X = Ifile[Xname][...]
                Nx = len(X)
                Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
                Xdim.Build(Ofile, X)
                NoX = False
    
            if (NoY):
                print("    Reading {0:s} ({1:s})".format(InFname, Yname))
                Y = Ifile[Yname][...]
                Ny = len(Y)
                Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
                Ydim.Build(Ofile, Y)
                NoY = False
    
            if (VcoordType == 'z'):
                if (NoZ):
                    print("    Reading {0:s} ({1:s})".format(InFname, Zname))
                    Z = Ifile[Zname][...]
                    Nz = len(Z)
                    Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
                    Zdim.Build(Ofile, Z)
                    NoZ = False
            elif (VcoordType == 'p'):
                if (NoP):
                    print("    Reading {0:s} ({1:s})".format(InFname, Zname))
                    P = Ifile[Zname][...]
                    Np = len(P)
                    Pdim = h5u.DimCoards(Pname, 1, [ Np ], 'p')
                    Pdim.Build(Ofile, P)
                    NoP = False
    
            if (NoT):
                print("    Reading {0:s} ({1:s})".format(InFname, Tname))
                T = Ifile[Tname][...]
                SimT = T / 3600.0 - 42
                Nt = len(T)
                Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
                Tdim.Build(Ofile, T)
    
                # Find indices of PSAP and SAP time periods
                SimT = T / 3600.0 - 42
    
                Select = np.where((SimT >= PsapStart) * (SimT <= PsapEnd))
                PS_T1 = Select[0][0]
                PS_T2 = Select[0][-1]
    
                Select = np.where((SimT >= SapStart) * (SimT <= SapEnd))
                S_T1 = Select[0][0]
                S_T2 = Select[0][-1]
    
                NoT = False
    
            Ifile.close()
    
            # Want to reduce time dimensions by taking mean value across
            # the time period.
            if (AvgPeriod == 'pre_sal'):
                AvgVar = np.squeeze(np.mean(Var[PS_T1:PS_T2,:,:], axis=0))
            elif (AvgPeriod == 'sal'):
                AvgVar = np.squeeze(np.mean(Var[S_T1:S_T2,:,:], axis=0))
    
            # Write out PSAP and SAP averged vars
            print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
            if (VcoordType == 'z'):
                VarDset = h5u.DsetCoards(OutVname, 2, [ Nz, Nx ])
                VarDset.Build(Ofile, AvgVar, Zdim, Xdim)
            elif (VcoordType == 'p'):
                VarDset = h5u.DsetCoards(OutVname, 2, [ Np, Nx ])
                VarDset.Build(Ofile, AvgVar, Pdim, Xdim)
            print('')
    
        Ofile.close()
