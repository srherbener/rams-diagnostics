#!/usr/bin/env python3
#
# Script to generate moist thermodynamic variables from temperature, pressure
# and vapor mixing ratio variables.
#

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
import numpy as np
import ConfigTsd as conf
import Hdf5Utils as h5u
import ThermoUtils as tu


Tstring = conf.SetTimeString()

SimList = [
    'TSD_SAL_DUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_NONSAL_NODUST'
    ]
Nsims = len(SimList)

# For the last argument (vertical coordinate type), use "AP" for pressure, "AS" for sigma-z,
# and "AC" for z.
VarList = [
    # sigma-z coords
    [ "press_lite", "/press", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "density_lite", "/density", "AS" ],
#    [ "press_lite", "/press", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "relhum_lite", "/relhum", "AS" ],
#    [ "press_lite", "/press", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "theta_lite", "/theta", "AS" ],
#    [ "press_lite", "/press", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "theta_e_lite", "/theta_e", "AS" ],
#    [ "press_lite", "/press", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "theta_es_lite", "/theta_es", "AS" ],
#    [ "press_lite", "/press", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "entropy_lite", "/entropy", "AS" ],
#    [ "press_lite", "/press", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "entropy_s_lite", "/entropy_s", "AS" ],
#
    [ "press_lite", "/press", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "density_nv_lite", "/density", "AS" ],
#    [ "press_lite", "/press", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "relhum_nv_lite", "/relhum", "AS" ],
#    [ "press_lite", "/press", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "theta_nv_lite", "/theta", "AS" ],
#    [ "press_lite", "/press", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "theta_e_nv_lite", "/theta_e", "AS" ],
#    [ "press_lite", "/press", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "theta_es_nv_lite", "/theta_es", "AS" ],
#    [ "press_lite", "/press", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "entropy_nv_lite", "/entropy", "AS" ],
#    [ "press_lite", "/press", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "entropy_s_nv_lite", "/entropy_s", "AS" ],

    # pressure coords
    # use the tempc file and "/z_coords" for the pressure values
    [ "tempc_lite", "/z_coords", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "density_lite", "/density", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "relhum_lite", "/relhum", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "theta_lite", "/theta", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "theta_e_lite", "/theta_e", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "theta_es_lite", "/theta_es", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "entropy_lite", "/entropy", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_lite", "/tempc", "vapor_lite", "/vapor", "entropy_s_lite", "/entropy_s", "AP" ],
#
    [ "tempc_lite", "/z_coords", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "density_nv_lite", "/density", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "relhum_nv_lite", "/relhum", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "theta_nv_lite", "/theta", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "theta_e_nv_lite", "/theta_e", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "theta_es_nv_lite", "/theta_es", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "entropy_nv_lite", "/entropy", "AP" ],
#    [ "tempc_lite", "/z_coords", "tempc_nv_lite", "/tempc_basic", "vapor_nv_lite", "/vapor_basic", "entropy_s_nv_lite", "/entropy_s", "AP" ],
    ]
Nvars = len(VarList)

FileTemplate = "HDF5/<SIM>/HDF5/<VAR>-<SIM>-<VCOORD>-2006-08-20-120000-g3.h5"

Xname = "/x_coords"
Yname = "/y_coords"
Zname = "/z_coords"
Tname = "/t_coords"

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Generating variables for simulations: {0:s}".format(Sim))
    print("")

    for ivar in range(Nvars):
        Pfvname = VarList[ivar][0]
        Pvname  = VarList[ivar][1]
        Tfvname = VarList[ivar][2]
        Tvname  = VarList[ivar][3]
        Vfvname = VarList[ivar][4]
        Vvname  = VarList[ivar][5]
        Ofvname = VarList[ivar][6]
        Ovname  = VarList[ivar][7]
        Vcoord  = VarList[ivar][8]

        Pfname = FileTemplate.replace("<SIM>", Sim).replace("<VCOORD>", Vcoord).replace("<VAR>", Pfvname)
        Tfname = FileTemplate.replace("<SIM>", Sim).replace("<VCOORD>", Vcoord).replace("<VAR>", Tfvname)
        Vfname = FileTemplate.replace("<SIM>", Sim).replace("<VCOORD>", Vcoord).replace("<VAR>", Vfvname)
        Ofname = FileTemplate.replace("<SIM>", Sim).replace("<VCOORD>", Vcoord).replace("<VAR>", Ofvname)

        # Open the files, note in the case of pressure coordinates Pfile and Tfile
        # will both be opening teh temperature file.
        Pfile = h5py.File(Pfname, mode='r')
        Tfile = h5py.File(Tfname, mode='r')
        Vfile = h5py.File(Vfname, mode='r')

        Ofile = h5py.File(Ofname, mode='w')

        # Build the coordinates
        print("  Reading: {0:s} ({1:s})".format(Pfname, Xname))
        print("  Writing: {0:s} ({1:s})".format(Ofname, Xname))
        X = Pfile[Xname][...]
        Nx = len(X)
        Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
        Xdim.Build(Ofile, X)
    
        print("  Reading: {0:s} ({1:s})".format(Pfname, Yname))
        print("  Writing: {0:s} ({1:s})".format(Ofname, Yname))
        Y = Pfile[Yname][...]
        Ny = len(Y)
        Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
        Ydim.Build(Ofile, Y)
    
        print("  Reading: {0:s} ({1:s})".format(Pfname, Zname))
        print("  Writing: {0:s} ({1:s})".format(Ofname, Zname))
        Z = Pfile[Zname][...]
        Nz = len(Z)
        if (Vcoord == "AP"):
            Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'p')
        else:
            Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
        Zdim.Build(Ofile, Z)
    
        print("  Reading: {0:s} ({1:s})".format(Pfname, Tname))
        print("  Writing: {0:s} ({1:s})".format(Ofname, Tname))
        T = Pfile[Tname][...]
        Nt = len(T)
        Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
        Tdim.Build(Ofile, T)
    
        print("")
    
        # Create the output dataset
        OutDset = h5u.DsetCoards(Ovname, 4, [ Nt, Nz, Ny, Nx ], chunks=( 1, Nz, Ny, Nx ))
        OutDset.Create(Ofile)
    
        # In order to keep memory requirements down, process one time step at a time
        print("  Reading: {0:s} ({1:s})".format(Pfname, Pvname))
        print("  Reading: {0:s} ({1:s})".format(Tfname, Tvname))
        print("  Reading: {0:s} ({1:s})".format(Vfname, Vvname))

        # If we have pressure coordinates, then read in the z_coords and expand
        # to a 3D field. This can be used for every time step since the pressure
        # levels remain constant.
        if (Vcoord == "AP"):
            P = Pfile[Pvname][...] * 100 # covert mb (hPa) to Pa
            P = P.reshape([ Nz, 1, 1 ])
            P = np.tile(P, [ 1, Ny, Nx ])

        for it in range(Nt):
            # Input fields are dimensioned as (t,z,y,x). So slice on the first dimension.

            # If we do not have pressure coordinates, then read in pressure values
            # from the input file. If we do have pressure coordinates, P has already
            # been constructed.
            if (Vcoord != "AP"):
                P = Pfile[Pvname][it,...] * 100 # convert mb (hPa) to Pa

            T = Tfile[Tvname][it,...] + 273.15 # convert deg C to K
            V = Vfile[Vvname][it,...] * 1e-3   # convert g/kg to kg/kg

            if (Ovname == "/density"):
                VAR = tu.DryAirDensity(T, P)
            elif (Ovname == "/relhum"):
                VAR = tu.RelHum(T, P, V)
            elif (Ovname == "/theta"):
                VAR = tu.PotTemp(T, P)
            elif (Ovname == "/theta_e"):
                VAR = tu.EquivPotTemp(T, P, V)
            elif (Ovname == "/theta_es"):
                VAR = tu.SatEquivPotTemp(T, P)
            elif (Ovname == "/entropy"):
                VAR = tu.Entropy(T, P, V)
            elif (Ovname == "/entropy_s"):
                VAR = tu.SatEntropy(T, P)
            else:
                print("WARNING: unrecognized output variable ({0:s}), filling output with zeros".format(Ovname))
                VAR = np.zeros([ Nz, Ny, Nx ])

            # Write this timestep into the output file.
            Ofile[Ovname][it,:,:,:] = VAR

            if ((it % 10) == 0):
               print("    Processing time step: {0:d}".format(it))
    
        print("")


        # Attach the dimensions to output fields
        print("  Writing: {0:s} ({1:s})".format(Ofname, Ovname))
        OutDset.AttachDims(Ofile, Tdim, Zdim, Ydim, Xdim)
        print("")

        Pfile.close()
        Tfile.close()
        Vfile.close()

        Ofile.close()
