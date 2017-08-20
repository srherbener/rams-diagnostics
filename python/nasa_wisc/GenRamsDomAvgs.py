#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/nasa_wisc".format(os.environ['HOME']))

import h5py
import glob
import numpy as np
import ConfigRce as conf
import Hdf5Utils as h5u


Tstring = conf.SetTimeString()

SimList = [
    'RCE_1km',
    'RCE_1km_SM',
    'RCE_1km_DP',
    'RCE_1km_DM',
    ]
Nsims = len(SimList)

VarList = [
    'total_cond',
    'precip_rate',
    'vapor',
    ]
Nvars = len(VarList)

InFilePatternTemplate = "RAMS/<SIM>/RAMS/<SIM>-L-*.h5"
OutFileTemplate = "DIAGS/rams_dom_avgs_<SIM>.h5"

for isim in range(Nsims):
    Sim = SimList[isim]
    print("************************************************************************************")
    print("Generating horizontal domain averages for simulation: {0:s}".format(Sim))
    print("")

    InFilePattern = InFilePatternTemplate.replace("<SIM>", Sim)
    FileList = glob.glob(InFilePattern)
    Nt = len(FileList)

    OutFname = OutFileTemplate.replace("<SIM>", Sim)
    OutFile = h5py.File(OutFname, mode='w')

    Tstep =  0
    for InFname in sorted(FileList):
        Tstep = Tstep + 1
        print("  Reading: {0:s}, Time step: {1:d}".format(InFname, Tstep))
        InFile = h5py.File(InFname, mode='r')

        for ivar in range(Nvars):
            Var = VarList[ivar]
            print("    Variable: {0:s}".format(Var))

            if (Var == 'total_cond'):
                RTP = InFile['RTP'][...]
                RV  = InFile['RV'][...]
                VAR = RTP - RV

            elif (Var == 'precip_rate'):
                VAR = InFile['PCPRR'][...]

            elif (Var == 'vapor'):
                VAR = InFile['RV'][...]

            else:
                print("      Warning: undefined variable ({0:s}), skipping this variable".format(Var))
                continue

            # cast a 2D variable into a 3D variable with a z-dimension lenght of 1
            if (len(VAR.shape) == 2):
                VAR = VAR.reshape((1, VAR.shape[0], VAR.shape[1]))

            # if this is the first time step, build the datasets in the output file for the variable
            if (Tstep == 1):
                Nz = VAR.shape[0]
                OutFile.create_dataset(Var, (Nt, Nz), compression="gzip", compression_opts=6, shuffle=True)

            # Exclude the lateral boundary points since RAMS will sometimes set these to zero
            # (they are technically outside the domain, ie boundary values).
            OutFile[Var][Tstep-1,...] = np.mean(np.mean(VAR[:,1:-1,1:-1],axis=2),axis=1) 

        InFile.close()
        print("")

    OutFile.close()
    print("")
