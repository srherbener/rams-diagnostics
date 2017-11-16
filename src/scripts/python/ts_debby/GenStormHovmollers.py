#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
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

HovmollerList = [
    #[ 'DIAGS/hist_meas_az_vapor_<SIM>.h5', '/all_core_vapor_ts',  '/all_core_vapor'  ],
    #[ 'DIAGS/hist_meas_az_vapor_<SIM>.h5', '/all_rb_vapor_ts',    '/all_rb_vapor'    ],
    #[ 'DIAGS/hist_meas_az_vapor_<SIM>.h5', '/lead_core_vapor_ts', '/lead_core_vapor' ],
    #[ 'DIAGS/hist_meas_az_vapor_<SIM>.h5', '/lead_rb_vapor_ts',   '/lead_rb_vapor'   ],
    #[ 'DIAGS/hist_meas_az_vapor_<SIM>.h5', '/lead_env_vapor_ts',  '/lead_env_vapor'  ],

    #[ 'DIAGS/hist_meas_az_theta_<SIM>.h5', '/all_core_theta_ts',  '/all_core_theta'  ],
    #[ 'DIAGS/hist_meas_az_theta_<SIM>.h5', '/all_rb_theta_ts',    '/all_rb_theta'    ],
    #[ 'DIAGS/hist_meas_az_theta_<SIM>.h5', '/lead_core_theta_ts', '/lead_core_theta' ],
    #[ 'DIAGS/hist_meas_az_theta_<SIM>.h5', '/lead_rb_theta_ts',   '/lead_rb_theta'   ],
    #[ 'DIAGS/hist_meas_az_theta_<SIM>.h5', '/lead_env_theta_ts',  '/lead_env_theta'  ],

    #[ 'DIAGS/hist_meas_az_theta_e_<SIM>.h5', '/all_core_theta_e_ts',  '/all_core_theta_e'  ],
    #[ 'DIAGS/hist_meas_az_theta_e_<SIM>.h5', '/all_rb_theta_e_ts',    '/all_rb_theta_e'    ],
    #[ 'DIAGS/hist_meas_az_theta_e_<SIM>.h5', '/lead_core_theta_e_ts', '/lead_core_theta_e' ],
    #[ 'DIAGS/hist_meas_az_theta_e_<SIM>.h5', '/lead_rb_theta_e_ts',   '/lead_rb_theta_e'   ],
    #[ 'DIAGS/hist_meas_az_theta_e_<SIM>.h5', '/lead_env_theta_e_ts',  '/lead_env_theta_e'  ],

    #[ 'DIAGS/hist_meas_az_cloud_<SIM>.h5', '/all_core_cloud_mass_ts',  '/all_core_cloud'  ],
    #[ 'DIAGS/hist_meas_az_cloud_<SIM>.h5', '/all_rb_cloud_mass_ts',    '/all_rb_cloud'    ],
    #[ 'DIAGS/hist_meas_az_cloud_<SIM>.h5', '/lead_core_cloud_mass_ts', '/lead_core_cloud' ],
    #[ 'DIAGS/hist_meas_az_cloud_<SIM>.h5', '/lead_rb_cloud_mass_ts',   '/lead_rb_cloud'   ],
    #[ 'DIAGS/hist_meas_az_cloud_<SIM>.h5', '/lead_env_cloud_mass_ts',  '/lead_env_cloud'  ],

    #[ 'DIAGS/hist_meas_az_rain_<SIM>.h5', '/all_core_rain_mass_ts',  '/all_core_rain'  ],
    #[ 'DIAGS/hist_meas_az_rain_<SIM>.h5', '/all_rb_rain_mass_ts',    '/all_rb_rain'    ],
    #[ 'DIAGS/hist_meas_az_rain_<SIM>.h5', '/lead_core_rain_mass_ts', '/lead_core_rain' ],
    #[ 'DIAGS/hist_meas_az_rain_<SIM>.h5', '/lead_rb_rain_mass_ts',   '/lead_rb_rain'   ],
    #[ 'DIAGS/hist_meas_az_rain_<SIM>.h5', '/lead_env_rain_mass_ts',  '/lead_env_rain'  ],

    #[ 'DIAGS/hist_meas_az_cloud_cond_<SIM>.h5', '/all_core_cloud_cond_ts',  '/all_core_cloud_cond'  ],
    #[ 'DIAGS/hist_meas_az_cloud_cond_<SIM>.h5', '/all_rb_cloud_cond_ts',    '/all_rb_cloud_cond'    ],
    #[ 'DIAGS/hist_meas_az_cloud_cond_<SIM>.h5', '/lead_core_cloud_cond_ts', '/lead_core_cloud_cond' ],
    #[ 'DIAGS/hist_meas_az_cloud_cond_<SIM>.h5', '/lead_rb_cloud_cond_ts',   '/lead_rb_cloud_cond'   ],
    #[ 'DIAGS/hist_meas_az_cloud_cond_<SIM>.h5', '/lead_env_cloud_cond_ts',  '/lead_env_cloud_cond'  ],

    #[ 'DIAGS/hist_meas_az_cloud_evap_<SIM>.h5', '/all_core_cloud_evap_ts',  '/all_core_cloud_evap'  ],
    #[ 'DIAGS/hist_meas_az_cloud_evap_<SIM>.h5', '/all_rb_cloud_evap_ts',    '/all_rb_cloud_evap'    ],
    #[ 'DIAGS/hist_meas_az_cloud_evap_<SIM>.h5', '/lead_core_cloud_evap_ts', '/lead_core_cloud_evap' ],
    #[ 'DIAGS/hist_meas_az_cloud_evap_<SIM>.h5', '/lead_rb_cloud_evap_ts',   '/lead_rb_cloud_evap'   ],
    #[ 'DIAGS/hist_meas_az_cloud_evap_<SIM>.h5', '/lead_env_cloud_evap_ts',  '/lead_env_cloud_evap'  ],

    #[ 'DIAGS/hist_meas_az_rain_cond_<SIM>.h5', '/all_core_rain_cond_ts',  '/all_core_rain_cond'  ],
    #[ 'DIAGS/hist_meas_az_rain_cond_<SIM>.h5', '/all_rb_rain_cond_ts',    '/all_rb_rain_cond'    ],
    #[ 'DIAGS/hist_meas_az_rain_cond_<SIM>.h5', '/lead_core_rain_cond_ts', '/lead_core_rain_cond' ],
    #[ 'DIAGS/hist_meas_az_rain_cond_<SIM>.h5', '/lead_rb_rain_cond_ts',   '/lead_rb_rain_cond'   ],
    #[ 'DIAGS/hist_meas_az_rain_cond_<SIM>.h5', '/lead_env_rain_cond_ts',  '/lead_env_rain_cond'  ],

    #[ 'DIAGS/hist_meas_az_rain_evap_<SIM>.h5', '/all_core_rain_evap_ts',  '/all_core_rain_evap'  ],
    #[ 'DIAGS/hist_meas_az_rain_evap_<SIM>.h5', '/all_rb_rain_evap_ts',    '/all_rb_rain_evap'    ],
    #[ 'DIAGS/hist_meas_az_rain_evap_<SIM>.h5', '/lead_core_rain_evap_ts', '/lead_core_rain_evap' ],
    #[ 'DIAGS/hist_meas_az_rain_evap_<SIM>.h5', '/lead_rb_rain_evap_ts',   '/lead_rb_rain_evap'   ],
    #[ 'DIAGS/hist_meas_az_rain_evap_<SIM>.h5', '/lead_env_rain_evap_ts',  '/lead_env_rain_evap'  ],

    #[ 'DIAGS/hist_meas_az_ice_dep_<SIM>.h5', '/all_core_ice_dep_ts', '/all_core_ice_dep' ],
    #[ 'DIAGS/hist_meas_az_ice_dep_<SIM>.h5', '/all_rb_ice_dep_ts',   '/all_rb_ice_dep'   ],

    #[ 'DIAGS/hist_meas_az_rain2ice_<SIM>.h5', '/all_core_rain2ice_ts', '/all_core_rain2ice' ],
    #[ 'DIAGS/hist_meas_az_rain2ice_<SIM>.h5', '/all_rb_rain2ice_ts',   '/all_rb_rain2ice'   ],

    #[ 'DIAGS/hist_meas_az_w_<SIM>.h5', '/all_core_updraft_ts', '/all_core_updraft' ],
    #[ 'DIAGS/hist_meas_az_w_<SIM>.h5', '/all_rb_updraft_ts',   '/all_rb_updraft'   ],
    #[ 'DIAGS/hist_meas_az_w_<SIM>.h5', '/all_core_dndraft_ts', '/all_core_dndraft' ],
    #[ 'DIAGS/hist_meas_az_w_<SIM>.h5', '/all_rb_dndraft_ts',   '/all_rb_dndraft'   ],

    [ 'DIAGS/hist_meas_az_p_vapor_lite_<SIM>.h5', '/all_core_vapor_ts',      '/all_p_core_vapor_lite'      ],
    [ 'DIAGS/hist_meas_az_p_vapor_lite_<SIM>.h5', '/all_rb_vapor_ts',        '/all_p_rb_vapor_lite'        ],
    [ 'DIAGS/hist_meas_az_p_vapor_lite_<SIM>.h5', '/all_nv_core_vapor_ts',   '/all_p_nv_core_vapor_lite'   ],
    [ 'DIAGS/hist_meas_az_p_vapor_lite_<SIM>.h5', '/all_nv_rb_vapor_ts',     '/all_p_nv_rb_vapor_lite'     ],

    [ 'DIAGS/hist_meas_az_p_theta_lite_<SIM>.h5', '/all_core_theta_ts',      '/all_p_core_theta_lite'      ],
    [ 'DIAGS/hist_meas_az_p_theta_lite_<SIM>.h5', '/all_rb_theta_ts',        '/all_p_rb_theta_lite'        ],
    [ 'DIAGS/hist_meas_az_p_theta_lite_<SIM>.h5', '/all_nv_core_theta_ts',   '/all_p_nv_core_theta_lite'   ],
    [ 'DIAGS/hist_meas_az_p_theta_lite_<SIM>.h5', '/all_nv_rb_theta_ts',     '/all_p_nv_rb_theta_lite'     ],

    [ 'DIAGS/hist_meas_az_p_theta_e_lite_<SIM>.h5', '/all_core_theta_e_ts',      '/all_p_core_theta_e_lite'      ],
    [ 'DIAGS/hist_meas_az_p_theta_e_lite_<SIM>.h5', '/all_rb_theta_e_ts',        '/all_p_rb_theta_e_lite'        ],
    [ 'DIAGS/hist_meas_az_p_theta_e_lite_<SIM>.h5', '/all_nv_core_theta_e_ts',   '/all_p_nv_core_theta_e_lite'   ],
    [ 'DIAGS/hist_meas_az_p_theta_e_lite_<SIM>.h5', '/all_nv_rb_theta_e_ts',     '/all_p_nv_rb_theta_e_lite'     ],

    [ 'DIAGS/hist_meas_az_p_relhum_lite_<SIM>.h5', '/all_core_relhum_ts',      '/all_p_core_relhum_lite'      ],
    [ 'DIAGS/hist_meas_az_p_relhum_lite_<SIM>.h5', '/all_rb_relhum_ts',        '/all_p_rb_relhum_lite'        ],
    [ 'DIAGS/hist_meas_az_p_relhum_lite_<SIM>.h5', '/all_nv_core_relhum_ts',   '/all_p_nv_core_relhum_lite'   ],
    [ 'DIAGS/hist_meas_az_p_relhum_lite_<SIM>.h5', '/all_nv_rb_relhum_ts',     '/all_p_nv_rb_relhum_lite'     ],

    [ 'DIAGS/hist_meas_az_p_entropy_lite_<SIM>.h5', '/all_core_entropy_ts',      '/all_p_core_entropy_lite'      ],
    [ 'DIAGS/hist_meas_az_p_entropy_lite_<SIM>.h5', '/all_rb_entropy_ts',        '/all_p_rb_entropy_lite'        ],
    [ 'DIAGS/hist_meas_az_p_entropy_lite_<SIM>.h5', '/all_nv_core_entropy_ts',   '/all_p_nv_core_entropy_lite'   ],
    [ 'DIAGS/hist_meas_az_p_entropy_lite_<SIM>.h5', '/all_nv_rb_entropy_ts',     '/all_p_nv_rb_entropy_lite'     ],

    [ 'DIAGS/hist_meas_az_p_ice_sub_lite_<SIM>.h5', '/all_core_ice_sub_ts',      '/all_p_core_ice_sub_lite'      ],
    [ 'DIAGS/hist_meas_az_p_ice_sub_lite_<SIM>.h5', '/all_rb_ice_sub_ts',        '/all_p_rb_ice_sub_lite'        ],
    [ 'DIAGS/hist_meas_az_p_ice_sub_lite_<SIM>.h5', '/all_nv_core_ice_sub_ts',   '/all_p_nv_core_ice_sub_lite'   ],
    [ 'DIAGS/hist_meas_az_p_ice_sub_lite_<SIM>.h5', '/all_nv_rb_ice_sub_ts',     '/all_p_nv_rb_ice_sub_lite'     ],

    [ 'DIAGS/hist_meas_az_p_ice_dep_lite_<SIM>.h5', '/all_core_ice_dep_ts',      '/all_p_core_ice_dep_lite'      ],
    [ 'DIAGS/hist_meas_az_p_ice_dep_lite_<SIM>.h5', '/all_rb_ice_dep_ts',        '/all_p_rb_ice_dep_lite'        ],
    [ 'DIAGS/hist_meas_az_p_ice_dep_lite_<SIM>.h5', '/all_nv_core_ice_dep_ts',   '/all_p_nv_core_ice_dep_lite'   ],
    [ 'DIAGS/hist_meas_az_p_ice_dep_lite_<SIM>.h5', '/all_nv_rb_ice_dep_ts',     '/all_p_nv_rb_ice_dep_lite'     ],

    [ 'DIAGS/hist_meas_az_p_liq_evap_lite_<SIM>.h5', '/all_core_liq_evap_ts',      '/all_p_core_liq_evap_lite'      ],
    [ 'DIAGS/hist_meas_az_p_liq_evap_lite_<SIM>.h5', '/all_rb_liq_evap_ts',        '/all_p_rb_liq_evap_lite'        ],
    [ 'DIAGS/hist_meas_az_p_liq_evap_lite_<SIM>.h5', '/all_nv_core_liq_evap_ts',   '/all_p_nv_core_liq_evap_lite'   ],
    [ 'DIAGS/hist_meas_az_p_liq_evap_lite_<SIM>.h5', '/all_nv_rb_liq_evap_ts',     '/all_p_nv_rb_liq_evap_lite'     ],

    [ 'DIAGS/hist_meas_az_p_liq_cond_lite_<SIM>.h5', '/all_core_liq_cond_ts',      '/all_p_core_liq_cond_lite'      ],
    [ 'DIAGS/hist_meas_az_p_liq_cond_lite_<SIM>.h5', '/all_rb_liq_cond_ts',        '/all_p_rb_liq_cond_lite'        ],
    [ 'DIAGS/hist_meas_az_p_liq_cond_lite_<SIM>.h5', '/all_nv_core_liq_cond_ts',   '/all_p_nv_core_liq_cond_lite'   ],
    [ 'DIAGS/hist_meas_az_p_liq_cond_lite_<SIM>.h5', '/all_nv_rb_liq_cond_ts',     '/all_p_nv_rb_liq_cond_lite'     ],

    [ 'DIAGS/hist_meas_az_p_cld2raint_lite_<SIM>.h5', '/all_core_cld2raint_ts',      '/all_p_core_cld2raint_lite'      ],
    [ 'DIAGS/hist_meas_az_p_cld2raint_lite_<SIM>.h5', '/all_rb_cld2raint_ts',        '/all_p_rb_cld2raint_lite'        ],
    [ 'DIAGS/hist_meas_az_p_cld2raint_lite_<SIM>.h5', '/all_nv_core_cld2raint_ts',   '/all_p_nv_core_cld2raint_lite'   ],
    [ 'DIAGS/hist_meas_az_p_cld2raint_lite_<SIM>.h5', '/all_nv_rb_cld2raint_ts',     '/all_p_nv_rb_cld2raint_lite'     ],

    [ 'DIAGS/hist_meas_az_p_ice2raint_lite_<SIM>.h5', '/all_core_ice2raint_ts',      '/all_p_core_ice2raint_lite'      ],
    [ 'DIAGS/hist_meas_az_p_ice2raint_lite_<SIM>.h5', '/all_rb_ice2raint_ts',        '/all_p_rb_ice2raint_lite'        ],
    [ 'DIAGS/hist_meas_az_p_ice2raint_lite_<SIM>.h5', '/all_nv_core_ice2raint_ts',   '/all_p_nv_core_ice2raint_lite'   ],
    [ 'DIAGS/hist_meas_az_p_ice2raint_lite_<SIM>.h5', '/all_nv_rb_ice2raint_ts',     '/all_p_nv_rb_ice2raint_lite'     ],

    [ 'DIAGS/hist_meas_az_p_melticet_lite_<SIM>.h5', '/all_core_melticet_ts',      '/all_p_core_melticet_lite'      ],
    [ 'DIAGS/hist_meas_az_p_melticet_lite_<SIM>.h5', '/all_rb_melticet_ts',        '/all_p_rb_melticet_lite'        ],
    [ 'DIAGS/hist_meas_az_p_melticet_lite_<SIM>.h5', '/all_nv_core_melticet_ts',   '/all_p_nv_core_melticet_lite'   ],
    [ 'DIAGS/hist_meas_az_p_melticet_lite_<SIM>.h5', '/all_nv_rb_melticet_ts',     '/all_p_nv_rb_melticet_lite'     ],

    [ 'DIAGS/hist_meas_az_p_rain2icet_lite_<SIM>.h5', '/all_core_rain2icet_ts',      '/all_p_core_rain2icet_lite'      ],
    [ 'DIAGS/hist_meas_az_p_rain2icet_lite_<SIM>.h5', '/all_rb_rain2icet_ts',        '/all_p_rb_rain2icet_lite'        ],
    [ 'DIAGS/hist_meas_az_p_rain2icet_lite_<SIM>.h5', '/all_nv_core_rain2icet_ts',   '/all_p_nv_core_rain2icet_lite'   ],
    [ 'DIAGS/hist_meas_az_p_rain2icet_lite_<SIM>.h5', '/all_nv_rb_rain2icet_ts',     '/all_p_nv_rb_rain2icet_lite'     ],

    [ 'DIAGS/hist_meas_az_p_rimecldt_lite_<SIM>.h5', '/all_core_rimecldt_ts',      '/all_p_core_rimecldt_lite'      ],
    [ 'DIAGS/hist_meas_az_p_rimecldt_lite_<SIM>.h5', '/all_rb_rimecldt_ts',        '/all_p_rb_rimecldt_lite'        ],
    [ 'DIAGS/hist_meas_az_p_rimecldt_lite_<SIM>.h5', '/all_nv_core_rimecldt_ts',   '/all_p_nv_core_rimecldt_lite'   ],
    [ 'DIAGS/hist_meas_az_p_rimecldt_lite_<SIM>.h5', '/all_nv_rb_rimecldt_ts',     '/all_p_nv_rb_rimecldt_lite'     ],

    ]
Nsets = len(HovmollerList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Creating Hovmoller data for simulation: {0:s}".format(Sim))
    print('')

    OutFname = "DIAGS/storm_hovs_{0:s}.h5".format(Sim)
    Ofile = h5py.File(OutFname, mode='w')

    for iset in range(Nsets):
        InFname  = HovmollerList[iset][0].replace("<SIM>", Sim)
        InVname  = HovmollerList[iset][1]
        OutVname = HovmollerList[iset][2]
    
        # Read input data. If this is the first set, read in the coordinates
        # and build the dimensions.
        print("  Reading {0:s} ({1:s})".format(InFname, InVname))
        Ifile = h5py.File(InFname, mode='r')
        Var = Ifile[InVname][...]
        Var = Var.transpose() # transpose since 2D contour plot routines expect
                              # data to be in the form (y,x), or for these
                              # hovmollers (z,t)
        if (iset == 0):
            X = Ifile[Xname][...]
            Y = Ifile[Yname][...]
            Z = Ifile[Zname][...]
            T = Ifile[Tname][...]

            Nx = len(X)
            Ny = len(Y)
            Nz = len(Z)
            Nt = len(T)

            Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
            Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
            Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
            Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)

            Xdim.Build(Ofile, X)
            Ydim.Build(Ofile, Y)
            Zdim.Build(Ofile, Z)
            Tdim.Build(Ofile, T)

        Ifile.close()

        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
        VarDset = h5u.DsetCoards(OutVname, 2, [ Nz, Nt ])
        VarDset.Build(Ofile, Var, Zdim, Tdim)
        print('')

    Ofile.close()
