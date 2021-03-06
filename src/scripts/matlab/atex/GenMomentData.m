function [ ] = GenMomentData(ConfigFile)
% GenMomentData generate w moment and w flux data

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;

    %Tfirst = 0; % calc moment for each time step, the take average over time steps.
    Tfirst = 1; % take average over time steps first, then calculate moments.

    % These lists are organized as 2D cell arrays. Each row is one complete spec for one variable.
    % Syntax for rows:
    %    { 'input file prefix' 'input file var name' 'input var term number' 'input var order number' 'output var name' }
    VarList = {
      % means
      { 'w_M3'                  'w-w-w'     1 1 'w'              }
      { 'w_M3_ud0p10'           'w-w-w'     1 1 'w_ud0p10'       }
      { 'w_M3_up0p10'           'w-w-w'     1 1 'w_up0p10'       }
      { 'w_M3_dn0p10'           'w-w-w'     1 1 'w_dn0p10'       }
      { 'w_M3_stall'            'w-w-w'     1 1 'w_stall'        }
      { 'w_M3_up0p10_stall'     'w-w-w'     1 1 'w_up0p10_stall' }
      { 'w_M3_dn0p10_stall'     'w-w-w'     1 1 'w_dn0p10_stall' }
      { 'w_M3_all_cld'          'w-w-w'     1 1 'w_all_cld'      }
      { 'w_M3_col_up0p10'       'w-w-w'     1 1 'w_col_up0p10'   }
      { 'w_M3_col_dn0p10'       'w-w-w'     1 1 'w_col_dn0p10'   }
      { 'w_M3_col_ud0p10'       'w-w-w'     1 1 'w_col_ud0p10'   }
      { 'w_M3_up0p10_all_cld'   'w-w-w'     1 1 'w_up0p10_all_cld' }
      { 'w_M3_dn0p10_all_cld'   'w-w-w'     1 1 'w_dn0p10_all_cld' }
      { 'w_M3_ud0p10_all_cld'   'w-w-w'     1 1 'w_ud0p10_all_cld' }
      { 'w_M3_up0p01_all_cld'   'w-w-w'     1 1 'w_up0p01_all_cld' }
      { 'w_M3_dn0p01_all_cld'   'w-w-w'     1 1 'w_dn0p01_all_cld' }
      { 'w_M3_ud0p01_all_cld'   'w-w-w'     1 1 'w_ud0p01_all_cld' }

      { 'w_theta_flux'              'w-theta'   2 1 'theta'              }
      { 'w_theta_flux_ud0p10'       'w-theta'   2 1 'theta_ud0p10'       }
      { 'w_theta_flux_up0p10'       'w-theta'   2 1 'theta_up0p10'       }
      { 'w_theta_flux_dn0p10'       'w-theta'   2 1 'theta_dn0p10'       }
      { 'w_theta_flux_stall'        'w-theta'   2 1 'theta_stall'        }
      { 'w_theta_flux_all_cld'      'w-theta'   2 1 'theta_all_cld'      }

      { 'theta_e_M1'            'theta_e'   1 1 'theta_e'        }
      { 'theta_e_M1_ud0p10'     'theta_e'   1 1 'theta_e_ud0p10' }
%      { 'theta_e_M1_up0p10'     'theta_e'   1 1 'theta_e_up0p10' }
%      { 'theta_e_M1_dn0p10'     'theta_e'   1 1 'theta_e_dn0p10' }

      { 'w_theta_v_flux'              'w-theta_v' 2 1 'theta_v'              }
      { 'w_theta_v_flux_ud0p10'       'w-theta_v' 2 1 'theta_v_ud0p10'       }
      { 'w_theta_v_flux_up0p10'       'w-theta_v' 2 1 'theta_v_up0p10'       }
      { 'w_theta_v_flux_dn0p10'       'w-theta_v' 2 1 'theta_v_dn0p10'       }
      { 'w_theta_v_flux_stall'        'w-theta_v' 2 1 'theta_v_stall'        }
      { 'w_theta_v_flux_up0p10_stall' 'w-theta_v' 2 1 'theta_v_up0p10_stall' }
      { 'w_theta_v_flux_dn0p10_stall' 'w-theta_v' 2 1 'theta_v_dn0p10_stall' }
      { 'w_theta_v_flux_all_cld'      'w-theta_v' 2 1 'theta_v_all_cld'      }
      { 'w_theta_v_flux_col_up0p10'   'w-theta_v' 2 1 'theta_v_col_up0p10'   }
      { 'w_theta_v_flux_col_dn0p10'   'w-theta_v' 2 1 'theta_v_col_dn0p10'   }
      { 'w_theta_v_flux_col_ud0p10'   'w-theta_v' 2 1 'theta_v_col_ud0p10'   }
      { 'w_theta_v_flux_up0p10_all_cld' 'w-theta_v' 2 1 'theta_v_up0p10_all_cld'   }
      { 'w_theta_v_flux_dn0p10_all_cld' 'w-theta_v' 2 1 'theta_v_dn0p10_all_cld'   }
      { 'w_theta_v_flux_ud0p10_all_cld' 'w-theta_v' 2 1 'theta_v_ud0p10_all_cld'   }
      { 'w_theta_v_flux_up0p01_all_cld' 'w-theta_v' 2 1 'theta_v_up0p01_all_cld'   }
      { 'w_theta_v_flux_dn0p01_all_cld' 'w-theta_v' 2 1 'theta_v_dn0p01_all_cld'   }
      { 'w_theta_v_flux_ud0p01_all_cld' 'w-theta_v' 2 1 'theta_v_ud0p01_all_cld'   }

      { 'w_vapor_flux'              'w-vapor'   2 1 'vapor'              }
      { 'w_vapor_flux_ud0p10'       'w-vapor'   2 1 'vapor_ud0p10'       }
      { 'w_vapor_flux_up0p10'       'w-vapor'   2 1 'vapor_up0p10'       }
      { 'w_vapor_flux_dn0p10'       'w-vapor'   2 1 'vapor_dn0p10'       }
      { 'w_vapor_flux_stall'        'w-vapor'   2 1 'vapor_stall'        }
      { 'w_vapor_flux_up0p10_stall' 'w-vapor'   2 1 'vapor_up0p10_stall' }
      { 'w_vapor_flux_dn0p10_stall' 'w-vapor'   2 1 'vapor_dn0p10_stall' }
      { 'w_vapor_flux_all_cld'      'w-vapor'   2 1 'vapor_all_cld'      }
      { 'w_vapor_flux_col_up0p10'   'w-vapor'   2 1 'vapor_col_up0p10'   }
      { 'w_vapor_flux_col_dn0p10'   'w-vapor'   2 1 'vapor_col_dn0p10'   }
      { 'w_vapor_flux_col_ud0p10'   'w-vapor'   2 1 'vapor_col_ud0p10'   }
      { 'w_vapor_flux_up0p10_all_cld' 'w-vapor'   2 1 'vapor_up0p10_all_cld'   }
      { 'w_vapor_flux_dn0p10_all_cld' 'w-vapor'   2 1 'vapor_dn0p10_all_cld'   }
      { 'w_vapor_flux_ud0p10_all_cld' 'w-vapor'   2 1 'vapor_ud0p10_all_cld'   }
      { 'w_vapor_flux_up0p01_all_cld' 'w-vapor'   2 1 'vapor_up0p01_all_cld'   }
      { 'w_vapor_flux_dn0p01_all_cld' 'w-vapor'   2 1 'vapor_dn0p01_all_cld'   }
      { 'w_vapor_flux_ud0p01_all_cld' 'w-vapor'   2 1 'vapor_ud0p01_all_cld'   }

      { 'w_speed_flux'          'w-speed'   2 1 'speed'          }
      { 'w_speed_flux_ud0p10'   'w-speed'   2 1 'speed_ud0p10'   }
      { 'w_speed_flux_up0p10'   'w-speed'   2 1 'speed_up0p10'   }
      { 'w_speed_flux_dn0p10'   'w-speed'   2 1 'speed_dn0p10'   }
      { 'w_speed_flux_stall'    'w-speed'   2 1 'speed_stall'    }
      { 'w_speed_flux_all_cld'  'w-speed'   2 1 'speed_all_cld'  }

      { 'cloud_M1'              'cloud'     1 1 'cloud'          }
      { 'cloud_M1_c0p01'        'cloud'     1 1 'cloud_c0p01'    }
      { 'cloud_M1_c0p10'        'cloud'     1 1 'cloud_c0p10'    }
      { 'cloud_M1_c0p01_stall'  'cloud'     1 1 'cloud_c0p01_stall'    }
      { 'cloud_M1_c0p10_stall'  'cloud'     1 1 'cloud_c0p10_stall'    }
      { 'cloud_M1_col_cld'      'cloud'     1 1 'cloud_col_cld'        }


      % fluxes (covariances)
      { 'w_theta_v_flux'              'w-theta_v' 1 2 'w-theta_v'              }
      { 'w_theta_v_flux_ud0p10'       'w-theta_v' 1 2 'w-theta_v_ud0p10'       }
      { 'w_theta_v_flux_up0p10'       'w-theta_v' 1 2 'w-theta_v_up0p10'       }
      { 'w_theta_v_flux_dn0p10'       'w-theta_v' 1 2 'w-theta_v_dn0p10'       }
      { 'w_theta_v_flux_stall'        'w-theta_v' 1 2 'w-theta_v_stall'        }
      { 'w_theta_v_flux_up0p10_stall' 'w-theta_v' 1 2 'w-theta_v_up0p10_stall' }
      { 'w_theta_v_flux_dn0p10_stall' 'w-theta_v' 1 2 'w-theta_v_dn0p10_stall' }
      { 'w_theta_v_flux_all_cld'      'w-theta_v' 1 2 'w-theta_v_all_cld'      }
      { 'w_theta_v_flux_col_up0p10'   'w-theta_v' 1 2 'w-theta_v_col_up0p10' }
      { 'w_theta_v_flux_col_dn0p10'   'w-theta_v' 1 2 'w-theta_v_col_dn0p10' }
      { 'w_theta_v_flux_col_ud0p10'   'w-theta_v' 1 2 'w-theta_v_col_ud0p10' }
      { 'w_theta_v_flux_up0p10_all_cld' 'w-theta_v' 1 2 'w-theta_v_up0p10_all_cld' }
      { 'w_theta_v_flux_dn0p10_all_cld' 'w-theta_v' 1 2 'w-theta_v_dn0p10_all_cld' }
      { 'w_theta_v_flux_ud0p10_all_cld' 'w-theta_v' 1 2 'w-theta_v_ud0p10_all_cld' }
      { 'w_theta_v_flux_up0p01_all_cld' 'w-theta_v' 1 2 'w-theta_v_up0p01_all_cld' }
      { 'w_theta_v_flux_dn0p01_all_cld' 'w-theta_v' 1 2 'w-theta_v_dn0p01_all_cld' }
      { 'w_theta_v_flux_ud0p01_all_cld' 'w-theta_v' 1 2 'w-theta_v_ud0p01_all_cld' }

      { 'w_theta_flux'           'w-theta'   1 2 'w-theta'            }
      { 'w_theta_flux_ud0p10'    'w-theta'   1 2 'w-theta_ud0p10'     }
      { 'w_theta_flux_up0p10'    'w-theta'   1 2 'w-theta_up0p10'     }
      { 'w_theta_flux_dn0p10'    'w-theta'   1 2 'w-theta_dn0p10'     }
      { 'w_theta_flux_stall'     'w-theta'   1 2 'w-theta_stall'      }
      { 'w_theta_flux_all_cld'   'w-theta'   1 2 'w-theta_all_cld'    }

      { 'w_vapor_flux'              'w-vapor'   1 2 'w-vapor'              }
      { 'w_vapor_flux_ud0p10'       'w-vapor'   1 2 'w-vapor_ud0p10'       }
      { 'w_vapor_flux_up0p10'       'w-vapor'   1 2 'w-vapor_up0p10'       }
      { 'w_vapor_flux_dn0p10'       'w-vapor'   1 2 'w-vapor_dn0p10'       }
      { 'w_vapor_flux_stall'        'w-vapor'   1 2 'w-vapor_stall'        }
      { 'w_vapor_flux_up0p10_stall' 'w-vapor'   1 2 'w-vapor_up0p10_stall' }
      { 'w_vapor_flux_dn0p10_stall' 'w-vapor'   1 2 'w-vapor_dn0p10_stall' }
      { 'w_vapor_flux_all_cld'      'w-vapor'   1 2 'w-vapor_all_cld'      }
      { 'w_vapor_flux_col_up0p10'   'w-vapor'   1 2 'w-vapor_col_up0p10'   }
      { 'w_vapor_flux_col_dn0p10'   'w-vapor'   1 2 'w-vapor_col_dn0p10'   }
      { 'w_vapor_flux_col_ud0p10'   'w-vapor'   1 2 'w-vapor_col_ud0p10'   }
      { 'w_vapor_flux_up0p10_all_cld' 'w-vapor'   1 2 'w-vapor_up0p10_all_cld'   }
      { 'w_vapor_flux_dn0p10_all_cld' 'w-vapor'   1 2 'w-vapor_dn0p10_all_cld'   }
      { 'w_vapor_flux_ud0p10_all_cld' 'w-vapor'   1 2 'w-vapor_ud0p10_all_cld'   }
      { 'w_vapor_flux_up0p01_all_cld' 'w-vapor'   1 2 'w-vapor_up0p01_all_cld'   }
      { 'w_vapor_flux_dn0p01_all_cld' 'w-vapor'   1 2 'w-vapor_dn0p01_all_cld'   }
      { 'w_vapor_flux_ud0p01_all_cld' 'w-vapor'   1 2 'w-vapor_ud0p01_all_cld'   }

      { 'w_speed_flux'           'w-speed'   1 2 'w-speed'            }
      { 'w_speed_flux_ud0p10'    'w-speed'   1 2 'w-speed_ud0p10'     }
      { 'w_speed_flux_up0p10'    'w-speed'   1 2 'w-speed_up0p10'     }
      { 'w_speed_flux_dn0p10'    'w-speed'   1 2 'w-speed_dn0p10'     }
      { 'w_speed_flux_stall'     'w-speed'   1 2 'w-speed_stall'      }
      { 'w_speed_flux_all_cld'   'w-speed'   1 2 'w-speed_all_cld'    }

      % variances
      { 'w_M3'                     'w-w-w'     1 2 'w-w'              }
      { 'w_M3_ud0p10'              'w-w-w'     1 2 'w-w_ud0p10'       }
      { 'w_M3_up0p10'              'w-w-w'     1 2 'w-w_up0p10'       }
      { 'w_M3_dn0p10'              'w-w-w'     1 2 'w-w_dn0p10'       }
      { 'w_M3_stall'               'w-w-w'     1 2 'w-w_stall'        }
      { 'w_M3_up0p10_stall'        'w-w-w'     1 2 'w-w_up0p10_stall' }
      { 'w_M3_dn0p10_stall'        'w-w-w'     1 2 'w-w_dn0p10_stall' }
      { 'w_M3_all_cld'             'w-w-w'     1 2 'w-w_all_cld'      }
      { 'w_M3_col_up0p10'          'w-w-w'     1 2 'w-w_col_up0p10'   }
      { 'w_M3_col_dn0p10'          'w-w-w'     1 2 'w-w_col_dn0p10'   }
      { 'w_M3_col_ud0p10'          'w-w-w'     1 2 'w-w_col_ud0p10'   }
      { 'w_M3_up0p10_all_cld'      'w-w-w'     1 2 'w-w_up0p10_all_cld'   }
      { 'w_M3_dn0p10_all_cld'      'w-w-w'     1 2 'w-w_dn0p10_all_cld'   }
      { 'w_M3_ud0p10_all_cld'      'w-w-w'     1 2 'w-w_ud0p10_all_cld'   }
      { 'w_M3_up0p01_all_cld'      'w-w-w'     1 2 'w-w_up0p01_all_cld'   }
      { 'w_M3_dn0p01_all_cld'      'w-w-w'     1 2 'w-w_dn0p01_all_cld'   }
      { 'w_M3_ud0p01_all_cld'      'w-w-w'     1 2 'w-w_ud0p01_all_cld'   }

      % skews
      { 'w_M3'                     'w-w-w'     1 3 'w-w-w'              }
      { 'w_M3_ud0p10'              'w-w-w'     1 3 'w-w-w_ud0p10'       }
      { 'w_M3_up0p10'              'w-w-w'     1 3 'w-w-w_up0p10'       }
      { 'w_M3_dn0p10'              'w-w-w'     1 3 'w-w-w_dn0p10'       }
      { 'w_M3_stall'               'w-w-w'     1 3 'w-w-w_stall'        }
      { 'w_M3_up0p10_stall'        'w-w-w'     1 3 'w-w-w_up0p10_stall' }
      { 'w_M3_all_cld'             'w-w-w'     1 3 'w-w-w_all_cld'      }
      { 'w_M3_col_up0p10'          'w-w-w'     1 3 'w-w-w_col_up0p10'   }
      { 'w_M3_col_dn0p10'          'w-w-w'     1 3 'w-w-w_col_dn0p10'   }
      { 'w_M3_col_ud0p10'          'w-w-w'     1 3 'w-w-w_col_ud0p10'   }
      { 'w_M3_up0p10_all_cld'      'w-w-w'     1 3 'w-w-w_up0p10_all_cld'   }
      { 'w_M3_dn0p10_all_cld'      'w-w-w'     1 3 'w-w-w_dn0p10_all_cld'   }
      { 'w_M3_ud0p10_all_cld'      'w-w-w'     1 3 'w-w-w_ud0p10_all_cld'   }
      { 'w_M3_up0p01_all_cld'      'w-w-w'     1 3 'w-w-w_up0p01_all_cld'   }
      { 'w_M3_dn0p01_all_cld'      'w-w-w'     1 3 'w-w-w_dn0p01_all_cld'   }
      { 'w_M3_ud0p01_all_cld'      'w-w-w'     1 3 'w-w-w_ud0p01_all_cld'   }

      };

    % do one hour time average at beginning, middle and end of sampling period
    % also do time average over entire sampling inerval
    TSstart = 12;
    TSend = 13;
    TMstart = 23.5;
    TMend = 24.5;
    TEstart = 35;
    TEend = 36;

    TAstart = 12;
    TAend = 36;

    TotalN = 158404; % excluding borders --> 398 * 398

    for icase = 1:length(Config.Cases)
      Case = Config.Cases(icase).Cname;
      OutFname = sprintf('%s/moments_%s.h5', Ddir, Case);

      fprintf('***************************************************************\n');
      fprintf('Generating moment/flux profiles:\n');
      fprintf('  Case: %s\n', Case);
      fprintf('  Output file: %s\n', OutFname);
      fprintf('\n');

      % The vars and coords are written into the file in append mode so write 
      % the file here without append mode in order to create a new file each time
      % this script is run.
      hdf5write(OutFname, 'header', Case);

      for ivar = 1:length(VarList)
        InFprefix = VarList{ivar}{1};
        InVname   = VarList{ivar}{2};
        InTerm    = VarList{ivar}{3};
        InOrder   = VarList{ivar}{4};
        OutVname  = VarList{ivar}{5};

        InFname = sprintf('%s/%s_%s.h5', Ddir, InFprefix, Case);

        fprintf('  Input file: %s\n', InFname);
        fprintf('    Var name: %s\n', InVname);
        fprintf('    Var term: %d\n', InTerm);
        fprintf('    Var order: %d\n', InOrder);
        fprintf('\n');

        % read in and 
        TERMS = squeeze(hdf5read(InFname, InVname));
        NPTS = squeeze(hdf5read(InFname, 'num_points'));
        Z = squeeze(hdf5read(InFname, 'z_coords'));
        T = squeeze(hdf5read(InFname, 't_coords')) / 3600;   % hr

        Nz = length(Z);
        Nt = length(T);

        % *_M1 files will lose their first two dimensions, so put them back
        if (strcmp(InVname, 'cloud') || strcmp(InVname, 'theta_e'))
          TERMS = reshape(TERMS, [ 1 1 Nz Nt ]);
        end

        % find the indices of the start, mid, end and "all" time intervals
        TS1 = find(T >= TSstart, 1, 'first');
        TS2 = find(T <= TSend, 1, 'last');

        TM1 = find(T >= TMstart, 1, 'first');
        TM2 = find(T <= TMend, 1, 'last');

        TE1 = find(T >= TEstart, 1, 'first');
        TE2 = find(T <= TEend, 1, 'last');

        TA1 = find(T >= TAstart, 1, 'first');
        TA2 = find(T <= TAend, 1, 'last');

        % MOMENTS is organized: (Nterms, Norder, Nz)
        [ MOMENTS_TS    OUT_NPTS ] = GenMoments(TERMS, NPTS, TS1, TS1, Tfirst);
        [ MOMENTS_TM    OUT_NPTS ] = GenMoments(TERMS, NPTS, TM1, TM2, Tfirst);
        [ MOMENTS_TE    OUT_NPTS ] = GenMoments(TERMS, NPTS, TE1, TE2, Tfirst);
        [ MOMENTS_TA OUT_NPTS ] = GenMoments(TERMS, NPTS, TA1, TA2, Tfirst);

        OUT_MOMENTS_TS = squeeze(MOMENTS_TS(InTerm, InOrder, :));
        OUT_MOMENTS_TM = squeeze(MOMENTS_TM(InTerm, InOrder, :));
        OUT_MOMENTS_TE = squeeze(MOMENTS_TE(InTerm, InOrder, :));
        OUT_MOMENTS_TA = squeeze(MOMENTS_TA(InTerm, InOrder, :));

        % GenMoments() will fill levels with zero count (NPTS == 0) with nans. For
        % cloud mass we want these moments to be zero.
        if (strcmp(InVname, 'cloud'))
          OUT_MOMENTS_TS(isnan(OUT_MOMENTS_TS)) = 0;
          OUT_MOMENTS_TM(isnan(OUT_MOMENTS_TM)) = 0;
          OUT_MOMENTS_TE(isnan(OUT_MOMENTS_TE)) = 0;
          OUT_MOMENTS_TA(isnan(OUT_MOMENTS_TA)) = 0;
        end

        % Generate a fraction statistic - the ratio of number of points selected 
        % to total number of points in domain
        PROF_FRAC = NPTS ./ TotalN;
        PROF_FRAC_S = squeeze(mean(PROF_FRAC(:,TS1:TS2), 2));
        PROF_FRAC_M = squeeze(mean(PROF_FRAC(:,TM1:TM2), 2));
        PROF_FRAC_E = squeeze(mean(PROF_FRAC(:,TE1:TE2), 2));
        PROF_FRAC_A = squeeze(mean(PROF_FRAC(:,TA1:TA2), 2));


        % Write out data - put in dummy x, y and t coordinates
        Xdummy = 1;
        Ydummy = 1;
        Tdummy = 1;

        OutVar = reshape(OUT_MOMENTS_TS, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_start', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(OUT_MOMENTS_TM, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_mid', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(OUT_MOMENTS_TE, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_end', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(OUT_MOMENTS_TA, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_all', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(PROF_FRAC_S, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_frac_start', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(PROF_FRAC_M, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_frac_mid', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(PROF_FRAC_E, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_frac_end', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(PROF_FRAC_A, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_frac_all', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        fprintf('\n');
      end

      % all vars are written out to the file, now write out the coordinates
      hdf5write(OutFname, 'x_coords', Xdummy, 'WriteMode', 'append');
      hdf5write(OutFname, 'y_coords', Ydummy, 'WriteMode', 'append');
      hdf5write(OutFname, 'z_coords', Z,      'WriteMode', 'append');
      hdf5write(OutFname, 't_coords', Tdummy, 'WriteMode', 'append');
    end
