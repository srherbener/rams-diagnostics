function [ ] = GenFsAverages()

  Ddir = 'DIAGS';

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
    'TSD_NONSAL_NODUST'
    };

  Ncases = length(CaseList);

  VarList = {
    %    Input File Name              Input Var Name     Output Var Name

    % Vertical coords -> pressure
    { 'DIAGS/storm_meas_<CASE>.h5'   '/max_wind_t'     'avg_max_wind_t'    }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/ps_max_wind_t'  'ps_avg_max_wind_t' }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/s_max_wind_t'   's_avg_max_wind_t'  }

    { 'DIAGS/storm_meas_<CASE>.h5'   '/avg_wind_t'     'avg_wind_t'    }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/ps_avg_wind_t'  'ps_avg_wind_t' }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/s_avg_wind_t'   's_avg_wind_t'  }

    { 'DIAGS/storm_meas_<CASE>.h5'   '/rmw_t'          'avg_rmw_t'         }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/ps_rmw_t'       'ps_avg_rmw_t'      }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/s_rmw_t'        's_avg_rmw_t'       }

    { 'DIAGS/storm_meas_<CASE>.h5'   '/min_slp'        'avg_min_press'     }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/ps_min_slp'     'ps_avg_min_press'  }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/s_min_slp'      's_avg_min_press'   }

    { 'DIAGS/storm_meas_<CASE>.h5'   '/ike'            'avg_ike'           }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/ps_ike'         'ps_avg_ike'        }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/s_ike'          's_avg_ike'         }
 
    { 'DIAGS/storm_meas_<CASE>.h5'   '/pcprate'        'avg_pcprate'       }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/ps_pcprate'     'ps_avg_pcprate'    }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/s_pcprate'      's_avg_pcprate'     }

    % Vertical coords -> pressure
    { 'DIAGS/storm_meas_<CASE>.h5'   '/max_wind_t_p'       'avg_max_wind_t_p'    }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/ps_max_wind_t_p'    'ps_avg_max_wind_t_p' }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/s_max_wind_t_p'     's_avg_max_wind_t_p'  }

    { 'DIAGS/storm_meas_<CASE>.h5'   '/avg_wind_t_p'       'avg_wind_t_p'        }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/ps_avg_wind_t_p'    'ps_avg_wind_t_p'     }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/s_avg_wind_t_p'     's_avg_wind_t_p'      }

    { 'DIAGS/storm_meas_<CASE>.h5'   '/rmw_t_p'            'avg_rmw_t_p'         }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/ps_rmw_t_p'         'ps_avg_rmw_t_p'      }
    { 'DIAGS/storm_meas_<CASE>.h5'   '/s_rmw_t_p'          's_avg_rmw_t_p'       }

    };

  Nvars = length(VarList);

  % Put all of the results into one output file.
  % If the file exists, remove it so that the HDF5 commands
  % can create a new datasets.
  OutFile = sprintf('%s/fs_averages.h5', Ddir);
  if (exist(OutFile, 'file') == 2)
    delete(OutFile);
  end

  % Create time averages of metrics that will be used for the factor separation.
  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating averages of metrics for factor separation analysis:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for ivar = 1:Nvars
      InFile   = VarList{ivar}{1};
      InVname  = VarList{ivar}{2};
      OutVname = VarList{ivar}{3};

      InFile = regexprep(InFile, '<CASE>' , Case);

      % read in var, form time average, write out result
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      VAR = squeeze(h5read(InFile, InVname));
 
      % calculate the time averages
      AVG = mean(VAR);
      Vsize = 1;

      % save the averages (raw data for the factor separation)
      OutVar = sprintf('/%s/%s', Case, OutVname);
      fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
      h5create(OutFile, OutVar,  Vsize);
      h5write( OutFile, OutVar,  AVG);

    end
    fprintf('\n');

  end
end 
