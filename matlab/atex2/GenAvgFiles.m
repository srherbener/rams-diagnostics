function [ ] = GenAvgFiles()
% GenAvgFiles generate averages from tsavg hda files

  Tdir = 'TsAveragedData';
  Ddir = 'DIAGS';

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  FileHeader = 'ATEX averaged data';

  CaseList = {
    'z.atex.ccn0050.sst293'
    'z.atex.ccn0100.sst293'
    'z.atex.ccn0200.sst293'
    'z.atex.ccn0400.sst293'
    'z.atex.ccn0800.sst293'
    'z.atex.ccn1600.sst293'

    'z.atex.ccn0050.sst298'
    'z.atex.ccn0100.sst298'
    'z.atex.ccn0200.sst298'
    'z.atex.ccn0400.sst298'
    'z.atex.ccn0800.sst298'
    'z.atex.ccn1600.sst298'
    };
Ncases = length(CaseList);

  VarSets = {
    { 'hda_cloud_ot'    { 'cot' 'cot_all_cld'  } 'avg_cot'         }
    { 'hda_cloud_mask'  { 'cloud_frac'  }        'avg_dom_cfrac'   }

    { 'hda_cloud_depth' { 'cdepth_all_cld'     } 'avg_cdepth'      }
    { 'hda_lwp2cdepth'  { 'lwp2cdepth_all_cld' } 'avg_lwp2cdepth'  }

    { 'hda_vapcldt'     { 'cloud_cond_all_cld' } 'avg_cloud_cond'  }
    { 'hda_vapcldt'     { 'cloud_evap_all_cld' } 'avg_cloud_evap'  }

    { 'hda_vapraint'    { 'rain_cond_all_cld'  } 'avg_rain_cond'   }
    { 'hda_vapraint'    { 'rain_evap_all_cld'  } 'avg_rain_evap'   }

    { 'hda_vapdrizt'    { 'driz_cond_all_cld'  } 'avg_driz_cond'   }
    { 'hda_vapdrizt'    { 'driz_evap_all_cld'  } 'avg_driz_evap'   }

    { 'hda_net_lw_flux' { 'lw_flux_all_cld'    } 'avg_net_lw_flux' }

    { 'hda_cloud'       { 'cloud_c0p01'        } 'avg_cloud'       }
    { 'hda_cloud_diam'  { 'cloud_diam_c0p01'   } 'avg_cloud_diam'  }
    { 'hda_cloud_num'   { 'cloud_num_c0p01'    } 'avg_cloud_num'   }
    
    { 'hda_rain'        { 'rain_r0p01'         } 'avg_rain'        }
    { 'hda_rain_diam'   { 'rain_diam_r0p01'    } 'avg_rain_diam'   }
    { 'hda_rain_num'    { 'rain_num_r0p01'     } 'avg_rain_num'    }

    };
  Nset = length(VarSets);

  % For time averaging
  Tstart = 24;
  Tend   = 48;
  Tname  = 'TALL';

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('***************************************************************\n');
    fprintf('Generating avg data:\n');
    fprintf('  Case: %s\n', Case);

    for iset = 1:Nset
      InVname     = VarSets{iset}{1};
      OutVarList  = VarSets{iset}{2};
      OutFprefix  = VarSets{iset}{3};

      fprintf('    Input files:\n');

      % Write a header into the file so that a write statement without append mode
      % will be run first which consequently will erase an existing file and replace
      % it with the contents about to be generated here.
      OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);
      hdf5write(OutFile, 'Header', FileHeader); 

      Nvar = length(OutVarList);
      for ivar = 1:Nvar
        OutVarName = OutVarList{ivar};

        % Special case for cfrac_strnp, etc. These have the strnp, cumul etc appended to
        % then end of the dataset name inside the HDF5 file.
        if (regexp(OutVarName, '^cfrac_'))
          InVarName = sprintf('%s_%s', InVname, regexprep(OutVarName, '^cfrac_', ''));
        else
          InVarName = InVname;
        end

        if (regexp(OutVarName, '^lcl'))  % grab hda_lcl file from Ddir instead of Tdir
          InFile = sprintf('%s/hda_%s_%s.h5', Ddir, OutVarName , Case);
        else
          InFile = sprintf('%s/hda_%s_%s.h5', Tdir, OutVarName , Case);
        end
        fprintf('      %s --> %s\n', InFile, InVarName);

        % grab the hda data and the t coordinates
        HDA  = squeeze(hdf5read(InFile, InVarName));
        T = squeeze(hdf5read(InFile, 't_coords'))/3600; % hours
        Z = squeeze(hdf5read(InFile, 'z_coords'));

        Nt = length(T); 
        Nz = length(Z);

        OutAvgName  = sprintf('%s_avg', OutVarName);
        OutNptsName = sprintf('%s_npts', OutVarName);

        % do time series of spatial averaging
        % HDA is either (2,z,t) or (2,t), where along the first dimension
        %   entry 1 are the sums
        %   entry 2 are the counts
        if (ndims(HDA) == 3)
          AVG  = squeeze(sum(HDA(1,:,:), 2));
          NPTS = squeeze(sum(HDA(2,:,:), 2));
        elseif (ndims(HDA) == 2)
          AVG  = squeeze(HDA(1,:));
          NPTS = squeeze(HDA(2,:));
        else
          % force AVG to be nan
          AVG = 0;
          NPTS = 0;
        end
        AVG = AVG ./ NPTS;
        
        % time series
        OutName = sprintf('%s_TSERIES', OutAvgName);
        hdf5write(OutFile, OutName, AVG, 'WriteMode', 'append'); 
        OutName = sprintf('%s_TSERIES', OutNptsName);
        hdf5write(OutFile, OutName, NPTS, 'WriteMode', 'append');

        % time averaging
        T1 = find(T >= Tstart, 1, 'first');
        T2 = find(T <= Tend,   1, 'last');

        % sum up bins across selected times, and convert HDA counts to an average
        [ AVG NPTS ] = CountsToAvg(HDA, T1, T2);

        % With the data selection it is possible to get no data points selected (all counts
        % equal to zero in hda file). When this happens, get nans in Avg since the sum
        % is zero. Ie, the AVG entry is 0/0 --> nan.  Change nans back to zeros to help
        % make plots look nicer.
        AVG(isnan(AVG)) = 0;

        OutName = sprintf('%s_%s', OutAvgName, Tname);
        hdf5write(OutFile, OutName, AVG, 'WriteMode', 'append'); 
        OutName = sprintf('%s_%s', OutNptsName, Tname);
        hdf5write(OutFile, OutName, NPTS, 'WriteMode', 'append');
      end
      fprintf('\n');

      % output coords so that ReadSelectXyzt can handle the output files
      Xdummy = 1;
      Ydummy = 1;

      fprintf('    Writing: %s\n', OutFile);

      hdf5write(OutFile, '/x_coords', Xdummy, 'WriteMode', 'append');
      hdf5write(OutFile, '/y_coords', Ydummy, 'WriteMode', 'append');
      hdf5write(OutFile, '/z_coords', Z,      'WriteMode', 'append');
      hdf5write(OutFile, '/t_coords', T,      'WriteMode', 'append');

      fprintf('\n');
    end
    fprintf('\n');
  end
end
