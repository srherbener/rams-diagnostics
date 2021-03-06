function [ ] = GenAvgFilesCtype(ConfigFile)
% GenAvgFilesCtype generate averages from tsavg hda files

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Tdir = Config.TsavgDir;
  Ddir = Config.DiagDir;

  FileHeader = 'ATEX PDF data';

  VarSets = {
%    { 'hda_pcprr'           { 'pcprr'  'pcprr_strnp'  'pcprr_strat'  'pcprr_cumul'  'pcprr_all_cld'  'pcprr_stall'  } 'avg_ctype_pcprr'  }
%    { 'hda_cloud_opt_thick' { 'cot'    'cot_strnp'    'cot_strat'    'cot_cumul'    'cot_all_cld'    'cot_stall'    } 'avg_ctype_cot'    }
%    { 'hda_cloud_depth'     { 'cdepth' 'cdepth_strnp' 'cdepth_strat' 'cdepth_cumul' 'cdepth_all_cld' 'cdepth_stall' } 'avg_ctype_cdepth' }
%    { 'hda_vint_cond'       { 'lwp'    'lwp_strnp'    'lwp_strat'    'lwp_cumul'    'lwp_all_cld'    'lwp_stall'    } 'avg_ctype_lwp'    }
%
%    { 'hda_lwp2cdepth'      { 'lwp2cdepth' 'lwp2cdepth_strnp' 'lwp2cdepth_strat' 'lwp2cdepth_cumul' 'lwp2cdepth_all_cld' 'lwp2cdepth_stall' } 'avg_ctype_lwp2cdepth'    }
%
%    { 'hda_cloud_mask'      { 'cfrac'  'cfrac_strnp'  'cfrac_strat'  'cfrac_cumul'  'cfrac_stmix' 'cfrac_scmix' 'cfrac_stnopr' 'cfrac_stdriz', 'cfrac_strain' } 'avg_ctype_cfrac'  }
%
%    { 'hda_lcl'             { 'lcl'    'lcl_stall'                                                                  } 'avg_ctype_lcl'  }




%    { 'hda_vapcldt'  { 'cloud_cond_stall'    'col_cloud_cond' 'cloud_cond_all_cld' 'cond_in_cloud' 'cond_in_cloud_nz'    } 'avg_ctype_cloud_cond'  }
%    { 'hda_vapcldt'  { 'cloud_evap_stall'    'col_cloud_evap' 'cloud_evap_all_cld' 'evap_in_cloud' 'evap_in_cloud_nz'    } 'avg_ctype_cloud_evap'  }
%    { 'hda_vapcldt'         { 'col_cloud_cond_evap'                                                                 } 'avg_ctype_cloud_cond_evap'  }

%    { 'hda_vapraint'  { 'rain_cond_all_cld' } 'avg_ctype_rain_cond'  }
%    { 'hda_vapraint'  { 'rain_evap_all_cld' } 'avg_ctype_rain_evap'  }

%    { 'hda_vapdrizt'  { 'driz_cond_all_cld' } 'avg_ctype_driz_cond'  }
%    { 'hda_vapdrizt'  { 'driz_evap_all_cld' } 'avg_ctype_driz_evap'  }

%    { 'hda_lh_vapt'  { 'col_lat_heat'      'lat_heat_all_cld' 'lat_heat_in_cloud' 'lat_heat_in_cloud_nz'     } 'avg_ctype_lat_heat'  }
%    { 'hda_lh_vapt'  { 'col_lat_cool'      'lat_cool_all_cld' 'lat_cool_in_cloud' 'lat_cool_in_cloud_nz'     } 'avg_ctype_lat_cool'  }
%    { 'hda_lh_vapt'         { 'col_lat_heat_cool'                                                                   } 'avg_ctype_lat_heat_cool'  }

    { 'hda_net_lw_flux'  { 'lw_flux_all_cld'     } 'avg_ctype_net_lw_flux'  }
    };
  Nset = length(VarSets);

  TimeSelects = {
    { 12    36    'TALL'   }      
    { 12    13    'TSTART' }
    { 23.5  24.5  'TMID'   }
    { 35    36    'TEND'   }
    };
  Ntsel = length(TimeSelects);

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;

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

        OutAvgName  = sprintf('%s_avg', OutVarName);
        OutNptsName = sprintf('%s_npts', OutVarName);

        for its = 1:Ntsel
          Tstart = TimeSelects{its}{1};
          Tend   = TimeSelects{its}{2};
          Tname  = TimeSelects{its}{3};

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
      end
      fprintf('\n');

      % output coords so that ReadSelectXyzt can handle the output files
      Xdummy = 1;
      Ydummy = 1;
      if (strcmp(OutFprefix, 'avg_ctype_cloud_cond') || ...
          strcmp(OutFprefix, 'avg_ctype_cloud_evap') || ...
          strcmp(OutFprefix, 'avg_ctype_cloud_cond_evap') || ...
          strcmp(OutFprefix, 'avg_ctype_driz_cond') || ...
          strcmp(OutFprefix, 'avg_ctype_driz_evap') || ...
          strcmp(OutFprefix, 'avg_ctype_rain_cond') || ...
          strcmp(OutFprefix, 'avg_ctype_rain_evap') || ...
          strcmp(OutFprefix, 'avg_ctype_lat_heat') || ...
          strcmp(OutFprefix, 'avg_ctype_lat_cool') || ...
          strcmp(OutFprefix, 'avg_ctype_net_lw_flux') || ...
          strcmp(OutFprefix, 'avg_ctype_lat_heat_cool'))
        Zout = Z;
      else
        Zout = 1;
      end
      Tdummy = 1;

      fprintf('    Writing: %s\n', OutFile);

      hdf5write(OutFile, '/x_coords', Xdummy, 'WriteMode', 'append');
      hdf5write(OutFile, '/y_coords', Ydummy, 'WriteMode', 'append');
      hdf5write(OutFile, '/z_coords', Zout,   'WriteMode', 'append');
      hdf5write(OutFile, '/t_coords', Tdummy, 'WriteMode', 'append');

      fprintf('\n');
    end
    fprintf('\n');
  end
end
