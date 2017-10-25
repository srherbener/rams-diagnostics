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

% air density (kg/m3) corresponding to ATEX grid levels
Rho = [ 1.2274
    1.2226
    1.2179
    1.2133
    1.2086
    1.2040
    1.1994
    1.1948
    1.1902
    1.1857
    1.1812
    1.1767
    1.1722
    1.1677
    1.1633
    1.1589
    1.1545
    1.1501
    1.1457
    1.1414
    1.1371
    1.1328
    1.1285
    1.1243
    1.1200
    1.1158
    1.1116
    1.1074
    1.1033
    1.0991
    1.0950
    1.0909
    1.0868
    1.0828
    1.0787
    1.0747
    1.0707
    1.0667
    1.0627
    1.0588
    1.0549
    1.0509
    1.0470
    1.0432
    1.0393
    1.0355
    1.0316
    1.0278
    1.0240
    1.0203
    1.0165
    1.0128
    1.0090
    1.0053
    1.0016
    0.9980
    0.9943
    0.9907
    0.9870
    0.9834
    0.9798
    0.9763
    0.9727
    0.9692
    0.9656
    0.9621
    0.9586
    0.9552
    0.9517
    0.9482
    0.9448
    0.9414
    0.9380
    0.9346
    0.9312
    0.9279
    0.9245
    0.9212
    0.9179
    0.9146
    0.9113
    0.9080
    0.9048
    0.9015
    0.8983
    0.8951
    0.8919
    0.8887
    0.8856
    0.8824
    0.8793
    0.8761
    0.8730
    0.8699
    0.8668
    0.8638
    0.8607
    0.8577
    0.8546
    0.8516
    0.8486
    ];

  VarSets = {
    { 'hda_cloud_ot'    { 'cot' 'cot_all_cld'  } 'avg_cot'         }
    { 'hda_cloud_frac'  { 'cloud_frac'         } 'avg_dom_cfrac'   }
    { 'hda_inv_height'  { 'inv_height'         } 'avg_inv_height'  }
    { 'hda_cdepth'      { 'cdepth_all_cld'     } 'avg_cdepth'      }

    { 'hda_lwp2cdepth'  { 'lwp2cdepth_all_cld' } 'avg_lwp2cdepth'  }


    { 'hda_cloud_cond'  { 'cloud_cond_all_cld' } 'avg_cloud_cond'  }
    { 'hda_cloud_evap'  { 'cloud_evap_all_cld' } 'avg_cloud_evap'  }

    { 'hda_rain_cond'   { 'rain_cond_all_cld'  } 'avg_rain_cond'   }
    { 'hda_rain_evap'   { 'rain_evap_all_cld'  } 'avg_rain_evap'   }

    { 'hda_driz_cond'   { 'driz_cond_all_cld'  } 'avg_driz_cond'   }
    { 'hda_driz_evap'   { 'driz_evap_all_cld'  } 'avg_driz_evap'   }

    { 'hda_lw_flux'     { 'lw_flux_all_cld'    } 'avg_net_lw_flux' }


    { 'hda_theta'       { 'theta'       } 'avg_theta'       }

    { 'hda_scbot'       { 'scbot'       } 'avg_scbot'       }
    { 'hda_sctop'       { 'sctop'       } 'avg_sctop'       }

    { 'hda_cloud'       { 'cloud'       } 'avg_cloud'       }
    { 'hda_cloud_diam'  { 'cloud_diam'  } 'avg_cloud_diam'  }
    { 'hda_cloud_num'   { 'cloud_num'   } 'avg_cloud_num'   }
   
    { 'hda_rain'        { 'rain'        } 'avg_rain'        }
    { 'hda_rain_diam'   { 'rain_diam'   } 'avg_rain_diam'   }
    { 'hda_rain_num'    { 'rain_num'    } 'avg_rain_num'    }

    };
  Nset = length(VarSets);

  % For time and height averaging
  Tstart = 24;
  Tend   = 48;
  Tname  = 'tavg';

  Zstart = 0;
  Zend   = 4000;
  Zname  = 'zavg';

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
        InVarName = sprintf('/%s', InVarName);

        if (regexp(OutVarName, '^lcl'))  % grab hda_lcl file from Ddir instead of Tdir
          InFile = sprintf('%s/hda_%s_%s.h5', Ddir, OutVarName , Case);
        else
          InFile = sprintf('%s/hda_%s_%s.h5', Tdir, OutVarName , Case);
        end
        fprintf('      %s (%s)\n', InFile, InVarName);

        % grab the hda data and the t coordinates
        HDA  = squeeze(hdf5read(InFile, InVarName));
        T = squeeze(hdf5read(InFile, 't_coords'))/3600; % hours
        Z = squeeze(hdf5read(InFile, 'z_coords'));

        Nt = length(T); 
        Nz = length(Z);

        % for height averaging
        Z1 = find(Z >= Zstart, 1, 'first');
        Z2 = find(Z <= Zend,   1, 'last');

        % for time averaging
        T1 = find(T >= Tstart, 1, 'first');
        T2 = find(T <= Tend,   1, 'last');

        % If doing a hydrometeor number concentration, the
        % HDA units are #/kg. Convert to #/cm3 using:
        %
        %   #/cm3 = #/kg * Rho (kg/m3) * 1e-6 (m3/cm3)
        %
        if (strcmp(InVarName, '/hda_cloud_num') || ...
            strcmp(InVarName, '/hda_rain_num'))
          if (ndims(HDA) == 3)
            % HDA is (2,z,t)
            HDA(1,:,:) = squeeze(HDA(1,:,:)) .* repmat(Rho,[1 Nt]) .* 1e-6;
          elseif (ndims(HDA) == 2)
            % HDA is (2,t)
            % Rho(2) is first level above the surface
            HDA(1,:) = squeeze(HDA(1,:)) .* Rho(2) .* 1e-6;
          end
        end

        % do time series of spatial averaging
        % HDA is either (2,z,t) or (2,t), where along the first dimension
        %   entry 1 are the sums
        %   entry 2 are the counts
        if (ndims(HDA) == 3)
          AVG  = squeeze(HDA(1,:,:));
          NPTS = squeeze(HDA(2,:,:));
          AVG_ZALL  = squeeze(sum(HDA(1,Z1:Z2,:), 2));
          NPTS_ZALL = squeeze(sum(HDA(2,Z1:Z2,:), 2));
        elseif (ndims(HDA) == 2)
          AVG  = squeeze(HDA(1,:));
          NPTS = squeeze(HDA(2,:));
          AVG_ZALL  = AVG;
          NPTS_ZALL = NPTS;
        else
          % force AVG to be nan
          AVG = 0;
          NPTS = 0;
          AVG_ZALL  = 0;
          NPTS_ZALL = 0;
        end

        % Form the averages
        AVG = AVG ./ NPTS;
        AVG_ZALL = AVG_ZALL ./ NPTS_ZALL;
        [ AVG_TALL NPTS_TALL ] = CountsToAvg(HDA, T1, T2);

        % levels and timesteps separated
        OutName = sprintf('%s_avg', OutVarName);
        hdf5write(OutFile, OutName, AVG, 'WriteMode', 'append'); 
        OutName = sprintf('%s_avg_npts', OutVarName);
        hdf5write(OutFile, OutName, NPTS, 'WriteMode', 'append');

        % levels averaged and timesteps separated
        OutName = sprintf('%s_%s', OutVarName, Zname);
        hdf5write(OutFile, OutName, AVG_ZALL, 'WriteMode', 'append'); 
        OutName = sprintf('%s_%s_npts', OutVarName, Zname);
        hdf5write(OutFile, OutName, NPTS_ZALL, 'WriteMode', 'append');

        % levels separated and timesteps averaged
        OutName = sprintf('%s_%s', OutVarName, Tname);
        hdf5write(OutFile, OutName, AVG_TALL, 'WriteMode', 'append'); 
        OutName = sprintf('%s_%s_npts', OutVarName, Tname);
        hdf5write(OutFile, OutName, NPTS_TALL, 'WriteMode', 'append');
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
