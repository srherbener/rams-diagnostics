function [ ] = GenStormMeas(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Adir = Config.AzavgDir;
  Tdir = Config.TsavgDir;
  Ddir = Config.DiagDir;

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
    'TSD_NONSAL_NODUST'
    };

  % Description of measurements
  MeasList = {
    { 'min_slp'  'azavg_hist' 'hist_press' 'press' 'farea' 0.10 }
    { 'max_wind' 'azavg_hist' 'hist_speed' 'speed' 'farea' 0.90 }
    };

  for icase = 1:length(Config.Cases)
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating storm measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for imeas = 1:length(MeasList)
      Mname    = MeasList{imeas}{1};
      Msource  = MeasList{imeas}{2};
      Mfprefix = MeasList{imeas}{3};
      Mvname   = MeasList{imeas}{4};
      Mmethod  = MeasList{imeas}{5};
      Mparam   = MeasList{imeas}{6};

      if (strncmp(Msource, 'azavg', 5))
        Mfile = sprintf('%s/%s_%s.h5', Adir, Mfprefix, Case);
      elseif (strncmp(Msource, 'tsavg', 5))
        Mfile = sprintf('%s/%s_%s.h5', Tdir, Mfprefix, Case);
      else
        Mfile = sprintf('%s_%s.h5', Mfprefix, Case);
      end

      Mvname = sprintf('/%s', Mvname);

      fprintf('  Measurement: %s\n', Mname);
      fprintf('    Reading: %s (%s)\n', Mfile, Mvname);
      fprintf('    Method: %s\n', Mmethod);
      if (strcmp(Mmethod, 'pct'))
        fprintf('      Percentile: %.4f\n', Mparam);
      end
      if (strcmp(Mmethod, 'farea'))
        fprintf('      Fractional area: %.4f\n', Mparam);
      end
      fprintf('\n');

      % Read in data, which is coming from either azavg or tsavg, meaning
      % that it will be 4D -> (x,y,z,t) regardless if all four dimensions
      % are actually used. Mname will indicate what size MDATA is and what
      % each of X, Y, Z, T represent.
      %
      MDATA = squeeze(h5read(Mfile, Mvname));
      X = squeeze(h5read(Mfile, '/x_coords'));
      Y = squeeze(h5read(Mfile, '/y_coords'));
      Z = squeeze(h5read(Mfile, '/z_coords'));
      T = squeeze(h5read(Mfile, '/t_coords'));

      [ Nx Ny Nz Nt ] = size(MDATA);

      if (strcmp(Msource, 'azavg_hist'))
        % have azimuthally selected histogram counts (histograms
        % corresponding to radial bands)
        %
        % MDATA --> (radius,counts,height,time)
        %    X --> radius values
        %    Y --> bin values
        %    Z --> height values
        %    T --> time values
        %
        % Need to call ReduceHists() to change bin counts to 
        % a single number.
        %
        RDATA = ReduceHists(MDATA, 2, Y, Mmethod, Mparam);

        if (strcmp(Mname, 'min_slp') || strcmp(Mname, 'max_wind'))
          % these measurements need to be taken from the k = 2 level
          % Mdata will be (r,z,t) at this point.
          RDATA = squeeze(RDATA(:,2,:));
       
          % RDATA is now (r,t)
          if (strncmp(Mname, 'min', 3))
            TSERIES = squeeze(min(RDATA, [], 1));
          elseif (strncmp(Mname, 'max', 3))
            TSERIES = squeeze(max(RDATA, [], 1));
          end
        end
      end

      % Write out measurement
      OutFile = sprintf('%s/%s_%s.h5', Ddir, Mname, Case);
      fprintf('    Writing: %s\n', OutFile)
      fprintf('\n');

      % create a new file
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end

      % Write out both RDATA and TSERIES for min_slp and max_wind
      if (strcmp(Mname, 'min_slp') || strcmp(Mname, 'max_wind'))
        h5create(OutFile, '/radial_time_series', size(RDATA));
        h5create(OutFile, '/time_series', size(TSERIES));
        h5create(OutFile, '/radius', size(X));
        h5create(OutFile, '/time', size(T));

        h5write(OutFile, '/radial_time_series', RDATA);
        h5write(OutFile, '/time_series', TSERIES);
        h5write(OutFile, '/radius', X);
        h5write(OutFile, '/time', T);
      end
    end

  end
end 
