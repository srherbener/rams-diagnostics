function [ ] = GenFactAvgs(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;
  Tdir = Config.TsavgDir;

  % Ignore the first 10 hours of the time series since this contains the
  % end of the rapid intensification phases of the simulated storms.
  T1 = 10;

  % Put all of the results into one output file.
  % If the file exists, remove it so that the HDF5 commands
  % can create a new datasets.
  OutFile = sprintf('%s/avg_factors.h5', Ddir);
  if (exist(OutFile, 'file') == 2)
    delete(OutFile);
  end

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;

    fprintf('*****************************************************************\n');
    fprintf('Generating averages of metrics for factor separation analysis:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Create averages for: max wind, min pressure, IKE and RMW. These are all
    % time series found in various files.
    WindFile = sprintf('%s/max_speed10m_%s.h5', Tdir, Case);
    WindVar = '/max_speed10m';
    PressFile = sprintf('%s/min_sea_press_%s.h5', Tdir, Case);
    PressVar = '/min_sea_press';
    IkeFile = sprintf('%s/horiz_ke_%s.h5', Tdir, Case);
    IkeVar = '/horiz_ke';
    RmwFile = sprintf('%s/ts_size_%s.h5', Ddir, Case);
    RmwVar = '/rmw';

    fprintf('  Reading: %s (%s)\n', WindFile, WindVar);
    fprintf('  Reading: %s (%s)\n', PressFile, PressVar);
    fprintf('  Reading: %s (%s)\n', IkeFile, IkeVar);
    fprintf('  Reading: %s (%s)\n', RmwFile, RmwVar);
    fprintf('\n');

    % All data are one dimensional time series
    WIND  = squeeze(h5read(WindFile, WindVar));
    PRESS = squeeze(h5read(PressFile, PressVar));
    IKE   = squeeze(h5read(IkeFile, IkeVar));
    RMW   = squeeze(h5read(RmwFile, RmwVar));

    % Find average
    WIND_AVG = nanmean(WIND(T1:end));
    PRESS_AVG = nanmean(PRESS(T1:end));
    IKE_AVG = nanmean(IKE(T1:end));
    RMW_AVG = nanmean(RMW(T1:end));

    % Write out measurement
    fprintf('  Writing: %s\n', OutFile)

    OutVar = sprintf('/%s/avg_wind', Case);
    fprintf ('    %s\n', OutVar);
    h5create(OutFile, OutVar,  size(WIND_AVG));
    h5write( OutFile, OutVar,  WIND_AVG);

    OutVar = sprintf('/%s/avg_press', Case);
    fprintf ('    %s\n', OutVar);
    h5create(OutFile, OutVar, size(PRESS_AVG));
    h5write( OutFile, OutVar, PRESS_AVG);

    OutVar = sprintf('/%s/avg_ike', Case);
    fprintf ('    %s\n', OutVar);
    h5create(OutFile, OutVar,   size(IKE_AVG));
    h5write( OutFile, OutVar,   IKE_AVG);

    OutVar = sprintf('/%s/avg_rmw', Case);
    fprintf ('    %s\n', OutVar);
    h5create(OutFile, OutVar,   size(RMW_AVG));
    h5write( OutFile, OutVar,   RMW_AVG);

    fprintf('\n');
  end
end 
