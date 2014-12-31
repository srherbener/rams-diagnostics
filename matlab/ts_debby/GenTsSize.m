function [ ] = GenTsSize(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;
  Adir = Config. AzavgDir;

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;

    fprintf('*****************************************************************\n');
    fprintf('Generating time series of size metrics:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Read in max wind data, and generate RMX and Radius of 34kt (17.5 m/s) winds
    %InFile = sprintf('%s/max_wind_%s.h5', Ddir, Case);
    %InVar = '/radial_time_series';
    InFile = sprintf('%s/speed10m_%s.h5', Adir, Case);
    InVar = '/speed10m';

    fprintf('  Reading: %s (%s)\n', InFile, InVar);
    MWIND = squeeze(h5read(InFile, InVar));
    %R     = squeeze(h5read(InFile, '/radius')) ./ 1000; % km
    %T     = squeeze(h5read(InFile, '/time')) ./ 3600; % hr
    R     = squeeze(h5read(InFile, '/x_coords')) ./ 1000; % km
    T     = squeeze(h5read(InFile, '/t_coords')) ./ 3600; % hr

    Nr = length(R);
    Nt = length(T);

    % MWIND will be: (r,t)
    % For each time step, find the last occurrence of the maximum wind speed and
    % of the >= 34kt wind speed. Then convert the position of these to the
    % corresponding radius value.
    RMW = zeros([ Nt 1 ]);
    R34 = zeros([ Nt 1 ]);

    for it = 1:Nt
      WIND = squeeze(MWIND(:,it));

      Wmax = max(WIND);
      R1 = find(WIND == Wmax, 1, 'last');
      R2 = find(WIND >= 17.5, 1, 'last');

      RMW(it) = R(R1);
      if (isempty(R2))
        R34(it) = nan;
      else
        R34(it) = R(R2);
      end
    end

    % Write out measurement
    OutFile = sprintf('%s/ts_size_%s.h5', Ddir, Case);
    fprintf('    Writing: %s\n', OutFile)
    fprintf('\n');

    % If the file exists, remove it so that the HDF5 commands
    % can create a new file.
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end
    h5create(OutFile, '/rmw',  size(RMW));
    h5create(OutFile, '/r34',  size(R34));
    h5create(OutFile, '/time', size(T));

    h5write(OutFile, '/rmw',  RMW);
    h5write(OutFile, '/r34',  R34);
    h5write(OutFile, '/time', T);
  end
end 
