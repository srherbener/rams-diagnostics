function [ ] = GenTimeAvg( ConfigFile )
% GenTimeAvg generate time average of a variable
%

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

DiagDir = Config.DiagDir;
UndefVal = Config.UndefVal;

% Make sure output directory exists
if (exist(DiagDir, 'dir') ~= 7)
    mkdir(DiagDir);
end

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  for i_tavg = 1: length(Config.Tavg)
    Name    = Config.Tavg(i_tavg).Name;
    InDir   = Config.Tavg(i_tavg).InDir;
    Fprefix = Config.Tavg(i_tavg).Fprefix;
    Vname   = Config.Tavg(i_tavg).Rvar;
    Tmin    = Config.Tavg(i_tavg).Tmin;
    Tmax    = Config.Tavg(i_tavg).Tmax;

    fprintf('***********************************************************************\n');
    fprintf('Generating Time Average:\n');
    fprintf('  Name: %s\n', Name);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('  Tmin: %.1f\n', Tmin);
    fprintf('  Tmax: %.1f\n', Tmax);
    fprintf('\n');

    InFile  = sprintf('%s/%s_%s.h5', InDir, Fprefix, Case);
    OutFile = sprintf('%s/%s_%s.h5', DiagDir, Name, Case);
    Hdset   = sprintf('/%s', Vname);

    % Read in the variable data - organized as (x,y,z,t) where
    fprintf('Reading file: %s\n', InFile);
    fprintf('\n');
    X = hdf5read(InFile, '/x_coords');
    Y = hdf5read(InFile, '/y_coords');
    Z = hdf5read(InFile, '/z_coords');
    T = hdf5read(InFile, '/t_coords');
    HDATA = hdf5read(InFile, Hdset);

    % Find the time range
    T1 = find(T >= Tmin, 1, 'first');
    T2 = find(T <= Tmax, 1, 'last');

    % Convert undef values to nans so they can be excluded from averaging
    % Note: nanmean(<all_nans) --> nan
    HDATA(HDATA == UndefVal) = nan;
    TAVG = nanmean(HDATA(:,:,:,T1:T2), 4);

    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');

    Hdset = sprintf('/%s', Vname);
    hdf5write(OutFile, Hdset, TAVG);

    hdf5write(OutFile, '/x_coords', X, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', Y, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Z, 'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
