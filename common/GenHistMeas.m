function [ ] = GenHistMeas( ConfigFile )
% GenHistMeas generate various measurments from histogram data
%

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

DiagDir = Config.DiagDir;

% Make sure output directory exists
if (exist(DiagDir, 'dir') ~= 7)
    mkdir(DiagDir);
end

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  for ihmeas = 1: length(Config.Hmeas)
    Name = Config.Hmeas(ihmeas).Name;
    InDir = Config.Hmeas(ihmeas).InDir;
    Fprefix = Config.Hmeas(ihmeas).Fprefix;
    Vname = Config.Hmeas(ihmeas).Rvar;
    Method = Config.Hmeas(ihmeas).Method;

    fprintf('***********************************************************************\n');
    fprintf('Generating Histogram Measurements:\n');
    fprintf('  Method: %s\n', Method);
    fprintf('  Name: %s\n', Name);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    InFile  = sprintf('%s/%s_%s.h5', InDir, Fprefix, Case);
    OutFile = sprintf('%s/%s_%s.h5', DiagDir, Name, Case);
    Hdset   = sprintf('/%s', Vname);

    % Read in the histogram data. HIST will be organized as (r,b,z,t) where
    %    r --> radius
    %    b --> histogram bins
    %    z --> heights
    %    t --> time
    fprintf('Reading file: %s\n', InFile);
    fprintf('\n');
    R = hdf5read(InFile, '/x_coords');
    B = hdf5read(InFile, '/y_coords');
    Z = hdf5read(InFile, '/z_coords');
    T = hdf5read(InFile, '/t_coords');
    HIST = hdf5read(InFile, Hdset);

    % ReduceHists( Hdata, Hdim, Bins, Method)
    %   Method: 1 - weight mean
    %           2 - max
    %           3 - center of mass
    switch Method
      case 'wtmean'
        [ HMEAS ] = ReduceHists(HIST, 2, B, 1);
      case 'max'
        [ HMEAS ]   = ReduceHists(HIST, 2, B, 2);
      case 'com'
        [ HMEAS ]   = ReduceHists(HIST, 2, B, 3);
      otherwise
        fprintf('WARNING: Unrecongnized measurement method: %s\n', Method);
        fprintf('WARNING:   skipping this case\n');
        continue;
    end

    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');

    Hdset = sprintf('/%s', Vname);
    hdf5write(OutFile, Hdset, HMEAS);

    hdf5write(OutFile, '/x_coords', R, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', B, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Z, 'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
