function [ ] = GenDeltapTseries( ConfigFile )
% GenDeltapTseries generate time series of a pressure difference
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
  for i_vint_ts = 1: length(Config.DeltapTs)
    Name    = Config.DeltapTs(i_vint_ts).Name;
    InDir   = Config.DeltapTs(i_vint_ts).InDir;
    Fprefix = Config.DeltapTs(i_vint_ts).Fprefix;
    Vname   = Config.DeltapTs(i_vint_ts).Rvar;
    Rmin    = Config.DeltapTs(i_vint_ts).Rmin;
    Rmax    = Config.DeltapTs(i_vint_ts).Rmax;

    fprintf('***********************************************************************\n');
    fprintf('Generating Delta Pressure Time Series:\n');
    fprintf('  Name: %s\n', Name);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('  Rmin: %d\n', Rmin);
    fprintf('  Rmax: %d\n', Rmax);
    fprintf('\n');

    InFile  = sprintf('%s/%s_%s.h5', InDir, Fprefix, Case);
    OutFile = sprintf('%s/%s_%s.h5', DiagDir, Name, Case);
    Hdset   = sprintf('/%s', Vname);

    % Read in the variable data - organized as (r,y,z,t) where y and z are dummy dimensions
    % Keep the data in this form so that LinePlot can be used to make the figures.
    fprintf('Reading file: %s\n', InFile);
    fprintf('\n');
    R = hdf5read(InFile, '/x_coords');
    T = hdf5read(InFile, '/t_coords');
    HDATA = hdf5read(InFile, Hdset);

    % Find the indices for Rmin and Rmax
    R1 = find(R >= Rmin, 1, 'first');
    R2 = find(R <= Rmax, 1, 'last');

    % Find the pressure difference
    DELTAP = HDATA(R2,:,:,:) - HDATA(R1,:,:,:);

    % For selection in line plots, we want X,Y,Z (dummy dimensions) to all
    % be size = 1, and value = 1.
    X(1) = 1;
    Y(1) = 1;
    Z(1) = 1; 

    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');

    Hdset = sprintf('/DeltapTs_%s', regexprep(Vname, 'ProfTs_', ''));
    hdf5write(OutFile, Hdset, DELTAP);

    hdf5write(OutFile, '/x_coords', X, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', Y, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Z, 'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
