function [ ] = GenVintTseries( ConfigFile )
% GenVintTseries generate time series of a vertically integrated quantity
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
  for i_vint_ts = 1: length(Config.VintTs)
    Name    = Config.VintTs(i_vint_ts).Name;
    InDir   = Config.VintTs(i_vint_ts).InDir;
    Fprefix = Config.VintTs(i_vint_ts).Fprefix;
    Vname   = Config.VintTs(i_vint_ts).Rvar;

    fprintf('***********************************************************************\n');
    fprintf('Generating Vertically Integrated Time Series:\n');
    fprintf('  Name: %s\n', Name);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    InFile  = sprintf('%s/%s_%s.h5', InDir, Fprefix, Case);
    OutFile = sprintf('%s/%s_%s.h5', DiagDir, Name, Case);
    Hdset   = sprintf('/%s', Vname);

    % Read in the variable data - organized as (z,t)
    fprintf('Reading file: %s\n', InFile);
    fprintf('\n');
    Z = hdf5read(InFile, '/Height');
    T = hdf5read(InFile, '/Time');
    HDATA = hdf5read(InFile, Hdset);

    % The first model level is below the surface so don't use it.
    % For each vertical integration
    %   sum up HDATA * DELTAZ
    %     go from Z(2) to Z(n) for summing
    %     use DELTAZ = Z(n) - Z(n-1)
    %
    % First create an array with the deltaz values replicated for
    % each column so that it can be used in an element-wise multiply
    % with the input data.
    [ Nz Nt ] = size(HDATA);
    DELTAZ = Z(2:end) - Z(1:end-1);
    DELTAZ = repmat(DELTAZ, [ 1 Nt ]);
    VINT = nansum((HDATA(2:end,:) .* DELTAZ), 1);

    % For GenLinePlots, reshape into 4D array -> (x,y,z,t)
    X(1) = 1;
    Y(1) = 1;
    OUT_Z(1) = 1;
    VINT = reshape(VINT, [ 1 1 1 Nt ]);

    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');

    Hdset = sprintf('/VintTs_%s', regexprep(Vname, 'ProfTs_', ''));
    hdf5write(OutFile, Hdset, VINT);

    hdf5write(OutFile, '/x_coords', X, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', Y, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', OUT_Z, 'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
