function [ ] = GenAvgCloudDepth(ConfigFile)
% GenAvgCloudDepth generate time series of the average cloud depth

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);

% Assume that you have histogram data for cloud depth
% Data are organized as (x,y,z,t) where
%   x - depth bins
%   y - dummy
%   z - dummy
%   t - time
    
Tdir = Config.TsavgDir;
InFilePrefix = 'hist_cdepth';
OutFilePrefix = 'horiz_avg_cdepth';
InVarName = 'hist_cloud_depth';
OutVarName = 'havg_cloud_depth';

for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;
    InFile = sprintf('%s/%s_%s.h5', Tdir, InFilePrefix, Case);
    OutFile = sprintf('%s/%s_%s.h5', Tdir, OutFilePrefix, Case);

    fprintf('***************************************************************\n');
    fprintf('Generating time series of average cloud depth:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('  Input file: %s\n', InFile);
    fprintf('  Output file: %s\n', OutFile);
    fprintf('  Input variable name: %s\n', InVarName);
    fprintf('  Output variable name: %s\n', OutVarName);
    fprintf('\n');

    % Get the cloud depth
    CD  = squeeze(hdf5read(InFile, InVarName));

    % copy the coordinate values to the output to keep GenTsPlots happy
    X = hdf5read(InFile, 'x_coords');
    Y = hdf5read(InFile, 'y_coords');
    Z = hdf5read(InFile, 'z_coords');
    T = hdf5read(InFile, 't_coords');

    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    Nt = length(T);

    % toss out the zeros - first bin --> want to do exclude clear sky
    CD = CD(2:end,:);
    X = X(2:end,:);

    % Do a weighted mean of the histogram data for the average value
    % Make the variable appear as if it came from the tsavg diagnostic
    % Replace nans with zeros - nans are the times when there is all
    % open sky in the domain so zero is appropriate.
    CD_AVG = ReduceHists(CD, 1, X, 'wtmean');
    CD_AVG(isnan(CD_AVG)) = 0;
    CD_AVG = reshape(CD_AVG, [ 1 1 1 Nt ]);

    hdf5write(OutFile, OutVarName, CD_AVG);
    hdf5write(OutFile, 'x_coords', [ 1 ], 'WriteMode', 'append');
    hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
    hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
    hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');

end
