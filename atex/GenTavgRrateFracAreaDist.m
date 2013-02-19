function [ ] = GenTavgRrateFracAreaDist(ConfigFile)
% GenTavgRrateFracAreaDist generate time average of the rain rate (pcprr) histogram data as fractional area of domain

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);

Tdir = Config.TsavgDir;
OutDir = Config.DiagDir;

RrateHistVar = 'hist_pcprr';

% Create the output directory if it doesn't exist
if (exist(OutDir,'dir') == 0)
  mkdir(OutDir);
end

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  RrateHistFile = sprintf('%s/%s_%s.h5', Tdir, RrateHistVar, Case);
  OutFile = sprintf('%s/tavg_%s_%s.h5', OutDir, RrateHistVar, Case);

  fprintf('***************************************************************\n');
  fprintf('Generating time averaged rain rate fractional area distribution:\n');
  fprintf('  Case: %s\n', Case);
  fprintf('  Input rain rain histogram file: %s\n', RrateHistFile);
  fprintf('  Output file: %s\n', OutFile);

  % The histogram is organized as (x,y,z,t) with x being the bins, y and z being dummy
  % coordinates (size = 1), and t the time steps. Use squeeze to compress into a 2D array
  % (x,t) that has bins and time steps.
  RRHIST = squeeze(hdf5read(RrateHistFile, RrateHistVar));
  BINS = hdf5read(RrateHistFile, 'x_coords');
  Nx = hdf5read(RrateHistFile, 'Nx');
  Ny = hdf5read(RrateHistFile, 'Ny');

  % Calculate the total number of points in the domain (2D --> Nx * Ny) multiplied
  % by the total number of time steps.
  [ Nb, Nt ] = size(RRHIST);
  Npts = Nx * Ny * Nt; % number of domain points after summing up over all time points
  fprintf('  Number of domain points: %d\n', Nx*Ny);
  fprintf('  Number of time points: %d\n', Nt);
  fprintf('\n');

  AVGDIST = sum(RRHIST,2) ./ Npts;

  hdf5write(OutFile, '/AvgDist', AVGDIST);
  hdf5write(OutFile, '/Bins', BINS, 'WriteMode', 'append');
  hdf5write(OutFile, '/Npts', Npts, 'WriteMode', 'append');
end
