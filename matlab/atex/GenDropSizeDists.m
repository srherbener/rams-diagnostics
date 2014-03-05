function [ ] = GenDropSizeDists(ConfigFile)
% GenDropSizeDists generate drop size distributions from RAMS Gamma pdfs

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

Cfprefix = 'cloud_xxx';
Cvar = 'cloud_xxx';
Rfprefix = 'rain_xxx';
Rvar = 'rain_xxx';

OutFprefix = 'dsd_xxx';

% for icase = 1:length(Config.Cases)
%     Case = Config.Cases(icase).Cname;
%     Cfile = sprintf('%s/%s_%s.h5', Tdir, Cfprefix, Case);
%     Rfile = sprintf('%s/%s_%s.h5', Tdir, Rfprefix, Case);
%     
%     OutFile = sprintf('%s/%s_%s.h5', Tdir, OutFprefix, Case);
% 
%     fprintf('***************************************************************\n');
%     fprintf('Generating time series of average cloud depth:\n');
%     fprintf('  Case: %s\n', Case);
%     fprintf('  Input cloud file: %s\n', Cfile);
%     fprintf('    Variable name: %s\n', Cvar);
%     fprintf('  Input rain file: %s\n', Rfile);
%     fprintf('    Variable name: %s\n', Rvar);
%     
%     fprintf('  Output file: %s\n', OutFile);
%     fprintf('\n');
%
%     % Get the cloud depth
%     CD  = squeeze(hdf5read(InFile, InVarName));
% 
%     % copy the coordinate values to the output to keep GenTsPlots happy
%     X = hdf5read(InFile, 'x_coords');
%     Y = hdf5read(InFile, 'y_coords');
%     Z = hdf5read(InFile, 'z_coords');
%     T = hdf5read(InFile, 't_coords');
% 
%     Nx = length(X);
%     Ny = length(Y);
%     Nz = length(Z);
%     Nt = length(T);
% 
%     % toss out the zeros - first bin --> want to do exclude clear sky
%     CD = CD(2:end,:);
%     X = X(2:end,:);
% 
%     % Do a weighted mean of the histogram data for the average value
%     % Make the variable appear as if it came from the tsavg diagnostic
%     % Replace nans with zeros - nans are the times when there is all
%     % open sky in the domain so zero is appropriate.
%     CD_AVG = ReduceHists(CD, 1, X, 'wtmean');
%     CD_AVG(isnan(CD_AVG)) = 0;
%     CD_AVG = reshape(CD_AVG, [ 1 1 1 Nt ]);
% 
%     hdf5write(OutFile, OutVarName, CD_AVG);
%     hdf5write(OutFile, 'x_coords', [ 1 ], 'WriteMode', 'append');
%     hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
%     hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
%     hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
% 
% end

LogSizes = -7:0.1:-2;
D = 10 .^ (LogSizes);

r = 0.1;      % kg/kg
Nt = 4.1e8;   % #/m3
Alpha = 524;  % kg/m3
Beta = 3;
Nu = 4; % Cloud
%Nu = 4; % Drizzle
%Nu = 2; % Rain

[ nD vD ] = RamsGammaDists(Alpha, Beta, Nu, r, Nt, D);


figure;
semilogx(D, nD);
figure;
semilogx(D, vD);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RamsGammaDists()
%
% This function returns the gamma based number and volume distributions given
% the two prognosed moments (mixing ratio and total number concentration).
%
% There are also two parameters that come from power law droplet mass:
% alpha (mass coefficient) and beta(power of diameter).
%
% Alpha is the mass coefficient (kg/m3)
% Beta is the mass power (dimensionless)
% Nu is shape parameter of gamma distribution (dimensionless)
% R is the mixing ratio (kg/kg)
% Nt is the total number concentration (#/m3)
% D is a vector of diameters on which to calculate the distribution (m)
%
% NumDist is the resutling vector holding the number distribution
% VolDist is the resutling vector holding the volume distribution
%
function [ NumDist VolDist  ] = RamsGammaDists(Alpha, Beta, Nu, R, Nt, D)

  % RAMS uses a characteristic diameter (Dn) that is diagnosed from the two moments
  % of the distribution.
  %
  %   Dn = [ (R * gamma(Nu) * density of air) / (Nt * Alpha * gamma(Nu + Beta)) ] ^(1/Beta)
  %
  % Once Dn is determined then NumDist can be formed by:
  %
  %   NumDist = (Nt / gamma(Nu)) * (D/Dn)^(Nu-1) * (1/Dn) * exp(-D/Dn)
  %
  % Once NumDist is determined, VolDist can be formed by:
  %
  %   VolDist = (pi/6) * D^3 * NumDist

  RhoAir = 1000; % 1000 kg/m3

  % characteristi diameter
  GammaNu = gamma(Nu);
  GammaNuBeta = gamma(Nu+Beta);
  Dn = ((R * GammaNu * RhoAir) / (Nt * Alpha * GammaNuBeta)) ^ (1/Beta);
  
  % gamma number distribution, units: 1/m * 1/m2
  Dnorm = D ./ Dn;
  NumDist = (Nt ./ GammaNu) .* ( Dnorm .^(Nu-1)) .* (1 ./ Dn) .* exp(-(Dnorm));

  % create volume distribution: m2 * 1/m2
  VolDist = (pi/6) .* (D.^3) .* NumDist;
    
end
