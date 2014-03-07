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
    
Ddir = Config.DiagDir;

% Generate two cases:
%   No filter
%   Filter based on all liquid mix ratios >= 0.01 g/kg
CloudList = {
  { 'cloud_M1'       'cloud' 'cloud_cm3_M1'       'cloud_concen_cm3' }
  { 'cloud_M1_c0p01' 'cloud' 'cloud_cm3_M1_c0p01' 'cloud_concen_cm3' }
%  { 'cloud_M1_c0p10' 'cloud' 'cloud_cm3_M1_c0p10' 'cloud_concen_cm3' }
  };

DrizzleList = {
  { 'drizzle_M1'         'drizzle'  'drizzle_cm3_M1'          'drizzle_concen_cm3'   }
  { 'drizzle_M1_d0p010'  'drizzle'  'drizzle_cm3_M1_d0p010'   'drizzle_concen_cm3'   }
%  { 'drizzle_M1_d0p001'  'drizzle'  'drizzle_cm3_M1_d0p001'   'drizzle_concen_cm3'   }
  };

RainList = {
  { 'rain_M1'        'rain'  'rain_m3_M1'         'rain_concen_m3'   }
  { 'rain_M1_r0p01'  'rain'  'rain_m3_M1_r0p01'   'rain_concen_m3'   }
%  { 'rain_M1_r0p10'  'rain'  'rain_m3_M1_r0p10'   'rain_concen_m3'   }
  };

OutFprefixList = {
  'dsd'
  'dsd_0p01'
%  'dsd_0p10'
  };

NumPtsVar = 'num_points';

% RAMS reference state (DN01D) from the ATEX runs, kg/m3
AirDensity = [
  1.186
  1.177
  1.167
  1.158
  1.148
  1.139
  1.129
  1.120
  1.110
  1.100
  1.090
  1.080
  1.070
  1.061
  1.051
  1.036
  1.017
  0.998
  0.988
  0.979
  0.970
  0.960
  0.951
  0.942
  0.933
  0.924
  0.915
  0.907
  0.898
  0.889
  0.881
  0.872
  0.864
  0.855
  0.847
  0.839
  0.831
  0.822
  0.814
  0.806
  0.798
  ];

% Equally spaced log(D), D in m
LogSizes = -7:0.1:-2;  % 0.1um to 10mm
D = 10 .^ (LogSizes); % m

% Liquid water, assume spherical drops:
Alpha = 524;  % kg/m3
Beta = 3;     % dimensionless

% Shape parameter for gamma distributions
CloudNu = 4;    % dimensionless
DrizzleNu = 4;   % dimensionless
RainNu = 2;    % dimensionless

  for icase = 1:length(Config.Cases)
    for ivar = 1:length(CloudList)
      Case = Config.Cases(icase).Cname;
  
      RcFile = sprintf('%s/%s_%s.h5', Ddir, CloudList{ivar}{1}, Case);
      RcVar  = CloudList{ivar}{2};
      NcFile = sprintf('%s/%s_%s.h5', Ddir, CloudList{ivar}{3}, Case);
      NcVar  = CloudList{ivar}{4};

      RdFile = sprintf('%s/%s_%s.h5', Ddir, DrizzleList{ivar}{1}, Case);
      RdVar  = DrizzleList{ivar}{2};
      NdFile = sprintf('%s/%s_%s.h5', Ddir, DrizzleList{ivar}{3}, Case);
      NdVar  = DrizzleList{ivar}{4};

      RrFile = sprintf('%s/%s_%s.h5', Ddir, RainList{ivar}{1}, Case);
      RrVar  = RainList{ivar}{2};
      NrFile = sprintf('%s/%s_%s.h5', Ddir, RainList{ivar}{3}, Case);
      NrVar  = RainList{ivar}{4};
      
      OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefixList{ivar}, Case);
  
      fprintf('***************************************************************\n');
      fprintf('Generating DSD:\n');
      fprintf('  Case: %s\n', Case);
      fprintf('  Input cloud mixing ratio file: %s\n', RcFile);
      fprintf('    Variable name: %s\n', RcVar);
      fprintf('  Input cloud number concentration file: %s\n', NcFile);
      fprintf('    Variable name: %s\n', NcVar);
      fprintf('  Input drizzle mixing ratio file: %s\n', RdFile);
      fprintf('    Variable name: %s\n', RdVar);
      fprintf('  Input drizzle number concentration file: %s\n', NdFile);
      fprintf('    Variable name: %s\n', NdVar);
      fprintf('  Input rain mixing ratio file: %s\n', RrFile);
      fprintf('    Variable name: %s\n', RrVar);
      fprintf('  Input rain number concentration file: %s\n', NrFile);
      fprintf('    Variable name: %s\n', NrVar);
      
      fprintf('  Output file: %s\n', OutFile);
      fprintf('\n');
  
      % Cloud mixing ratio is in g/kg, convert to kg/kg -> multiply by 1e-3
      % Cloud number concentration is in #/cm3, convert to #/m3 -> multiply by 1e6
      %
      % Drizzle mixing ratio is in g/kg, convert to kg/kg -> multiply by 1e-3
      % Drizzle number concentration is in #/cm3, convert to #/m3 -> multiply by 1e6
      %
      % Rain mixing ratio is in g/kg, convert to kg/kg -> multiply by 1e-3
      % Rain number concentration is in #/m3
      %
      % In places where the number of points is zero, the sum will also be zero.
      % Dividing sum by number of points in these cases will do 0/0 which will
      % produce a nan. Change the nans back to zeros since want zero result for
      % zero sums.
      %
      RC_SUMS = squeeze(hdf5read(RcFile, RcVar)) .* 1e-3;
      RC_NPTS = squeeze(hdf5read(RcFile, NumPtsVar));
      RC = RC_SUMS ./ RC_NPTS;
      RC(isnan(RC)) = 0;
  
      NC_SUMS = squeeze(hdf5read(NcFile, NcVar)) .* 1e6;
      NC_NPTS = squeeze(hdf5read(NcFile, NumPtsVar));
      NC = NC_SUMS ./ NC_NPTS;
      NC(isnan(NC)) = 0;
  
      RD_SUMS = squeeze(hdf5read(RdFile, RdVar)) .* 1e-3;
      RD_NPTS = squeeze(hdf5read(RdFile, NumPtsVar));
      RD = RD_SUMS ./ RD_NPTS;
      RD(isnan(RD)) = 0;
  
      ND_SUMS = squeeze(hdf5read(NdFile, NdVar)) .* 1e6;
      ND_NPTS = squeeze(hdf5read(NdFile, NumPtsVar));
      ND = ND_SUMS ./ ND_NPTS;
      ND(isnan(ND)) = 0;
  
      RR_SUMS = squeeze(hdf5read(RrFile, RrVar)) .* 1e-3;
      RR_NPTS = squeeze(hdf5read(RrFile, NumPtsVar));
      RR = RR_SUMS ./ RR_NPTS;
      RR(isnan(RR)) = 0;
  
      NR_SUMS = squeeze(hdf5read(NrFile, NrVar)) .* 1e6;
      NR_NPTS = squeeze(hdf5read(NrFile, NumPtsVar));
      NR = NR_SUMS ./ NR_NPTS;
      NR(isnan(NR)) = 0;

      % coordinates for output
      X = [ 1 ];
      Y = D .* 1e6;  % output diameters in microns
      Z = hdf5read(RcFile, 'z_coords');
      T = hdf5read(RcFile, 't_coords');

      % Mixing ratios and number concentrations are organized as (z,t)
      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % Mixing ratios are in #/m3, but want #/kg for the RamsGammaNdist() function. Need to
      % multiply by the density of air.
      RhoAir = repmat(AirDensity, [ 1 Nt ]);
      RC = RC .* RhoAir;
      RD = RD .* RhoAir;
      RR = RR .* RhoAir;

      % Form individual number distributions, then add them together to form a single distribution
      % with all components. Then used the combined number distribution to form the volume
      % distribution.
      [ CloudNumDist CloudCharDiam ] = RamsGammaNdist(Alpha, Beta, CloudNu, RC, NC, D);
      [ DrizzleNumDist DrizzleCharDiam ] = RamsGammaNdist(Alpha, Beta, DrizzleNu, RD, ND, D);
      [ RainNumDist RainCharDiam ]   = RamsGammaNdist(Alpha, Beta, RainNu,  RR, NR, D);

      % Output distribution from RamsGammaNdist will have units m-1 * m-3, covert to um-1 * cm-3 by multiplying by 1e-12
      CloudNumDist   = CloudNumDist .* 1e-12;
      DrizzleNumDist = DrizzleNumDist .* 1e-12;
      RainNumDist    = RainNumDist .* 1e-12;
      NumDist = (CloudNumDist + DrizzleNumDist + RainNumDist);
  
      % Form volume, make sure to use diameters in microns (Y).
      % Y is organized as (d) and we want to repeat it so that it is organized as (d,z,t)
      % in order to match up with NumDist. Cannot use a single repmat call since that will
      % simply tile Y into a 2D array. Instead call successive repmats to form a 3D array.
      % The inner repmat creates (z,d) array, the outer repmat creates a (z,d,t) array, and
      % the permute rearranges into a (d,z,t) array. This depends on Y being a row vector.
      DiamsExp = permute(repmat(repmat(Y, [ Nz 1 ]), [ 1 1 Nt ]), [ 2 1 3 ]);
      VolDist = (pi/6) .* (DiamsExp.^3) .* NumDist; % units are new um2 * cm-3

      % Save the diameters (in Y) and characteristic diameters (convert m to um)
      OutVar = reshape(NumDist, [ Nx Ny Nz Nt ]);
      hdf5write(OutFile, 'NumDist', OutVar);

      OutVar = reshape(VolDist, [ Nx Ny Nz Nt ]);
      hdf5write(OutFile, 'VolDist', OutVar, 'WriteMode', 'append');

      OutVar = reshape(CloudNumDist, [ Nx Ny Nz Nt ]);
      hdf5write(OutFile, 'CloudNumDist', OutVar, 'WriteMode', 'append');
      OutVar = reshape(DrizzleNumDist, [ Nx Ny Nz Nt ]);
      hdf5write(OutFile, 'DrizzleNumDist', OutVar, 'WriteMode', 'append');
      OutVar = reshape(RainNumDist, [ Nx Ny Nz Nt ]);
      hdf5write(OutFile, 'RainNumDist', OutVar, 'WriteMode', 'append');

      OutVar = reshape(CloudCharDiam.*1e6, [ 1 1 Nz Nt ]);
      hdf5write(OutFile, 'CloudCharDiam', OutVar, 'WriteMode', 'append');
      OutVar = reshape(DrizzleCharDiam.*1e6, [ 1 1 Nz Nt ]);
      hdf5write(OutFile, 'DrizzleCharDiam', OutVar, 'WriteMode', 'append');
      OutVar = reshape(RainCharDiam.*1e6, [ 1 1 Nz Nt ]);
      hdf5write(OutFile, 'RainCharDiam', OutVar, 'WriteMode', 'append');

      hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
      hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
      hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
      hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
      fprintf('\n');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RamsGammaNdist()
%
% This function returns the gamma based number concentration distribution
% given the two prognosed moments (mixing ratio and total number
% concentration) plus the gamma shape parameter nu.
%
% There are also two parameters that come from power law droplet mass:
% alpha (mass coefficient) and beta (power of diameter).
%
% Alpha is the mass coefficient (kg/m3)
% Beta is the mass power (dimensionless)
% Nu is shape parameter of gamma distribution (dimensionless)
%
% MixRat is the mixing ratio (kg/kg)
% TotNum is the total number concentration (#/kg)
% Diams is a vector of diameters on which to calculate the distribution (m)
%
% NumDist is the resutling vector holding the number distribution
% VolDist is the resutling vector holding the volume distribution
%
function [ NumDist CharDiam ] = RamsGammaNdist(Alpha, Beta, Nu, MixRat, TotNum, Diams)

  % MixRat and TotNum have dimensions (z,t)
  % Want to add a third dimension which is the diameters (Diams) used for the distribution.
  % Put the new Diams dimension first.

  % RAMS uses a characteristic diameter (CharDiam) that is diagnosed from the two moments
  % of the distribution.
  %
  %   CharDiam = [ (MixRat * gamma(Nu)) / (TotNum * Alpha * gamma(Nu + Beta)) ] ^(1/Beta)
  %
  % Once CharDiam is determined then NumDist can be formed by:
  %
  %   NumDist = (TotNum / gamma(Nu)) * (Diams/CharDiam)^(Nu-1) * (1/CharDiam) * exp(-Diams/CharDiam)
  %

  % characteristic diameter, units m
  GammaNu = gamma(Nu);
  GammaNuBeta = gamma(Nu+Beta);
  CharDiam = ((MixRat .* GammaNu) ./ (TotNum .* Alpha .* GammaNuBeta)) .^ (1./Beta);
  CharDiam(isnan(CharDiam)) = 0;

  % At this point CharDiam and TotNum are (z,t); Diams is (d); GammaNu is scalar.
  % Repeat matrices so that new dimensions are (d,z,t). Cannot do repmat(X, [ Nd 1 1 ]) since
  % this will just tile X into a new 2D array. Instead place the Diams dimension last -> (z,t,d),
  % form the distributions (in the last dimension), the permute the results so that the
  % output has the form (d,z,t).
  Nd = length(Diams);
  [ Nz Nt ] = size(TotNum);

  CharDiamExp = repmat(CharDiam, [ 1 1 Nd ]);
  TotNumExp   = repmat(TotNum, [ 1 1 Nd ]);
  % Diams is organized as (d), need to add z and t dimensions one at a time using repmat since repmat
  % will just tile when z and t are added simutaneously. The following code depends on Diams being
  % a row vector. The inner repmat creates a (z,d) array, the outer repmat forms a (z,d,t) array and
  % the permute rearranges the arrat to (z,t,d) so that it matches up with CharDiamExp and TotNumExp.
  DiamsExp = permute(repmat(repmat(Diams,[ Nz 1 ]), [ 1 1 Nt ]), [ 1 3 2 ]);

  % gamma number distribution, units: m-1 * m-3
  Dnorm = DiamsExp ./ CharDiamExp;
  Dnorm(isnan(Dnorm)) = 0;
  NumDist = (TotNumExp ./ GammaNu) .* ( Dnorm .^(Nu-1)) .* (1 ./ CharDiamExp) .* exp(-(Dnorm)); % organized as (z,t,d)
  NumDist(isnan(NumDist)) = 0;

  NumDist = permute(NumDist, [ 3 1 2 ]); % organize as (d,z,t)
end
