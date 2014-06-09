function [ U, V, T, Zg, RH ] = GenInitDpFiles(Case)
% GenInitDpFiles generate initial conditions in a dp file format for RAMS

  % Initialization is based on techniques described in
  %   Boutle et al., 2010
  %   Boutle et al., 2011
  %   Polvani and Esler, 2007
  %
  % Case selects between the two ETC cases described by Thorncroft et al., 1993.
  %   1 --> LC1  (anticyclonic lifecycle)
  %   2 --> LC2  (cyclonic lifecycle)
  %
  % The basic idea is to create a temperature gradient (north-south) with
  % a westerly jet in order to mimic the baroclinic zone along the Northern
  % Hemispher mid-latitude jet stream. Data is homogeneous in the zonal direction.
  %
  % The file format is described in CSU-RAMS manual:
  %
  %   "CSU-RAMS: Standard Input Formats for Pressure Coordinate and Observational Data (Data Prep Stage)" 
  %
  % written by Stephen M. Saleeby, Feb 12, 2014
  %
  %
  % There are 5 atmospheric variables:
  %    u velocity  (zonal, m/s)
  %    v velocity  (meridonal, m/s)
  %    temperature (K)
  %    geopotential height (m)
  %    relative humidity (fraction)
  %
  % and 6 soil variables:
  %    top soil level volumetric soil moisture (m^3/m^3)
  %    next level down volumetric soil moisture (m^3/m^3)
  %    top soil level soil temperature (K)
  %    next level down soil temperature (K)
  %    snow water equivalent (kg/m^2)
  %    snow depth (m)

  % Emulating a winter-time storm: use Jan 1, 2014, 12Z as the start time
  Year = 2014;
  Month = 1;
  Day = 1;
  Hour = 1200;

  % Grid:
  %   one degree in lat and lon
  %   lat is -90 to 90
  %   lon is -180 to 179
  %   vertical is pressure, 26 levels, 1000mb through 10mb

  LAT = -90:90;
  LON = -180:179;
  PRESS = [ 1000 975 950 925 900 850 800 750 700 650 600 550 500 450 400 350 300 250 200 150 100 70 50 30 20 10 ];

  Nx = length(LON);
  Ny = length(LAT);
  Nz = length(PRESS);

  % Generate winds from analytic formulas
  [ U V ] = GenInitWindField(LON, LAT, PRESS, Case);

  T  = zeros([ Nx Ny Nz ]);
  Zg = zeros([ Nx Ny Nz ]);
  RH = zeros([ Nx Ny Nz ]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ U V ] = GenInitWindField(Lon, Lat, Press, Case)
% GenWindField generate initial wind field

  % Parameters from Polvani and Esler, 2007
  U0 = 45;      % max jet speed (m/s)
  Zt = 13000;   % temperature scale height (m)
  H  = 7500;    % atmospheric scale height (m)
  P0 = 1000;    % reference pressure (mb)

  Us = 45;      % m/s
  Phi_s = 35;   % degrees
  Delta_s = 20; % degrees
  Zs = 10000;   % m

  % Use formulas from Polvani and Esler, 2007 for the 2 Thorncroft et al., 1993 ETC cases
  %   Case == 1 --> LC1
  %   Case == 2 --> LC2
  %
  % Formulas use log-pressure height, z = H * ln(P0/Press)
  %
  % LC1
  %   U1 = U0 * F * Z/Zt * exp (- ((Z/Zt)^2 -1)/2)
  %         F = sin^3[pi * sin^2(lat)] for lat > 0
  %           = 0 for lat <= 0
  %   V1 = 0;
  %
  % LC2
  %   U2 = U1 + US
  %   where US = Us * exp(-Z/Zs) * sin^2((lat-lat_s)/delta_s) * exp (- ((lat-lat_s)/delta_s)^2)
  %     US is a meridonal shear term
  %   V2 = 0;

  Nx = length(Lon);
  Ny = length(Lat);
  Nz = length(Press);

  U = zeros([ Nx Ny Nz ]);
  V = zeros([ Nx Ny Nz ]);

  % Create log-pressure heights -> Z, and scaled heights -> Zscale = Z/Zt
  Z = H .* log(P0 ./ Press);
  Zscale = repmat((Z ./ Zt), [ Ny 1 ]);

  % Create latitude in radians
  LatRad = (pi / 180) .* Lat;

  % Create lat dependent factor of Thorncroft et al. formula for U1
  F = zeros([ Ny 1 ]);
  for i = 1:Ny
    if (Lat(i) <= 0)
      F(i) = 0;
    else
      F(i) = sin(pi * sin(LatRad(i))^2)^3;
    end
  end
  F = repmat(F, [ 1 Nz ]);

  % form 2D wind field in lat and z directions
  % then repeat for all lon
  U_2D = U0 .* F .* Zscale .* exp(-((Zscale).^2 - 1)./2);
  U = repmat(reshape(U_2D, [ 1 Ny Nz ]), [ Nx 1 1 ]);

  % If doing LC2, then add in the meridonal shear term
  if (Case == 2)
    % Z scale
    Zscale = repmat((Z ./ Zs), [ Ny 1 ]);

    % Lat dependent factor of Thorncroft et al. formula for US
    LatScale = (Lat - Phi_s) ./ Delta_s;  % use lat in degrees
    F = sin(2.*LatRad).^2 .* LatScale .* exp(-(LatScale.^2));
    F = repmat(F', [ 1 Nz ]); % Note transpose on F inside repmat()
    
    % Create shear term in 2D, then replicate for all lon
    US_2D = -Us .* exp(-Zscale) .* F;
    US = repmat(reshape(US_2D, [ 1 Ny Nz]), [ Nx 1 1 ]);

    U = U + US;
  end
end
