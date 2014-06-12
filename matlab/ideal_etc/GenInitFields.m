function [ U, V, T, Zg, RH, LON, LAT, PRESS ] = GenInitFields(Case)
% GenInitFields generate initial conditions bases on Thorncroft et al, 1993 (LC1 and LC2)

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
  %    v velocity  (meridional, m/s)
  %    temperature (K)
  %    geopotential height (m)
  %    relative humidity (fraction)
  %

  % Parameters from Polvani and Esler, 2007
  U0 = 45;            % max jet speed (m/s)
  Zt = 13000;         % temperature scale height (m)
  Us = 45;            % m/s
  Phi_s = 35;         % degrees
  Delta_s = 20;       % degrees
  Zs = 10000;         % m
  T0 = 280;           % K, surface temp at 45N
  Gamma0 = -0.0065;   % K/m, roughly moist adiabatic lapse rate
  Alpha = 10;         % dimensionless
  R = 287;            % J/kg/K
  Omega = 7.297e-5;   % 1/s
  RadEarth = 6.371e6; % m
  G0 = 9.8;           % m/s^2

  H = 7500;  % m
  P0 = 1000; % mb
  

  % Grid:
  %   one degree in lat and lon
  %   lat is -90 to 90
  %   lon is -180 to 179
  %   vertical is pressure, 26 levels, 1000mb through 10mb

  LAT = -90:90;
  LON = -180:179;
  PRESS = [ 1000 975 950 925 900 850 800 750 700 650 600 550 500 450 400 350 300 250 200 150 100 70 50 30 20 10 ];

  
  % Create log-pressure heights
  Z = H .* log(P0 ./ PRESS);
  
  Nx = length(LON);
  Ny = length(LAT);
  Nz = length(PRESS);

  % First generate 2D fields (lat,p) for all variables, then replicate for all longitudes.

  % Generate zonal wind from analytic formulas. 
  [ U V ] = GenInitWindField(LAT, Z, Case, U0, Zt, Us, Phi_s, Delta_s, Zs);
  
  % Generate temperature field that is in thermal balance with U
  [ T ] = GenTempBalanced(LAT, Z, PRESS, U, T0, Gamma0, Zt, Alpha, R, Omega, RadEarth);

  % Generate initial moisture field
  [ RH ] = GenInitRelHumField(LAT, Z, Zt);

  % Generate geopotential heights
  [ Zg ] = GenGeoPotField(LAT, PRESS, T, R, G0);

  % Create homogeneous fields in the zonal direction
  U  = repmat(reshape(U,  [ 1 Ny Nz ]), [ Nx 1 1]);
  V  = repmat(reshape(V,  [ 1 Ny Nz ]), [ Nx 1 1]);
  T  = repmat(reshape(T,  [ 1 Ny Nz ]), [ Nx 1 1]);
  Zg = repmat(reshape(Zg, [ 1 Ny Nz ]), [ Nx 1 1]);
  RH = repmat(reshape(RH, [ 1 Ny Nz ]), [ Nx 1 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ U V ] = GenInitWindField(Lat, Z, Case, U0, Zt, Us, Phi_s, Delta_s, Zs)
% GenWindField generate initial wind field in 2d (lat,z)

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

  Ny = length(Lat);
  Nz = length(Z);

  % Create scaled heights -> Zscale = Z/Zt
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
  U = U0 .* F .* Zscale .* exp(-((Zscale).^2 - 1)./2);

  % If doing LC2, then add in the meridonal shear term
  if (Case == 2)
    % Z scale
    Zscale = repmat((Z ./ Zs), [ Ny 1 ]);

    % Lat dependent factor of Thorncroft et al. formula for US
    LatScale = (Lat - Phi_s) ./ Delta_s;  % use lat in degrees
    F = sin(2.*LatRad).^2 .* LatScale .* exp(-(LatScale.^2));
    F(LatRad < 0) = 0;
    F = repmat(F', [ 1 Nz ]); % Note transpose on F inside repmat()
    
    % Create shear term in 2D, then replicate for all lon
    US = -Us .* exp(-Zscale) .* F;

    U = U + US;
  end
  
  V = zeros([ Ny Nz ]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ T ] = GenTempBalanced(Lat, Z, Press, U, T0, Gamma0, Zt, Alpha, R, Omega, RadEarth)
% GenTempBalance - generate a temperature field that is thermally balanced with U
  
  Ny = length(Lat);
  Nz = length(Z);
  
  % Methods of creating balanced temperature field is from Polvani
  % and Esler, 1993, except doing this in pressure coordinates instead
  % of sigma (p/p0) coordinates.

  % Create latitude in radians
  LatRad = (pi / 180) .* Lat;

  
  % T = Tref + Tbal

  % Reference temperature profile, Tref
  Tref = T0 + (Gamma0 ./ ((Zt^-Alpha + Z.^-Alpha).^(1/Alpha)));
  Tref = repmat(Tref, [ Ny 1 ]);
  
  % Tbal is an integral formulated from balancing the initial wind
  % field with temperature using the primitive equations on a sphere.
  %
  %  Tbal = integral[ (p/R) * (a*2*omega*sin(lat) + 2*u*tan(lat)) * (du/dp) ] dlat
  %
  PR = repmat(Press, [ Ny 1 ]) ./ R;

  SinLat = repmat(reshape(sin(LatRad), [ Ny 1 ]), [ 1 Nz ]);
  TanLat = repmat(reshape(tan(LatRad), [ Ny 1 ]), [ 1 Nz ]);
  Geom = (RadEarth .* 2 .* Omega .* SinLat) + (2 .* U .* TanLat);
  
  % du/dp
  %  1) repeat top and bottom values in Press and U
  %  2) take averages between vertical entries in Press and U
  %  3) form du/dp by differences of averages from 2)
  %
  Ptemp = [ Press(1) Press Press(end) ];
  Utemp = [ U(:,1) U U(:,end) ];
  Pavg = (Ptemp(1:end-1) + Ptemp(2:end)) .* 0.5;
  Uavg = (Utemp(:,1:end-1) + Utemp(:,2:end)) .* 0.5;
  DuDp = (Uavg(:,2:end) - Uavg(:,1:end-1)) ./ repmat((Pavg(2:end) - Pavg(1:end-1)), [ Ny 1 ]);

  % Dlat
  Dlat = LatRad(2:end) - LatRad(1:end-1);
  Dlat = [ Dlat Dlat(end) ];
  Dlat = repmat(Dlat', [ 1 Nz ]);

  % integrand
  Integrand = PR .* Geom .* DuDp .* Dlat;

  % integral goes from 0 to Lat so do a cumulative sum along the
  % columns of Integrand to get the integral value.
  Tint = cumsum(Integrand, 1);

  T = Tref + Tint;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ RH ] = GenInitRelHumField(Lat, Z, Zt)
% GenInitRelHumField - generate a relative humidity field

  Ny = length(Lat);
  Nz = length(Z);

  % Use method in Boutle et al., 2010

  % latitudinal variation
  %   R = 1                          Lat < 25N
  %   R = 0.5                        Lat > 65N
  %   R = 1 - 0.5 * ((Lat-25)/40)    25 <= Lat <= 65
  %
  R = 1 - (0.5 .* ((Lat-25)./40));
  R(Lat < 25) = 1;
  R(Lat > 65) = 0.5;

  R = repmat(reshape(R, [ Ny 1]), [ 1 Nz ]);

  % RH calculation (as a fraction instead of percent)
  %   RH = 0.80 * (1-0.9R*(Z/Zt)^(5/4))        Z <  Zt
  %   RH = 0.05                                Z >= Zt
  %
  Z_2D = repmat(Z, [ Ny 1 ]);
  Zfact = (Z_2D ./ Zt).^(1.25);
  RH = 0.8 .* (1 - (0.9 .* R .* Zfact));
  RH(Z_2D >= Zt) = 0.05;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Zg ] = GenGeoPotField(Lat, Press, T, R, G0)
% GenGeoPotField - generate a geopotential height field

  Ny = length(Lat);
  Nz = length(Press);

  % Assume near surface that the geopotential and actual height are equal.
  % Use CRC formula for height as a function of pressure to determine
  % the lowest level geopotential height.
  %
  % Then use hypsometric equation to determine all other geopotential
  % heights.

  Zg = zeros([ Ny Nz ]);

  % First level - CRC formula:
  %   z = 44331.5 - 4946.62 * p^0.190263
  %   p in pascals, z in meters
  %
  % Note Press is in mb (hPa)
  %
  Zg(:,1) = 44331.5 - (4946.62 .* ((Press(1) .* ones([ Ny 1]) .* 100).^0.190263));

  % 2 though Nz levels - Hypsometric eqn:
  %   z2 - z1 = R/G0 * Tavg * ln(p1/p2)
  %               OR
  %   z2 = z1 + R/G0 * Tavg * ln(p1/p2)
  Tavg  = (T(:,1:end-1) + T(:,2:end)) .* 0.5;
  Pfact = repmat(log(Press(1:end-1) ./ Press(2:end)), [ Ny 1 ]);
  for i = 2:Nz
    Zg(:,i) = Zg(:,i-1) + ((R/G0) .* Tavg(:,i-1) .* Pfact(:,i-1));
  end
end
