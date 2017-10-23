function [ Tprime ] = GenInitTperturb(LON, LAT, PRESS, That, LatHat, WaveNum)
% GenInitTperturb generate initial temperature perturbation

  % That is magnitude of the perturbation
  % LatHat is used to set the center of the perturbations
  % WaveNum is zonal wavenumber

  Nx = length(LON);
  Ny = length(LAT);
  Nz = length(PRESS);

  % Perturbation is 3D field
  % Set up LON, LAT values into 3D arrays, then apply
  % Formula from Polvani and Esler, 2007.
  % Lat and Lon need to be in radians
  %
  LON_3D = repmat((LON.*pi./180)', [ 1 Ny Nz ]);
  LAT_3D = repmat((LAT.*pi./180),  [ Nx 1 Nz ]);
  LatHatRad = LatHat .* pi ./180;

  Tprime = That .* cos(WaveNum .* LON_3D) .* (sech(WaveNum .* (LAT_3D - LatHatRad))).^2;
end
