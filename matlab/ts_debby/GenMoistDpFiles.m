function [ ] = GenMoistDpFiles()
% GenMoistDpFiles function to add moisture to Dprep files of the SAL

  DinDir = 'DP_FILES/ORIG';
  DoutDir = 'DP_FILES/MOIST';

  % Limits for finding SAL
  RhLimit  = 0.45;
  Pmax     = 900;   % 900mb -> 1km
  Pmin     = 550;   % 550mb -> 5km

  % Use the built-in topography map to define continents so that data over Africa can be
  % left alone.
  WORLD = load('topo', 'topo');
  CONT = WORLD.topo >= 0; % continents exist where elevation >= sea level (zero)

  % The DP file data and CONT array are all on 1 degree grids, but the DP longitude goes from
  % -180 to 179, whereas the CONT array goes from 0 to 359; DP latitude goes from -90 to 90,
  % whereas the CONT array goes from -89 to 90.
  %
  % Change CONT to match what will be read in from the DP files.
  CONT = [ CONT(:,181:360) CONT(:,1:180) ]; % shift to place 0 deg long (prime meridian) in the center
  CONT = [ CONT; CONT(1,:) ];               % repeat 89S to fill in 90S
  CONT = CONT';                             % transpose to match definition of DP file arrays

  % Limit the region that the interpolation will take place to the Atlantic Ocean off the west
  % coase of Africa
  SwLat = 10;
  SwLon = -60;
  NeLat = 40;
  NeLon = -5;

  % Walk through all of the input Dprep files and add moisture to the SAL region
  DinPattern = sprintf('%s/dp-*', DinDir);
  DP_FILES = dir(DinPattern);
  Nfiles = length(DP_FILES);

%  for i = 1:Nfiles
  for i = 1:1
    DpFile = DP_FILES(i).name;
    InFile = sprintf('%s/%s', DinDir, DpFile);
    OutFile = sprintf('%s/%s', DoutDir, DpFile);

    fprintf('*********************************************************************************\n');
    fprintf('Adding moisture to Dprep file:\n');

    fprintf('    Reading:%s\n', InFile);
    [ U V T Zg RH Lon Lat Press Year Month Day Time ] = ReadDpFile(InFile);
    [ Nx Ny Np ] = size(RH);

    % If this is the first file, form a mask that will define the area where interpolation
    % can take place.
    if (i == 1)
      I1 = find(Lon >= SwLon, 1, 'first');
      I2 = find(Lon <= NeLon, 1, 'last');
      J1 = find(Lat >= SwLat, 1, 'first');
      J2 = find(Lat <= NeLat, 1, 'last');

      INTERP_AREA = ~CONT; % exclude continental land
      INTERP_AREA(1:I1,:) = false;
      INTERP_AREA(I2:end,:) = false;
      INTERP_AREA(:,1:J1) = false;
      INTERP_AREA(:,J2:end) = false;
    end

    % Create a mask of where the SAL is located based on RH being less than a threshold
    % Do this level by level since the region of low RH will vary from level to level
    P1 = find(Press <= Pmax, 1, 'first');
    P2 = find(Press >= Pmin, 1, 'last');

    fprintf('    Adding moisture\n');
    MOIST_RH = RH;
    for k = P1:P2
      MASK = squeeze(RH(:,:,k)) < RhLimit; % 0's where RH >= RhLimit, 1's where RH < RhLimit
      MASK = MASK & INTERP_AREA;    % Limit MASK to Atlantic region west of Africa

      HDAT = squeeze(MOIST_RH(:,:,k));
      HDAT(MASK) = nan;
      HDAT = inpaintn(HDAT);
      HDAT(HDAT<0) = 0;
      MOIST_RH(:,:,k) = HDAT;

    end

    fprintf('    Writing:%s\n', OutFile);
    WriteDpFile(U, V, T, Zg, MOIST_RH, Lon, Lat, Press, Year, Month, Day, Time, OutFile);

    fprintf('\n');
  end
end
