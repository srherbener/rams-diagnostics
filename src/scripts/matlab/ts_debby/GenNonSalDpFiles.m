function [ ] = GenNonSalDpFiles()
% GenNonSalDpFiles function to add moisture to Dprep files of the SAL

  DinDir = 'DP_FILES/ORIG';
  DoutDir = 'DP_FILES/NONSAL';

  % Limits for finding SAL
  Pmax     = 900;   % 900mb -> 1km
  Pmin     = 550;   % 550mb -> 5km

  % From Dunion & Marron (2008). Values for the non-SAL mean RH profile for August (their Fig. 6)
  % Dunion & Marron (2008) Fig. 6 only goes up to 200 mb, so only fill in these values.
  %   Should be okay since the SAL doesn't even get close to the 200 mb height.
  % Keep this in sync with pressure levels in Dprep file:
  %       i   Press(i)    Non-SAL Aug profile(i) (%)
  %       1    1000                 83
  %       2     975                 83
  %       3     950                 83
  %       4     925                 83
  %       5     900                 82
  %       6     850                 76
  %       7     800                 73
  %       8     750                 67
  %       9     700                 63
  %      10     650                 58
  %      11     600                 53
  %      12     550                 50
  %      13     500                 46
  %      14     450                 43
  %      15     400                 40
  %      16     350                 40
  %      17     300                 40
  %      18     250                 39
  %      19     200                 38
  %      20     150                  -
  %      21     100                  -
  %      22      70                  -
  %      23      50                  -
  %      24      30                  -
  %      25      20                  -
  %      26      10                  -
  %
  RhLimitList = [ 83 83 83 83 82 76 73 67 63 58 53 50 46 43 40 40 40 39 38 ] .* 0.01; % Dprep uses fraction instead of %

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

  for i = 1:Nfiles
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
    % Find the SAL location using RH
    P1 = find(Press <= Pmax, 1, 'first');
    P2 = find(Press >= Pmin, 1, 'last');

    fprintf('    Adding moisture\n');
    MOIST_RH = RH;
    COOL_T = T;
    VOLUME_MASK = zeros([ Nx Ny Np ]);
    for k = P1:P2
      RhLimit = RhLimitList(k);
      MASK = squeeze(RH(:,:,k)) < RhLimit; % 0's where RH >= RhLimit, 1's where RH < RhLimit
      MASK = MASK & INTERP_AREA;    % Limit MASK to Atlantic region west of Africa

      % Adjust RH and T
      HDAT = squeeze(MOIST_RH(:,:,k));
      TDAT = squeeze(COOL_T(:,:,k));

      HDAT(MASK) = nan;
      TDAT(MASK) = nan;

      HDAT = inpaintn(HDAT);
      TDAT = inpaintn(TDAT);

      MOIST_RH(:,:,k) = HDAT;
      COOL_T(:,:,k)   = TDAT;

      % Build up volume mask
      VOLUME_MASK(:,:,k) = MASK;
    end

    fprintf('    Writing:%s\n', OutFile);
    WriteDpFile(U, V, COOL_T, Zg, MOIST_RH, Lon, Lat, Press, Year, Month, Day, Time, OutFile);

    OutMaskFile = sprintf('%s-MASK.h5', OutFile);
    fprintf('    Writing mask: %s\n', OutMaskFile);
    hdf5write(OutMaskFile, 'MASK', VOLUME_MASK);
    hdf5write(OutMaskFile, 'LAT', Lat, 'WriteMode', 'append');
    hdf5write(OutMaskFile, 'LON', Lon, 'WriteMode', 'append');
    hdf5write(OutMaskFile, 'PRESS', Press, 'WriteMode', 'append');

    fprintf('\n');
  end
end
