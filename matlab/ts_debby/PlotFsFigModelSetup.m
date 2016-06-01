function [ ] = PlotFsFigModelSetup()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

  CaseList = {
    { 'TSD_SAL_DUST'       'SD'   'black' }
    { 'TSD_NONSAL_DUST'    'NSD'  'red'   }
    { 'TSD_SAL_NODUST'     'SND'  'blue'  }
    { 'TSD_NONSAL_NODUST'  'NSND' 'green' }
    };
  Nc = length(CaseList);

  % Storm center locations, max wind and min pressure
  for ic = 1:Nc
    Case  = CaseList{ic}{1};
    Label = CaseList{ic}{2};
    Color = CaseList{ic}{3};

    % Storm track
    InFname = sprintf('HDF5/%s/HDF5/storm_center-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    InVname = '/press_cent_xloc';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    SimTrackLons(:,ic) = squeeze(h5read(InFname, InVname));

    InVname = '/press_cent_yloc';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    SimTrackLats(:,ic) = squeeze(h5read(InFname, InVname));

    LegText{ic} = Label;
    LineColors{ic} = Color;
  end


  % Get grid locations and create map lat,lon
  GlocFname = 'GridLocs.txt';
  Grid1Loc = ReadGridLoc(GlocFname, 'Grid1:');
  Grid2Loc = ReadGridLoc(GlocFname, 'Grid2:');
  Grid3Loc = ReadGridLoc(GlocFname, 'Grid3:');

  % Make the map extend a little ways beyond grid1 extent.
  MapInc = 10;
  BigMapLoc   = Grid1Loc + [ -MapInc MapInc -MapInc MapInc ];
  MapInc = 4;
  SmallMapLoc = Grid3Loc + [ -MapInc MapInc -MapInc MapInc ];

  % The lat,lon in the first entry of every storm track will be the same
  % since all sims started from the same history file.
  InitLat =  13.5;
  InitLon = -23.4;

  %*************************************************************************
  % Read in and assemble the aerosol profiles for the DUST case
  % After the text scan, InData will contain:
  %    InData{1}  ->  Z
  %    InData{2}  ->  CCN # conc
  %    InData{3}  ->  Dust 1 (small) # conc
  %    InData{4}  ->  Dust 2 (large) # conc
  %    InData{5}  ->  IN # conc
  %
  %    all # conc in #/cc
  %
  DustInFile = 'z.NAMMA_RAMS_LEVELS_aerosols.txt';
  fprintf('  Reading: %s\n', DustInFile);
  InFileId = fopen(DustInFile);
  InData = textscan(InFileId, '%f %f %f %f %f');
  fclose(InFileId);

  Znamma   = InData{1} ./ 1000; % km
  DUST_CONC = [ InData{3} InData{4} ];

  % Vertically integrated dust mass
  VintDustFile = 'HDF5/TSD_SAL_DUST/HDF5/vint_dust-TSD_SAL_DUST-AS-2006-08-20-120000-g1.h5';
  VintDustVname = '/vertint_dust';

  % Read in and create a snapshot of the initial (06Z, 22Aug)
  % vertically integrated dust field Aug22, 06Z
  fprintf('Reading: %s (%s)\n', VintDustFile, VintDustVname);
  X    = squeeze(h5read(VintDustFile, '/x_coords'));
  Y    = squeeze(h5read(VintDustFile, '/y_coords'));

  Nx = length(X);
  Ny = length(Y);

  % Take the first time step (initialized dust field)
  Dstart = [ 1 1 1 ];
  Dcount = [ Nx Ny 1 ];
  DUST = squeeze(h5read(VintDustFile, VintDustVname, Dstart, Dcount));

  % plot
  Fig = figure;

  Paxes = subplot(2,2,1);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.12;
  Ploc(3) = Ploc(3) + 0.2;
  Ploc(2) = Ploc(2) - 0.1;
  Ploc(4) = Ploc(4) + 0.2;
  set(Paxes, 'Position', Ploc);
  PlotGridLocs(Paxes, BigMapLoc, Grid1Loc, Grid2Loc, Grid3Loc, 'a', '', Fsize);
  
  % Initial dust profiles
  Paxes = subplot(2,2,2);
  PlotDpFigInitDust(Paxes, DUST_CONC, Znamma, 'b', '', Fsize);

  % Vertically integrated dust mass
  Paxes = subplot(2,2,3);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.10;
  Ploc(3) = Ploc(3) + 0.12;
  Ploc(2) = Ploc(2) - 0.10;
  Ploc(4) = Ploc(4) + 0.10;
  set(Paxes, 'Position', Ploc);
  PlotVintDust(Paxes, SmallMapLoc, Grid3Loc, X, Y, DUST, InitLat, InitLon, 'c', '06Z, 22Aug', Fsize, 1, 1);

  Paxes = subplot(2,2,4);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.05;
  Ploc(3) = Ploc(3) + 0.12;
  Ploc(2) = Ploc(2) - 0.10;
  Ploc(4) = Ploc(4) + 0.10;
  set(Paxes, 'Position', Ploc);
  PlotDpFigTrack(Paxes, SmallMapLoc, Grid3Loc, SimTrackLons, SimTrackLats, LegText, LineColors, 'd', '', Fsize);

  OutFile = sprintf('%s/FsFigModelSetup.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = PlotGridLocs(Paxes, MapLoc, Grid1Loc, Grid2Loc, Grid3Loc, Pmarker, Ptitle, Fsize)

  axes(Paxes);

  LegendFsize = 12;
  LineW = 2;

  MapLat1 = MapLoc(1);
  MapLat2 = MapLoc(2);
  MapLon1 = MapLoc(3);
  MapLon2 = MapLoc(4);

  G1_Lat1 = Grid1Loc(1);
  G1_Lat2 = Grid1Loc(2);
  G1_Lon1 = Grid1Loc(3);
  G1_Lon2 = Grid1Loc(4);

  G2_Lat1 = Grid2Loc(1);
  G2_Lat2 = Grid2Loc(2);
  G2_Lon1 = Grid2Loc(3);
  G2_Lon2 = Grid2Loc(4);

  G3_Lat1 = Grid3Loc(1);
  G3_Lat2 = Grid3Loc(2);
  G3_Lon1 = Grid3Loc(3);
  G3_Lon2 = Grid3Loc(4);

  MapLats = [ MapLat1 MapLat2 ];
  MapLons = [ MapLon1 MapLon2 ];

  % Draw a box for each grid
  Grid1Lats = [ G1_Lat1 G1_Lat2 G1_Lat2 G1_Lat1 G1_Lat1 ];
  Grid1Lons = [ G1_Lon1 G1_Lon1 G1_Lon2 G1_Lon2 G1_Lon1 ];

  Grid2Lats = [ G2_Lat1 G2_Lat2 G2_Lat2 G2_Lat1 G2_Lat1 ];
  Grid2Lons = [ G2_Lon1 G2_Lon1 G2_Lon2 G2_Lon2 G2_Lon1 ];

  Grid3Lats = [ G3_Lat1 G3_Lat2 G3_Lat2 G3_Lat1 G3_Lat1 ];
  Grid3Lons = [ G3_Lon1 G3_Lon1 G3_Lon2 G3_Lon2 G3_Lon1 ];

  set(Paxes, 'FontSize', Fsize);
  worldmap(MapLats, MapLons);
  setm(Paxes, 'MapProjection', 'miller');
  load coast;  % creates vars "lat" and "long" that contain coast lines
  plotm(lat, long, 'Color', 'k', 'LineWidth', LineW);

  % Plot grid boxes
  G1color = str2rgb('red');
  G2color = str2rgb('green');
  G3color = str2rgb('blue');

  G1line = linem(Grid1Lats, Grid1Lons, 'LineWidth', LineW, 'LineStyle', '-', 'Color', G1color);
  G2line = linem(Grid2Lats, Grid2Lons, 'LineWidth', LineW, 'LineStyle', '-', 'Color', G2color);
  G3line = linem(Grid3Lats, Grid3Lons, 'LineWidth', LineW, 'LineStyle', '-', 'Color', G3color);

%  LegLines = [ G1line G2line G3line ];
%  LegText  = { 'Grid1' 'Grid2' 'Grid3' };
%  legend(LegLines, LegText, 'Location', 'SouthOutside', 'FontSize', LegendFsize);
  textm(-17, -68, 'Grid1', 'Color', G1color, 'FontSize', 10);
  textm( 27, -40, 'Grid2', 'Color', G2color, 'FontSize',  8);
  textm( 18, -35, 'Grid3', 'Color', G3color, 'FontSize',  7);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = PlotDpFigInitDust(Paxes, DUST, Z, Pmarker, Ptitle, Fsize)
  axes(Paxes);

  Zmin =  0; %km
  Zmax = 13; %km

  % Make the plot
  Xlab = 'N_d (cm^-^3)';
  Ylab = 'Height (km)';

  Xlimits = [ -10 250 ];
  Xticks = [ 0 100 200 ];

  LegFsize = 10;
  LineWidth = 2;

  set(Paxes, 'FontSize', Fsize);
  plot(DUST, Z, 'LineWidth', LineWidth);
  set (gca, 'FontSize', Fsize);
  xlim(Xlimits);
  ylim([ Zmin Zmax ]);
  set (gca, 'XTick', Xticks);

  xlabel(Xlab);
  ylabel(Ylab);

  legend({ 'Small Dust' 'Large Dust' }, 'Location', 'NorthEast', 'FontSize', LegFsize);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlotVintDust(Paxes, MapLoc, Grid3Loc, X, Y, Z, InitLat, InitLon, Pmarker, Ptitle, Fsize, ShowX, ShowY)

  axes(Paxes);

  MapLat1 = MapLoc(1);
  MapLat2 = MapLoc(2);
  MapLon1 = MapLoc(3);
  MapLon2 = MapLoc(4);

  MapLats = [ MapLat1 MapLat2 ];
  MapLons = [ MapLon1 MapLon2 ];

  G3_Lat1 = Grid3Loc(1);
  G3_Lat2 = Grid3Loc(2);
  G3_Lon1 = Grid3Loc(3);
  G3_Lon2 = Grid3Loc(4);

  Grid3Lats = [ G3_Lat1 G3_Lat2 G3_Lat2 G3_Lat1 G3_Lat1 ];
  Grid3Lons = [ G3_Lon1 G3_Lon1 G3_Lon2 G3_Lon2 G3_Lon1 ];

  LineW = 2;

  set(Paxes, 'FontSize', Fsize);
  worldmap(MapLats, MapLons);
  setm(Paxes, 'MapProjection', 'miller');

  % Want this plot to show where the initial SAL regin exists within
  % grid 3. The initial vertically integrated dust field is going to
  % contain two values: 
  %    Inside SAL:    1.8e6
  %    Outside SAL:   zero
  %
  % Replace the Outside SAL values with nans so that they don't show
  % up on the plot.
  ZSAL = Z;
  ZSAL(ZSAL < 10) = nan;

  % The values of dust in the SAL are the same so cut down the size
  % of ZSAL in order to speed up the contourfm performance.
  Dinc = 1;
  XSAL = X(1:Dinc:end);
  YSAL = Y(1:Dinc:end);
  ZSAL = ZSAL(1:Dinc:end,1:Dinc:end);

  % X is LON, Y is LAT
  % Want to show log values of Z, but there is not a built-in
  % method for this.
  %   1. send contourfm the log values
  %      this means you have to work in log values for the caxis
  %      and tick marks, etc.
  %      e.g., caxis [ -2 6 ] translates to 1e-2 to 1e6
  %   2. contourfm/contourcmap create a colorbar with
  %      the selected colors for the contour levels displayed
  %      like tiles and a scaled tick mark value that places
  %      integer values at the center of each color. This makes
  %      the edges of the tiles equal to tick values like 0.5 1.5 2.5 ...
  %      where the centers of the tiles are equal to tick values
  %      1 2 3 ... n where n is the number of contour levels.
  % 
  % So, in this example the range for Z is 1e-2 to 1e6. This
  % calls for:
  %   caxis: -2 to 6
  %   contour levels: -2:0.5:6 (sixteen levels)
  %   color bar axes: center of tiles have values 1 2 3 ... 16
  %                   edges have values 0.5 1.5 2.5 ... 16.5
  %
  %   This means that there is a linear scale between the colorbar
  %   tick mark values and the actual values in the plot:
  %
  %      tick value        plot value    origninal value
  %         0.5               -2            0.01
  %         1.5               -1.5
  %         2.5               -1            0.1
  %         3.5               -0.5
  %         4.5                0            1
  %         ...                ...
  %        15.5                5.5
  %        16.6                6            1000000
  %                  
  %
  
  Clevs = [ -2:0.5:6 ];
  Cticks = [ 2.5 6.5 10.5 14.5 ];
  CtickLabels = { '10^-^1' '10^1' '10^3' '10^5' };
  contourfm(double(YSAL), double(XSAL), double(ZSAL)', Clevs, 'LineStyle', 'none'); 
  caxis(Paxes, [ -2 6 ]);

  colormap('copper');

  % Draw the coast after the contourfm call so that coast lines will appear on top
  % of the contour plot
  Coast = load('coast');
  plotm(Coast.lat, Coast.long, 'Color', 'k', 'LineWidth', LineW);
  textm(20, -35, 'SAL', 'FontSize', 10, 'Color', 'k');

  % Read in and place hurricane symbol where TS Debby is located
  H_SYMBOL = imread('IMAGES/HurricaneSymbolBlue.png');
  HsymInc = 2; % make the symbol 4 degrees in diameter
  HsymLat1 = InitLat - HsymInc;
  HsymLat2 = InitLat + HsymInc;
  HsymLon1 = InitLon - HsymInc;
  HsymLon2 = InitLon + HsymInc;

  [ Nlat Nlon Ncolors ] = size(H_SYMBOL);
  HsymLatInc = (HsymLat2-HsymLat1) / (Nlat-1);
  HsymLonInc = (HsymLon2-HsymLon1) / (Nlon-1);
  HLAT = fliplr(HsymLat1:HsymLatInc:HsymLat2);
  HLON = HsymLon1:HsymLonInc:HsymLon2;
  [ HLON_GRID HLAT_GRID ] = meshgrid(HLON, HLAT);

  geoshow(HLAT_GRID, HLON_GRID, H_SYMBOL);

  % Mark Africa
  textm(14, -12, 'Africa', 'FontSize', 10, 'Color', 'k', 'Rotation', 90);

  % Mark grid 3
  Grid3Color = str2rgb('blue');
  linem(Grid3Lats, Grid3Lons, 'LineWidth', LineW, 'LineStyle', '-', 'Color', Grid3Color);
  textm(8, -39, 'Grid3', 'FontSize', 12, 'Color', Grid3Color);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = PlotDpFigTrack(Paxes, MapLoc, Grid3Loc, SimTrackLons, SimTrackLats, SimTrackLabels, SimTrackColors, Pmarker, Ptitle, Fsize)

  axes(Paxes);

  LegText = { 'NHC' SimTrackLabels{1:end} };
  
  LegendFsize = 8;
  LineW = 2;

  [ Npts Nlines ] = size(SimTrackLats);
  
  MapLat1 = MapLoc(1);
  MapLat2 = MapLoc(2);
  MapLon1 = MapLoc(3);
  MapLon2 = MapLoc(4);

  MapLats = [ MapLat1 MapLat2 ];
  MapLons = [ MapLon1 MapLon2 ];

  G3_Lat1 = Grid3Loc(1);
  G3_Lat2 = Grid3Loc(2);
  G3_Lon1 = Grid3Loc(3);
  G3_Lon2 = Grid3Loc(4);

  Grid3Lats = [ G3_Lat1 G3_Lat2 G3_Lat2 G3_Lat1 G3_Lat1 ];
  Grid3Lons = [ G3_Lon1 G3_Lon1 G3_Lon2 G3_Lon2 G3_Lon1 ];

  % The simulation runs from Aug 22, 6Z to Aug 24, 18Z.
  %
  % Locations from NHC report: Aug 22, 6Z through Aug 24, 18Z.
  NhcTrackLats = [ 12.6 13.4 14.2 14.9 15.7 16.7 17.6 18.4 19.2 20.1 20.9 ];
  NhcTrackLons = [ 23.9 25.3 26.7 28.1 29.5 31.0 32.4 33.9 35.5 37.1 38.7 ] * -1; 

  set(Paxes, 'FontSize', Fsize);
  worldmap(MapLats, MapLons);
  setm(Paxes, 'MapProjection', 'miller');
  load coast;  % creates vars "lat" and "long" that contain coast lines
  plotm(lat, long, 'Color', 'k', 'LineWidth', LineW);

  NhcTrack = linem(NhcTrackLats, NhcTrackLons, 'LineWidth', LineW, 'Color', 'k', 'LineStyle', 'none', 'Marker', '+');
  for i = 1:Nlines
    Lats = double(squeeze(SimTrackLats(:,i)));
    Lons = double(squeeze(SimTrackLons(:,i)));
    Color = str2rgb(SimTrackColors{i});
    SimTrack(i) = linem(Lats, Lons, 'Color', Color, 'LineWidth', LineW);
  end

  % Mark Africa
  textm(14, -12, 'Africa', 'FontSize', 10, 'Color', 'k', 'Rotation', 90);

  % Mark grid 3
  Grid3Color = str2rgb('blue');
  linem(Grid3Lats, Grid3Lons, 'LineWidth', LineW, 'LineStyle', '-', 'Color', Grid3Color);
  textm(9, -39, 'Grid3', 'FontSize', 12, 'Color', Grid3Color);

  legend( [ NhcTrack SimTrack ], LegText,'Location', 'EastOutside', 'FontSize', LegendFsize);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end

