function [ ] = PlotDpFigGenRes()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

  % Get grid 3 location
  GlocFname = 'GridLocs.txt';
  Grid3Loc = ReadGridLoc(GlocFname, 'Grid3:');

  % Map locations:
  %  grid3 extended a small amount
  MapLocSmall = Grid3Loc + [ -2 +3 -2 +6 ];
  MapLocLarge = [ 0 40 -55 5 ];

  % Image and vint dust don't quite fit right within grid3 when plotted. Move grid3
  % up about 0.5 degree to make up for this. Do this after setting the map confines
  % so that the alignment of the SAL image remains lined up with the coast line of
  % Africa.
  Grid3Loc = Grid3Loc + [ +0.25 +0.25 0 0 ];

  % Storm center locations
  Hfile = 'HDF5/TSD_SAL_DUST/HDF5/storm_center-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  HdsetLon = '/press_cent_xloc';
  HdsetLat = '/press_cent_yloc';
  fprintf('Reading: %s\n', Hfile);
  fprintf('  Track longitude: %s\n', HdsetLon);
  fprintf('  Track latitude: %s\n', HdsetLat);

  % Slons and Slats are (t)
  SimTrackLons = squeeze(h5read(Hfile, HdsetLon));
  SimTrackLats = squeeze(h5read(Hfile, HdsetLat));

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

  % SAL strength image
  SalPic = imread('IMAGES/SAL_Aug23_12Z_DP_Fig1.png');

  % Vertically integrated dust mass
  VintDustFile = 'HDF5/TSD_SAL_DUST/HDF5/vint_dust-TSD_SAL_DUST-AS-2006-08-20-120000-g1.h5';
  VintDustVname = '/vertint_dust';

  % Read in and create a snapshot of the vertically integrated dust field at Aug23, 12Z
  fprintf('Reading: %s (%s)\n', VintDustFile, VintDustVname);
  X    = squeeze(h5read(VintDustFile, '/x_coords'));
  Y    = squeeze(h5read(VintDustFile, '/y_coords'));
  Nx = length(X);
  Ny = length(Y);
  % VINT_DUST is (x,y,t)
  %   t = 61 corresponds to Aug23, 12Z
  %   convert to g/m^2
  VINT_DUST = squeeze(h5read(VintDustFile, VintDustVname,[ 1 1 61 ], [ Nx Ny 1 ])) .* 1e-6;

  % plot
  Fig = figure;

  % Map doesn't quite fit in the vertical center of the subplot region. Needs to get
  % bumped upward just a little bit.
  Paxes = subplot(2,2,1);
  PlotDpFigTrack(Paxes, MapLocSmall, Grid3Loc, SimTrackLons, SimTrackLats, 'a', '', Fsize);
  
  % Initial dust profiles
  % Doesn't quite fit in the vertical center of the subplot region. Needs to get
  % bumped upward just a little bit.
  Paxes = subplot(2,2,2);
  PlotDpFigInitDust(Paxes, DUST_CONC, Znamma, 'b', '', Fsize);

  % SAL strength image
  Paxes = subplot(2,2,3);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.02;
  Ploc(3) = Ploc(3) * 1.04;
  Ploc(2) = Ploc(2) - 0.02;
  Ploc(4) = Ploc(4) * 1.04;
  set(Paxes, 'Position', Ploc);
  PlaceSalImage(Paxes, MapLocLarge, Grid3Loc, SalPic, 'c', '12Z, 23Aug', Fsize);

  % Vertically integrated dust mass
  Paxes = subplot(2,2,4);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.05;
  Ploc(3) = Ploc(3) * 1.10;
  Ploc(2) = Ploc(2) - 0.05;
  Ploc(4) = Ploc(4) * 1.10;
  set(Paxes, 'Position', Ploc);
  PlotVintDust(Paxes, MapLocLarge, Grid3Loc, X, Y, VINT_DUST, 'd', '12Z, 23Aug', Fsize, 1, 1);

  OutFile = sprintf('%s/DpFig1_GenRes.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = PlotDpFigTrack(Paxes, MapLoc, Grid3Loc, SimTrackLons, SimTrackLats, Pmarker, Ptitle, Fsize)

  axes(Paxes);

  LegText = { 'NHC' 'Model' };
  
  LegendFsize = 8;
  LineW = 2;
  
  LatBounds = [ MapLoc(1) MapLoc(2) ];
  LonBounds = [ MapLoc(3) MapLoc(4) ];

  Grid3Lats = [ Grid3Loc(1) Grid3Loc(2) Grid3Loc(2) Grid3Loc(1) Grid3Loc(1) ];
  Grid3Lons = [ Grid3Loc(3) Grid3Loc(3) Grid3Loc(4) Grid3Loc(4) Grid3Loc(3) ];

  G3color = str2rgb('blue');
  G3linew = 2.5;
  
  % The simulation runs from Aug 22, 6Z to Aug 24, 18Z.
  %
  % Locations from NHC report: Aug 22, 6Z through Aug 24, 18Z.
  NhcTrackLats = [ 12.6 13.4 14.2 14.9 15.7 16.7 17.6 18.4 19.2 20.1 20.9 ];
  NhcTrackLons = [ 23.9 25.3 26.7 28.1 29.5 31.0 32.4 33.9 35.5 37.1 38.7 ] * -1; 

  % SAL_AR gets clipped at the top of grid3
  SalArTop = min( [ 24.0 Grid3Loc(2) ]);
  SalRegionLats = [ 12.0 SalArTop SalArTop 12.0 12.0 ];
  SalRegionLons = [ 36.5 36.5 25.0 25.0 36.5 ] * -1; 

  set(Paxes, 'FontSize', Fsize);
  worldmap(LatBounds, LonBounds);
  setm(Paxes, 'MapProjection', 'miller');
  load coast;  % creates vars "lat" and "long" that contain coast lines
  plotm(lat, long, 'Color', 'k', 'LineWidth', LineW);

  NhcTrack = linem(NhcTrackLats, NhcTrackLons, 'LineWidth', LineW, 'Color', 'k', 'LineStyle', 'none', 'Marker', '+');
  SimTrack = linem(double(SimTrackLats), double(SimTrackLons), 'LineWidth', LineW);

  % Add SAL analysis region
  linem(SalRegionLats, SalRegionLons, 'LineWidth', LineW, 'Color', 'r');
  textm(21, -34, 'SAL\_AR', 'FontSize', 10, 'Color', 'r');

  % Mark Africa
  textm(14, -12, 'Africa', 'FontSize', 10, 'Color', 'k', 'Rotation', 90);

  % Mark Grid3
  linem(Grid3Lats, Grid3Lons, 'LineWidth', G3linew, 'LineStyle', '-', 'Color', G3color);
  textm( 24.1, -23, 'Grid3', 'Color', G3color, 'FontSize',  10);

  legend( [ NhcTrack SimTrack' ], LegText,'Location', 'SouthWest', 'FontSize', LegendFsize);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlaceSalImage(Paxes, MapLoc, Grid3Loc, SalImage, Pmarker, Ptitle, Fsize)

  axes(Paxes);

  LineW = 2;
  G3linew = 2.5;

  LatBounds = [ MapLoc(1) MapLoc(2) ];
  LonBounds = [ MapLoc(3) MapLoc(4) ];

  Grid3Lats = [ Grid3Loc(1) Grid3Loc(2) Grid3Loc(2) Grid3Loc(1) Grid3Loc(1) ];
  Grid3Lons = [ Grid3Loc(3) Grid3Loc(3) Grid3Loc(4) Grid3Loc(4) Grid3Loc(3) ];

  G3color = str2rgb('white');

  % create a gap so that image lines up with vint_dust image
  Gap = 0.05;
  Ploc = get(Paxes, 'Position');

  % move the plot up a bit
  Ploc(2) = Ploc(2) + Gap;
  Ploc(4) = Ploc(4) - Gap;
  set(Paxes, 'Position', Ploc);

  set(Paxes, 'FontSize', Fsize);
  worldmap(LatBounds, LonBounds);
  setm(Paxes, 'MapProjection', 'miller');


  % Create gridded versions of the latitude and longitude value
  % that the image covers.
  [ Nlat Nlon Ncolors ] = size(SalImage);
  Lat1 = 3.1;
  Lat2 = 35.4;
  Lon1 = -55;
  Lon2 = 5;

  LatInc = (Lat2-Lat1)/(Nlat-1);
  LonInc = (Lon2-Lon1)/(Nlon-1);

  LAT = fliplr(Lat1:LatInc:Lat2);
  LON = Lon1:LonInc:Lon2;
  [ LonGrid LatGrid ] = meshgrid(LON, LAT);

  geoshow(LatGrid, LonGrid, SalImage);

  % Draw the coast after the contourfm call so that coast lines will appear on top
  % of the contour plot
  Coast = load('coast');
  plotm(Coast.lat, Coast.long, 'Color', 'k', 'LineWidth', LineW);

  % Mark Africa
  textm(25, -6, 'Africa', 'FontSize', 10, 'Color', 'w');

  % Mark Grid3
  linem(Grid3Lats, Grid3Lons, 'LineWidth', G3linew, 'LineStyle', '-', 'Color', G3color);
  textm(25, -25, 'Grid3', 'Color', G3color, 'FontSize',  10);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlotVintDust(Paxes, MapLoc, Grid3Loc, X, Y, Z, Pmarker, Ptitle, Fsize, ShowX, ShowY)

  axes(Paxes);

  LatBounds = [ MapLoc(1) MapLoc(2) ];
  LonBounds = [ MapLoc(3) MapLoc(4) ];

  Grid3Lats = [ Grid3Loc(1) Grid3Loc(2) Grid3Loc(2) Grid3Loc(1) Grid3Loc(1) ];
  Grid3Lons = [ Grid3Loc(3) Grid3Loc(3) Grid3Loc(4) Grid3Loc(4) Grid3Loc(3) ];

  G3color = str2rgb('blue');

  LineW = 2;
  G3linew = 2.5;

  % create a gap for the colorbar
  Gap = 0.05;
  Ploc = get(Paxes, 'Position');
  CbarLoc = Ploc;

  % move the plot up a bit
  Ploc(2) = Ploc(2) + Gap;
  Ploc(4) = Ploc(4) - Gap;
  set(Paxes, 'Position', Ploc);

  % define the region for the colorbar in the gap
  CbarLoc(4) = Gap * 0.5;

  set(Paxes, 'FontSize', Fsize);
  worldmap(LatBounds, LonBounds);
  setm(Paxes, 'MapProjection', 'miller');

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
  
%  % Plot on log scale
%  Clevs = [ -2:0.5:6 ];
%  Cticks = [ 2.5 6.5 10.5 14.5 ];
%  CtickLabels = { '10^-^1' '10^1' '10^3' '10^5' };
%  contourfm(double(Y), double(X), double(log10(Z))', Clevs, 'LineStyle', 'none'); 
%  caxis(Paxes, [ -2 6 ]);

  % Plot on linear scale
  Clevs = [ 0:0.1:3 ];
  Cticks = [ 0.5 5.5 10.5 15.5 20.5 25.5 30.5 ];
  CtickLabels = { '0.0' '0.5' '1.0' '1.5' '2.0' '2.5' '3.0' };
  contourfm(double(Y), double(X), double(Z)', Clevs, 'LineStyle', 'none'); 
  caxis([ 0 3 ]);

% The subsequent pause statements are placed in the code to work around an
% issue where the colormap and tick marks on the colorbar would not get set
% properly. For some reason, it seems to fix the problem when there are
% pauses inbetween the graphics commands.

  % Set colormap and create a colorbar
  %   Do this before setting the tick marks on the colorbar since the colormap
  %   call resets the tick marks to one per contour level which is too many.
  Cbar = contourcmap('hot', 'ColorBar', 'on', 'Location', 'horizontal' );
pause(1)
  Cmap = colormap('hot');
  Cmap = Cmap(end:-1:1,:); % reverse the entries
  colormap(Paxes,Cmap);

pause(1)
  % Move the colorbar down a bit so that it won't overlap with the x-axis labels
  set(Cbar, 'Position', CbarLoc);

pause(1)
  % Reduce the amount of tick marks so that the labels won't overlap
  set(Cbar, 'XTick', Cticks);
  set(Cbar, 'XTickLabel', CtickLabels);

pause(1)
  % Draw the coast after the contourfm call so that coast lines will appear on top
  % of the contour plot
  Coast = load('coast');
  plotm(Coast.lat, Coast.long, 'Color', 'k', 'LineWidth', LineW);

pause(1)
  % Mark Africa
  textm(12, -7, 'Africa', 'FontSize', 10, 'Color', 'k');

pause(1)
  % Mark Grid3
  linem(Grid3Lats, Grid3Lons, 'LineWidth', G3linew, 'LineStyle', '-', 'Color', G3color);
  textm( 4.5, -40, 'Grid3', 'Color', G3color, 'FontSize',  10);

pause(1)
  % Mark TS Debby location
  plotm( [ 17 ], [ -30 ], 'Color', 'r', 'LineStyle', 'none', 'Marker', 'x', 'LineWidth',  LineW, 'MarkerSize', 10);

pause(1)
  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end
