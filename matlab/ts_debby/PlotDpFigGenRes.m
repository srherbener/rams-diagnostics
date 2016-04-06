function [ ] = PlotDpFigGenRes()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

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
  %SalPic = 'IMAGES/SAL_Aug23_12Z.png';
  SalPic = 'IMAGES/SAL_Aug23_12Z_grid3.png';

  % Vertically integrated dust mass
  VintDustFile = 'HDF5/TSD_SAL_DUST/HDF5/vint_dust-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  VintDustVname = '/vertint_dust';

  % Read in and create a snapshot of the vertically integrated dust field at Aug23, 12Z
  fprintf('Reading: %s (%s)\n', VintDustFile, VintDustVname);
  DUST = squeeze(h5read(VintDustFile, VintDustVname));
  X    = squeeze(h5read(VintDustFile, '/x_coords'));
  Y    = squeeze(h5read(VintDustFile, '/y_coords'));

  % DUST is (x,y,t)
  VINT_DUST = squeeze(DUST(:,:,61)); % t = 61 corresponds to Aug23, 12Z

  % plot
  Fig = figure;

  % Map doesn't quite fit in the vertical center of the subplot region. Needs to get
  % bumped upward just a little bit.
  Paxes = subplot(2,2,1);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.05;
  Ploc(2) = Ploc(2) + 0.03;
  set(Paxes, 'Position', Ploc);
  PlotDpFigTrack(Paxes, SimTrackLons, SimTrackLats, 'a', '', Fsize);
  
  % Initial dust profiles
  % Doesn't quite fit in the vertical center of the subplot region. Needs to get
  % bumped upward just a little bit.
  Paxes = subplot(2,2,2);
  Ploc = get(Paxes, 'Position');
  Ploc(2) = Ploc(2) + 0.03;
  set(Paxes, 'Position', Ploc);
  PlotDpFigInitDust(Paxes, DUST_CONC, Znamma, 'b', '', Fsize);

  % SAL strength image
  Paxes = subplot(2,2,3);
  PlaceSalImage(Paxes, SalPic, 'c', '12Z, Aug23', Fsize);

  % Vertically integrated dust mass
  Paxes = subplot(2,2,4);
  PlotVintDust(Paxes, X, Y, VINT_DUST, 'd', '12Z, Aug23', Fsize, 1, 1);

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
function [] = PlotDpFigTrack(Paxes, SimTrackLons, SimTrackLats, Pmarker, Ptitle, Fsize)

  axes(Paxes);

  LegText = { 'NHC' 'Model' };
  
  LegendFsize = 8;
  LineW = 2;
  
  LatBounds = [ 5 26 ];
  LonBounds = [ -42 -8 ];
  
  % The simulation runs from Aug 22, 6Z to Aug 24, 18Z.
  %
  % Locations from NHC report: Aug 22, 6Z through Aug 24, 18Z.
  NhcTrackLats = [ 12.6 13.4 14.2 14.9 15.7 16.7 17.6 18.4 19.2 20.1 20.9 ];
  NhcTrackLons = [ 23.9 25.3 26.7 28.1 29.5 31.0 32.4 33.9 35.5 37.1 38.7 ] * -1; 

  SalRegionLats = [ 12.0 24.0 24.0 12.0 12.0 ];
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
  textm(22, -35, 'SAL\_AR', 'FontSize', 10, 'Color', 'r');

  % Mark Africa
  textm(14, -12, 'Africa', 'FontSize', 10, 'Color', 'k', 'Rotation', 90);

  legend( [ NhcTrack SimTrack' ], LegText,'Location', 'SouthWest', 'FontSize', LegendFsize);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlaceSalImage(Paxes, ImFile, Pmarker, Ptitle, Fsize)

  axes(Paxes);

  set(Paxes, 'FontSize', Fsize);

%  LatBounds = [ 5 26 ];
%  LonBounds = [ -42 -8 ];
%  worldmap(LatBounds, LonBounds);
%  setm(Paxes, 'MapProjection', 'miller');

  % Axes position is the location of the axes within the entire figure.
  % 0 -> left side, or bottom
  % 1 -> right side, or top
  %  position -> [ left_x bottom_y width_x height_y ]
  Apos = get(Paxes, 'Position');
  ImgScale = 1.05; % make image 12% larger
  Apos(1) = 0.05;                % shift to the left
  Apos(2) = Apos(2) - 0.02;      % shift downward
  Apos(3) = Apos(3) * ImgScale;
  Apos(4) = Apos(4) * ImgScale;

  set(Paxes, 'Position', Apos);

  imshow(ImFile);

%  % create a gap so that image lines up with vint_dust image
%  Gap = 0.05;
%  Ploc = get(Paxes, 'Position');
%
%  % move the plot up a bit
%  Ploc(2) = Ploc(2) + Gap;
%  Ploc(4) = Ploc(4) - Gap;
%  set(Paxes, 'Position', Ploc);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlotVintDust(Paxes, X, Y, Z, Pmarker, Ptitle, Fsize, ShowX, ShowY)

  LatBounds = [ 5 26 ];
  LonBounds = [ -42 -8 ];

  LineW = 2;

  axes(Paxes);

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
  
  Clevs = [ -2:0.5:6 ];
  Cticks = [ 2.5 6.5 10.5 14.5 ];
  CtickLabels = { '10^-^1' '10^1' '10^3' '10^5' };
  contourfm(double(Y), double(X), double(log10(Z))', Clevs, 'LineStyle', 'none'); 
  caxis(Paxes, [ -2 6 ]);

  % "parula" is the default colormap
  Cbar = contourcmap('parula', 'ColorBar', 'on', 'Location', 'horizontal' );
  set(Cbar, 'Position', CbarLoc);
  set(Cbar, 'Xtick', Cticks);
  set(Cbar, 'XTickLabel', CtickLabels);

  % Draw the coast after the contourfm call so that coast lines will appear on top
  % of the contour plot
  Coast = load('coast');
  plotm(Coast.lat, Coast.long, 'Color', 'k', 'LineWidth', LineW);

  % Mark Africa
  textm(14, -12, 'Africa', 'FontSize', 10, 'Color', 'k', 'Rotation', 90);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end
