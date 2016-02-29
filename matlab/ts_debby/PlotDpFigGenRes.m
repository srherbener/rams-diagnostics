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

  % Hovmoller data for dust
  InFile = 'DIAGS/storm_hovmollers_TSD_SAL_DUST.h5';
  InVname = '/sal_aero_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  HOV_DUST = squeeze(h5read(InFile, InVname));
  Z    = squeeze(h5read(InFile, '/z_coords'))./1000;         % convert to km
  T    = squeeze(h5read(InFile, '/t_coords'))./3600 - 42;    % convert to sim time in hours

  % For plotting track, use simulation track data from Aug 22, 18Z
  % through Aug 24, 18Z. This corresponds to the best track locations
  % that reside inside the SAL analysis region. These times correspond
  % to simulation times 12 h and 60 h.
  T1 = find(T >= 12, 1, 'first');
  T2 = find(T <= 60, 1, 'last');

  % Time series of Md
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_aero_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MD = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MD = TS_SAL_MD .* 1e-12;  % convert to Tg


  % plot
  Fig = figure;

  % Map doesn't quite fit in the vertical center of the subplot region. Needs to get
  % bumped upward just a little bit.
  Paxes = subplot(3,2,1);
  Ploc = get(Paxes, 'Position');
  Ploc(2) = Ploc(2) + 0.03;
  set(Paxes, 'Position', Ploc);
  PlotDpFigTrack(Paxes, SimTrackLons(T1:T2), SimTrackLats(T1:T2), 'a', '', Fsize);
  
  % Initial dust profiles
  % Doesn't quite fit in the vertical center of the subplot region. Needs to get
  % bumped upward just a little bit.
  Paxes = subplot(3,2,2);
  Ploc = get(Paxes, 'Position');
  Ploc(2) = Ploc(2) + 0.03;
  set(Paxes, 'Position', Ploc);
  PlotDpFigInitDust(Paxes, DUST_CONC, Znamma, 'b', '', Fsize);

  % aero mass hovmoller
  Paxes = subplot(3,2,[3 4]);
  PlotDpFigDustHov(Paxes, T, Z, HOV_DUST, 'c', 'SAL\_AR: M_d', Fsize, 0, 1, 0, 2);

  % Total aero mass time series
  Paxes = subplot(3,2,[5 6]);
  PlotDpFigTseries(Paxes, T, TS_SAL_MD, 'd', 'SAL\_AR', 'M_d (Tg)', Fsize, 1, 1, { }, 'none');
  
  OutFile = sprintf('%s/DpFigGenRes.jpg', Pdir);
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
  
  % East side of SAL region was picked to be 24W. The first Best Track inside the SAL
  % region corresponds to Aug 22, 18Z. The simulation runs from Aug 22, 6Z to Aug 24, 18Z.
  %
  % Plot out the common time between the Best Track data and the simulation data. This would
  % start at Aug 22, 18Z and end at Aug 24, 18Z.
  
  % Locations from NHC report: Aug 22, 18Z through Aug 24, 18Z.
  NhcTrackLats = [ 14.2 14.9 15.7 16.7 17.6 18.4 19.2 20.1 20.9 ];
  NhcTrackLons = [ 26.7 28.1 29.5 31.0 32.4 33.9 35.5 37.1 38.7 ] * -1; 

  SalRegionLats = [ 12 24 24 12 12 ];
  SalRegionLons = [ 41 41 24 24 41 ] * -1; 

  set(Paxes, 'FontSize', Fsize);
%  m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
%  m_coast('color', 'k'); % k --> black
%  m_grid('linestyle','none','box','fancy','tickdir','out');
  worldmap(LatBounds, LonBounds);
  load coast;  % creates vars "lat" and "long" that contain coast lines
  plotm(lat, long, 'Color', 'k', 'LineWidth', LineW);

%  NhcTrack = m_line(NhcTrackLons, NhcTrackLats, 'linewi', LineW, 'color', 'k', 'linestyle', 'none', 'marker', '+');
%  SimTrack = m_line(SimTrackLons, SimTrackLats, 'linewi', LineW);
  NhcTrack = linem(NhcTrackLats, NhcTrackLons, 'LineWidth', LineW, 'Color', 'k', 'LineStyle', 'none', 'Marker', '+');
  SimTrack = linem(double(SimTrackLats), double(SimTrackLons), 'LineWidth', LineW);

  % Add SAL analysis region
%  m_line(SalRegionLons, SalRegionLats, 'linewi', LineW, 'color', 'r');
  linem(SalRegionLats, SalRegionLons, 'LineWidth', LineW, 'Color', 'r');
  textm(22, -35, 'SAL\_AR', 'FontSize', 10, 'Color', 'r');

  % Mark Africa
  textm(14, -12, 'Africa', 'FontSize', 10, 'Color', 'k', 'Rotation', 90);

  legend( [ NhcTrack SimTrack' ], LegText,'Location', 'SouthEast', 'FontSize', LegendFsize);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end
