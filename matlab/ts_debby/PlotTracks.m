function [ ] = PlotTracks(ConfigFile, Case)

  [ Config ] = ReadConfig(ConfigFile);
  
  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  %Case = 'TSD_3GRIDS';
  %Case = 'TSD_DRY_NODUST';
  Case = 'TSD_DRY_DUST';
  %Case = 'TSD_MOIST_NODUST';
  %Case = 'TSD_MOIST_DUST';

  %SimLabel = 'DN';
  SimLabel = 'DD';
  %SimLabel = 'MN';
  %SimLabel = 'MD';
  
  % for the TS Debby simulations:
  LatBounds = [ 5 26 ];
  LonBounds = [ -42 -8 ];
  
  % Aug 22, 2006, 6Z through Aug 24, 2006, 18Z (11 points)
  NhcTrackLats = [ 12.6 13.4 14.2 14.9 15.7 16.7 17.6 18.4 19.2 20.1 20.9 ];
  NhcTrackLons = [ 23.9 25.3 26.7 28.1 29.5 30.0 32.4 33.9 35.5 37.1 38.7 ] * -1; 
  
  Hfile = sprintf('HDF5/storm_center-%s-AS-2006-08-20-120000-g3.h5', Case);
  HdsetLon = 'press_cent_xloc';
  HdsetLat = 'press_cent_yloc';
  
  fprintf('Reading: %s\n', Hfile);
  fprintf('  Track longitude: %s\n', HdsetLon);
  fprintf('  Track latitude: %s\n', HdsetLat);
  SimTrackLons = squeeze(hdf5read(Hfile, HdsetLon));
  SimTrackLats = squeeze(hdf5read(Hfile, HdsetLat));
  
  % plot
  FigTracks = figure;
  
  Fsize = 22;
  LegendFsize = 15;

  LineW = 2;
  
  set(gca, 'FontSize', Fsize);
  m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
  m_coast('color', 'k'); % k --> black
  m_grid('linestyle','none','box','fancy','tickdir','out');
  NhcTrack = m_line(NhcTrackLons, NhcTrackLats, 'linewi', LineW, 'color', 'k', 'linestyle', 'none', 'marker', '+');
  SimTrack = m_line(SimTrackLons, SimTrackLats, 'linewi', LineW); 
  legend([ NhcTrack SimTrack ], 'NHC Best Track', SimLabel, 'Location', 'NorthWest');
  
%  T = title('(a)');
%  set(T, 'Units', 'Normalized');
%  set(T, 'HorizontalAlignment', 'Left');
%  Tpos = get(T, 'Position');
%  Tpos(1) = 0; % line up with left edge of plot area
%  set(T, 'Position', Tpos);
  title ('Storm Tracks');
  
  OutFile = sprintf('%s/TsDebbyTracks_%s.jpg', Pdir, Case);
  fprintf('Writing: %s\n', OutFile);
  saveas(FigTracks, OutFile);
  close(FigTracks);

end
