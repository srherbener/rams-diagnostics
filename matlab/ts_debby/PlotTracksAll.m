function [ ] = PlotTracksAll(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);
  
  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Hdir = 'HDF5';
  
  Cases = {
   'TSD_DRY_DUST'
   'TSD_DRY_NODUST'
   'TSD_MOIST_DUST'
   'TSD_MOIST_NODUST'
   };
  
  LegText = {
   'NHC Best Track'
   'DD'
   'DN'
   'MD'
   'MN'
   };
  
  % for the TS Debby simulations:
  LatBounds = [ 5 26 ];
  LonBounds = [ -42 -8 ];
  
  NhcTrackLats = [ 12.6 13.4 14.2 14.9 15.7 16.7 17.6 18.4 19.2 20.1 20.9 ];
  NhcTrackLons = [ 23.9 25.3 26.7 28.1 29.5 31.0 32.4 33.9 35.5 37.1 38.7 ] * -1; 
  
  % read in the Vt data
  Nc = length(Cases);
  for icase = 1:Nc
    Hfile = sprintf('%s/storm_center-%s-AS-2006-08-20-120000-g3.h5', Hdir, Cases{icase});
    HdsetLon = 'press_cent_xloc';
    HdsetLat = 'press_cent_yloc';
    fprintf('Reading: %s\n', Hfile);
    fprintf('  Track longitude: %s\n', HdsetLon);
    fprintf('  Track latitude: %s\n', HdsetLat);
  
    % Slons and Slats are (t)
    Slons = squeeze(hdf5read(Hfile, HdsetLon));
    Slats = squeeze(hdf5read(Hfile, HdsetLat));
    if (icase == 1)
      Nt = size(Slons,1);
      SimTrackLons = zeros([ Nt Nc ]);
      SimTrackLats = zeros([ Nt Nc ]);
    end
  
    % lines go into columns
    SimTrackLons(:,icase) = Slons;
    SimTrackLats(:,icase) = Slats;
  end
  
  % plot
  FigTracks = figure;
  
  Fsize = 22;
  LegendFsize = 15;

  LineW = 2;
  
  set(gca, 'FontSize', 25);
  m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
  m_coast('color', 'k'); % k --> black
  m_grid('linestyle','none','box','fancy','tickdir','out');
  NhcTrack = m_line(NhcTrackLons, NhcTrackLats, 'linewi', LineW, 'color', 'k', 'linestyle', 'none', 'marker', '+');
  SimTrack = m_line(SimTrackLons, SimTrackLats, 'linewi', LineW);
  %title('Storm Tracks');
  legend( [ NhcTrack SimTrack' ], LegText,'Location', 'SouthWest', 'FontSize', LegendFsize);
  
  OutFile = sprintf('%s/TsDebbyTracksAll.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(FigTracks, OutFile);
  close(FigTracks);

end
