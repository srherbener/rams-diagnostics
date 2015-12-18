function [ ] = PlotTracksAll()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Cases = {
%   'TSD_DRY_DUST'
%   'TSD_DRY_NODUST'
%   'TSD_MOIST_DUST'
%   'TSD_MOIST_NODUST'

   'TSD_SAL_DUST'
%   'TSD_SAL_NODUST'
%   'TSD_NONSAL_DUST'
%   'TSD_NONSAL_NODUST'
   };
  
  LegText = {
   'NHC Best Track'

%   'DRY\_DUST'
%   'DRY\_NODUST'
%   'MOIST\_DUST'
%   'MOIST\_NODUST'
%   'DRY\_DUST'
%   'DRY'
%   'DUST'
%   'BASE'

   'Model'
%   'SAL\_DUST'
%   'SAL\_NODUST'
%   'NONSAL\_DUST'
%   'NONSAL\_NODUST'
   };
   Nc = length(Cases);
  
  % for the TS Debby simulations:
  LatBounds = [ 5 26 ];
  LonBounds = [ -42 -8 ];
  
  NhcTrackLats = [ 12.6 13.4 14.2 14.9 15.7 16.7 17.6 18.4 19.2 20.1 20.9 ];
  NhcTrackLons = [ 23.9 25.3 26.7 28.1 29.5 31.0 32.4 33.9 35.5 37.1 38.7 ] * -1; 

  % first find the max number of time steps
  MaxNt = 0;
  for icase = 1:Nc
    Case = Cases{icase};
    Hfile = sprintf('HDF5/%s/HDF5/storm_center-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    Hinfo = h5info(Hfile, '/t_coords');
    Tsize = Hinfo.Dataspace.Size;
    if (Tsize > MaxNt)
      MaxNt = Tsize;
    end
  end

  % allocate array to hold the plot data
  % initialize to nan so that empty spaces won't produce lines on the plot
  SimTrackLons = nan([ MaxNt Nc ]);
  
  % read in the Vt data
  for icase = 1:Nc
    Case = Cases{icase};
    Hfile = sprintf('HDF5/%s/HDF5/storm_center-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    HdsetLon = '/press_cent_xloc';
    HdsetLat = '/press_cent_yloc';
    fprintf('Reading: %s\n', Hfile);
    fprintf('  Track longitude: %s\n', HdsetLon);
    fprintf('  Track latitude: %s\n', HdsetLat);
  
    % Slons and Slats are (t)
    Slons = squeeze(h5read(Hfile, HdsetLon));
    Slats = squeeze(h5read(Hfile, HdsetLat));
    Nt = size(Slons,1);
  
    % lines go into columns
    SimTrackLons(1:Nt,icase) = Slons;
    SimTrackLats(1:Nt,icase) = Slats;
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
