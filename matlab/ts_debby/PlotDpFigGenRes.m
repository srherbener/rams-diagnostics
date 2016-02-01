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
  Nt = length(SimTrackLons);

  % Hovmoller data for dust
  InFile = 'DIAGS/storm_hovmollers_TSD_SAL_DUST.h5';
  InVname = '/sal_aero_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  HOV_DUST = squeeze(h5read(InFile, InVname));
  Z    = squeeze(h5read(InFile, '/z_coords'))./1000;         % convert to km
  T    = squeeze(h5read(InFile, '/t_coords'))./3600 - 42;    % convert to sim time in hours

  % Time series of Md
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_aero_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MD = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MD = TS_SAL_MD .* 1e-12;  % convert to billions of kg


  % plot
  Fig = figure;

  Paxes = subplot(3,1,1);
  PlotDpFigTrack(Paxes, SimTrackLons, SimTrackLats, 'a', '', Fsize);
  
  Paxes = subplot(3,1,2);
  PlotDpFigDustHov(Paxes, T, Z, HOV_DUST, 'b', 'SAL', Fsize, 0, 1, 0, 1);

  Paxes = subplot(3,1,3);
  PlotDpFigTseries(Paxes, T, TS_SAL_MD, 'c', 'SAL', 'M_d (10^9 kg)', Fsize, 1, 1);
  
  OutFile = sprintf('%s/DpFigGenRes.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = PlotDpFigTrack(Paxes, SimTrackLons, SimTrackLats, Pmarker, Ptitle, Fsize)

  LegText = { 'NHC Best Track' 'Model' };
  
  LegendFsize = 12;
  LineW = 2;
  
  LatBounds = [ 5 26 ];
  LonBounds = [ -42 -8 ];
  
  NhcTrackLats = [ 12.6 13.4 14.2 14.9 15.7 16.7 17.6 18.4 19.2 20.1 20.9 ];
  NhcTrackLons = [ 23.9 25.3 26.7 28.1 29.5 31.0 32.4 33.9 35.5 37.1 38.7 ] * -1; 

  SalRegionLats = [ 12 24 24 12 12 ];
  SalRegionLons = [ 41 41 24 24 41 ] * -1; 

  set(Paxes, 'FontSize', Fsize);
  m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
  m_coast('color', 'k'); % k --> black
  m_grid('linestyle','none','box','fancy','tickdir','out');
  NhcTrack = m_line(NhcTrackLons, NhcTrackLats, 'linewi', LineW, 'color', 'k', 'linestyle', 'none', 'marker', '+');
  SimTrack = m_line(SimTrackLons, SimTrackLats, 'linewi', LineW);

  % Add SAL analysis region
  m_line(SalRegionLons, SalRegionLats, 'linewi', LineW, 'color', 'r');

  legend( [ NhcTrack SimTrack' ], LegText,'Location', 'EastOutside', 'FontSize', LegendFsize);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end
