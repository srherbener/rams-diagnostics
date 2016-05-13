function [ ] = PlotFsFigAnalysisSetup()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 25;

  % SD: Storm center locations
  Hfile = 'HDF5/TSD_SAL_DUST/HDF5/storm_center-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  HdsetLon = '/press_cent_xloc';
  HdsetLat = '/press_cent_yloc';

  fprintf('Reading: %s\n', Hfile);
  fprintf('  Track longitude: %s\n', HdsetLon);
  fprintf('  Track latitude: %s\n', HdsetLat);

  SimTrackLons = squeeze(h5read(Hfile, HdsetLon));
  SimTrackLats = squeeze(h5read(Hfile, HdsetLat));

  % plot
  Fig = figure;

  PlotTrackXsection(gca, SimTrackLons, SimTrackLats, 'a', '', Fsize);

  OutFile = sprintf('%s/FsFigAnalysisSetup.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = PlotTrackXsection(Paxes, SimTrackLons, SimTrackLats, Pmarker, Ptitle, Fsize)

  axes(Paxes);

  LineW = 2;
  LegendFsize = 8;
  
  LatBounds = [   5  26 ];
  LonBounds = [ -42  -8 ];

  set(Paxes, 'FontSize', Fsize);
  worldmap(LatBounds, LonBounds)
  setm(Paxes, 'FontSize', Fsize);
  setm(Paxes, 'MapProjection', 'miller');
  load coast;
  plotm(lat, long, 'Color', 'k', 'LineWidth', LineW);

  SimTrack = linem(double(SimTrackLats), double(SimTrackLons), 'Color', 'k', 'LineStyle', '-', 'LineWidth', LineW);

  % Mark Africa
  textm(14, -12, 'Africa', 'FontSize', 10, 'Color', 'k', 'Rotation', 90);

  % Mark cross section lines
  % strack - along storm track
  % ptrack - perpendicular to storm track
  Slat1 =  13.0;
  Slon1 = -21.3;
  Slat2 =  20.1;
  Slon2 = -39.6;

  Plat1 =   7.4;
  Plon1 = -30.0;
  Plat2 =  23.0;
  Plon2 = -24.0;

  StrackLats = [  13.0  20.1 ];
  StrackLons = [ -21.3 -39.6 ];
  PtrackLats = [   7.4  23.0 ];
  PtrackLons = [ -30.0 -24.0 ];

  Scolor = 'r';
  Pcolor = 'b';

  Strack = linem([ Slat1 Slat2 ], [ Slon1 Slon2 ], 'Color', Scolor, 'LineStyle', '--', 'LineWidth', LineW);
  textm(Slat1, Slon1+0.1, 'A', 'Color', Scolor, 'FontSize', Fsize);
  textm(Slat2, Slon2-1.0, 'B', 'Color', Scolor, 'FontSize', Fsize);
  Ptrack = linem([ Plat1 Plat2 ], [ Plon1 Plon2 ], 'Color', Pcolor, 'LineStyle', '--', 'LineWidth', LineW);
  textm(Plat1-1.3, Plon1-0.9, 'C', 'Color', Pcolor, 'FontSize', Fsize);
  textm(Plat2+0.3, Plon2-0.6, 'D', 'Color', Pcolor, 'FontSize', Fsize);

  legend( [ SimTrack' Strack' Ptrack' ], { 'SD' 'STRACK' 'PTRACK' },'Location', 'SouthWest', 'FontSize', LegendFsize);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end
