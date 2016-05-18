function [ ] = PlotFsFigAnalysisSetup()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 12;

  CaseList = {
      { 'TSD_SAL_DUST'       'SD'   'black' }
      { 'TSD_NONSAL_DUST'    'NSD'  'red'   }
      { 'TSD_SAL_NODUST'     'SND'  'blue'  }
      { 'TSD_NONSAL_NODUST'  'NSND' 'green' }
    };
  Nc = length(CaseList);

  % SD: Storm center locations
  InFname = 'HDF5/TSD_SAL_DUST/HDF5/storm_center-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  InVnameLon = '/press_cent_xloc';
  InVnameLat = '/press_cent_yloc';

  fprintf('Reading: %s (%s, %s)\n', InFname, InVnameLon, InVnameLat);

  SimTrackLons = squeeze(h5read(InFname, InVnameLon));
  SimTrackLats = squeeze(h5read(InFname, InVnameLat));

  % PreSAL and SAL wind structure
  for i = 1:Nc
    Case       = CaseList{i}{1};
    LegText{i} = CaseList{i}{2};
    Colors{i}  = CaseList{i}{3};

    InFname = sprintf('DIAGS/hist_meas_az_speed_%s.h5', Case);
    InVname = '/all_ps_speed_maxlev';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    PS_MAXW(:,i) = squeeze(h5read(InFname, InVname));

    if (i == 1)
      X = squeeze(h5read(InFname, '/x_coords'))./1000; % Radius in km
    end

    InVname = '/all_s_speed_maxlev';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    S_MAXW(:,i) = squeeze(h5read(InFname, InVname));
  end

  % plot
  Fig = figure;

  Paxes = subplot(3,1,1);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) + 0.15;
  Ploc(2) = Ploc(2) - 0.02;
  Ploc(4) = Ploc(4) + 0.05;
  set(Paxes, 'Position', Ploc);
  PlotTrackXsection(Paxes, SimTrackLons, SimTrackLats, 'a', '', Fsize);

  Xlim = [ 0 500 ];
  Ylim = [ 5  20 ];
  Paxes = subplot(3,1,2);
  PlotFsFigLine(Paxes, X, PS_MAXW, 'b', 'Pre-SAL', 'Radius (km)', Xlim, 0, 'Max Wind (m s^-^1)', Ylim, 1, Fsize, LegText, 'none', Colors);

  Paxes = subplot(3,1,3);
  PlotFsFigLine(Paxes, X, S_MAXW, 'c', 'SAL', 'Radius (km)', Xlim, 1, 'Max Wind (m s^-^1)', Ylim, 1, Fsize, LegText, 'NorthEast', Colors);

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

  Scolor = 'r';
  Pcolor = 'b';

  Strack = linem([ Slat1 Slat2 ], [ Slon1 Slon2 ], 'Color', Scolor, 'LineStyle', ':', 'LineWidth', LineW);
  textm(Slat1, Slon1+0.4, 'A', 'Color', Scolor, 'FontSize', Fsize-2);
  textm(Slat2, Slon2-2.1, 'B', 'Color', Scolor, 'FontSize', Fsize-2);
  Ptrack = linem([ Plat1 Plat2 ], [ Plon1 Plon2 ], 'Color', Pcolor, 'LineStyle', ':', 'LineWidth', LineW);
  textm(Plat1-1.3, Plon1-2.7, 'C', 'Color', Pcolor, 'FontSize', Fsize-2);
  textm(Plat2+1.3, Plon2-0.6, 'D', 'Color', Pcolor, 'FontSize', Fsize-2);

  legend( [ SimTrack' Strack' Ptrack' ], { 'SD' 'STRACK' 'PTRACK' },'Location', 'NorthEastOutside', 'FontSize', LegendFsize);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end
