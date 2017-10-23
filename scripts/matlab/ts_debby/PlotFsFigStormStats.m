function [ ] = PlotFsFigStats()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 11;

  CaseList = {
    { 'TSD_SAL_DUST'       'SD'   'crimson' }
%    { 'TSD_SAL_DUST'       'SD'   'black' }
%    { 'TSD_NONSAL_DUST'    'NSD'  'red'   }
%    { 'TSD_SAL_NODUST'     'SND'  'blue'  }
%    { 'TSD_NONSAL_NODUST'  'NSND' 'green' }
    };
  Nc = length(CaseList);
    
  GlocFname = 'GridLocs.txt';
  Grid3Loc = ReadGridLoc(GlocFname, 'Grid3:');

  MapInc = 4;
  MapLoc = Grid3Loc + [ -MapInc MapInc -MapInc MapInc ];

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

    % Max wind
    InFname = sprintf('DIAGS/storm_meas_%s.h5', Case);
    InVname = '/max_wind';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    MaxWind(:,ic) = squeeze(h5read(InFname, InVname));
    SimTime = squeeze(h5read(InFname, '/t_coords'))./3600 - 42; % Sim time starting with zero

    % Min pressure
    InVname = '/min_slp';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    MinPress(:,ic) = squeeze(h5read(InFname, InVname));

    LegText{ic} = Label;
    LineColors{ic} = Color;
  end

  % plot
  Fig = figure;

  Paxes = subplot(2,2,[ 1 2 ]);
  Ploc = get(Paxes, 'Position');
  Ploc(2) = Ploc(2) - 0.07;
  Ploc(4) = Ploc(4) + 0.10;
  set(Paxes, 'Position', Ploc);
  PlotDpFigTrack(Paxes, MapLoc, Grid3Loc, SimTrackLons, SimTrackLats, LegText, LineColors, 'a', '', Fsize);

  Paxes = subplot(2,2,3);
  PlotFsFigMaxWind(Paxes, SimTime, MaxWind, 'b', '', 1, 'Max Wind (m s^-^1)', [ 0 25 ], 1, Fsize, LegText, 'SouthEast', LineColors);

  Paxes = subplot(2,2,4);
  PlotFsFigMinPress(Paxes, SimTime, MinPress, 'c', '', 1, 'Min Press (mb)', [ 995 1010 ], 1, Fsize, LegText, 'none', LineColors);
  
  OutFile = sprintf('%s/FsFigStormStats.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
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
  textm(8, -39, 'Grid3', 'FontSize', 12, 'Color', Grid3Color);

  legend( [ NhcTrack SimTrack ], LegText,'Location', 'EastOutside', 'FontSize', LegendFsize);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(Paxes, T);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = PlotFsFigMaxWind(Paxes, T, MaxWind, Pmarker, Ptitle, Xshow, Ylabel, Ylim, Yshow, Fsize, LegText, LegLoc, Colors)

  axes(Paxes);

  LineW = 2;
  LegendFsize = 8;

  [ Npts Nlines ] = size(MaxWind);

  % NHC Best Track (every six hours) data
  % time step 1 from the simulation is where the NHC data starts
  % NHC wind is in mph, multiply by 0.447 to convert to m/s
  %
  % TAFB - Tropical Analysis and Forecast Branch
  % TAFB wind is in kts, multiply by 0.514 to convert to m/s
  %
  % time == 0 is Aug 22, 2006 at 06Z
  %
  NhcWind  = [ 35 35 35 40 45 50 50 50 50 50 50 ] * 0.447;
  TafbWind  = [ 30 30 30 35 35 35 35 35 35 35 35 ] * 0.514;
  NhcTimes = (0:6:60);

  set(Paxes, 'FontSize', Fsize);
  set(Paxes, 'LineWidth', LineW);
  set(Paxes, 'TickLength', [ 0.025 0.025 ]);
  for i = 1:Nlines
    Color = str2rgb(Colors{i});
    SimLines(i) = line(T, squeeze(MaxWind(:,i)), 'Color', Color, 'LineWidth', LineW);
  end

  NhcLine  = line(NhcTimes, NhcWind, 'Color', 'k', 'LineStyle', 'none', 'Marker', '+', 'LineWidth', LineW);
  TafbLine = line(NhcTimes, TafbWind, 'Color', 'k', 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', LineW);

  set(Paxes, 'Xtick', [ 6 18 30 42 54 ]);
  if (Xshow > 0)
    set(Paxes, 'XtickLabel', { ' 12Z\newline22Aug' '  0Z\newline23Aug' ' 12Z\newline23Aug' '  0Z\newline24Aug', ' 12Z\newline24Aug' });
  else
    set(Paxes, 'XtickLabel', {});
  end

  ylim(Ylim);
  if (Yshow > 0)
    ylabel(Ylabel)
  else
    set(Paxes, 'YtickLabel', {});
  end

  if (~strcmp(LegLoc, 'none'))
    LegLines = [ NhcLine TafbLine SimLines ];
    LegLabels = { 'NHC' 'TAFB' LegText{1:end} };
    legend(LegLines, LegLabels, 'Location', LegLoc, 'FontSize', LegendFsize);
  end

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(Paxes, T);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = PlotFsFigMinPress(Paxes, T, MinPress, Pmarker, Ptitle, Xshow, Ylabel, Ylim, Yshow, Fsize, LegText, LegLoc, Colors)

  axes(Paxes);

  LineW = 2;
  LegendFsize = 8;

  [ Npts Nlines ] = size(MinPress);

  % NHC Best Track (every six hours) data
  % time step 1 from the simulation is where the NHC data starts
  %
  % TFAB - Tropical Analysis and Forecast Branch
  %
  % time == 0 is Aug 22, 2006 at 06Z
  %
  NhcPress = [ 1007 1007 1005 1003 1002 1001 1001 1000  999 1000 1000 ];
  TafbPress = [ 1009 1009 1009 1005 1005 1005 1005 1005 1005 1005 1005 ];
  NhcTimes = (0:6:60);

  set(Paxes, 'FontSize', Fsize);
  set(Paxes, 'LineWidth', LineW);
  set(Paxes, 'TickLength', [ 0.025 0.025 ]);
  for i = 1:Nlines
    Color = str2rgb(Colors{i});
    SimLines(i) = line(T, squeeze(MinPress(:,i)), 'Color', Color, 'LineWidth', LineW);
  end

  NhcLine  = line(NhcTimes, NhcPress, 'Color', 'k', 'LineStyle', 'none', 'Marker', '+', 'LineWidth', LineW);
  TafbLine = line(NhcTimes, TafbPress, 'Color', 'k', 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', LineW);

  set(Paxes, 'Xtick', [ 6 18 30 42 54 ]);
  if (Xshow > 0)
    set(Paxes, 'XtickLabel', { ' 12Z\newline22Aug' '  0Z\newline23Aug' ' 12Z\newline23Aug' '  0Z\newline24Aug', ' 12Z\newline24Aug' });
  else
    set(Paxes, 'XtickLabel', {});
  end

  ylim(Ylim);
  if (Yshow > 0)
    ylabel(Ylabel)
  else
    set(Paxes, 'YtickLabel', {});
  end

  if (~strcmp(LegLoc, 'none'))
    LegLines = [ NhcLine TafbLine SimLines ];
    LegLabels = { 'NHC' 'TAFB' LegText{1:end} };
    legend(LegLines, LegLabels, 'Location', LegLoc, 'FontSize', LegendFsize);
  end

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(Paxes, T);
  end
end
