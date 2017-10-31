function [ ] = PlotFsFigWind()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  Fsize = 13;

  WmaxMeas = 'max';  % 'avg' -> use avg Vt measurement
                     % 'max' -> use max Vt measurement


  % Read in factor separation data
  % Arrays in files have the time series in columns
  % In the sims var:
  %     column    sim
  %        1      NSND  (S0)
  %        2      SND   (S1)
  %        3      NSD   (S2)
  %        4      SD    (S12)
  %
  % In the factors var:
  %     column   factor
  %        1       F0
  %        2       F1
  %        3       F2
  %        4       F12

  SimColors = { 'green' 'red' 'blue' 'black' };
  SimLegText = { 'NSND' 'SND' 'NSD' 'SD' };
  LegLoc = 'NorthEastOutside';

  InFname = 'DIAGS/fs_factor_time_series.h5';

  InVname = sprintf('/%s_wind_t_sims', WmaxMeas);
  fprintf('Reading %s (%s)\n', InFname, InVname);
  SIMS = squeeze(h5read(InFname, InVname));

  InVname = sprintf('/%s_wind_t_factors', WmaxMeas);
  fprintf('Reading %s (%s)\n', InFname, InVname);
  F = squeeze(h5read(InFname, InVname));
  FACTORS = [ squeeze(SIMS(:,1)) squeeze(F(:,2)) squeeze(F(:,3)) squeeze(F(:,4)) squeeze(SIMS(:,4)) ];

  T = squeeze(h5read(InFname, '/t_coords'))./3600 - 42; % sim time in hours

  % Plot: 2 panels (2x1)
  Fig = figure;

  % pull up the bottom of the top panel and push down the top
  % of the bottom panel in order to make room for labeling
  PlocInc = 0.01;
  switch(WmaxMeas)
    case 'avg'
      OutFile = sprintf('%s/FsFigWindAvgTseries.jpg', Pdir);
      SimYlim = [ 0 15 ];
    case 'max'
      OutFile = sprintf('%s/FsFigWindMaxTseries.jpg', Pdir);
      SimYlim = [ 0 25 ];
  end
  Ptitle = sprintf('V_t (%s): Sims', WmaxMeas);

  Paxes = subplot(2, 1, 1);
  Ploc = get(Paxes, 'Position');
  Ploc(2) = Ploc(2) + PlocInc;
  Ploc(3) = Ploc(3) * 0.97;    % get the ends of the plots to line up
  Ploc(4) = Ploc(4) - PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigTseries(Paxes, T, SIMS, SimColors, 'a', Ptitle, 'Speed (ms^-^1)', 'linear', SimYlim, Fsize, 1, 1, SimLegText, LegLoc);

  % Mark the PRESAL (10-30 h) and SAL (40-60 h) time periods
  line([ 10 30 ], [ 5 5 ], 'Color', 'k', 'LineWidth', 2);
  line([ 10 10 ], [ 4 6 ], 'Color', 'k', 'LineWidth', 2);
  line([ 30 30 ], [ 4 6 ], 'Color', 'k', 'LineWidth', 2);
  text(15, 6.5, 'PSAP', 'FontSize', Fsize);

  line([ 40 60 ], [ 5 5 ], 'Color', 'k', 'LineWidth', 2);
  line([ 40 40 ], [ 4 6 ], 'Color', 'k', 'LineWidth', 2);
  line([ 60 60 ], [ 4 6 ], 'Color', 'k', 'LineWidth', 2);
  text(48, 6.5, 'SAP', 'FontSize', Fsize);

  % Factors time series
  Ptitle = sprintf('V_t (%s): Factors', WmaxMeas);

  FacColors = { 'green' 'cyan' 'orange' 'magenta' 'black' };
  FacLegText = { 'NSND' 'F1' 'F2' 'F12' 'SD' };
  switch(WmaxMeas)
    case 'avg'
      FacYlim = [ -5 15 ];
    case 'max'
      FacYlim = [ -5 25 ];
  end

  Paxes = subplot(2, 1, 2);
  Ploc = get(Paxes, 'Position');
  Ploc(2) = Ploc(2) + PlocInc;
  Ploc(3) = Ploc(3) * 0.97;    % get the ends of the plots to line up
  Ploc(4) = Ploc(4) - PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigTseries(Paxes, T, FACTORS, FacColors, 'b', Ptitle, 'Speed (ms^-^1)', 'linear', FacYlim, Fsize, 1, 1, FacLegText, LegLoc);

  % add ref for zero
  line([ T(1) T(end) ], [ 0 0 ], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
