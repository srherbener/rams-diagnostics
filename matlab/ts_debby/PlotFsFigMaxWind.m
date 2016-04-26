function [ ] = PlotFsFigMaxWind()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  CaseList = {
    { 'TSD_SAL_DUST'      'SD'   'black' }
    { 'TSD_SAL_NODUST'    'SND'  'blue'  }
    { 'TSD_NONSAL_DUST'   'NSD'  'red'   }
    { 'TSD_NONSAL_NODUST' 'NSND' 'green' }
    };
  Ncases = length(CaseList);
  
  Fsize = 13;

  % Read in the max wind time series
  for icase = 1:Ncases
    Case  = CaseList{icase}{1};
    Label = CaseList{icase}{2};
    Color = CaseList{icase}{3};

    InFname = sprintf('DIAGS/storm_meas_%s.h5', Case);
    InVname = '/max_wind';
    fprintf('Reading %s (%s)\n', InFname, InVname);
    if (icase == 1)
      T = squeeze(h5read(InFname, '/t_coords'))./3600 - 42; % convert to sim time in hours starting with zero
      Nt = length(T);
      WMAX_TS = zeros([Nt Ncases]);
    end
    WMAX_TS(:,icase)   = squeeze(h5read(InFname, InVname));
    WmaxLegText{icase} = Label;
    WmaxColors{icase}  = Color;
  end
  WmaxLegLoc  = 'NorthEastOutside';
  fprintf('\n');

  % Read in factor separation data
  InFname = 'DIAGS/fs_factors.h5';

  InVname = '/s_max_wind_bar_frac_change';
  fprintf('Reading %s (%s)\n', InFname, InVname);
  SAL_FS_WMAX = squeeze(h5read(InFname, InVname));

  InVname = '/s_max_wind_bar_int';
  fprintf('Reading %s (%s)\n', InFname, InVname);
  SAL_FS_WMAX_INT = squeeze(h5read(InFname, InVname));

  InVname = '/ps_max_wind_bar_frac_change';
  fprintf('Reading %s (%s)\n', InFname, InVname);
  PRESAL_FS_WMAX = squeeze(h5read(InFname, InVname));

  InVname = '/ps_max_wind_bar_int';
  fprintf('Reading %s (%s)\n', InFname, InVname);
  PRESAL_FS_WMAX_INT = squeeze(h5read(InFname, InVname));

  % In the fs_factors.h5 file, the four numbers of interest
  % are split up into two vectors so that plotting results
  % in one set of bars in which two different colors are used.
  %
  % This is accomplished by plotting the bars as stacked and 
  % placing nans (which get plotted, but are not noticeable) in
  % the appropriate places so that you don't see any stacking.
  %
  %     *_frac_change --> [ S1  S2  S12 nan ]  (simulation averages)
  %     *_int         --> [ nan nan nan F12 ]  (factors)
  %
  % For this plot combine the two vectors and place both presal and sal data
  % on grouped bar graph.

  BDATA = [ PRESAL_FS_WMAX(1:3)' PRESAL_FS_WMAX_INT(4); SAL_FS_WMAX(1:3)' SAL_FS_WMAX_INT(4) ]' .* 100;
  BarColors = { 'dodgerblue' 'cyan' };
  BarLabels = { 'NS' 'ND' 'NSD' 'INT' };
  BarLegText = { 'Pre-SAL' 'SAL' };
  BarLegLoc = 'NorthEastOutside';
  fprintf('\n');

  % Plot: 2 panels (2x1)
  OutFile = sprintf('%s/FsFigMaxWind.jpg', Pdir);
  Fig = figure;

  % pull up the bottom of the top panel and push down the top
  % of the bottom panel in order to make room for labeling
  PlocInc = 0.01;

  Paxes = subplot(2, 1, 1);
  Ploc = get(Paxes, 'Position');
  Ploc(2) = Ploc(2) + PlocInc;
  Ploc(3) = Ploc(3) * 0.97;    % get the ends of the plots to line up
  Ploc(4) = Ploc(4) - PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigTseries(Paxes, T, WMAX_TS, WmaxColors, 'a', '', 'Speed (ms^-^1)', 'linear', [ 10 25  ], Fsize, 1, 1, WmaxLegText, WmaxLegLoc);

  Paxes = subplot(2, 1, 2);
  Ploc = get(Paxes, 'Position');
  Ploc(4) = Ploc(4) - PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigBgraph(Paxes, BDATA, BarColors, 'b', '', 'Factor', BarLabels, 'Difference (%)', [ -20 20 ], Fsize, 1, 1, BarLegText, BarLegLoc );

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
