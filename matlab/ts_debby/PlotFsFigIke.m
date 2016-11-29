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

  % Read in the ike time series
  for icase = 1:Ncases
    Case  = CaseList{icase}{1};
    Label = CaseList{icase}{2};
    Color = CaseList{icase}{3};

    InFname = sprintf('DIAGS/storm_meas_%s.h5', Case);
    InVname = '/ike';
    fprintf('Reading %s (%s)\n', InFname, InVname);
    if (icase == 1)
      T = squeeze(h5read(InFname, '/t_coords'))./3600 - 42; % convert to sim time in hours starting with zero
      Nt = length(T);
      IKE_TS = zeros([Nt Ncases]);
    end
    IKE_TS(:,icase)   = squeeze(h5read(InFname, InVname)).*1e-7; % scale for plot
    IkeLegText{icase} = Label;
    IkeColors{icase}  = Color;
  end
  IkeLegLoc  = 'NorthEastOutside';
  fprintf('\n');

  % Read in factor separation data
  InFname = 'DIAGS/fs_factors.h5';

  InVname = '/s_ike_bar_frac_change';
  fprintf('Reading %s (%s)\n', InFname, InVname);
  SAL_FS_IKE = squeeze(h5read(InFname, InVname));

  InVname = '/s_ike_bar_int';
  fprintf('Reading %s (%s)\n', InFname, InVname);
  SAL_FS_IKE_INT = squeeze(h5read(InFname, InVname));

  InVname = '/ps_ike_bar_frac_change';
  fprintf('Reading %s (%s)\n', InFname, InVname);
  PRESAL_FS_IKE = squeeze(h5read(InFname, InVname));

  InVname = '/ps_ike_bar_int';
  fprintf('Reading %s (%s)\n', InFname, InVname);
  PRESAL_FS_IKE_INT = squeeze(h5read(InFname, InVname));

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

  BDATA = [ PRESAL_FS_IKE(1:3)' PRESAL_FS_IKE_INT(4); SAL_FS_IKE(1:3)' SAL_FS_IKE_INT(4) ]' .* 100;
  BarColors = { 'dodgerblue' 'cyan' };
  BarLabels = { 'NS' 'ND' 'NSND' 'INT' };
  BarLegText = { 'PSAP' 'SAP' };
  BarLegLoc = 'NorthEastOutside';
  fprintf('\n');

  % Plot: 2 panels (2x1)
  OutFile = sprintf('%s/FsFigIke.jpg', Pdir);
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
  PlotFsFigTseries(Paxes, T, IKE_TS, IkeColors, 'a', '', 'IKE (J X 10^7)', 'linear', [ 0 5 ], Fsize, 1, 1, IkeLegText, IkeLegLoc);

  % Mark the PRESAL (10-30 h) and SAL (40-60 h) time periods
  line([ 10 30 ], [ 1    1    ], 'Color', 'k', 'LineWidth', 2);
  line([ 10 10 ], [ 0.75 1.25 ], 'Color', 'k', 'LineWidth', 2);
  line([ 30 30 ], [ 0.75 1.25 ], 'Color', 'k', 'LineWidth', 2);
  text(15, 1.25, 'PSAP', 'FontSize', Fsize);

  line([ 40 60 ], [ 1    1    ], 'Color', 'k', 'LineWidth', 2);
  line([ 40 40 ], [ 0.75 1.25 ], 'Color', 'k', 'LineWidth', 2);
  line([ 60 60 ], [ 0.75 1.25 ], 'Color', 'k', 'LineWidth', 2);
  text(48, 1.25, 'SAP', 'FontSize', Fsize);

  Paxes = subplot(2, 1, 2);
  Ploc = get(Paxes, 'Position');
  Ploc(4) = Ploc(4) - PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigBgraph(Paxes, BDATA, BarColors, 'b', '', 'Factor', BarLabels, 'Difference (%)', [ -40 0 ], Fsize, 1, 1, BarLegText, BarLegLoc );

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
