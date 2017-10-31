function [ ] = PlotFsFigIke()

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

  InVname = '/s_ike_bar_factors';
  fprintf('Reading %s (%s)\n', InFname, InVname);
  SAL_FS_IKE = squeeze(h5read(InFname, InVname)) .* 1e-7; % scale for plot

  InVname = '/ps_ike_bar_factors';
  fprintf('Reading %s (%s)\n', InFname, InVname);
  PRESAL_FS_IKE = squeeze(h5read(InFname, InVname)) .* 1e-7; % scale for plot

  % Plot: 2 panels (2x1)
  OutFile = sprintf('%s/FsFigIke.jpg', Pdir);
  Fig = figure;

  % pull up the bottom of the top panel and push down the top
  % of the bottom panel in order to make room for labeling
  PlocInc = 0.01;
  IkeYlim = [ 0 4 ];

  Paxes = subplot(2, 2, [ 1 2 ]);
  Ploc = get(Paxes, 'Position');
  Ploc(2) = Ploc(2) + PlocInc;
  Ploc(3) = Ploc(3) * 0.97;    % get the ends of the plots to line up
  Ploc(4) = Ploc(4) - PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigTseries(Paxes, T, IKE_TS, IkeColors, 'a', '', 'IKE (J X 10^7)', 'linear', IkeYlim, Fsize, 1, 1, IkeLegText, IkeLegLoc);

  % Mark the PRESAL (10-30 h) and SAL (40-60 h) time periods
  line([ 10 30 ], [ 1    1    ], 'Color', 'k', 'LineWidth', 2);
  line([ 10 10 ], [ 0.75 1.25 ], 'Color', 'k', 'LineWidth', 2);
  line([ 30 30 ], [ 0.75 1.25 ], 'Color', 'k', 'LineWidth', 2);
  text(15, 1.25, 'PSAP', 'FontSize', Fsize);

  line([ 40 60 ], [ 1    1    ], 'Color', 'k', 'LineWidth', 2);
  line([ 40 40 ], [ 0.75 1.25 ], 'Color', 'k', 'LineWidth', 2);
  line([ 60 60 ], [ 0.75 1.25 ], 'Color', 'k', 'LineWidth', 2);
  text(48, 1.25, 'SAP', 'FontSize', Fsize);


  BarColors = { 'dodgerblue' 'cyan' 'orange' 'yellow' 'red' };
  BarLabels = { 'NSND' 'F1' 'F2' 'F12' 'SD' };
  BarLegText = 'none';
  BarLegLoc = '';
  Bylim = [ 2 3.5 ];

  Paxes = subplot(2, 2, 3);
  Ploc = get(Paxes, 'Position');
  Ploc(4) = Ploc(4) - PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigBgraph(Paxes, PRESAL_FS_IKE, BarColors, 'b', '', 'Factor', BarLabels, 'Magnitude (X 10^7)', Bylim, Fsize, 1, 1, BarLegText, BarLegLoc );

  Paxes = subplot(2, 2, 4);
  Ploc = get(Paxes, 'Position');
  Ploc(4) = Ploc(4) - PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigBgraph(Paxes, SAL_FS_IKE, BarColors, 'c', '', 'Factor', BarLabels, 'Magnitude (X 10^7)', Bylim, Fsize, 1, 1, BarLegText, BarLegLoc );

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
