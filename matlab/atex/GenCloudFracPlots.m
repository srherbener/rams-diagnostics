function [ ] = GenCloudFracPlots(ConfigFile)
% GenCloudFracPlots function to plot cloud fraction time series

  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;
  Pdir = Config.PlotDir;

  InFile = sprintf('%s/%s', Ddir, 'ensemble_ts.h5');

  % Plot specs
  %   Name
  %   MinVar
  %   MaxVar
  %   AvgVar
  %   Times   --> { start end }  in hours
  %   Title   --> { PanelMarker Label }
  %   Xaxis   --> { [Xmin Xmax] Xscale Xticks Xshow Xlabel }
  %   Yaxis   --> { [Ymin Ymax] Yscale Yticks Yshow Ylabel }
  %   Lcolors --> { min max avg }
  %   Lstyles --> { min max avg }
  %   Legend  --> { { min max avg } Location }
  %   OutFile

  PlotList = {
    % 
    {
    'S293'
    'cloud_frac_293_ens_min'
    'cloud_frac_293_ens_max'
    'cloud_frac_293_ens_avg'
    { 12 36 }
    { 'd' 'S293' }
    { [0 24]  'linear' [0 5 10 15 20]        1 'Analysis Time (h)' }
    { [0 1.1] 'linear' [0 0.2 0.4 0.6 0.8 1] 1 'CF'                }
    { 'blue' 'orangered' 'black' }
    { '-' '-' '-' }
    { { 'Min' 'Max' 'Avg' }  'SouthWest' }
    'ens_avg_cfrac_TALL_CO_S293.jpg'
    }

    {
    'S298'
    'cloud_frac_298_ens_min'
    'cloud_frac_298_ens_max'
    'cloud_frac_298_ens_avg'
    { 12 36 }
    { 'd' 'S298' }
    { [0 24]  'linear' [0 5 10 15 20]        1 'Analysis Time (h)' }
    { [0 1.1] 'linear' [0 0.2 0.4 0.6 0.8 1] 1 'CF'                }
    { 'blue' 'orangered' 'black' }
    { '-' '-' '-' }
    { { 'Min' 'Max' 'Avg' }  'SouthWest' }
    'ens_avg_cfrac_TALL_CO_S298.jpg'
    }

    {
    'S303'
    'cloud_frac_303_ens_min'
    'cloud_frac_303_ens_max'
    'cloud_frac_303_ens_avg'
    { 12 36 }
    { 'd' 'S303' }
    { [0 24]  'linear' [0 5 10 15 20]        1 'Analysis Time (h)' }
    { [0 1.1] 'linear' [0 0.2 0.4 0.6 0.8 1] 1 'CF'                }
    { 'blue' 'orangered' 'black' }
    { '-' '-' '-' }
    { { 'Min' 'Max' 'Avg' }  'SouthWest' }
    'ens_avg_cfrac_TALL_CO_S303.jpg'
    }

    };

  Nplots = length(PlotList);

  Lgscales = [];
  AddMeas = 'none';


  % make the plots
  for iplot = 1:Nplots
    Pname   = PlotList{iplot}{1};

    MinVar  = PlotList{iplot}{2};
    MaxVar  = PlotList{iplot}{3};
    AvgVar  = PlotList{iplot}{4};

    Tstart  = PlotList{iplot}{5}{1};
    Tend    = PlotList{iplot}{5}{2};

    Pmarker = { PlotList{iplot}{6}{1} };
    Plabel  = PlotList{iplot}{6}{2};

    Xlim   = PlotList{iplot}{7}{1};
    Xscale = PlotList{iplot}{7}{2};
    Xticks = PlotList{iplot}{7}{3};
    Xshow  = PlotList{iplot}{7}{4};
    Xlabel = PlotList{iplot}{7}{5};

    Ylim   = PlotList{iplot}{8}{1};
    Yscale = PlotList{iplot}{8}{2};
    Yticks = PlotList{iplot}{8}{3};
    Yshow  = PlotList{iplot}{8}{4};
    Ylabel = PlotList{iplot}{8}{5};

    Lcolors = PlotList{iplot}{9};
    Lstyles = PlotList{iplot}{10};

    LegText = PlotList{iplot}{11}{1};
    LegLoc  = PlotList{iplot}{11}{2};

    OutFile = sprintf('%s/%s', Pdir, PlotList{iplot}{12});

    fprintf ('Generating plot: %s\n', Pname);
    fprintf ('  Reading: %s (%s %s %s)\n', InFile, MinVar, MaxVar, AvgVar);
    fprintf('\n');

    % read in and organize data into arrays for Plot2dSet
    MIN = squeeze(hdf5read(InFile, MinVar));
    MAX = squeeze(hdf5read(InFile, MaxVar));
    AVG = squeeze(hdf5read(InFile, AvgVar));
    T   = squeeze(hdf5read(InFile, 't_coords'))./3600;

    T1 = find(T >= Tstart, 1, 'first');
    T2 = find(T <= Tend,   1, 'last');

    Xall = [ T(T1:T2)-12   T(T1:T2)-12   T(T1:T2)-12   ]'; % note transpose, subtract off the initial 12 hour period
    Yall = [ MIN(T1:T2) MAX(T1:T2) AVG(T1:T2) ]';

    % make the plot
    clear AxisProps;

    i_ap = 1;
    AxisProps(i_ap).Name = 'FontSize';
    AxisProps(i_ap).Val = 25;
    i_ap = i_ap + 1;

    AxisProps(i_ap).Name = 'Xlim';
    AxisProps(i_ap).Val = Xlim;
    i_ap = i_ap + 1;

    AxisProps(i_ap).Name = 'XTick';
    AxisProps(i_ap).Val = Xticks;
    i_ap = i_ap + 1;

    AxisProps(i_ap).Name = 'Ylim';
    AxisProps(i_ap).Val = Ylim;
    i_ap = i_ap + 1;

    AxisProps(i_ap).Name = 'YTick';
    AxisProps(i_ap).Val = Yticks;
    i_ap = i_ap + 1;

    fprintf ('  Writing: %s\n', OutFile);
    Fig = figure;
    set(Fig, 'Visible', 'off');
    Plot2dSet( Xall, Yall, Plabel, Pmarker, Xlabel, Ylabel, Lcolors, Lstyles, Lgscales, LegText, LegLoc, AxisProps, AddMeas, Fig );
    saveas(Fig, OutFile);
    close(Fig); 
    fprintf('\n');
  end

end
