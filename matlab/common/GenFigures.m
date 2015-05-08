function [ Config AxesSpecs DataSpecs LegSpecs ] = GenFigures(ConfigFile)
% GenFigures function to generate figures

  [ Config ] = ReadFigureConfig(ConfigFile);
  UndefVal = Config.UndefVal;
  Case = 'RCE_BASE';

  fprintf('Generating Figures: %s\n', ConfigFile);

  for i_panel = 1:length(Config.FigPanels)
    % clear out legend structure
    clear LegendSpecs;

    fprintf('  Panel: %s\n', Config.FigPanels(i_panel).Name);

    % Check that plot data sets and axes got associated
    i_paxis = Config.FigPanels(i_panel).PAnum;
    i_pset = Config.FigPanels(i_panel).PSnum;
    if (i_paxis == 0)
      fprintf('WARNING: skipping FigPanel number %d due to no associated PlotAxes\n', i_panel)
      continue;
    end
    if (i_pset == 0)
      fprintf('WARNING: skipping FigPanel number %d due to no associated PlotSet\n', i_panel)
      continue;
    end

    % Fill in the AxesSpecs structure
    [ AxesSpecs ] = GenAxesSpecs(Config, i_panel, i_paxis);

    % Fill in the DataSpecs structure
    % Last arg is indent spacing for formatting messages.
    [ DataSpecs LegText DSokay ] = GenDataSpecs(Config, Case, i_panel, i_pset, '    ');
    fprintf('\n');
    if (DSokay == 0)
      continue;
    end

    % Legend specs
    LegSpecs.Text  = LegText;
    LegSpecs.Loc   = Config.FigPanels(i_panel).LegLoc;
    LegSpecs.Fsize = Config.FigPanels(i_panel).LegFsize;

  end
  fprintf('\n');

end
