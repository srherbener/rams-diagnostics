function [ ] = GenFigures(ConfigFile)
% GenFigures function to generate figures

  [ Config ] = ReadFigureConfig(ConfigFile);

  fprintf('Generating Figures: %s\n', ConfigFile);

  for i_fig = 1:length(Config.Figures)
    fprintf('  Figure: %s\n', Config.Figures(i_fig).Name);

    i_cset = Config.Figures(i_fig).CSnum;
    if (i_cset == 0)
      fprintf('WARNING: skipping Figure number %d due to no associated CaseSet\n', i_fig)
      continue;
    end

    Npanels = Config.Figures(i_fig).Npanels;
    Prows   = Config.Figures(i_fig).Psize(1);
    Pcols   = Config.Figures(i_fig).Psize(2);

    % Create one figure for each case in the case list
    for i_case = 1:Config.CaseSets(i_cset).Ncases
      FigCase = Config.CaseSets(i_cset).Cases(i_case).Cname;
      fprintf('    Case: %s\n', FigCase); 
      fprintf('    Panel Rows: %d\n', Prows); 
      fprintf('    Panel Columns: %d\n', Pcols); 

      % Create the figure
      Fig = figure;

      % Construct the plots for each panel
      for i_fig_panel = 1:Npanels
        % Grab the current panel from the figure spec
        i_panel = Config.Figures(i_fig).Panels(i_fig_panel).FPnum;
        if (i_panel == 0)
          fprintf('WARNING: skipping Figure number %d, panel number %d due to no associated Panels\n', i_fig, i_fig_panel);
          continue;
        end

        Ploc = Config.Figures(i_fig).Panels(i_fig_panel).Ploc;
        fprintf('      Panel: %s\n', Config.FigPanels(i_panel).Name);

        % clear out legend structure
        clear LegendSpecs;
    
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
        [ AxesSpecs ] = GenAxesSpecs(Config, i_panel, i_paxis, FigCase);
    
        % Fill in the DataSpecs structure
        % Last arg is indent spacing for formatting messages.
        [ DataSpecs, LegText, DSokay ] = GenDataSpecs(Config, FigCase, i_panel, i_pset, '        ');
        fprintf('\n');
        if (DSokay == 0)
          continue;
        end
    
        % Legend specs
        LegendSpecs.Text  = LegText;
        LegendSpecs.Loc   = Config.FigPanels(i_panel).LegLoc;
        LegendSpecs.Fsize = Config.FigPanels(i_panel).LegFsize;

        % Create the panel
        Axes = subplot(Prows, Pcols, Ploc);
        DrawPlotData(Axes, DataSpecs, Config.PlotSets(i_pset).Type);
        SetPlotAxes(Axes, AxesSpecs);
        
        % Place the legend
        if (~strcmp(LegendSpecs.Loc, 'none'))
          legend(LegendSpecs.Text, 'Location', LegendSpecs.Loc, 'FontSize', LegendSpecs.Fsize);
        end
      end
      if (strcmp(FigCase, 'none'))
        OutFile = Config.Figures(i_fig).OutFile;
      else
        OutFile = regexprep(Config.Figures(i_fig).OutFile, '<CASE>', FigCase);
      end
      fprintf('      Writing: %s\n', OutFile);

      saveas(Fig, OutFile);
      close(Fig);

      fprintf('\n');
    end
  end
end
