function [ DataSpecs LegText DSokay ] = GenDataSpecs(Config, Case, i_panel, i_pset, Indent)
% GenDataSpecs function to generate data specfifications for plots

  Smooth  = Config.FigPanels(i_panel).Smooth;
  Flength = Config.FigPanels(i_panel).Flength;
  Nlines  = Config.PlotSets(i_pset).Nlines;

  DSokay = 1;

  for i_line = 1:Nlines
    i_pdata = Config.PlotSets(i_pset).Lines(i_line).PDnum;
    if (i_pdata == 0)
      fprintf('WARNING: skipping FigPanel number %d due to no associated PlotData on PlotSet number %d, Line number %d\n', i_panel, i_pset, i_line)
      DSokay = 0;
      continue;
    end
 
    % Line properties
    DataSpecs(i_line).Lwidth  = Config.PlotSets(i_pset).Lines(i_line).Lwidth;
    DataSpecs(i_line).Lcolor  = Config.PlotSets(i_pset).Lines(i_line).Lcolor;
    DataSpecs(i_line).Lstyle  = Config.PlotSets(i_pset).Lines(i_line).Lstyle;
    DataSpecs(i_line).Lgscale = Config.PlotSets(i_pset).Lines(i_line).Lgscale;

    % Legend text for this line
    LegText{i_line} = Config.PlotSets(i_pset).Lines(i_line).Legend;

    % Read specs for X and Y data, substitute for case name in file names
    Xvar    = Config.PlotData(i_pdata).Xvar;
    Xfile   = regexprep(Config.PlotData(i_pdata).Xfile, '<CASE>', Case);
    Xselect = Config.PlotData(i_pdata).Xselect;
    Xscale  = Config.PlotData(i_pdata).Xscale;
    Xoffset = Config.PlotData(i_pdata).Xoffset;
    Xzoom   = Config.PlotSets(i_pset).Lines(i_line).Xzoom;

    Yvar    = Config.PlotData(i_pdata).Yvar;
    Yfile   = regexprep(Config.PlotData(i_pdata).Yfile, '<CASE>', Case);
    Yselect = Config.PlotData(i_pdata).Yselect;
    Yscale  = Config.PlotData(i_pdata).Yscale;
    Yoffset = Config.PlotData(i_pdata).Yoffset;
    Yzoom   = Config.PlotSets(i_pset).Lines(i_line).Yzoom;

    % Read in and process the data
    % X data
    fprintf('%sReading: %s (%s)\n', Indent, Xfile, Xvar);
    Xdata = ReadSelectHdf5(Xfile, Xvar, Xselect);
    Xdata(Xdata == Config.UndefVal) = nan;
    if (strcmp(Smooth, 'x') || strcmp(Smooth, 'xy'))
      Xdata = SmoothFillTseries(Xdata, length(Xdata), Flength);
    end
    DataSpecs(i_line).Xdata = (Xdata .* Xscale) + Xoffset;

    % Y data
    fprintf('%sReading: %s (%s)\n', Indent, Yfile, Yvar);
    Ydata = ReadSelectHdf5(Yfile, Yvar, Yselect);
    Ydata(Ydata == Config.UndefVal) = nan;
    if (strcmp(Smooth, 'y') || strcmp(Smooth, 'xy'))
      Ydata = SmoothFillTseries(Ydata, length(Ydata), Flength);
    end
    DataSpecs(i_line).Ydata = (Ydata .* Yscale) + Yoffset;

  end
end
