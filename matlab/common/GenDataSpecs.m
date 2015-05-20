function [ DataSpecs, LegText, DSokay ] = GenDataSpecs(Config, FigCase, i_panel, i_pset, Indent)
% GenDataSpecs function to generate data specfifications for plots

  Smooth  = Config.FigPanels(i_panel).Smooth;
  Flength = Config.FigPanels(i_panel).Flength;
  Ndsets  = Config.PlotSets(i_pset).Ndsets;

  DSokay = 1;

  for i_dset = 1:Ndsets
    i_pdata = Config.PlotSets(i_pset).DataSets(i_dset).PDnum;
    if (i_pdata == 0)
      fprintf('WARNING: skipping FigPanel number %d due to no associated PlotData on PlotSet number %d, Line number %d\n', i_panel, i_pset, i_dset)
      DSokay = 0;
      continue;
    end
 
    % Line properties
    DataSpecs(i_dset).Lwidth  = Config.PlotSets(i_pset).DataSets(i_dset).Lwidth;
    DataSpecs(i_dset).Lcolor  = Config.PlotSets(i_pset).DataSets(i_dset).Lcolor;
    DataSpecs(i_dset).Lstyle  = Config.PlotSets(i_pset).DataSets(i_dset).Lstyle;
    DataSpecs(i_dset).Lgscale = Config.PlotSets(i_pset).DataSets(i_dset).Lgscale;

    % Legend text for this line
    LegText{i_dset} = Config.PlotSets(i_pset).DataSets(i_dset).Legend;

    % Line specific case
    % The line specific case takes precedence over the figure case
    LineCase = Config.PlotSets(i_pset).DataSets(i_dset).Case;
    if (strcmp(LineCase, 'none'))
      Case = FigCase;
    else
      Case = LineCase;
    end

    % Set the type of data (line or bar, for now)
    DataSpecs(i_dset).Ptype = Config.PlotData(i_pdata).Type;
    
    if (strcmp(DataSpecs(i_dset).Ptype, 'line') || ...
        strcmp(DataSpecs(i_dset).Ptype, 'bar'))
    
      % Get specs for x data, substitute for case name in file name
      Xvar    = Config.PlotData(i_pdata).Xvar;
      Xfile   = regexprep(Config.PlotData(i_pdata).Xfile, '<CASE>', Case);
      Xselect = Config.PlotData(i_pdata).Xselect;
      Xscale  = Config.PlotData(i_pdata).Xscale;
      Xoffset = Config.PlotData(i_pdata).Xoffset;
      Xzoom   = Config.PlotSets(i_pset).DataSets(i_dset).Xzoom;
      
      % Read in and process X data
      fprintf('%sReading: %s (%s)\n', Indent, Xfile, Xvar);
      Xdata = ReadSelectHdf5(Xfile, Xvar, Xselect);
      Xdata(Xdata == Config.UndefVal) = nan;
      if (strcmp(Smooth, 'x') || strcmp(Smooth, 'xy'))
        Xdata = SmoothFillTseries(Xdata, length(Xdata), Flength);
      end
      DataSpecs(i_dset).Xdata = (Xdata .* Xscale) + Xoffset;
    end

    if (strcmp(DataSpecs(i_dset).Ptype, 'line'))
        
      % Get specs for y data, substitute for case name in file name  
      Yvar    = Config.PlotData(i_pdata).Yvar;
      Yfile   = regexprep(Config.PlotData(i_pdata).Yfile, '<CASE>', Case);
      Yselect = Config.PlotData(i_pdata).Yselect;
      Yscale  = Config.PlotData(i_pdata).Yscale;
      Yoffset = Config.PlotData(i_pdata).Yoffset;
      Yzoom   = Config.PlotSets(i_pset).DataSets(i_dset).Yzoom;
      
      % Read in and process y data
      fprintf('%sReading: %s (%s)\n', Indent, Yfile, Yvar);
      Ydata = ReadSelectHdf5(Yfile, Yvar, Yselect);
      Ydata(Ydata == Config.UndefVal) = nan;
      if (strcmp(Smooth, 'y') || strcmp(Smooth, 'xy'))
        Ydata = SmoothFillTseries(Ydata, length(Ydata), Flength);
      end
      DataSpecs(i_dset).Ydata = (Ydata .* Yscale) + Yoffset;
    end
    
  end
end
