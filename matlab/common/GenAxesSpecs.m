function [ AxesSpecs ] = GenAxesSpecs(Config, i_panel, i_paxis)
% GenAxesSpecs function to generate axes specfifications for plots

  Ptitle  = Config.FigPanels(i_panel).Title.Main;
  if (isempty(Config.FigPanels(i_panel).Title.Pmarkers))
    Pmarker = '';
  else
    Pmarker = Config.FigPanels(i_panel).Title.Pmarkers{1};
  end

  XAshow = Config.FigPanels(i_panel).XAshow;
  YAshow = Config.FigPanels(i_panel).YAshow;

  i_ap = 0;

  % Fill in the AxesSpecs structure
  Xlabel      = Config.PlotAxes(i_paxis).Xlabel;
  Xunits      = Config.PlotAxes(i_paxis).Xunits;
  Xmin        = Config.PlotAxes(i_paxis).Xmin;
  Xmax        = Config.PlotAxes(i_paxis).Xmax;
  Xscale      = Config.PlotAxes(i_paxis).Xscale;
  Xticks      = Config.PlotAxes(i_paxis).Xticks;
  XtickLabels = Config.PlotAxes(i_paxis).XtickLabels;

  Ylabel      = Config.PlotAxes(i_paxis).Ylabel;
  Yunits      = Config.PlotAxes(i_paxis).Yunits;
  Ymin        = Config.PlotAxes(i_paxis).Ymin;
  Ymax        = Config.PlotAxes(i_paxis).Ymax;
  Yscale      = Config.PlotAxes(i_paxis).Yscale;
  Yticks      = Config.PlotAxes(i_paxis).Yticks;
  YtickLabels = Config.PlotAxes(i_paxis).YtickLabels;

  % axes font size
  i_ap = i_ap + 1;
  AxesSpecs.Props(i_ap).Name = 'FontSize';
  AxesSpecs.Props(i_ap).Val  = Config.PlotAxes(i_paxis).Fsize;

  % Line and tick weight (line widht, tick mark length)
  i_ap = i_ap + 1;
  AxesSpecs.Props(i_ap).Name = 'LineWidth';
  AxesSpecs.Props(i_ap).Val  = Config.PlotAxes(i_paxis).Lwidth;
  i_ap = i_ap + 1;
  AxesSpecs.Props(i_ap).Name = 'TickLength';
  AxesSpecs.Props(i_ap).Val  = Config.PlotAxes(i_paxis).Tlength;

  % Title
  AxesSpecs.Panel = ~isempty(Pmarker);
  if (AxesSpecs.Panel)
    AxesSpecs.Title = sprintf('(%s) %s', Pmarker, Ptitle);
  else
    AxesSpecs.Title = Ptitle;
  end

  %%%%%%%% X Axis %%%%%%%%%

  % x axis label, tick marks
  if (XAshow > 0)
    % show x-axis
    if (strcmp(Xunits, ' '))
      AxesSpecs.Xlabel   = sprintf('%s', Xlabel);
    else
      AxesSpecs.Xlabel   = sprintf('%s (%s)', Xlabel, Xunits);
    end

    % axis tick marks
    if (~isempty(Xticks))
      i_ap = i_ap + 1;
      AxesSpecs.Props(i_ap).Name = 'XTick';
      AxesSpecs.Props(i_ap).Val  = Xticks;
    end

    % axis tick labels
    if (~isempty(XtickLabels))
      i_ap = i_ap + 1;
      AxesSpecs.Props(i_ap).Name = 'XTickLabel';
      AxesSpecs.Props(i_ap).Val  = XtickLabels;
    end
  else
    % "hide" x-axis: no label, erase tick labels, keep tick marks
    AxesSpecs.Xlabel = '';
    i_ap = i_ap + 1;
    AxesSpecs.Props(i_ap).Name = 'XTickLabel';
    AxesSpecs.Props(i_ap).Val  = {};
  end

  % direction and limits of x axis
  if (Xmin > Xmax)
    i_ap = i_ap + 1;
    AxesSpecs.Props(i_ap).Name = 'Xlim';
    AxesSpecs.Props(i_ap).Val = [ Xmax Xmin ]; 

    i_ap = i_ap + 1;
    AxesSpecs.Props(i_ap).Name = 'XDir';
    AxesSpecs.Props(i_ap).Val = 'Reverse';
  else
    i_ap = i_ap + 1;
    AxesSpecs.Props(i_ap).Name = 'Xlim';
    AxesSpecs.Props(i_ap).Val = [ Xmin Xmax ]; 
  end

  % x axis scale
  i_ap = i_ap + 1;
  AxesSpecs.Props(i_ap).Name = 'Xscale';
  AxesSpecs.Props(i_ap).Val  = Xscale;

  %%%%%%%% Y Axis %%%%%%%%%

  % y axis label, tick marks
  if (YAshow > 0)
    % show y-axis
    if (strcmp(Yunits, ' '))
      AxesSpecs.Ylabel   = sprintf('%s', Ylabel);
    else
      AxesSpecs.Ylabel   = sprintf('%s (%s)', Ylabel, Yunits);
    end

    % y axis tick marks
    if (~isempty(Yticks))
      i_ap = i_ap + 1;
      AxesSpecs.Props(i_ap).Name = 'YTick';
      AxesSpecs.Props(i_ap).Val  = Yticks;
    end

    % axis tick labels
    if (~isempty(YtickLabels))
      i_ap = i_ap + 1;
      AxesSpecs.Props(i_ap).Name = 'YTickLabel';
      AxesSpecs.Props(i_ap).Val  = YtickLabels;
    end
  else
    % "hide" y-axis: no label, erase tick labels, keep tick marks
    i_ap = i_ap + 1;
    AxesSpecs.Ylabel = '';
    AxesSpecs.Props(i_ap).Name = 'YTickLabel';
    AxesSpecs.Props(i_ap).Val  = {};
  end

  % direction and limits of y axis
  if (Ymin > Ymax)
    i_ap = i_ap + 1;
    AxesSpecs.Props(i_ap).Name = 'Ylim';
    AxesSpecs.Props(i_ap).Val = [ Ymax Ymin ]; 

    i_ap = i_ap + 1;
    AxesSpecs.Props(i_ap).Name = 'YDir';
    AxesSpecs.Props(i_ap).Val = 'Reverse';
  else
    i_ap = i_ap + 1;
    AxesSpecs.Props(i_ap).Name = 'Ylim';
    AxesSpecs.Props(i_ap).Val = [ Ymin Ymax ]; 
  end

  % y axis scale
  i_ap = i_ap + 1;
  AxesSpecs.Props(i_ap).Name = 'Yscale';
  AxesSpecs.Props(i_ap).Val  = Yscale;

end
