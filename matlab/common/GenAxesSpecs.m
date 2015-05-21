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

  AxesSpecs.Clims = [ ]; % need empty list for check in SetPlotAxes()
  i_ap = 0;

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

  % Assume that you have Naxes equal to 1, 2 or 3
  %  Axis number 1 --> x axis specs
  %  Axis number 2 --> y axis specs
  %  Axis number 3 --> z axis specs
  %
  Naxes = Config.PlotAxes(i_paxis).Naxes;

  if (Naxes >= 1)
    % have a spec for the x axis 

    Xlabel      = Config.PlotAxes(i_paxis).Axes(1).Label;
    Xunits      = Config.PlotAxes(i_paxis).Axes(1).Units;
    Xmin        = Config.PlotAxes(i_paxis).Axes(1).Min;
    Xmax        = Config.PlotAxes(i_paxis).Axes(1).Max;
    Xscale      = Config.PlotAxes(i_paxis).Axes(1).Scale;
    Xticks      = Config.PlotAxes(i_paxis).Axes(1).Ticks;
    XtickLabels = Config.PlotAxes(i_paxis).Axes(1).TickLabels;

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
  end

  if (Naxes >= 2)
    % have a spec for the y axis 

    Ylabel      = Config.PlotAxes(i_paxis).Axes(2).Label;
    Yunits      = Config.PlotAxes(i_paxis).Axes(2).Units;
    Ymin        = Config.PlotAxes(i_paxis).Axes(2).Min;
    Ymax        = Config.PlotAxes(i_paxis).Axes(2).Max;
    Yscale      = Config.PlotAxes(i_paxis).Axes(2).Scale;
    Yticks      = Config.PlotAxes(i_paxis).Axes(2).Ticks;
    YtickLabels = Config.PlotAxes(i_paxis).Axes(2).TickLabels;

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

  if (Naxes >= 3)
    % have a spec for the z axis 

    Zlabel      = Config.PlotAxes(i_paxis).Axes(3).Label;
    Zunits      = Config.PlotAxes(i_paxis).Axes(3).Units;
    Zmin        = Config.PlotAxes(i_paxis).Axes(3).Min;
    Zmax        = Config.PlotAxes(i_paxis).Axes(3).Max;
    Zscale      = Config.PlotAxes(i_paxis).Axes(3).Scale;
    Zticks      = Config.PlotAxes(i_paxis).Axes(3).Ticks;
    ZtickLabels = Config.PlotAxes(i_paxis).Axes(3).TickLabels;

    % direction and limits of z axis
    if (Zmin > Zmax)
      AxesSpecs.Clims = [ Zmax Zmin ]; 
  
      i_ap = i_ap + 1;
      AxesSpecs.Props(i_ap).Name = 'ZDir';
      AxesSpecs.Props(i_ap).Val = 'Reverse';
    else
      AxesSpecs.Clims = [ Zmin Zmax ]; 
    end
  
    % z axis scale
    i_ap = i_ap + 1;
    AxesSpecs.Props(i_ap).Name = 'Zscale';
    AxesSpecs.Props(i_ap).Val  = Zscale;
  end
end
