function [ ] = GenLinePlot( Fig, AxisSpecs, DataSpecs, LegendSpecs )
%GenLinePlot Plot a set of lines the same panel
%   Fig - handle to a figure that was opened by the caller
%
%   AxisProps - structure containing a list of axis property names and
%               associated values that are desired to be set
%
%   DataSpecs - structure containing lists of X and Y values describing
%               lines that are to be drawn on the plot
%
%   LegendSpecs - struction containing instructions for making a legend 
%

  % select the passed in figure
  figure(Fig);

  % Create the axes and apply all the passed in properties
  Axes = axes;

  AxisProps = AxisSpecs.Props;
  Nprops = length(AxisProps);
  for i = 1:Nprops
    set(Axes, AxisProps(i).Name, AxisProps(i).Val);
  end

  title(AxisSpecs.Title);
  xlabel(AxisSpecs.Xlabel);
  ylabel(AxisSpecs.Ylabel);

  % Draw the lines
  Nlines = length(DataSpecs);
  for i = 1:Nlines
    X = DataSpecs(i).Xdata;
    Y = DataSpecs(i).Ydata;
    Lcolor = str2rgb(DataSpecs(i).Lcolor);
    Lwidth = DataSpecs(i).Lwidth;
    Lstyle = DataSpecs(i).Lstyle;

    line(X, Y, 'LineStyle', Lstyle, 'LineWidth', Lwidth, 'Color', Lcolor);
  end

  % Place the legend
  if (~strcmp(LegendSpecs.Loc, 'none'))
    legend(LegendSpecs.Text, 'Location', LegendSpecs.Loc, 'FontSize', LegendSpecs.Fsize);
  end

end
