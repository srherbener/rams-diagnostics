function [ ] = MakeLinePlot( Fig, AxisSpecs, DataSpecs )
%MakeLinePlot Plot a set of lines the same panel
%   Fig - handle to a figure that was opened by the caller
%
%   AxisProps - structure containing a list of axis property names and
%               associated values that are desired to be set
%
%   DataSpecs - structure containgin lists of X and Y values describing
%               lines that are to be drawn on the plot
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
  Nlines = length(DataSpecs.X);
  for i = 1:Nlines
    X = DataSpecs.X{i};
    Y = DataSpecs.Y{i};
    Lcolor = str2rgb(DataSpecs.Lcolor{i});
    Lwidth = DataSpecs.Lwidth{i};
    Lstyle = DataSpecs.Lstyle{i};

    line(X, Y, 'LineStyle', Lstyle, 'LineWidth', Lwidth, 'Color', Lcolor);
  end

end
