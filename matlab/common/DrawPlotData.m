function [ ] = DrawPlotData( Axes, DataSpecs, Ptype )
%DrawPlotData draw a set of lines/bars/etc on a given set of axes
%   Axes - handle to a set of axes created by the caller
%
%   DataSpecs - structure containing lists of X and Y values describing
%               lines that are to be drawn on the plot
%
%   Ptype - Plot type, 'line' or 'bar'
%

  % select the passed in axes
  axes(Axes);

  % Draw the lines
  Ndsets = length(DataSpecs);
  for i = 1:Ndsets
    if (strcmp(Ptype, 'line'))
      X = DataSpecs(i).Xdata;
      Y = DataSpecs(i).Ydata;
      Lcolor = str2rgb(DataSpecs(i).Lcolor);
      Lwidth = DataSpecs(i).Lwidth;
      Lstyle = DataSpecs(i).Lstyle;

      line(X, Y, 'LineStyle', Lstyle, 'LineWidth', Lwidth, 'Color', Lcolor);
    elseif (strcmp(Ptype, 'bar'))
      X = DataSpecs(i).Xdata;
      Y = DataSpecs(i).Ydata;
      bar(X,Y);
    end
  end
  
end
