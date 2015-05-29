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

  % Draw the plot
  Ndsets = length(DataSpecs);
  hold on;
  for i = 1:Ndsets
    if (strcmp(Ptype, 'line'))
      X = DataSpecs(i).Xdata;
      Y = DataSpecs(i).Ydata;
      Lcolor = str2rgb(DataSpecs(i).Lcolor);
      Lwidth = DataSpecs(i).Lwidth;
      Lstyle = DataSpecs(i).Lstyle;

      line(X, Y, 'LineStyle', Lstyle, 'LineWidth', Lwidth, 'Color', Lcolor);
    elseif (strncmp(Ptype, 'bar', 3))
      Bcolor = str2rgb(DataSpecs(i).Lcolor);

      Bspecs = strsplit(Ptype, ':');
      Btype = 'grouped';
      if (length(Bspecs) > 1)
        Btype = Bspecs{2};
      end

      X = DataSpecs(i).Xdata;
      Y = DataSpecs(i).Ydata;
      bar(X,Y,Btype, 'FaceColor', Bcolor);
    elseif (strcmp(Ptype, 'contourf'))
      X = DataSpecs(i).Xdata;
      Y = DataSpecs(i).Ydata;
      Z = DataSpecs(i).Zdata;

      % check dimension sizes, and transpose Z if necessary
      if ((length(X) == size(Z,1)) && (length(Y) == size(Z,2)))
        Z = Z';
      end

      contourf(X, Y, Z, 'LineStyle', 'none');
      colorbar;
    end
  end
  hold off;
  
end
