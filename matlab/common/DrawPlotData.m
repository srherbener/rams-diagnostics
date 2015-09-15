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

  if (strcmp(Ptype, 'line'))
    % line plot
    %
    % Loop through each dataset (line) using the line
    % command to draw the curve. This will accommodate
    % different sized datasets as well as datasets containing
    % differing x values.
    %
    for i = 1:Ndsets
      X = DataSpecs(i).Xdata;
      Y = DataSpecs(i).Ydata;
      Lcolor = str2rgb(DataSpecs(i).Lcolor);
      Lwidth = DataSpecs(i).Lwidth;
      Lstyle = DataSpecs(i).Lstyle;

      line(X, Y, 'LineStyle', Lstyle, 'LineWidth', Lwidth, 'Color', Lcolor);
    end
  elseif (strncmp(Ptype, 'bar', 3))
    % bar graphs
    %
    % The bar command will only work with multiple datasets by
    % placing the data in arrays and calling the bar command once.
    % Note that running:
    %
    %   hold on
    %   bar()
    %   bar()
    %   ...
    %   hold off
    %
    % will only superimpose the bar graphs on top of each other.

    % collect the dataset values into arrays
    % assume that dataset lengths are all the same
    Nvals = length(DataSpecs(1).Xdata);
    X = zeros([ Nvals Ndsets ]);
    Y = zeros([ Nvals Ndsets ]);
    for i = 1:Ndsets
      X(:,i) = DataSpecs(i).Xdata;
      Y(:,i) = DataSpecs(i).Ydata;
    end

    % Determine the type of bar graph (grouped, stacked, etc)
    Bspecs = strsplit(Ptype, ':');
    Btype = 'grouped';
    if (length(Bspecs) > 1)
      Btype = Bspecs{2};
    end

    % create the graph
    bgraph = bar(X,Y,Btype);

    % Set the color for each set of bars. This must be done
    % by setting structure elements referenced by the bar
    % graph handle (bgraph).
    for i = 1:Ndsets
      Bcolor = str2rgb(DataSpecs(i).Lcolor);
      bgraph(i).FaceColor = Bcolor;
    end
  elseif (regexp(Ptype, 'contour'))
    % contour plot
    Nclevs = 20;

    % assume only one dataset
    X = DataSpecs(1).Xdata;
    Y = DataSpecs(1).Ydata;
    Z = DataSpecs(1).Zdata;

    % check dimension sizes, and transpose Z if necessary
    if ((length(X) == size(Z,1)) && (length(Y) == size(Z,2)))
      Z = Z';
    end

    if (regexp(Ptype, 'contourf'))
      contourf(X, Y, Z, Nclevs, 'LineStyle', 'none');
    else
      contour(X, Y, Z, Nclevs, 'LineStyle', 'none');
    end

    % used red/blue for difference plots
    if (strcmp(Ptype, 'diff_contourf'))
      colormap(Axes, redblue);
    else
      colormap(Axes, 'default');
    end
    cb = colorbar;
  end

end
