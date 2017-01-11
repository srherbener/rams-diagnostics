function [] = PlotFsFigLine(Paxes, X, Y, Pmarker, Ptitle, Xlabel, Xlim, Xscale, Xshow, Ylabel, Ylim, Yscale, Yshow, Fsize, LegText, LegLoc, Colors)

  % X and Y need to contain line data in each column. Therefore, the length of the columns
  % (ie, the number of rows) of X and Y need to match. X and Y can either be vectors or
  % matrices. If both are matrices, then the number of columns (ie, number of lines) need
  % to match.

  axes(Paxes);

  LegendFsize = 8;
  LineW = 2;

  [ Xnpts Xnlines ] = size(X);
  [ Ynpts Ynlines ] = size(Y);

  if (Xnpts ~= Ynpts)
    fprintf('ERROR: PlotFsFigLine: number of rows of X and Y need to match\n');
    fprintf('ERROR: PlotFsFigLine: skipping this plot\n');
    fprintf('\n');
    return
  end

  if ((Xnlines == 1) && (Ynlines > 1))
    Nlines = Ynlines;
    Xdata = repmat(X, [ 1 Ynlines ]);
    Ydata = Y;
  elseif ((Xnlines >1) && (Ynlines == 1))
    Nlines = Xnlines;
    Xdata = X;
    Ydata = repmat(Y, [ 1 Xnlines ]);
  elseif (Xnlines ~= Ynlines)
    % If we got to this test, Xnlines and Ynlines are either both 1 or both > 1.
    % If they don't match, then we have a problem and will skip this plot.
    fprintf('ERROR: PlotFsFigLine: number of columns in X and Y need to match when they are both > 1\n');
    fprintf('ERROR: PlotFsFigLine: skipping this plot\n');
    fprintf('\n');
    return
  else
    % If we got to here, Xnlines and Ynlines match.
    Nlines = Xnlines;
    Xdata = X;
    Ydata = Y;
  end

  set(Paxes, 'FontSize', Fsize);
  for i = 1:Nlines
    Xline = squeeze(Xdata(:,i));
    Yline = squeeze(Ydata(:,i));
    Color = str2rgb(Colors{i});
    Plines(i) = line(Xline, Yline, 'Color', Color, 'LineWidth', LineW);
  end

  set(Paxes, 'Xscale', Xscale);
  set(Paxes, 'Yscale', Yscale);

  xlim(Xlim);
  ylim(Ylim);

  if (Xshow > 0)
    xlabel(Xlabel);
  else
    set(Paxes, 'XtickLabel', {});
  end
  if (Yshow > 0)
    ylabel(Ylabel);
  else
    set(Paxes, 'YtickLabel', {});
  end

  if (~strcmp(LegLoc, 'none'))
    legend(Plines, LegText,'Location', LegLoc, 'FontSize', LegendFsize);
  end

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(Paxes, T);
  end
end
