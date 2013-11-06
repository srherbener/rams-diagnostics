function [ ] = Plot2dMap( Fig, X, Y, Z, Clevs, Cbounds, Xlab, Ylab, Ptitle, SelectData )
%Plot2dMap Create a countour plot of the 2D data given in vector Z.
%   This function will extract the 2D data given in selected columns of
%   D, convert the column data to 2D data and creat contour plots of that
%   data as subplots in a single figure (specified by Fig).
%
%   Vector X contains the values for the x-axis, and vector Y contains the
%   values for the y-axis.
%
%   It is assumned that the data in Z is organized as:
%     first Ncols (length of X) items --> row = 1, columns 1 through Ncols
%     next Ncols items                --> row = 2, columns 1 through Ncols
%     etc.
%
%   SelectData contains [ x1 x2 y1 y2 ] which are used to select a
%   rectangular section out of the entire map.

Fsize = 38;

PanelTitle = false;
if (regexp(Ptitle, '^PANEL:'))
    Ptitle = regexprep(Ptitle, '^PANEL:', '');
    PanelTitle = true;
end

figure(Fig);

Nrows = length(Y);
Ncols = length(X);
Map = Create2dMap(Z,Nrows,Ncols);

% Trim out the rectangular region defined by SelectData
x1 = SelectData(1);
x2 = SelectData(2);
y1 = SelectData(3);
y2 = SelectData(4);
XP = X(x1:x2);
YP = Y(y1:y2);
MapP = Map(y1:y2,x1:x2);

%
Xticks = [ 70 140 210 ];
Xticklabels = { '70' '140' '210' };

% Find the largest absolute value of the entries in Dmap for setting the
% colormap axis. Want to center this about zero so that blue represents
% negative values and red represents positive values.
contourf(XP, YP, MapP, Clevs);
set(gca, 'FontSize', Fsize);
set(gca, 'LineWidth', 2);
set(gca, 'TickLength', [ 0.025 0.025 ]);
set(gca, 'XTick', Xticks);
set(gca, 'XTickLabel', Xticklabels);
shading flat;
caxis(Cbounds);
colormap('redblue');
cbar = colorbar;
set(cbar, 'FontSize', Fsize);
hold on;

% Draw the zero contour
contour(XP, YP, MapP, [ 0 0 ], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 5);

if (~strcmp(Ptitle, ' '))
  if (PanelTitle)
      % The title is in a box that adjusts to the amount of characters in
      % the title. Ie, it doesn't do any good to do Left/Center/Right
      % alignment. But, the entire box can be moved to the left side of the
      % plot.
      T = title(Ptitle);
      set(T, 'Units', 'Normalized');
      set(T, 'HorizontalAlignment', 'Left');
      Tpos = get(T, 'Position');
      Tpos(1) = 0; % line up with left edge of plot area
      set(T, 'Position', Tpos);
  else
      title(Ptitle);
  end
end
xlabel(Xlab);
ylabel(Ylab);


% Fix up the positioning
Ppos = get(gca, 'Position'); % position of plot area
Ppos(1) = Ppos(1) * 1.00;
Ppos(2) = Ppos(2) * 1.00;
Ppos(3) = Ppos(3) * 0.85;
Ppos(4) = Ppos(4) * 0.85;
set(gca, 'Position', Ppos);

end

