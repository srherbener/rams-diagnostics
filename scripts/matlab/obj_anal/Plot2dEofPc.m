function [ ] = Plot2dEofPc(EOF, PC, X, Y, Clevs, Cbounds, E_Title, E_Xlabel, E_Ylabel, P_Title, P_Xlabel, P_Ylabel, SelectData, EofFile, PcFile)
%Plot2dEofPc plot out a specified EOF and PC of a set of 2D observations
%   This function will plot out the Nth EOF and PC of set of 2D
%   observations. The EOFs and PCs are in the columns of their respective
%   input arrays. X and Y are the coordinate values for the 2D field.
%
%   SelectData contains [ x1 x2 y1 y2 ] which are used to select a rectangular
%   region out of the entire 2D map.

Fsize = 35;

PanelTitle = false;
if (regexp(P_Title, '^PANEL:'))
    P_Title = regexprep(P_Title, '^PANEL:', '');
    PanelTitle = true;
end

Fig = figure;
Plot2dMap(Fig,X,Y,EOF, Clevs, Cbounds, E_Xlabel, E_Ylabel, E_Title, SelectData);
saveas(Fig, EofFile);
close(Fig);

Fig = figure;
plot(PC, 'Color', 'k', 'LineWidth', 4);
set(gca, 'FontSize', Fsize);
set(gca, 'LineWidth', 2);
set(gca, 'TickLength', [ 0.025 0.025 ]);

if (~strcmp(P_Title, ' '))
  if (PanelTitle)
      % The title is in a box that adjusts to the amount of characters in
      % the title. Ie, it doesn't do any good to do Left/Center/Right
      % alignment. But, the entire box can be moved to the left side of the
      % plot.
      T = title(P_Title);
      set(T, 'Units', 'Normalized');
      set(T, 'HorizontalAlignment', 'Left');
      Tpos = get(T, 'Position');
      Tpos(1) = 0; % line up with left edge of plot area
      set(T, 'Position', Tpos);
  else
      title(P_Title);
  end
end
xlabel(P_Xlabel);
ylabel(P_Ylabel);

xlim([ 0 150 ]);

% Fix up the positioning
Ppos = get(gca, 'Position'); % position of plot area
Ppos(1) = Ppos(1) * 1.05;
Ppos(2) = Ppos(2) * 1.05;
Ppos(3) = Ppos(3) * 0.90;
Ppos(4) = Ppos(4) * 0.90;
set(gca, 'Position', Ppos);

saveas(Fig, PcFile);
close(Fig);

end
