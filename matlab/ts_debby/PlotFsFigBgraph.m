function [] = PlotFsFigBgraph(Paxes, Y, Bcolors, Pmarker, Ptitle, Xlabel, Blabels, Ylabel, Ylim, Fsize, ShowX, ShowY, LegText, LegLoc)

  axes(Paxes);

  % Want to show accummulation of factors using a plot like a bar graph.
  % Each factor is shown as a bar with it's particular height (positive
  % or negative) going from the running sum of the previous factors. Like
  % a stacked bar graph but the pieces of the bars are fanned out beside
  % each other so that a negative piece doesn't obstruct a neighboring
  % positive piece.

  % The Y matrix is two rows, with the starting height in the first
  % row and the ending height in the second row. Ie, the number of
  % bars we want is equal to the number of columns in Y.
  % 
  % Proportion the x-axis so that the bar widths are 0.8 and
  % distances between bar centers is 1.
  Nbars = size(Y,2);
  Xlim = [ 0 Nbars+1 ];
  Xticks = [ 1:Nbars ];

  % Use patch to draw each bar at its specified height
  hold on;
  for i = 1:Nbars
    Bcolor = str2rgb(Bcolors{i});
    X1 = i - 0.4; 
    X2 = i + 0.4; 
    Y1 = Y(1,i);
    Y2 = Y(2,i);

    Xbar = [ X1 X1 X2 X2 ];
    Ybar = [ Y1 Y2 Y2 Y1 ];

    % Draw a bar for the factor magnitude
    patch(Xbar,Ybar,Bcolor);

%    % On just the factors F1, F2, F12, superimpose
%    % an arrow to denote the sign of the factor.
%    if ((i ~= 1) && (i ~= Nbars))
%      Xarrow1 = [  i  i ];
%      Yarrow1 = [ Y1 Y2 ];
%      Ainc = 0.08;
%      % lines for the arrowhead
%      if (Y1 <= Y2)
%        % positive value, place arrowhead on top
%        Xarrow2 = [  i-Ainc  i ];
%        Yarrow2 = [ Y2-Ainc Y2 ];
%
%        Xarrow3 = [  i+Ainc  i ];
%        Yarrow3 = [ Y2-Ainc Y2 ];
%      else
%        % negative value, place arrowhead on bottom
%        Xarrow2 = [  i-Ainc  i ];
%        Yarrow2 = [ Y2+Ainc Y2 ];
%
%       Xarrow3 = [  i+Ainc  i ];
%       Yarrow3 = [ Y2+Ainc Y2 ];
%     end
% 
%     line(Xarrow1, Yarrow1, 'Color', 'k');
%     line(Xarrow2, Yarrow2, 'Color', 'k');
%     line(Xarrow3, Yarrow3, 'Color', 'k');
%    end
  end

  set(Paxes, 'Xtick', Xticks);
  set(Paxes, 'XTickLabel', Blabels);
 
  set(Paxes, 'FontSize', Fsize);
  set(Paxes, 'LineWidth', 2);
  set(Paxes, 'TickLength', [ 0.025 0.025 ]);

  xlabel(Xlabel);
  xlim(Xlim);

  ylabel(Ylabel);
  ylim(Ylim);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

  if (~strcmp(LegText, 'none'))
    legend(LegText, 'Location', LegLoc);
  end
end
