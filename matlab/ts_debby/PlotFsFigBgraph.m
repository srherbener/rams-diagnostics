function [] = PlotFsFigBgraph(Paxes, Y, Bcolors, Pmarker, Ptitle, Xlabel, Blabels, Ylabel, Ylim, Fsize, ShowX, ShowY, LegText, LegLoc)

  axes(Paxes);

  Ysize = size(Y);
  Nbars = Ysize(1);
  Nsets = Ysize(2);

  % Create a dummy X vector
  X = (1:Nbars);

  bgraph = bar(X, Y, 'grouped');
  set(Paxes, 'FontSize', Fsize);
  for i = 1:Nsets
    Bcolor = str2rgb(Bcolors{i});
    bgraph(i).FaceColor = Bcolor;
  end

  set(Paxes, 'XTickLabel', Blabels);
  xlabel(Xlabel);

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
