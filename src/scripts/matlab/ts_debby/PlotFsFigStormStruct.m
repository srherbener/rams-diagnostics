function [ ] = PlotFsFigStormStruct(Paxes, X, Y, VT, UP, Pmarker, Ptitle, Xlabel, Xlim, Ylabel, Ylim, VtCmap, VtClim, VtClevs, UpClevs, Fsize, ShowX, ShowY)

  % Plot Vt in shaded contours
  PlotFsFigXsection(Paxes, X, Y, VT, Pmarker, Ptitle, Xlabel, Xlim, Ylabel, Ylim, VtCmap, VtClim, VtClevs, Fsize, ShowX, ShowY);
  
  % place updraft speeds in line contours over the shaded contours
  axes(Paxes)

  hold on;
  [ Cup Hup ] = contour(X, Y, UP, UpClevs, 'LineColor', 'w');
  clabel(Cup, 'Color', 'w');
  hold off;

end
