function [] = PlotFsFigImageWmap(Paxes, IrImage, IrCmap, LatBounds, LonBounds, Pmarker, Ptitle, Fsize, IrExtent)

  axes(Paxes);

  LineW = 2;

  set(Paxes, 'FontSize', Fsize);
  worldmap(LatBounds, LonBounds);
  setm(Paxes, 'MapProjection', 'miller', 'FontSize', 8);

  % Create gridded versions of the latitude and longitude value
  % that the image covers.
  [ Nlat Nlon Ncolors ] = size(IrImage);
  Lat1 = IrExtent(1);
  Lat2 = IrExtent(2);
  Lon1 = IrExtent(3);
  Lon2 = IrExtent(4);

  LatInc = (Lat2-Lat1)/(Nlat-1);
  LonInc = (Lon2-Lon1)/(Nlon-1);

  LAT = fliplr(Lat1:LatInc:Lat2);
  LON = Lon1:LonInc:Lon2;
  [ LonGrid LatGrid ] = meshgrid(LON, LAT);

  geoshow(LatGrid, LonGrid, IrImage, IrCmap);

  % Draw the coast after the contourfm call so that coast lines will appear on top
  % of the contour plot
  Coast = load('coast');
  plotm(Coast.lat, Coast.long, 'Color', 'k', 'LineWidth', LineW);

  % Mark Africa
  textm(14, -12, 'Africa', 'FontSize', 10, 'Color', 'k', 'Rotation', 90);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

end
