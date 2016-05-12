function [ ] = PlotFsFigObsSimStorm()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

  % TS IR images, gif files need to read in the image as well as the color map, then
  % use both to display.
  [ IrStart IrStartCmap ] = imread('IMAGES/2006AL04_4KMIRIMG_200608220600.GIF');
  [ IrMid   IrMidCmap   ] = imread('IMAGES/2006AL04_4KMIRIMG_200608231200.GIF');
  [ IrEnd   IrEndCmap   ] = imread('IMAGES/2006AL04_4KMIRIMG_200608241800.GIF');

  % Record the Lat,Lon extent of each image so that the image gets placed correctly
  % on the worldmap grid.
  % Format: [ Lat-of-south-edge Lat-of-north-edge Lon-of-west-edge Lon-of-east-edge ]
  IrStartExt = [  4.0 20.2 -34.1 -11.0 ];
  IrMidExt   = [  6.4 23.7 -40.9 -17.9 ];
  IrEndExt   = [ 12.5 27.5 -48.5 -25.5 ];

  % Vertically integrated condensate
  VcondFname = 'HDF5/TSD_SAL_DUST/HDF5/vint_cond-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  VcondVname = '/vertint_cond_nd';

  % Read in and create a snapshot of the vertically integrated dust field at Aug23, 12Z
  fprintf('Reading: %s (%s)\n', VcondFname, VcondVname);
  X = squeeze(h5read(VcondFname, '/x_coords'));  % lon
  Y = squeeze(h5read(VcondFname, '/y_coords'));  % lat
  T = squeeze(h5read(VcondFname, '/t_coords'))./3600 - 42; % sim hours, start with zero

  Nx = length(X);
  Ny = length(Y);
  Nt = length(T);

  % COND is (x,y,t)
  %   Start is time =  0 h
  %   Mid   is time = 30 h
  %   End   is time = 60 h
  Ccount = [ Nx Ny  1 ];
  T1 = find(T >= 0, 1, 'first');
  Cstart = [  1  1 T1 ];
  COND_START = squeeze(h5read(VcondFname, VcondVname, Cstart, Ccount));

  T1 = find(T >= 30, 1, 'first');
  Cstart = [  1  1 T1 ];
  COND_MID   = squeeze(h5read(VcondFname, VcondVname, Cstart, Ccount));

  T1 = find(T >= 60, 1, 'first');
  Cstart = [  1  1 T1 ];
  COND_END   = squeeze(h5read(VcondFname, VcondVname, Cstart, Ccount));

  % Cut down the size of the data for plotting
  % Also, contourfm wants x,y,z values to be double
  Inc = 10;
  X = double(X(1:Inc:end));
  Y = double(Y(1:Inc:end));
  COND_START = double(COND_START(1:Inc:end, 1:Inc:end));
  COND_MID   = double(COND_MID(1:Inc:end, 1:Inc:end));
  COND_END   = double(COND_END(1:Inc:end, 1:Inc:end));

  % plot
  Fig = figure;

  % IR images in left column
  Paxes = subplot(3,2,1);
  PlaceIrImage(Paxes, IrStart, IrStartCmap, 'a', '06Z, 22Aug (SD)', Fsize, IrStartExt)

  Paxes = subplot(3,2,3);
  PlaceIrImage(Paxes, IrMid, IrMidCmap, 'c', '12Z, 23Aug', Fsize, IrMidExt)

  Paxes = subplot(3,2,5);
  PlaceIrImage(Paxes, IrEnd, IrEndCmap, 'e', '18Z, 24Aug', Fsize, IrEndExt)

  % Vertically integrated condensate in right column
  Paxes = subplot(3,2,2);
  PlotVintCond(Paxes, X, Y, COND_START, 'b', '06Z, 22Aug', Fsize)

  Paxes = subplot(3,2,4);
  PlotVintCond(Paxes, X, Y, COND_MID, 'd', '12Z, 23Aug', Fsize)

  Paxes = subplot(3,2,6);
  PlotVintCond(Paxes, X, Y, COND_END, 'f', '18Z, 24Aug', Fsize)

  OutFile = sprintf('%s/FsFigObsSimStorm.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlaceIrImage(Paxes, IrImage, IrCmap, Pmarker, Ptitle, Fsize, IrExtent)

  axes(Paxes);

  LatBounds = [ 5 26 ];
  LonBounds = [ -42 -8 ];
  LineW = 2;

  set(Paxes, 'FontSize', Fsize);
  worldmap(LatBounds, LonBounds);
  setm(Paxes, 'MapProjection', 'miller');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlotVintCond(Paxes, X, Y, Z, Pmarker, Ptitle, Fsize)

  axes(Paxes);

  LatBounds = [ 5 26 ];
  LonBounds = [ -42 -8 ];

  LineW = 2;

  set(Paxes, 'FontSize', Fsize);
  worldmap(LatBounds, LonBounds);
  setm(Paxes, 'MapProjection', 'miller');

  % Color levels and limits
  Cmin = 0.4;
  Cmax = 4.0;
  Cinc = 0.2;
  Clim = [ Cmin Cmax ];
  Clevs = Cmin:Cinc:Cmax;

  % replace small values with nans so that they don't appear
  % as dark regions in the contourfm plots
  Z(Z < Cmin) = nan;

  contourfm(Y, X, Z', Clevs, 'LineStyle', 'none'); 
  caxis(Paxes, Clim);

  % Put in a colorbar. The colorbar labeling has too many labels and comes out
  % cluttered. Just read in the Ytick and YtickLabel values that contourcmap
  % created, select every nth tick and label, and then reset the colorbar.
  Cbar = contourcmap('jet', 'ColorBar', 'on', 'Location', 'vertical');

  Cticks = get(Cbar, 'Ytick');
  CtickLabels = cellstr(get(Cbar, 'YtickLabel'))';

  CtickInc = 2;
  Cticks = Cticks(1:CtickInc:end);
  CtickLabels = { CtickLabels{1:CtickInc:end} };

  set(Cbar, 'Ytick', Cticks);
  set(Cbar, 'YtickLabel', CtickLabels);

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

