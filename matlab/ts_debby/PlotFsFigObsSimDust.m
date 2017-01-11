function [ ] = PlotFsFigObsSimDust()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

  % TS Debby AOD images. PNG files have true color, so Cmap vars will be empty (which
  % is exactly what you want for geoshow.
  % AOD images are composites over one day
  [ AodStart AodStartCmap ] = imread('IMAGES/MYBGAOD_C10K.A2006234.006.0.1.EastAtlantic.png');
  [ AodMid   AodMidCmap   ] = imread('IMAGES/MYBGAOD_C10K.A2006235.006.0.1.EastAtlantic.png');
  [ AodEnd   AodEndCmap   ] = imread('IMAGES/MYBGAOD_C10K.A2006236.006.0.1.EastAtlantic.png');

  % Vertically integrated dust
  VdustFname = 'HDF5/TSD_SAL_DUST/HDF5/vint_dust-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  VdustVname = '/vertint_dust';

  % Read in and create a snapshot of the vertically integrated dust field at Aug23, 12Z
  fprintf('Reading: %s (%s)\n', VdustFname, VdustVname);
  X = squeeze(h5read(VdustFname, '/x_coords'));  % lon
  Y = squeeze(h5read(VdustFname, '/y_coords'));  % lat
  T = squeeze(h5read(VdustFname, '/t_coords'))./3600 - 42; % sim hours, start with zero

  Nx = length(X);
  Ny = length(Y);
  Nt = length(T);

  % DUST is (x,y,t)
  %   Start is time =  0 h
  %   Mid   is time = 30 h
  %   End   is time = 60 h
  Ccount = [ Nx Ny  1 ];
  T1 = find(T >= 0, 1, 'first');
  Cstart = [  1  1 T1 ];
  DUST_START = squeeze(h5read(VdustFname, VdustVname, Cstart, Ccount));

  T1 = find(T >= 30, 1, 'first');
  Cstart = [  1  1 T1 ];
  DUST_MID   = squeeze(h5read(VdustFname, VdustVname, Cstart, Ccount));

  T1 = find(T >= 60, 1, 'first');
  Cstart = [  1  1 T1 ];
  DUST_END   = squeeze(h5read(VdustFname, VdustVname, Cstart, Ccount));

  % Cut down the size of the data for plotting
  % Also, contourfm wants x,y,z values to be double
  Inc = 10;
  X = double(X(1:Inc:end));
  Y = double(Y(1:Inc:end));
  DUST_START = double(DUST_START(1:Inc:end, 1:Inc:end));
  DUST_MID   = double(DUST_MID(1:Inc:end, 1:Inc:end));
  DUST_END   = double(DUST_END(1:Inc:end, 1:Inc:end));

  % plot
  Fig = figure;

  WmapLatBounds = [ 5 26 ];
  WmapLonBounds = [ -42 -8 ];

  % AOD satellite images in left column
  Paxes = subplot(3,2,1);
  PlaceAodImage(Paxes, AodStart, AodStartCmap, 'a', '22Aug', Fsize)

  Paxes = subplot(3,2,3);
  PlaceAodImage(Paxes, AodMid, AodMidCmap, 'c', '23Aug', Fsize)

  Paxes = subplot(3,2,5);
  PlaceAodImage(Paxes, AodEnd, AodEndCmap, 'e', '24Aug', Fsize)

  % Vertically integrated condensate in right column
  Paxes = subplot(3,2,2);
  PlotVintDust(Paxes, X, Y, DUST_START, WmapLatBounds, WmapLonBounds, 'b', '06Z, 22Aug', Fsize)

  Paxes = subplot(3,2,4);
  PlotVintDust(Paxes, X, Y, DUST_MID, WmapLatBounds, WmapLonBounds, 'd', '12Z, 23Aug', Fsize)

  Paxes = subplot(3,2,6);
  PlotVintDust(Paxes, X, Y, DUST_END, WmapLatBounds, WmapLonBounds, 'f', '18Z, 24Aug', Fsize)

  OutFile = sprintf('%s/FsFigObsSimDust.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = PlaceAodImage(Paxes, Image, Cmap, Pmarker, Ptitle, Fsize)

  axes(Paxes);

  set(Paxes, 'FontSize', Fsize);

  imshow(Image, Cmap);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(Paxes, T);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlotVintDust(Paxes, X, Y, Z, LatBounds, LonBounds, Pmarker, Ptitle, Fsize)

  axes(Paxes);

  LineW = 2;

  set(Paxes, 'FontSize', Fsize);
  worldmap(LatBounds, LonBounds);
  setm(Paxes, 'MapProjection', 'miller');

  % X is LON, Y is LAT
  % Want to show log values of Z, but there is not a built-in
  % method for this.
  %   1. send contourfm the log values
  %      this means you have to work in log values for the caxis
  %      and tick marks, etc.
  %      e.g., caxis [ -2 6 ] translates to 1e-2 to 1e6
  %   2. contourfm/contourcmap create a colorbar with
  %      the selected colors for the contour levels displayed
  %      like tiles and a scaled tick mark value that places
  %      integer values at the center of each color. This makes
  %      the edges of the tiles equal to tick values like 0.5 1.5 2.5 ...
  %      where the centers of the tiles are equal to tick values
  %      1 2 3 ... n where n is the number of contour levels.
  % 
  % So, in this example the range for Z is 1e-2 to 1e6. This
  % calls for:
  %   caxis: -2 to 6
  %   contour levels: -2:0.5:6 (sixteen levels)
  %   color bar axes: center of tiles have values 1 2 3 ... 16
  %                   edges have values 0.5 1.5 2.5 ... 16.5
  %
  %   This means that there is a linear scale between the colorbar
  %   tick mark values and the actual values in the plot:
  %
  %      tick value        plot value    origninal value
  %         0.5               -2            0.01
  %         1.5               -1.5
  %         2.5               -1            0.1
  %         3.5               -0.5
  %         4.5                0            1
  %         ...                ...
  %        15.5                5.5
  %        16.6                6            1000000
  %                  
  %
  
  Clevs = [ -2:0.5:6 ];
  contourfm(double(Y), double(X), double(log10(Z))', Clevs, 'LineStyle', 'none'); 
  caxis(Paxes, [ Clevs(1) Clevs(end) ]);

  % "parula" is the default colormap
  Cbar = contourcmap('parula', 'ColorBar', 'on', 'Location', 'vertical' );
  Cticks = [ 2.5 6.5 10.5 14.5 ];
  CtickLabels = { '10^-^1' '10^1' '10^3' '10^5' };
  set(Cbar, 'Xtick', Cticks);
  set(Cbar, 'XTickLabel', CtickLabels);

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
    LeftJustTitle(Paxes, T);
  end
end

