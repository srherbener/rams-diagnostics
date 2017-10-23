function [ ] = PlotFsFigVisSimStorm()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

  % TS Debby Visible images. PNG files have true color, so Cmap vars will be empty (which
  % is exactly what you want for geoshow.
  % Visible images are composites over one day
  [ VisStart VisStartCmap ] = imread('IMAGES/MYBGLSR_C05K.A2006234.006.1.3.png');
  [ VisMid   VisMidCmap   ] = imread('IMAGES/MYBGLSR_C05K.A2006235.006.1.3.png');
  [ VisEnd   VisEndCmap   ] = imread('IMAGES/MYBGLSR_C05K.A2006236.006.1.3.png');

  % Record the Lat,Lon extent of each image so that the image gets placed correctly
  % on the worldmap grid.
  % Format: [ Lat-of-south-edge Lat-of-north-edge Lon-of-west-edge Lon-of-east-edge ]
  VisStartExt = [  0.0 40.9 -41.7   0.0 ];
  VisMidExt   = [  0.0 40.9 -41.7   0.0 ];
  VisEndExt   = [  0.0 40.9 -41.7   0.0 ];

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

  % Visible satellite images in left column
  WmapLatBounds = [ 5 26 ];
  WmapLonBounds = [ -42 -8 ];
  Paxes = subplot(3,2,1);
  PlotFsFigImageWmap(Paxes, VisStart, VisStartCmap, WmapLatBounds, WmapLonBounds, 'a', '22Aug', Fsize, VisStartExt)

  Paxes = subplot(3,2,3);
  PlotFsFigImageWmap(Paxes, VisMid, VisMidCmap, WmapLatBounds, WmapLonBounds, 'c', '23Aug', Fsize, VisMidExt)

  Paxes = subplot(3,2,5);
  PlotFsFigImageWmap(Paxes, VisEnd, VisEndCmap, WmapLatBounds, WmapLonBounds, 'e', '24Aug', Fsize, VisEndExt)

  % Vertically integrated condensate in right column
  Clevs = 0.4:0.2:4.0;
  Cmap = 'jet';
  ClabInc = 2;

  Paxes = subplot(3,2,2);
  PlotFsFigDiagWmap(Paxes, X, Y, COND_START, WmapLatBounds, WmapLonBounds, Clevs, Cmap, ClabInc, 'b', '06Z, 22Aug', Fsize)

  Paxes = subplot(3,2,4);
  PlotFsFigDiagWmap(Paxes, X, Y, COND_MID, WmapLatBounds, WmapLonBounds, Clevs, Cmap, ClabInc, 'd', '12Z, 23Aug', Fsize)

  Paxes = subplot(3,2,6);
  PlotFsFigDiagWmap(Paxes, X, Y, COND_END, WmapLatBounds, WmapLonBounds, Clevs, Cmap, ClabInc, 'f', '18Z, 24Aug', Fsize)

  OutFile = sprintf('%s/FsFigVisSimStorm.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
