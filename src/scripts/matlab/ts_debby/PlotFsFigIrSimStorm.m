function [ ] = PlotFsFigIrSimStorm()

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
  WmapLatBounds = [ 5 26 ];
  WmapLonBounds = [ -42 -8 ];

  Paxes = subplot(3,2,1);
  PlotFsFigImageWmap(Paxes, IrStart, IrStartCmap, WmapLatBounds, WmapLonBounds, 'a', '06Z, 22Aug', Fsize, IrStartExt)

  Paxes = subplot(3,2,3);
  PlotFsFigImageWmap(Paxes, IrMid, IrMidCmap, WmapLatBounds, WmapLonBounds, 'c', '12Z, 23Aug', Fsize, IrMidExt)

  Paxes = subplot(3,2,5);
  PlotFsFigImageWmap(Paxes, IrEnd, IrEndCmap, WmapLatBounds, WmapLonBounds, 'e', '18Z, 24Aug', Fsize, IrEndExt)

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

  OutFile = sprintf('%s/FsFigIrSimStorm.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
