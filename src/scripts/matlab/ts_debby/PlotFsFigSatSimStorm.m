function [ ] = PlotFsFigVisSimStorm()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 11;

  % get grid3 extent and make map a little bigger
  GlocFname = 'GridLocs.txt';
  Grid3Loc = ReadGridLoc(GlocFname, 'Grid3:');

  MapInc = 4;
  MapLoc = Grid3Loc + [ -MapInc MapInc -MapInc MapInc ];
  MapLats = [ MapLoc(1:2) ];
  MapLons = [ MapLoc(3:4) ];

  % TS Debby Visible images. PNG files have true color, so Cmap vars will be empty (which
  % is exactly what you want for geoshow.
  % Visible images are composites over one day
  [ VisStart VisStartCmap ] = imread('IMAGES/MYBGLSR_C05K.A2006234.006.1.3.png');
  [ VisMid   VisMidCmap   ] = imread('IMAGES/MYBGLSR_C05K.A2006235.006.1.3.png');
  [ VisEnd   VisEndCmap   ] = imread('IMAGES/MYBGLSR_C05K.A2006236.006.1.3.png');

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
  Paxes = subplot(3,3,1);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.05;
  Ploc(2) = Ploc(2) - 0.01;
  Ploc(4) = Ploc(4) + 0.05;
  set(Paxes, 'Position', Ploc);
  PlotFsFigImageWmap(Paxes, VisStart, VisStartCmap, MapLats, MapLons, 'a', '22Aug', Fsize, VisStartExt)

  Paxes = subplot(3,3,4);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.05;
  Ploc(2) = Ploc(2) - 0.04;
  Ploc(4) = Ploc(4) + 0.05;
  set(Paxes, 'Position', Ploc);
  PlotFsFigImageWmap(Paxes, VisMid, VisMidCmap, MapLats, MapLons, 'd', '23Aug', Fsize, VisMidExt)

  Paxes = subplot(3,3,7);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.05;
  Ploc(2) = Ploc(2) - 0.07;
  Ploc(4) = Ploc(4) + 0.05;
  set(Paxes, 'Position', Ploc);
  PlotFsFigImageWmap(Paxes, VisEnd, VisEndCmap, MapLats, MapLons, 'g', '24Aug', Fsize, VisEndExt)

  % Infrared satellite images in middle column
  Paxes = subplot(3,3,2);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.05;
  Ploc(2) = Ploc(2) - 0.01;
  Ploc(4) = Ploc(4) + 0.05;
  set(Paxes, 'Position', Ploc);
  PlotFsFigImageWmap(Paxes, IrStart, IrStartCmap, MapLats, MapLons, 'b', '06Z, 22Aug', Fsize, IrStartExt)

  Paxes = subplot(3,3,5);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.05;
  Ploc(2) = Ploc(2) - 0.04;
  Ploc(4) = Ploc(4) + 0.05;
  set(Paxes, 'Position', Ploc);
  PlotFsFigImageWmap(Paxes, IrMid, IrMidCmap, MapLats, MapLons, 'e', '12Z, 23Aug', Fsize, IrMidExt)

  Paxes = subplot(3,3,8);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.05;
  Ploc(2) = Ploc(2) - 0.07;
  Ploc(4) = Ploc(4) + 0.05;
  set(Paxes, 'Position', Ploc);
  PlotFsFigImageWmap(Paxes, IrEnd, IrEndCmap, MapLats, MapLons, 'h', '18Z, 24Aug', Fsize, IrEndExt)

  % Vertically integrated condensate in right column
  Clevs = 0.4:0.2:4.0;
  Cmap = 'jet';
  ClabInc = 2;

  Paxes = subplot(3,3,3);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.04;
  Ploc(2) = Ploc(2) - 0.01;
  Ploc(4) = Ploc(4) + 0.05;
  set(Paxes, 'Position', Ploc);
  PlotFsFigDiagWmap(Paxes, X, Y, COND_START, MapLats, MapLons, Clevs, Cmap, ClabInc, 'c', '06Z, 22Aug (SD)', Fsize)

  Paxes = subplot(3,3,6);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.04;
  Ploc(2) = Ploc(2) - 0.04;
  Ploc(4) = Ploc(4) + 0.05;
  set(Paxes, 'Position', Ploc);
  PlotFsFigDiagWmap(Paxes, X, Y, COND_MID, MapLats, MapLons, Clevs, Cmap, ClabInc, 'f', '12Z, 23Aug', Fsize)

  Paxes = subplot(3,3,9);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - 0.04;
  Ploc(2) = Ploc(2) - 0.07;
  Ploc(4) = Ploc(4) + 0.05;
  set(Paxes, 'Position', Ploc);
  PlotFsFigDiagWmap(Paxes, X, Y, COND_END, MapLats, MapLons, Clevs, Cmap, ClabInc, 'i', '18Z, 24Aug', Fsize)

  OutFile = sprintf('%s/FsFigSatSimStorm.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
