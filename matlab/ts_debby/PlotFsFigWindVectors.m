function [ ] = PlotFsFigWindVectors()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  % Read in sample wind fields from pre sal period (t = 20 hrs) and
  % at z = 3500 m. This is where the enhanced jet appears.
  Tsamp = 20; % sim time hrs
  Zsamp = 3.5; % km

  LatBounds = [ 7 23 ];
  LonBounds = [ -40 -14 ];

  Inc = 3;
  Fsize = 15;
  Qscale = 1.5;
  LineW = 2;

  fprintf('*******************************************************\n');
  fprintf('Plotting Wind Vectors:\n');
  fprintf('\n');

  
  Case = 'TSD_NONSAL_NODUST';
  InFile = sprintf('HDF5/%s/HDF5/u_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
  
  % Get the x and y dimension sizes
  Xinfo = h5info(InFile, '/x_coords');
  Yinfo = h5info(InFile, '/y_coords');
  Nx = (Xinfo.Dataspace.Size + (Inc-1))/Inc;
  Ny = (Yinfo.Dataspace.Size + (Inc-1))/Inc;
  
  % Read in x, y, get the dimension sizes from these vars
  % and set the z and t indices.
  X = squeeze(h5read(InFile, '/x_coords', [ 1 ], [ Nx ], [ Inc ]));
  Y = squeeze(h5read(InFile, '/y_coords', [ 1 ], [ Ny ], [ Inc ]));
  Z = squeeze(h5read(InFile, '/z_coords'))./1000; % km
  T = squeeze(h5read(InFile, '/t_coords'))./3600 - 42; % sim time, hrs

  Z1 = find(Z >= Zsamp, 1, 'first');
  T1 = find(T >= Tsamp, 1, 'first');

  % create a mesh for m_quiver
  [ MY MX ] = meshgrid(Y, X);

  InVname = '/u';
  fprintf('  Reading: %s (%s)\n', InFile, InVname);
  NSND_U = squeeze(h5read(InFile, InVname, [ 1 1 Z1 T1 ], [ Nx Ny 1 1 ], [ Inc Inc 1 1 ]));

  InFile = sprintf('HDF5/%s/HDF5/v_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
  InVname = '/v';
  fprintf('  Reading: %s (%s)\n', InFile, InVname);
  NSND_V = squeeze(h5read(InFile, InVname, [ 1 1 Z1 T1 ], [ Nx Ny 1 1 ], [ Inc Inc 1 1 ]));
  fprintf('\n');

  % Read in the other cases.
  Case = 'TSD_SAL_NODUST';
  InFile = sprintf('HDF5/%s/HDF5/u_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
  InVname = '/u';
  fprintf('  Reading: %s (%s)\n', InFile, InVname);
  SND_U = squeeze(h5read(InFile, InVname, [ 1 1 Z1 T1 ], [ Nx Ny 1 1 ], [ Inc Inc 1 1 ]));

  InFile = sprintf('HDF5/%s/HDF5/v_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
  InVname = '/v';
  fprintf('  Reading: %s (%s)\n', InFile, InVname);
  SND_V = squeeze(h5read(InFile, InVname, [ 1 1 Z1 T1 ], [ Nx Ny 1 1 ], [ Inc Inc 1 1 ]));
  fprintf('\n');

  Case = 'TSD_NONSAL_DUST';
  InFile = sprintf('HDF5/%s/HDF5/u_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
  InVname = '/u';
  fprintf('  Reading: %s (%s)\n', InFile, InVname);
  NSD_U = squeeze(h5read(InFile, InVname, [ 1 1 Z1 T1 ], [ Nx Ny 1 1 ], [ Inc Inc 1 1 ]));

  InFile = sprintf('HDF5/%s/HDF5/v_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
  InVname = '/v';
  fprintf('  Reading: %s (%s)\n', InFile, InVname);
  NSD_V = squeeze(h5read(InFile, InVname, [ 1 1 Z1 T1 ], [ Nx Ny 1 1 ], [ Inc Inc 1 1 ]));
  fprintf('\n');

  Case = 'TSD_SAL_DUST';
  InFile = sprintf('HDF5/%s/HDF5/u_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
  InVname = '/u';
  fprintf('  Reading: %s (%s)\n', InFile, InVname);
  SD_U = squeeze(h5read(InFile, InVname, [ 1 1 Z1 T1 ], [ Nx Ny 1 1 ], [ Inc Inc 1 1 ]));

  InFile = sprintf('HDF5/%s/HDF5/v_lite-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
  InVname = '/v';
  fprintf('  Reading: %s (%s)\n', InFile, InVname);
  SD_V = squeeze(h5read(InFile, InVname, [ 1 1 Z1 T1 ], [ Nx Ny 1 1 ], [ Inc Inc 1 1 ]));
  fprintf('\n');

  % Create the factors
  F1_U = SND_U - NSND_U;
  F1_V = SND_V - NSND_V;

  F2_U = NSD_U - NSND_U;
  F2_V = NSD_V - NSND_V;

  F12_U = SD_U - (SND_U + NSD_U) + NSND_U;
  F12_V = SD_V - (SND_V + NSD_V) + NSND_V;

  % Make a vector plot using quiver command
  OutFile = sprintf('%s/FsFigWindVectors.jpg', Pdir);
  fprintf('  Writing: %s\n', OutFile);
  fprintf('\n');

  Fig = figure;

%  worldmap(LatBounds, LonBounds);
%  setm(gca, 'MapProjection', 'miller');
%  load coast; % creates vars "lat" and "long" that contain coast lines
%  plotm(lat, long, 'Color', 'k', 'LineWidth', LineW);
%  quiverm(double(MY), double(MX), double(V), double(U), Qscale);
  Paxes = subplot(2,2,1);
  PlotFsFigVectorWmap(gca, Y, X, NSND_V, NSND_U, LatBounds, LonBounds, 'a', 'NSND', Fsize, Qscale);

  Paxes = subplot(2,2,2);
  PlotFsFigVectorWmap(gca, Y, X, F1_V, F1_U, LatBounds, LonBounds, 'b', 'F1', Fsize, Qscale);

  Paxes = subplot(2,2,3);
  PlotFsFigVectorWmap(gca, Y, X, F2_V, F2_U, LatBounds, LonBounds, 'c', 'F2', Fsize, Qscale);

  Paxes = subplot(2,2,4);
  PlotFsFigVectorWmap(gca, Y, X, F12_V, F12_U, LatBounds, LonBounds, 'd', 'F12', Fsize, Qscale);

  tightfig(Fig);

  saveas(Fig,OutFile);
  close(Fig);
end
