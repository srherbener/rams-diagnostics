function [ ] = PlotSpeedtThetae(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

Hdir = 'HDF5';
Adir = Config.AzavgDir;

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

SampleTimes = [ 7 31 55 ];
SampleNames = {
 'Aug 22, 12Z'
 'Aug 23, 12Z'
 'Aug 24, 12Z'
 };

% Get set up to read in horizontal winds
Ufile = sprintf('%s/u-TSD_3GRIDS-AS-2006-08-20-120000-g3.h5', Hdir);
Uvar = 'u';
fprintf('Reading: %s, Variable: %s\n', Ufile, Uvar);
Vfile = sprintf('%s/v-TSD_3GRIDS-AS-2006-08-20-120000-g3.h5', Hdir);
Vvar = 'v';
fprintf('Reading: %s, Variable: %s\n', Vfile, Vvar);
VTfile = sprintf('%s/speed_t_TSD_3GRIDS.h5', Adir);
VTvar = 'speed_t';
fprintf('Reading: %s, Variable: %s\n', VTfile, VTvar);

Udset = ncgeodataset(Ufile);
Vdset = ncgeodataset(Vfile);
VTdset = ncgeodataset(VTfile);
Uwind = Udset.geovariable(Uvar);
Vwind = Vdset.geovariable(Vvar);
VTwind = VTdset.geovariable(VTvar);

% Read in coord values (same in both U and V files)
X = Udset.data('x_coords');      % degrees lon
Y = Udset.data('y_coords');      % degrees lat
Z = Udset.data('z_coords')/1000; % height km
T = Udset.data('t_coords')/3600; % time hr

% Read in Vt coords
R = VTdset.data('x_coords')/1000; % radius km

% Pick out a good level for the wind barb plots - use ~500m for now since that is where
% the tangential winds are the strongest
WB_Z1 = find(Z <= 0.5, 1, 'last');

% Levels for the Vt plot
VT_R1 = find(R >= 0, 1, 'first');
VT_R2 = find(R <= 300, 1, 'last');
VT_Z1 = find(Z >= 0, 1, 'first');
VT_Z2 = find(Z <= 5, 1, 'last');

% Create wind barb plots: one for each time sample
% m_quiver wants coordinates in a meshgrid style
[ LON, LAT ] = meshgrid(X, Y);
LON = double(LON);
LAT = double(LAT);

LonBounds = [ -40 -13.5 ];
LatBounds = [   7  24   ];

% m_quiver tries to draw an arrow for every point - way too dense
% so need these factors to thin out wind data
Xinc = 40;
Yinc = 40;
Scale = 1;

% For the Vt plots
VT_R = R(VT_R1:VT_R2);
VT_Z = Z(VT_Z1:VT_Z2);
Clims = [ 0 20 ];

for it = 1:length(SampleTimes)
  T1 = SampleTimes(it);
  Time = T(T1);
  fprintf('Generating stream line plot: T = %d hr\n', Time);

  Ptitle = sprintf('TS Debby (2006), %s', SampleNames{it});

  %%%%%%%%%%%%%%%%%%%%%% WIND BARB PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read in horizontal winds at the given time and level.
  % The nctoolbox routines will preserve the dimension order in the
  % hdf5 file (as opposed to hdf5read which reverses the dimensions).
  % This works out since we need lat (y) in the rows and lon (x) in
  % the columns.
  % Need double for m_quiver.
  U = double(squeeze(Uwind.data(T1,WB_Z1,:,:)));
  V = double(squeeze(Vwind.data(T1,WB_Z1,:,:)));

  Fig = figure;
  set(gca, 'FontSize', 18);

  m_proj('miller', 'lat', LatBounds, 'lon', LonBounds);
  m_coast('color', 'k', 'linewidth', 3); % k --> black
  m_grid('linestyle','none','box','fancy','tickdir','out');

  hold on;
  m_quiver(LON(1:Xinc:end,1:Yinc:end), LAT(1:Xinc:end,1:Yinc:end), ...
             U(1:Xinc:end,1:Yinc:end),   V(1:Xinc:end,1:Yinc:end), ...
          Scale, 'linewidth', 1.5);

  title(Ptitle);

  OutFile = sprintf('%s/WindsStream_T%d.jpg', Pdir, Time);
  fprintf('  Writing file: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

  
  %%%%%%%%%%%%%%%%%%%%%% Vt PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  VT = double(squeeze(VTwind.data(T1,VT_Z1:VT_Z2,:,VT_R1:VT_R2)));

  Fig = figure;
  set(gca, 'FontSize', 18);

  contourf(VT_R, VT_Z, VT);
  shading flat;
  colorbar;
  caxis(Clims);

  title(Ptitle);
  xlabel('Radius (km)');
  ylabel('Height (km)');

  OutFile = sprintf('%s/WindsVt_T%d.jpg', Pdir, Time);
  fprintf('  Writing file: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end


end

