function [ ] = PlotWinds(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

Hdir = 'HDF5';
Adir = Config.AzavgDir;

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

Fsize = 20;

SampleTimes = [ 7 31 55 ];
SampleNames = {
 'Aug 22:12Z'
 'Aug 23:12Z'
 'Aug 24:12Z'
 };

Slabels = {
 '(a)'
 '(d)'
 '(g)'
 };

VTlabels = {
 '(c)'
 '(f)'
 '(i)'
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
Wfile = sprintf('%s/w_TSD_3GRIDS.h5', Adir);
Wvar = 'w';
fprintf('Reading: %s, Variable: %s\n', Wfile, Wvar);

Udset = ncgeodataset(Ufile);
Vdset = ncgeodataset(Vfile);
VTdset = ncgeodataset(VTfile);
Wdset = ncgeodataset(Wfile);
Uwind = Udset.geovariable(Uvar);
Vwind = Vdset.geovariable(Vvar);
VTwind = VTdset.geovariable(VTvar);
Wwind = Wdset.geovariable(Wvar);

% Read in coord values (same in both U and V files)
X = Udset.data('x_coords');      % degrees lon
Y = Udset.data('y_coords');      % degrees lat

% Read in Vt coords
R = VTdset.data('x_coords')/1000; % radius km
Z = VTdset.data('z_coords')/1000; % height km
T = VTdset.data('t_coords')/3600; % time hr

% Pick out a good level for the wind barb plots - use ~500m for now since that is where
% the tangential winds are the strongest
WB_Z1 = find(Z <= 0.5, 1, 'last');

% Levels for the Vt plot
VT_R1 = find(R >= 0, 1, 'first');
VT_R2 = find(R <= 400, 1, 'last');
VT_Z1 = find(Z >= 0, 1, 'first');
VT_Z2 = find(Z <= 5, 1, 'last');
W_R1 = find(R >= 0, 1, 'first');
W_R2 = find(R <= 400, 1, 'last');
W_Z1 = find(Z >= 0, 1, 'first');
W_Z2 = find(Z <= 15, 1, 'last');

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
VT_Clims = [ 0 15 ];
W_R = R(W_R1:W_R2);
W_Z = Z(W_Z1:W_Z2);
W_Clims = [ -0.2 0.2 ];

for it = 1:length(SampleTimes)
  T1 = SampleTimes(it);
  Time = T(T1);
  fprintf('Generating stream line plot: T = %d hr\n', Time);


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
  set(gca, 'FontSize', Fsize);

  m_proj('miller', 'lat', LatBounds, 'lon', LonBounds);
  m_coast('color', 'k', 'linewidth', 3); % k --> black
  m_grid('linestyle','none','box','fancy','tickdir','out');

  hold on;
  m_quiver(LON(1:Xinc:end,1:Yinc:end), LAT(1:Xinc:end,1:Yinc:end), ...
             U(1:Xinc:end,1:Yinc:end),   V(1:Xinc:end,1:Yinc:end), ...
          Scale, 'linewidth', 1.5);

  Ptitle = sprintf('%s %s', Slabels{it}, SampleNames{it});
  Thand = title(Ptitle);
  set(Thand, 'Units', 'Normalized');
  set(Thand, 'HorizontalAlignment', 'Left');
  Tpos = get(Thand, 'Position');
  Tpos(1) = 0; % line up with left edge of plot area
  set(Thand, 'Position', Tpos);

  OutFile = sprintf('%s/WindsStream_T%d.jpg', Pdir, it);
  fprintf('  Writing file: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

  
  %%%%%%%%%%%%%%%%%%%%%% Vt PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  VT = double(squeeze(VTwind.data(T1,VT_Z1:VT_Z2,:,VT_R1:VT_R2)));

  Fig = figure;

  contourf(VT_R, VT_Z, VT);
  shading flat;
  cbar = colorbar;
  set(cbar, 'FontSize', Fsize);
  caxis(VT_Clims);
  set(gca, 'FontSize', Fsize);
  set(gca, 'LineWidth', 2); 
  set(gca, 'TickLength', [ 0.025 0.025 ]); 

  Ptitle = sprintf('%s %s', VTlabels{it}, SampleNames{it});
  Thand = title(Ptitle);
  set(Thand, 'Units', 'Normalized');
  set(Thand, 'HorizontalAlignment', 'Left');
  Tpos = get(Thand, 'Position');
  Tpos(1) = 0; % line up with left edge of plot area
  set(Thand, 'Position', Tpos);

  xlabel('Radius (km)');
  ylabel('Height (km)');

  OutFile = sprintf('%s/WindsVt_T%d.jpg', Pdir, it);
  fprintf('  Writing file: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
  
  %%%%%%%%%%%%%%%%%%%%%% W PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  W = double(squeeze(Wwind.data(T1,W_Z1:W_Z2,:,W_R1:W_R2)));

  Fig = figure;

  contourf(W_R, W_Z, W);
  shading flat;
  cbar = colorbar;
  set(cbar, 'FontSize', Fsize);
  colormap(redblue);
  caxis(W_Clims);
  set(gca, 'FontSize', Fsize);
  set(gca, 'LineWidth', 2); 
  set(gca, 'TickLength', [ 0.025 0.025 ]); 

  Ptitle = sprintf('TS Debby (2006), %s: W (m/s)', SampleNames{it});
  title(Ptitle);
  xlabel('Radius (km)');
  ylabel('Height (km)');

  OutFile = sprintf('%s/WindsW_T%d.jpg', Pdir, it);
  fprintf('  Writing file: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end


end

