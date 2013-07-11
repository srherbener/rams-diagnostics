function [ ] = PlotHydro(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

Hdir = 'HDF5';
Adir = Config.AzavgDir;

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

Fsize = 20;

SampleTimes = [ 1 25 49 ];
SampleNames = {
 'Aug 22, 06Z'
 'Aug 23, 06Z'
 'Aug 24, 06Z'
 };

CLlabels = {
 '(a)'
 '(c)'
 '(e)'
 };

RNlabels = {
 '(b)'
 '(d)'
 '(f)'
 };

AGlabels = {
 '(b)'
 '(d)'
 '(f)'
 };

GRlabels = {
 '(b)'
 '(d)'
 '(f)'
 };

PRlabels = {
 '(b)'
 '(e)'
 '(h)'
 };

% Get set up to read in horizontal winds
CLfile = sprintf('%s/cloud_TSD_3GRIDS.h5', Adir);
CLvar = 'cloud';
fprintf('Reading: %s, Variable: %s\n', CLfile, CLvar);
RNfile = sprintf('%s/rain_TSD_3GRIDS.h5', Adir);
RNvar = 'rain';
fprintf('Reading: %s, Variable: %s\n', RNfile, RNvar);
AGfile = sprintf('%s/aggr_TSD_3GRIDS.h5', Adir);
AGvar = 'aggr';
fprintf('Reading: %s, Variable: %s\n', AGfile, AGvar);
GRfile = sprintf('%s/graup_TSD_3GRIDS.h5', Adir);
GRvar = 'graup';
fprintf('Reading: %s, Variable: %s\n', GRfile, GRvar);
PRfile = sprintf('%s/pcprate-TSD_3GRIDS-AS-2006-08-20-120000-g3.h5', Hdir);
PRvar = 'pcprate';
fprintf('Reading: %s, Variable: %s\n', PRfile, PRvar);
VCfile = sprintf('%s/vint_cond_TSD_3GRIDS.h5', Adir);
VCvar = 'vint_cond';
fprintf('Reading: %s, Variable: %s\n', VCfile, VCvar);

CLdset = ncgeodataset(CLfile);
CLhydro = CLdset.geovariable(CLvar);
RNdset = ncgeodataset(RNfile);
RNhydro = RNdset.geovariable(RNvar);
AGdset = ncgeodataset(AGfile);
AGhydro = AGdset.geovariable(AGvar);
GRdset = ncgeodataset(GRfile);
GRhydro = GRdset.geovariable(GRvar);
PRdset = ncgeodataset(PRfile);
PRhydro = PRdset.geovariable(PRvar);
VCdset = ncgeodataset(VCfile);
VChydro = VCdset.geovariable(VCvar);

% coords
R   = CLdset.data('x_coords')/1000; % radius km
Z   = CLdset.data('z_coords')/1000; % height km
T   = PRdset.data('t_coords')/3600; % time hr
LON = PRdset.data('x_coords');      % degrees lon
LAT = PRdset.data('y_coords');      % degrees lat
LonBounds = [ -40 -13.5 ];
LatBounds = [   7  24   ];

% Data Selection
CL_R1 = find(R >= 0, 1, 'first');
CL_R2 = find(R <= 400, 1, 'last');
CL_Z1 = find(Z >= 0, 1, 'first');
CL_Z2 = find(Z <= 15, 1, 'last');
RN_R1 = find(R >= 0, 1, 'first');
RN_R2 = find(R <= 400, 1, 'last');
RN_Z1 = find(Z >= 0, 1, 'first');
RN_Z2 = find(Z <= 15, 1, 'last');
AG_R1 = find(R >= 0, 1, 'first');
AG_R2 = find(R <= 400, 1, 'last');
AG_Z1 = find(Z >= 0, 1, 'first');
AG_Z2 = find(Z <= 15, 1, 'last');
GR_R1 = find(R >= 0, 1, 'first');
GR_R2 = find(R <= 400, 1, 'last');
GR_Z1 = find(Z >= 0, 1, 'first');
GR_Z2 = find(Z <= 15, 1, 'last');
VC_R1 = find(R >= 0, 1, 'first');
VC_R2 = find(R <= 400, 1, 'last');

% For plots
CL_R = R(CL_R1:CL_R2);
CL_Z = Z(CL_Z1:CL_Z2);
CL_Clims = [ 0 0.3 ];
RN_R = R(RN_R1:RN_R2);
RN_Z = Z(RN_Z1:RN_Z2);
RN_Clims = [ 0 0.4 ];
AG_R = R(AG_R1:AG_R2);
AG_Z = Z(AG_Z1:AG_Z2);
AG_Clims = [ 0 0.5 ];
GR_R = R(GR_R1:GR_R2);
GR_Z = Z(GR_Z1:GR_Z2);
GR_Clims = [ 0 0.05 ];
VC_R = R(VC_R1:VC_R2);
VC_Clims = [ 0 3 ];

PR_Clims = [ 0 20 ];

for it = 1:length(SampleTimes)
  T1 = SampleTimes(it);
  Time = T(T1);
  fprintf('Generating hydrometeor plots: T = %d hr\n', Time);

  %%%%%%%%%%%%%%%%%%%%%% Precip Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % The nctoolbox routines will preserve the dimension order in the
  % hdf5 file (as opposed to hdf5read which reverses the dimensions).
  % This works out since we need lat (y) in the rows and lon (x) in
  % the columns.
  PR = squeeze(PRhydro.data(T1,:,:,:));
  PR(PR < 0.3) = nan;

  Fig = figure;
  set(gca, 'FontSize', Fsize);

  m_proj('miller', 'lat', LatBounds, 'lon', LonBounds);
  m_coast('color', 'k', 'linewidth', 3); % k --> black
  m_grid('linestyle','none','box','fancy','tickdir','out');

  hold on;
  m_contourf(LON, LAT, PR);
  cbar = colorbar;
  set(cbar, 'FontSize', Fsize);
  caxis(PR_Clims);

  Ptitle = sprintf('%s %s', PRlabels{it}, SampleNames{it});
  Thand = title(Ptitle);
  set(Thand, 'Units', 'Normalized');
  set(Thand, 'HorizontalAlignment', 'Left');
  Tpos = get(Thand, 'Position');
  Tpos(1) = 0; % line up with left edge of plot area
  set(Thand, 'Position', Tpos);

  OutFile = sprintf('%s/PrecipRate_T%d.jpg', Pdir, Time);
  fprintf('  Writing file: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

  %%%%%%%%%%%%%%%%%%%%%% Cloud PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  CL = double(squeeze(CLhydro.data(T1,CL_Z1:CL_Z2,:,CL_R1:CL_R2)));

  Fig = figure;
  set(gca, 'FontSize', Fsize);

  contourf(CL_R, CL_Z, CL);
  shading flat;
  cbar = colorbar;
  set(cbar, 'FontSize', Fsize);
  caxis(CL_Clims);

  Ptitle = sprintf('%s %s', CLlabels{it}, SampleNames{it});
  Thand = title(Ptitle);
  set(Thand, 'Units', 'Normalized');
  set(Thand, 'HorizontalAlignment', 'Left');
  Tpos = get(Thand, 'Position');
  Tpos(1) = 0; % line up with left edge of plot area
  set(Thand, 'Position', Tpos);

  xlabel('Radius (km)');
  ylabel('Height (km)');

  OutFile = sprintf('%s/HydroCloud_T%d.jpg', Pdir, Time);
  fprintf('  Writing file: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
  
  %%%%%%%%%%%%%%%%%%%%%% Rain PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  RN = double(squeeze(RNhydro.data(T1,RN_Z1:RN_Z2,:,RN_R1:RN_R2)));

  Fig = figure;
  set(gca, 'FontSize', Fsize);

  contourf(RN_R, RN_Z, RN);
  shading flat;
  cbar = colorbar;
  set(cbar, 'FontSize', Fsize);
  caxis(RN_Clims);

  Ptitle = sprintf('%s %s', RNlabels{it}, SampleNames{it});
  Thand = title(Ptitle);
  set(Thand, 'Units', 'Normalized');
  set(Thand, 'HorizontalAlignment', 'Left');
  Tpos = get(Thand, 'Position');
  Tpos(1) = 0; % line up with left edge of plot area
  set(Thand, 'Position', Tpos);

  xlabel('Radius (km)');
  ylabel('Height (km)');

  OutFile = sprintf('%s/HydroRain_T%d.jpg', Pdir, Time);
  fprintf('  Writing file: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

  %%%%%%%%%%%%%%%%%%%%%% Aggregate PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  AG = double(squeeze(AGhydro.data(T1,AG_Z1:AG_Z2,:,AG_R1:AG_R2)));

  Fig = figure;
  set(gca, 'FontSize', Fsize);

  contourf(AG_R, AG_Z, AG);
  shading flat;
  cbar = colorbar;
  set(cbar, 'FontSize', Fsize);
  caxis(AG_Clims);

  Ptitle = sprintf('%s %s', AGlabels{it}, SampleNames{it});
  Thand = title(Ptitle);
  set(Thand, 'Units', 'Normalized');
  set(Thand, 'HorizontalAlignment', 'Left');
  Tpos = get(Thand, 'Position');
  Tpos(1) = 0; % line up with left edge of plot area
  set(Thand, 'Position', Tpos);

  xlabel('Radius (km)');
  ylabel('Height (km)');

  OutFile = sprintf('%s/HydroAggr_T%d.jpg', Pdir, Time);
  fprintf('  Writing file: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

  %%%%%%%%%%%%%%%%%%%%%% Graupel PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  GR = double(squeeze(GRhydro.data(T1,GR_Z1:GR_Z2,:,GR_R1:GR_R2)));

  Fig = figure;
  set(gca, 'FontSize', Fsize);

  contourf(GR_R, GR_Z, GR);
  shading flat;
  cbar = colorbar;
  set(cbar, 'FontSize', Fsize);
  caxis(GR_Clims);

  Ptitle = sprintf('%s %s', GRlabels{it}, SampleNames{it});
  Thand = title(Ptitle);
  set(Thand, 'Units', 'Normalized');
  set(Thand, 'HorizontalAlignment', 'Left');
  Tpos = get(Thand, 'Position');
  Tpos(1) = 0; % line up with left edge of plot area
  set(Thand, 'Position', Tpos);

  xlabel('Radius (km)');
  ylabel('Height (km)');

  OutFile = sprintf('%s/HydroGraup_T%d.jpg', Pdir, Time);
  fprintf('  Writing file: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end

% create hovmoller plot of vertically integrated condensate
VC = double(squeeze(VChydro.data(:,:,:,VC_R1:VC_R2)));

Fig = figure;
set(gca, 'FontSize', Fsize);

contourf(VC_R, T, VC);
shading flat;
cbar = colorbar;
set(cbar, 'FontSize', Fsize);
caxis(VC_Clims);

Ptitle = '(c)';
Thand = title(Ptitle);
set(Thand, 'Units', 'Normalized');
set(Thand, 'HorizontalAlignment', 'Left');
Tpos = get(Thand, 'Position');
Tpos(1) = 0; % line up with left edge of plot area
set(Thand, 'Position', Tpos);

xlabel('Radius (km)');
ylabel('Time (hr)');

OutFile = sprintf('%s/HovVintCond.jpg', Pdir);
fprintf('  Writing file: %s\n', OutFile);
saveas(Fig, OutFile);
close(Fig);


end

