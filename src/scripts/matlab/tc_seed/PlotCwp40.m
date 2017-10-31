%script to plot cloud water path at simtime 40 (1st time step of the T2
%file: vert_int-TCS_CNTL...

clear;

path(path,'./m_map');

% Grab the cloud water path directly from REVU output
CWP_TS41_73 = hdf5read('../TCS_CNTL/HDF5/T2/vint_cloud-TCS_CNTL-AS-1998-08-22-120000-g3.h5','/vertint_cloud');
CWP40 = squeeze(CWP_TS41_73(:,:,1))';
%CWP40(CWP40 < 2) = nan;

Lon = hdf5read('../TCS_CNTL/HDF5/T2/vint_cloud-TCS_CNTL-AS-1998-08-22-120000-g3.h5','/x_coords')';
Lat = hdf5read('../TCS_CNTL/HDF5/T2/vint_cloud-TCS_CNTL-AS-1998-08-22-120000-g3.h5','/y_coords')';

% Grab the 3D picture of cloud water from VAPOR
CLOUD_3D = imread('../VAPOR/cloud-TCS_CNTL-simtime40hr.png');

% Grab the azimuthally averaged tangential wind speed from azavg output
AZAVG_TWIND_TS41_73 = hdf5read('../AzAveragedData/speed_t_TCS_CNTL.h5','/speed_t');
AZAVG_TWIND40 = squeeze(AZAVG_TWIND_TS41_73(:,:,1))';

Levels = hdf5read('../TCS_CNTL/HDF5/T2/u-TCS_CNTL-AS-1998-08-22-120000-g3.h5','/z_coords')' / 1000;
Radii = (0:4:196);

% Grab the azimuthally averaged Theta-E from azavg ouptput
AZAVG_THETAE_TS41_73 = hdf5read('../AzAveragedData/theta_e_TCS_CNTL.h5','/theta_e');
AZAVG_THETAE40 = squeeze(AZAVG_THETAE_TS41_73(:,:,1))';

% plot - want some room at the top to put text for a collective title for
% all four panels. The standard position (from running subplot(2,2,1) then
% get(gca,'position') on all four panels is:
%   subplot     llx      lly   width   height
%    2,2,1     0.1300  0.5838  0.3347  0.3412
%    2,2,2     0.5703  0.5838  0.3347  0.3412
%    2,2,3     0.1300  0.1100  0.3347  0.3412
%    2,2,4     0.5703  0.1100  0.3347  0.3412
%
% Want to put a "0.1" space at the top so reduce all the y values
% by 0.95 (each of two rows gets reduced by 0.05).
p1 = [ 0.1300 0.5546 0.3347 0.3241 ];
p2 = [ 0.5703 0.5546 0.3347 0.3241 ];
p3 = [ 0.1300 0.1100 0.3347 0.3241 ];
p4 = [ 0.5703 0.1100 0.3347 0.3241 ];

FigTC = figure;
Fsize = 11;
TFsize = 14;

subplot('position', p1);
set(gca, 'FontSize', Fsize);
Clevs = (0:0.1:1);
contourf(Lon,Lat,CWP40,Clevs);
shading flat;
colorbar;
title('Cloud Water Path (mm)','FontSize',TFsize,'FontWeight', 'bold');
xlabel('Longitude');
ylabel('Latitude');

subplot('position', p2);
set(gca, 'FontSize', Fsize);
image(CLOUD_3D);
set(gca, 'Xtick', []);
set(gca, 'XtickLabel', []);
set(gca, 'Ytick', []);
set(gca, 'YtickLabel', []);
title('Cloud Water Mixing Ratio (g/kg)','FontSize',TFsize,'FontWeight', 'bold');
xlabel('3D Surface shown: 0.06 g/kg');

subplot('position',p3);
x1 = 1;
x2 = 30;
z1 = 2;
z2 = 30;
set(gca, 'FontSize', Fsize);
contourf(Radii(x1:x2), Levels(z1:z2), AZAVG_TWIND40(z1:z2,x1:x2),15);
shading flat;
colorbar;
title('Avg Tan Wind (m/s)','FontSize',TFsize,'FontWeight', 'bold');
xlabel('Radius (km)');
ylabel('Height (km)');

subplot('position',p4);
x1 = 1;
x2 = 30;
z1 = 2;
z2 = 14;
set(gca, 'FontSize', Fsize);
contourf(Radii(x1:x2), Levels(z1:z2), AZAVG_THETAE40(z1:z2,x1:x2),15);
shading flat;
colorbar;
title('Avg Theta-E (K)','FontSize',TFsize,'FontWeight', 'bold');
xlabel('Radius (km)');
ylabel('Height (km)');

% text title
axes('position',[0,0,1,1],'visible','off');
text(0.2400, 0.9600, 'Idealized TC Simulation: t = 40hrs ','FontSize', 16, 'FontWeight', 'bold');

saveas(FigTC, 'TCplots.jpg');
