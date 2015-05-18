function [ ] = PlotSpeedtThetae(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

Adir = Config.AzavgDir;

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

% read in the speed and theta_e
Hfile = sprintf('%s/speed_t_TSD_3GRIDS.h5', Adir);
Hdset = 'speed_t';
fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);
SPEED_T = squeeze(hdf5read(Hfile, Hdset));
R = hdf5read(Hfile, 'x_coords')/1000; % radius in km
Z = hdf5read(Hfile, 'z_coords')/1000; % height in km

Hfile = sprintf('%s/theta_e_TSD_3GRIDS.h5', Adir);
Hdset = 'theta_e';
fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);
THETA_E = squeeze(hdf5read(Hfile, Hdset));

% Take samples of the wind speed diagnostic. Each sample will be a vertical
% slice showing the azimuthally averaged tangetial wind speed. Time steps:
%    13 --> Aug 21, 0Z
%    37 --> Aug 22, 0Z
%    61 --> Aug 23, 0Z
Psamples = [ 13 37 61 ];
ST = SPEED_T(:,:,Psamples);
TE = THETA_E(:,:,Psamples);

% Tangential wind
FigSpeedT = figure;
Heights = (2:30);
Clevs = (0:2:20);
TextSize = 11;

ax = subplot(3,1,1);
set(gca, 'FontSize', TextSize);
contourf(R, Z(Heights),ST(:,Heights,1)',Clevs);
shading flat;
title( { 'TsDebby: Azimuthally Averaged Tangential Wind Speed', 'Aug 21, 0Z' } );
%xlabel('Radius (km)');
ylabel('Height (km)');
pos = get(ax, 'Position');
set (ax, 'Position', [ pos(1) pos(2) 0.85*pos(3) pos(4) ]); % move the right edge over to the left

ax = subplot(3,1,2);
set(gca, 'FontSize', TextSize);
contourf(R, Z(Heights),ST(:,Heights,2)', Clevs);
shading flat;
title( { 'Aug 22, 0Z' } );
%xlabel('Radius (km)');
ylabel('Height (km)');
pos = get(ax, 'Position');
set (ax, 'Position', [ pos(1) pos(2) 0.85*pos(3) pos(4) ]);

ax = subplot(3,1,3);
set(gca, 'FontSize', TextSize);
contourf(R, Z(Heights),ST(:,Heights,3)', Clevs);
shading flat;
title( { 'Aug 23, 0Z' } );
xlabel('Radius (km)');
ylabel('Height (km)');
pos = get(ax, 'Position');
set (ax, 'Position', [ pos(1) pos(2) 0.85*pos(3) pos(4) ]);

cbar = colorbar;
set (cbar, 'Position', [ 0.8314 0.11 0.0581 0.8150 ] );
set (cbar, 'FontSize', TextSize);

OutFile = sprintf('%s/TsDebbySpeedT.jpg', Pdir);
fprintf('Writing: %s\n', OutFile);
saveas(FigSpeedT, OutFile);
close(FigSpeedT);

% Theta-E
FigThetaE = figure;
set(gca, 'FontSize', 18);
Heights = (2:30);
Clevs = (340:2:360);
subplot(3,1,1);
contourf(R, Z(Heights),TE(:,Heights,1)', Clevs);
shading flat;
colorbar;
title( { 'TsDebby: Azimuthally Averaged Theta-E', 'Aug 21, 0Z' } );
%xlabel('Radius (km)');
ylabel('Height (km)');
subplot(3,1,2);
contourf(R, Z(Heights),TE(:,Heights,2)', Clevs);
shading flat;
colorbar;
title( { 'Aug 22, 0Z' } );
%xlabel('Radius (km)');
ylabel('Height (km)');
subplot(3,1,3);
contourf(R, Z(Heights),TE(:,Heights,3)', Clevs);
shading flat;
colorbar;
title( { 'Aug 23, 0Z' } );
xlabel('Radius (km)');
ylabel('Height (km)');

OutFile = sprintf('%s/TsDebbyThetaE.jpg', Pdir);
fprintf('Writing: %s\n', OutFile);
saveas(FigThetaE, OutFile);
close(FigThetaE);

end

