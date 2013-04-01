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
SPEED_T = hdf5read(Hfile, Hdset);

Hfile = sprintf('%s/theta_e_TSD_3GRIDS.h5', Adir);
Hdset = 'theta_e';
fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);
THETA_E = hdf5read(Hfile, Hdset);

% Take samples of the wind speed diagnostic. Each sample will be a vertical
% slice showing the azimuthally averaged tangetial wind speed. Time steps:
%    49 --> Aug 22, 0Z
%    73 --> Aug 23, 0Z
%    97 --> Aug 24, 0Z
Psamples = [ 49 73 97 ];
ST = SPEED_T(:,:,Psamples);
TE = THETA_E(:,:,Psamples);

% x and y coordinate values

Xcoords = (4:4:200);

Ycoords = [ -23.000  24.519    76.481   132.599   193.207   258.664   329.357   405.705   488.162   577.215   673.392   777.263   889.444  1010.600  1141.448  1282.764  1435.385  1600.216  1778.233  1970.492  2178.131  2402.382  2644.573  2906.138  3188.629  3493.719  3823.217  4179.074  4563.400  4978.473  5426.750  5910.890  6433.761  6998.461  7608.338  8267.005  8978.364  9746.633 10576.363 11472.807 12442.610 13447.062 14447.062 15447.062 16447.063 17447.063 18447.063 19447.063 20447.063 21447.063 22447.063 23447.063 24447.063 25447.063 26447.063 27447.063 ];


% Tangential wind
FigSpeedT = figure;
Heights = (2:25);
Clevs = (0:2:20);
TextSize = 11;

ax = subplot(3,1,1);
set(gca, 'FontSize', TextSize);
contourf(Xcoords, Ycoords(Heights),ST(:,Heights,1)',Clevs);
title( { 'TsDebby: Azimuthally Averaged Tangential Wind Speed', 'Aug 22, 0Z' } );
%xlabel('Radius (km)');
ylabel('Height (m)');
pos = get(ax, 'Position');
set (ax, 'Position', [ pos(1) pos(2) 0.85*pos(3) pos(4) ]); % move the right edge over to the left

ax = subplot(3,1,2);
set(gca, 'FontSize', TextSize);
contourf(Xcoords, Ycoords(Heights),ST(:,Heights,2)', Clevs);
title( { 'Aug 23, 0Z' } );
%xlabel('Radius (km)');
ylabel('Height (m)');
pos = get(ax, 'Position');
set (ax, 'Position', [ pos(1) pos(2) 0.85*pos(3) pos(4) ]);

ax = subplot(3,1,3);
set(gca, 'FontSize', TextSize);
contourf(Xcoords, Ycoords(Heights),ST(:,Heights,3)', Clevs);
title( { 'Aug 24, 0Z' } );
xlabel('Radius (km)');
ylabel('Height (m)');
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
Heights = (2:25);
Clevs = (340:2:360);
subplot(3,1,1);
contourf(Xcoords, Ycoords(Heights),TE(:,Heights,1)', Clevs);
colorbar;
ylabel('Height (m)');
subplot(3,1,2);
contourf(Xcoords, Ycoords(Heights),TE(:,Heights,2)', Clevs);
colorbar;
ylabel('Height (m)');
subplot(3,1,3);
contourf(Xcoords, Ycoords(Heights),TE(:,Heights,3)', Clevs);
colorbar;
xlabel('Radius (km)');
ylabel('Height (m)');

OutFile = sprintf('%s/TsDebbyThetaE.jpg', Pdir);
fprintf('Writing: %s\n', OutFile);
saveas(FigThetaE, OutFile);
close(FigThetaE);


%title('TS Debby Minumum SLP');
%xlabel('Time');
%set(gca,'xtick', (13:24:61));
%set(gca,'xticklabel', { 'Aug22:0Z', 'Aug23:0Z', 'Aug24:0Z' });
%ylabel('Pressure (mb)');
%ylim([ 995 1010 ]);
%legend([ NhcSLP SimSLP ], 'NHC Best Track', 'Simulated Storm', 'Location', 'NorthEast');



