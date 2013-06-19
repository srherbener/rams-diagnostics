function [ ] = PlotMaxWind(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

Adir = Config.AzavgDir;

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

% read in the sea level pressure
Hfile = sprintf('%s/speed_t_TSD_3GRIDS.h5', Adir);
Hdset = 'speed_t';
fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);
SPEED_T = squeeze(hdf5read(Hfile, Hdset));
TIMES = 1:size(SPEED_T,3);

% ST is organized as r,z,t so take the max along the first two dimensions
% Generate the time series of the maximum of the near surface Vt. Use Z = 2
% since Z = 1 is below surface in RAMS.
% Also, skip the first radial entry. There is a miscalculation in the
% speed_t diagnostic due to TS Debby travelling over the Verde Islands.
R1 = 2;
R2 = size(SPEED_T,1);
ST_MAX = max(squeeze(SPEED_T(R1:R2,2,:)), [] ,1);

% NHC Best Track (every six hours) data
% time step 1 from the simulation is where the NHC data starts
% NHC wind is in mph, multiply by 0.447 to convert to m/s
% Aug 22, 2006, 6Z through Aug 24, 2006, 18Z (11 points)
NHC_WIND  = [ 35   35   35   40   45   50   50   50   50   50   50   ] * 0.447;
NHC_TIMES = (1:6:61);

% plot
OutFile = sprintf('%s/TsDebbyWind.jpg', Pdir);

FigWind = figure;
set(gca, 'FontSize', 18);
SimST = plot(TIMES, ST_MAX, '-r', 'linewi', 2);
hold on;
NhcST = plot(NHC_TIMES, NHC_WIND, '+b', 'linewi', 2);
title('TS Debby Maximum Tangential Wind Speed');
xlabel('Time');
set(gca,'xtick', (13:24:61));
set(gca,'xticklabel', { 'Aug22:0Z', 'Aug23:0Z', 'Aug24:0Z' });
ylabel('Wind Speed (m/s)');
ylim([ 0 40 ]);
legend([ NhcST SimST ], 'NHC Best Track', 'Simulated Storm', 'Location', 'NorthEast');

fprintf('Writing: %s\n', OutFile);
saveas(FigWind, OutFile);
close(FigWind);
