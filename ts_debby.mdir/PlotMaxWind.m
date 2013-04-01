function [ ] = PlotMaxWind(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

Tdir = Config.TsavgDir;

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

% read in the sea level pressure
Hfile = sprintf('%s/speed_t_TSD_3GRIDS.h5', Tdir);
Hdset = 'speed_t';
fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);
SPEED_T = hdf5read(Hfile, Hdset);

% Throw away first 36 hrs of the simulation since the circulation hadn't
% formed yet.
ST = SPEED_T(:,:,37:end);
whos
TIMES = 1:size(ST,3);

% generate the time series of the maximum Vt of the horizontal domain
% ST is organized as x,y,t so take the max along the first two dimensions
ST_MAX = squeeze(max(max(ST,[],1),[],2));

% NHC Best Track (every six hours) data
% time step 43 from the simulation is where the NHC data starts
% NHC wind is in mph, multiply by 0.447 to convert to m/s
NHC_WIND  = [ 35   35   35   35   35   40   45   50   50   50   50   50   50   ] * 0.447;
NHC_TIMES = (7:6:79);

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
