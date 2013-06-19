function [ ] = PlotSlp(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

% read in the sea level pressure
Hfile = sprintf('AzAveragedData/sea_press_TSD_3GRIDS.h5');
Hdset = 'sea_press';
fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);
SPRESS = squeeze(hdf5read(Hfile, Hdset));
TIMES = 1:size(SPRESS,2);

% SPRESS is organized as r,t so take the minimum along the first dimension
% generate the time series of the minimum SLP of the horizontal domain
R1 = 1;
R2 = size(SPRESS,1);
SLP_MIN = min(SPRESS(R1:R2,:),[],1);

% NHC Best Track (every six hours) data
% time step 1 from the simulation is where the NHC data starts
% Aug 22, 2006, 6Z through Aug 24, 2006, 18Z (11 points)
NHC_SLP = [ 1007 1007 1005 1003 1002 1001 1001 1000  999 1000 1000 ];
NHC_TIMES = (1:6:61);


% plot
OutFile = sprintf('%s/TsDebbyPress.jpg', Pdir);

FigPress = figure;
set(gca, 'FontSize', 18);
SimSLP = plot(TIMES, SLP_MIN, '-r', 'linewi', 2);
hold on;
NhcSLP = plot(NHC_TIMES, NHC_SLP, '+b', 'linewi', 2);
title('TS Debby Minumum SLP');
xlabel('Time');
set(gca,'xtick', (13:24:61));
set(gca,'xticklabel', { 'Aug22:0Z', 'Aug23:0Z', 'Aug24:0Z' });
ylabel('Pressure (mb)');
ylim([ 990 1010 ]);
legend([ NhcSLP SimSLP ], 'NHC Best Track', 'Simulated Storm', 'Location', 'NorthEast');

fprintf('Writing: %s\n', OutFile);
saveas(FigPress, OutFile);
close(FigPress);
