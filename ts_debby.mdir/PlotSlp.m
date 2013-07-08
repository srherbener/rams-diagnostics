function [ ] = PlotSlp(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end


% read in the sea level pressure
Hfile = sprintf('AzAveragedData/hist_press_TSD_3GRIDS.h5');
Hdset = 'press';
fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);
HDATA = hdf5read(Hfile, Hdset);
R = hdf5read(Hfile, 'x_coords')/1000;  % radius (km)
P = hdf5read(Hfile, 'y_coords');  % pressure
Z = hdf5read(Hfile, 'z_coords');  % height
T = hdf5read(Hfile, 't_coords');  % times

% Select data
%  Radius from zero to 50km
%  All pressure bins
%  First model level above ground
%  All times
R1 = find(R >= 0, 1, 'first');
R2 = find(R <= 20, 1, 'last');
P1 = 1;
P2 = length(P);
Z1 = find(Z >= 0, 1, 'first');
T1 = 1;
T2 = length(T);

SPRESS = squeeze(HDATA(R1:R2,P1:P2,Z1,T1:T2));
TIMES = T1:T2;

% SPRESS is organized as r,p,t where each entry contains a count of the
% occurrences of pressure value, p, at radius, r. To get the miniumum
% pressure, at each time step, find the smallest pressure bin that has a
% non-zero count across all radii and use the pressure value for that bin.
SPRESS = squeeze(sum(SPRESS,1));

% SPRESS is now organized as (p,t). For each time step, find the smallest
% pressure bin that has a non-zero count, and use that bin to look up the
% pressure value associtated with it. This will be the minimum pressure
% for that time step.
Nt = length(TIMES);
SLP_MIN = zeros([1 Nt]);
for it = 1:Nt
  Pind = find(SPRESS(:,it) > 0, 1, 'first');
  SLP_MIN(it) = P(Pind);
end

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
ylim([ 995 1010 ]);
legend([ NhcSLP SimSLP ], 'NHC Best Track', 'Simulated Storm', 'Location', 'NorthEast');

fprintf('Writing: %s\n', OutFile);
saveas(FigPress, OutFile);
close(FigPress);
