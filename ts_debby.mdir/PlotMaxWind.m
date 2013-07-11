function [ ] = PlotMaxWind(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

Adir = Config.AzavgDir;

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

% read in the tangential wind
Hfile = sprintf('AzAveragedData/hist_speed_t_TSD_3GRIDS.h5');
Hdset = 'speed_t';
fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);
HDATA = hdf5read(Hfile, Hdset);
R = hdf5read(Hfile, 'x_coords')/1000;  % radius (km)
S = hdf5read(Hfile, 'y_coords');  % wind speed
Z = hdf5read(Hfile, 'z_coords')/1000;  % height (km)
T = hdf5read(Hfile, 't_coords');  % times

% Select data
%  All radii
%  All pressure bins
%  First model level above ground
%  All times
R1 = 1;
R2 = length(R);
S1 = 1;
S2 = length(S);
Z1 = find(Z >= 0, 1, 'first');
Z2 = find(Z <= 0.05, 1, 'last');
T1 = 1;
T2 = length(T);

SPEED = squeeze(HDATA(R1:R2,S1:S2,Z1:Z2,T1:T2));
TIMES = T1:T2;

% SPEED is organized as r,s,z,t where each entry contains a count of the
% occurrences of speed value, s, at radius, r. To get the maximum
% speed, at each time step, find the largest speed bin that has a
% non-zero count across all radii and use the speed value for that bin.
%
% UNCOMMENT THE FOLLOWING IF Z2 > Z1
%SPEED = sum(SPEED, 3); % sum across z
SPEED = squeeze(sum(SPEED,1)); % sum across r

% SPEED is now organized as (s,t). For each time step, find the largest
% speed bin that has a non-zero count, and use that bin to look up the
% speed value associtated with it. This will be the maximum speed
% for that time step.
Nt = length(TIMES);
ST_MAX = zeros([1 Nt]);
for it = 1:Nt
  Sind = find(SPEED(:,it) > 0, 1, 'last');
  ST_MAX(it) = S(Sind);
end

% NHC Best Track (every six hours) data
% time step 1 from the simulation is where the NHC data starts
% NHC wind is in mph, multiply by 0.447 to convert to m/s
% Aug 22, 2006, 6Z through Aug 24, 2006, 18Z (11 points)
NHC_WIND  = [ 35   35   35   40   45   50   50   50   50   50   50   ] * 0.447;
NHC_TIMES = (1:6:61);

% plot
OutFile = sprintf('%s/TsDebbyWind.jpg', Pdir);

FigWind = figure;
set(gca, 'FontSize', 20);
SimST = plot(TIMES, ST_MAX, '-r', 'linewi', 2);
hold on;
NhcST = plot(NHC_TIMES, NHC_WIND, '+b', 'linewi', 2);

T = title('(c)');
set(T, 'Units', 'Normalized');
set(T, 'HorizontalAlignment', 'Left');
Tpos = get(T, 'Position');
Tpos(1) = 0; % line up with left edge of plot area
set(T, 'Position', Tpos);

xlabel('Time');
set(gca,'xtick', (13:24:61));
set(gca,'xticklabel', { 'Aug22:0Z', 'Aug23:0Z', 'Aug24:0Z' });
ylabel('Wind Speed (m/s)');
ylim([ 0 30 ]);
legend([ NhcST SimST ], 'NHC Best Track', 'Simulated Storm', 'Location', 'NorthEast');

fprintf('Writing: %s\n', OutFile);
saveas(FigWind, OutFile);
close(FigWind);
