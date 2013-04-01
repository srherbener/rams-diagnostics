function [ ] = PlotTracks(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

% for the TS Debby simulations:
LatBounds = [ 5 26 ];
LonBounds = [ -42 -8 ];

NhcTrackLats = [ 11.6 12.0 12.6 13.4 14.2 14.9 15.7 16.7 17.6 18.4 19.2 20.1 20.9 ];
NhcTrackLons = [ 21.7 22.7 23.9 25.3 26.7 28.1 29.5 30.0 32.4 33.9 35.5 37.1 38.7 ] * -1; 

SimTrackLats = [ 12.5 13.0 13.5 14.0 14.5 15.5 16.0 16.5 17.5 18.0 19.0 19.0 19.5 ];
SimTrackLons = [ 20.0 21.5 23.0 24.0 25.0 26.5 27.5 29.0 30.5 31.5 33.0 34.0 36.0 ] * -1;

% plot
FigTracks = figure;
set(gca, 'FontSize', 18);
m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
m_coast('color', 'k'); % k --> black
m_grid('linestyle','none','box','fancy','tickdir','out');
NhcTrack = m_line(NhcTrackLons, NhcTrackLats, 'linewi', 3, 'color', 'b');
SimTrack = m_line(SimTrackLons, SimTrackLats, 'linewi', 3, 'color', 'r');
title('TS Debby Tracks');
legend([ NhcTrack SimTrack ], 'NHC Best Track', 'Simulated Track', 'Location', 'NorthWest');

OutFile = sprintf('%s/TsDebbyTracks.jpg', Pdir);
fprintf('Writing: %s\n', OutFile);
saveas(FigTracks, OutFile);
close(FigTracks);
