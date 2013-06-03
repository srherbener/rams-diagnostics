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

Hfile = 'FILTERS/all_TSD_3GRIDS.h5';
HdsetLon = 'min_press_xloc';
HdsetLat = 'min_press_yloc';

fprintf('Reading: %s\n', Hfile);
fprintf('  Track longitude: %s\n', HdsetLon);
fprintf('  Track latitude: %s\n', HdsetLat);
SimTrackLons = squeeze(hdf5read(Hfile, HdsetLon));
SimTrackLats = squeeze(hdf5read(Hfile, HdsetLat));

% plot
FigTracks = figure;
set(gca, 'FontSize', 18);
m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
m_coast('color', 'k'); % k --> black
m_grid('linestyle','none','box','fancy','tickdir','out');
NhcTrack = m_line(NhcTrackLons, NhcTrackLats, 'linewi', 3, 'color', 'b', 'linestyle', 'none', 'marker', '+');
SimTrack = m_line(SimTrackLons, SimTrackLats, 'linewi', 3, 'color', 'r'); 
title('TS Debby Tracks');
legend([ NhcTrack SimTrack ], 'NHC Best Track', 'Simulated Track', 'Location', 'NorthWest');

OutFile = sprintf('%s/TsDebbyTracks.jpg', Pdir);
fprintf('Writing: %s\n', OutFile);
saveas(FigTracks, OutFile);
close(FigTracks);
