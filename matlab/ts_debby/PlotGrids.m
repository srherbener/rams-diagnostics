function [ ] = PlotGrids()
% PlotGrids function to plot the grids used in TS Debby simulations

Pdir = 'Plots';
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

Fsize = 22;
LegendFsize = 15;

% for the TS Debby simulations:
LatBounds = [ -40 60 ];
LonBounds = [ -120 60 ];

% go around grid: SW, NW, NE, SE, SW
% repeat SW so that m_line will draw a closed figure

% RAMS reported numbers
%G1_lats = [ -12.0 30.1 30.1 -12.0 -12.0 ];
%G1_lons = [ -70.9 -80.4 36.4 26.9 -70.9 ];

%G2_lats = [ 6.6 23.0 23.7 7.0 6.6 ];
%G2_lons = [ -39.5 -40.9 -12.9 -13.6 -39.5  ];

%G3_lats = [ 6.9 22.7 23.4 7.3 6.9 ];
%G3_lons = [ -39.0 -40.4 -13.4 -14.0 -39.0 ];

G1_lats = [ -12 30 30 -12 -12 ];
G1_lons = [ -75 -75 30 30 -75 ];

G2_lats = [ 6 24 24 6 6 ];
G2_lons = [ -40 -40 -13 -13 -40 ];

G3_lats = [ 7 23 23 7 7 ];
G3_lons = [ -39 -39 -14 -14 -39 ];


% plot
OutFile = sprintf('%s/Fig2_ModelGrids.jpg', Pdir);

FigTracks = figure;
set(gca, 'FontSize', Fsize);
m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
m_coast('color', 'k'); % k --> black
m_grid('linestyle','none','box','fancy','tickdir','out');
Grid1 = m_line(G1_lons, G1_lats, 'linewi', 3, 'color', 'b');
Grid2 = m_line(G2_lons, G2_lats, 'linewi', 3, 'color', 'r');
Grid3 = m_line(G3_lons, G3_lats, 'linewi', 3, 'color', 'g');

% %title('TS Debby Simulation Grids');
% T = title('(a)');
% set(T, 'Units', 'Normalized');
% set(T, 'HorizontalAlignment', 'Left');
% Tpos = get(T, 'Position');
% Tpos(1) = 0; % line up with left edge of plot area
% set(T, 'Position', Tpos);

legend([ Grid1 Grid2 Grid3 ], { 'Grid1' 'Grid2' 'Grid3' }, 'Location', 'NorthWest', 'FontSize', LegendFsize);

fprintf('Writing: %s\n', OutFile);
saveas(FigTracks, OutFile);
close(FigTracks);
