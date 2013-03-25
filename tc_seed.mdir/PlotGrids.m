% Script to plot out the extents of the grids used in the RAMS simulations
%
% Use the package "m_map" for drawing the continental outlines, etc.

clear;

Pdir = 'Plots';
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

OutFile = sprintf('%s/RamsGrids.jpg', Pdir);

Fsize = 25;
Lwidth = 2;

% Show enough to see some continents
Lat = (0:0.5:30);
Lon = (300:0.5:345);
LatRange = [ min(Lat) max(Lat) ];
LonRange = [ min(Lon) max(Lon) ];

% Grids
Grid1Lon = [ 311.5 311.5 328.5 328.5 311.5 ];
Grid1Lat = [   6.5  23.0  23.0   6.5   6.5 ];

Grid2Lon = [ 316.5 316.5 323.8 323.8 316.5 ];
Grid2Lat = [  11.5  18.6  18.6  11.5  11.5 ];
    
Grid3Lon = [ 317.3 317.3 322.8 322.8 317.3 ];
Grid3Lat = [  12.3  17.7  17.7  12.3  12.3 ];
    
Fig = figure;

% Back ground
m_proj('miller', 'lat', LatRange, 'long', LonRange);
m_coast('patch', [ 0.7 0.7 0.7 ] ); % gray
m_grid('linestyle','none','box','fancy','tickdir','out', 'fontsize', Fsize, 'linewidth', Lwidth);

% Grid boxes
m_line(Grid1Lon, Grid1Lat, 'linewidth', Lwidth, 'color', 'k', 'linestyle', '-.');
m_line(Grid2Lon, Grid2Lat, 'linewidth', Lwidth, 'color', 'k', 'linestyle', '--');
m_line(Grid3Lon, Grid3Lat, 'linewidth', Lwidth, 'color', 'k', 'linestyle', '-');;

% Grid labels
m_text(Grid1Lon(2), Grid1Lat(2), 'G1', 'fontsize', Fsize, 'VerticalAlignment', 'bottom');
m_text(Grid2Lon(2), Grid2Lat(2), 'G2', 'fontsize', Fsize, 'VerticalAlignment', 'bottom');
m_text(320, 15, 'G3', 'fontsize', Fsize, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

% Continent labels
m_text(304, 3, 'S.A.', 'fontsize', Fsize, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

fprintf ('Writing: %s\n', OutFile);
saveas (Fig, OutFile);
close(Fig);
