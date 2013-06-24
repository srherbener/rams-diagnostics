% Script to plot out the extents of the grids used in the RAMS simulations
%
% Use the package "m_map" for drawing the continental outlines, etc.

clear;

Pdir = 'Plots';
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

OutFile = sprintf('%s/AerosolLoc.jpg', Pdir);

Fsize = 35;
CenterFsize = 25;
Lwidth = 2;

% Make the outline be the grid 3 boundaries
Grid3Lon = [ 317.3 317.3 322.8 322.8 317.3 ];
Grid3Lat = [  12.3  17.7  17.7  12.3  12.3 ];

G3_LonRange = [ min(Grid3Lon) max(Grid3Lon) ];
G3_LatRange = [ min(Grid3Lat) max(Grid3Lat) ];

G3_LabelX = G3_LonRange(1) + 0.2;
G3_LabelY = G3_LatRange(1) + 0.1;

G3_CenterX = sum(G3_LonRange) / 2;
G3_CenterY = sum(G3_LatRange) / 2;

% Don't have Lon/Lat for aerosol location, but do have i,j index values.
% Both x and y have 302 grid points which means 301 intervals.
AeroLon = G3_LonRange(1) + (([  10  10 290 290  10 ] ./ 301) .* (G3_LonRange(2) - G3_LonRange(1)));
AeroLat = G3_LatRange(1) + (([ 225 230 230 225 225 ] ./ 301) .* (G3_LatRange(2) - G3_LatRange(1)));
    
Fig = figure;

% Back ground
m_proj('miller', 'lat', G3_LatRange, 'long', G3_LonRange);
m_grid('xtick', [ 318 320 322 ], 'ytick', [ 13 15 17 ], 'linestyle','none','box','fancy','tickdir','out', 'fontsize', Fsize, 'linewidth', Lwidth);

% Grid boxes
m_line(Grid3Lon, Grid3Lat, 'linewidth', Lwidth, 'color', 'k', 'linestyle', '-');

% Grid labels
m_text(G3_LabelX, G3_LabelY, 'G3', 'fontsize', Fsize, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Storm center
m_text(G3_CenterX, G3_CenterY, { 'Storm' '+' 'Center' }, 'fontsize', CenterFsize, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

% Aerosol source location (plan view)
m_line(AeroLon, AeroLat, 'linewidth', Lwidth, 'color', 'k', 'linestyle', '-');

% The title is in a box that adjusts to the amount of characters in
% the title. Ie, it doesn't do any good to do Left/Center/Right
% alignment. But, the entire box can be moved to the left side of the
% plot.
T = title('(b)', 'FontSize', Fsize);
set(T, 'Units', 'Normalized');
set(T, 'HorizontalAlignment', 'Left');
Tpos = get(T, 'Position');
Tpos(1) = 0; % line up with left edge of plot area
set(T, 'Position', Tpos);

% fix position
Ppos = get(gca, 'Position');
Ppos(1) = Ppos(1) * 1.20;
Ppos(2) = Ppos(2) * 1.20;
Ppos(3) = Ppos(3) * 0.90;
Ppos(4) = Ppos(4) * 0.90;
set(gca, 'Position', Ppos);

fprintf ('Writing: %s\n', OutFile);
saveas (Fig, OutFile);
close(Fig);
