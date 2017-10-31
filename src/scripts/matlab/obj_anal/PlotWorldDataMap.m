function [ ] = PlotWorldDataMap( Figure, Wdata, Lat, Lon, Ptitle )
%PlotWorldDataMap Create a 2D contour plot superimposed over a world map
%
%   This function will take the row vector Wdata, which holds gridded data
%   points located across the world, and plot these data superimposed over
%   the continental outlines.
%
%   Figure is the handle to the figure in which the plot is to be drawn.
%
%   Lat, and Lon are vectors holding the latitude and longitude points. It
%   is up to the caller to define these properly so that the points in
%   Wdata line up with the latitude, longitude points in Lat and Lon. The
%   2D map will be created with Lat corresponding to rows and Lon
%   corresponding to columns.

% Using the package "m_map" for drawing the continental outlines, etc.

% Convert the row vector Wdata into a 2D array with rows corresponding to
% latitudes and columns corresponding to longitudes.

NumLat = length(Lat);
NumLon = length(Lon);
%Dmap = zeros(NumLat,NumLon);
%for ic = 1:length(Wdata)
%    % Figure out the 2D indices for the maps from the 1D index
%    % for the world data
%    i = floor((ic-1)/NumLon) + 1;
%    j = mod((ic-1),NumLon) + 1;
%    
%    Dmap(i,j) = Wdata(ic);
%end
% reshape fills in columns first, and we want to fill in rows first
% so use reshape to create the transpose of the result we want.
Dmap = reshape(Wdata,NumLon,NumLat)';

% Replace NaNs with zeros to get the "shading flat" mode working.
Dmap(isnan(Dmap)) = 0;

% Grab the min and max values of Lat and Lon in order to set the axis on
% the world map (m_proj call)
LatBounds = [ min(Lat) max(Lat) ];
LonBounds = [ min(Lon) max(Lon) ];

% Find the largest absolute value of the entries in Dmap for setting the
% colormap axis. Want to center this about zero so that blue represents
% negative values and red represents positive values.
Clim = max(max(abs(Dmap))) * 1.1; % need the range to be slightly
                                  % larger than the actual data
Cbounds = [ -Clim Clim ];

% plot the 2D map
figure(Figure);
m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
m_contourf(Lon,Lat,Dmap,100);
shading flat;
m_coast('color', 'k'); % k --> black
m_grid('linestyle','none','box','fancy','tickdir','out');
title(Ptitle);
caxis(Cbounds);
colormap(redblue);
colorbar;

end

