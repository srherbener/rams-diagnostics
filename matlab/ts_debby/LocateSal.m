function [ ] = LocateSal(InFile, InVname, RhLimit, LonRange, LatRange, Zrange, ColCount, OutFile)
% LocateSal function to determine location of SAL within the horizontal domain
%
% This function will look at RH in layers between 1 km and 5 km AGL, determine whether or not the
% dry layer characteristic of the SAL exists, and output 1's and 0's in an array to mark the location
% of the SAL (1 - SAL exists at this location, 0 - SAL does not exist at this location).
%
% RhLimit is the relative humidity value (in percent) that specifies what constitutes the SAL
% dry layer (RH < RhLimit).
%
% Zrange ([Zmin, Zmax]) mark the heights (in km) to look for the dry SAL layer, ie. only check RH from
% Zmin to Zmax heights.
%
% ColCount is the number of occurrences of RH < RhLimit required before that column is counted as
% part of the SAL.
%
% OutFile is the name of the output HDF5 file.

LonMin = LonRange(1);
LonMax = LonRange(2);

LatMin = LatRange(1);
LatMax = LatRange(2);

Zmin = Zrange(1);
Zmax = Zrange(2);

% RH will be organized as (x,y,z)
fprintf('Reading: %s\n', InFile);
fprintf('  Dataset: %s\n', InVname);
fprintf('\n');
RH = squeeze(hdf5read(InFile, InVname));
X = squeeze(hdf5read(InFile, '/x_coords'));      % deg lon
Y = squeeze(hdf5read(InFile, '/y_coords'));      % deg lat
Z = squeeze(hdf5read(InFile, '/z_coords'))/1000; % km AGL

Nx = length(X);
Ny = length(Y);

% We want the output to cover the entire lat/lon extent, but want
% to confine the 1's in the output to the box defined by LonRange
% and LatRange. Form a mask that has 1's inside the box and 0's outside.
LON_MASK = repmat((X >= LonMin & X <= LonMax), [ 1 Ny ]);
LAT_MASK = repmat((Y >= LatMin & Y <= LatMax)', [ Nx 1 ]);
AREA_MASK = LON_MASK & LAT_MASK;

% Create another mask (1's and 0's) that shows where the SAL is based on anywhere
% in the column between 1km and 5km AGL that has RH < RhLimit.
Z1 = find(Z >= Zmin, 1, 'first');
Z2 = find(Z <= Zmax, 1, 'last');

MASK = RH(:,:,Z1:Z2) < RhLimit;  % 0's where RH >= RhLimit, 1's where RH < RhLimit
MASK = squeeze(sum(MASK,3));     % 0 where RH >= RhLimit, n where RH < RhLimit n times in the column

MASK = double((MASK >= ColCount) & AREA_MASK);  % have to have at least ColCount occurrences of RH < RhLimit in the column

% Write out the mask
fprintf('Writing: %s\n', OutFile);
hdf5write(OutFile, 'SAL_LOC', MASK);
hdf5write(OutFile, '/x_coords', X, 'WriteMode', 'append');
hdf5write(OutFile, '/y_coords', Y, 'WriteMode', 'append');
fprintf('\n');

% graphic check
contourf(X, Y, MASK');
colorbar;

