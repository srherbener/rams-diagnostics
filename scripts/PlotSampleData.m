% Script to plot out sample data
%

clear;

% Sample data file has the CCN concentration values and SST values in two
% arrays. Read these in to form the retrieval of the data fields. Each data
% field is stored in:
%
%    /CCN_<ccn_concentration>/SST_<sst_value>/<variable>
%
% The CCN and SST arrays are 1D and corresponding entries show the
% combinations of CCN and SST values for each variable in the file.
%
% The 2D fields for the variables have x as the first dimension (rows) and
% y as the second dimension (columns). x corresponds to longitude values
% and y corresponds to latitude values.
%
% Create a multi-panel plot of one variable with all combinations of CCN
% and SST values.
h5_fin = 'DIAG/SampleData.h5';

CCN = hdf5read(h5_fin,'/CcnConcen');
SST = hdf5read(h5_fin,'/Sst');
Ns = length(CCN);

% Latitude and Longitude are 2D arrays that show the latitude values at all
% (x,y) locations in the variable fields, and the longitude values at all
% (x,y) locations. Note since rows are different longitude values and
% columns are different latitude values, then the rows of the 2D latitude
% each show the range of latitudes from the sim; and the columns of the 2D
% longitude each show the range of the longitudes from the sim. The domain
% is small in this sim so just use the first row/column to grab the
% latitude/longitude ranges.
Lat2D = hdf5read(h5_fin,'Lat');
Lon2D = hdf5read(h5_fin,'Lon');
Lat = Lat2D(1,:)';
Lon = Lon2D(:,1);
Nx = length(Lon);
Ny = length(Lat);


PCPRR = zeros(Ns,Nx,Ny);
CLOUDTOP_TEMPC = zeros(Ns,Nx,Ny);
VERTINT_COND = zeros(Ns,Nx,Ny);

% Read in variables
for i = 1:Ns
    C = CCN(i);
    S = SST(i);
    
    h5_dset = sprintf('/CCN_%d/SST_%d/PCPRR',C,S);
    fprintf('Reading: %s\n', h5_dset);
    % put the different longitudes in the columns and the different
    % latitudes in the rows (ie, save the transpose of the 2D field from
    % the hdf5 file dataset.
    TEMP = hdf5read(h5_fin, h5_dset);
    PCPRR(i,:,:) = TEMP';
    
    h5_dset = sprintf('/CCN_%d/SST_%d/CLOUDTOP_TEMPC',C,S);
    fprintf('Reading: %s\n', h5_dset);
    TEMP = hdf5read(h5_fin, h5_dset);
    CLOUDTOP_TEMPC(i,:,:) = TEMP';
    
    h5_dset = sprintf('/CCN_%d/SST_%d/VERTINT_COND',C,S);
    fprintf('Reading: %s\n', h5_dset);
    TEMP = hdf5read(h5_fin, h5_dset);
    VERTINT_COND(i,:,:) = TEMP';
end

% Plot the sample data - 4 panel plots each one showing one variable with
% all CCN, SST combinations.

PlotMultiPanelSample(PCPRR,Lon,Lat,CCN,SST,'DIAG/PCPRR.sample.jpg');
PlotMultiPanelSample(CLOUDTOP_TEMPC,Lon,Lat,CCN,SST,'DIAG/CLOUDTOP_TEMPC.sample.jpg');
PlotMultiPanelSample(VERTINT_COND,Lon,Lat,CCN,SST,'DIAG/VERTINT_COND.sample.jpg');