% Script to plot out sample data
%

clear;

DiagDir = 'DIAGS';
PlotDir = 'plots';

% Keep these in sync with ExtractSample.m
CCN =  [ 50     50 1600   1600 ];
GCCN = [  1 0.0001    1 0.0001 ];

CCN_NAMES  = { 'ccn0050' 'ccn0050' 'ccn1600' 'ccn1600' };
GCCN_NAMES = { 'gcn10m0' 'gcn10m4' 'gcn10m0' 'gcn10m4' };

% Sample data file has the CCN concentration values and GCCN values in two
% arrays. Read these in to form the retrieval of the data fields. Each data
% field is stored in:
%
%    /<ccn_concentration>/<gccn_concentration>/<variable>
%
% The CCN and GCCN arrays are 1D and corresponding entries show the
% combinations of CCN and SST values for each variable in the file.
%
% The 2D fields for the variables have x as the first dimension (rows) and
% y as the second dimension (columns). x corresponds to longitude values
% and y corresponds to latitude values.
%
% Create a multi-panel plot of one variable with all combinations of CCN
% and GCCN values.
h5_fin = sprintf('%s/SampleData.h5', DiagDir);

Tstep = hdf5read(h5_fin,'/Tstep');

Ns = length(CCN_NAMES);

% Latitude and Longitude are 1D arrays that show the latitude values at all
% y locations in the variable fields, and the longitude values at all
% x locations.
%
Lon = hdf5read(h5_fin,'Lon');
Lat = hdf5read(h5_fin,'Lat');
Nx = length(Lon);
Ny = length(Lat);

SimTime = (Tstep - 1) / 12; % 5 min per time step, so divide by 12 to get hours
StimeStr = sprintf('Simulation Time = %d hrs', SimTime);

PCPRR = zeros(Ns,Nx,Ny);
CLOUDTOP_TEMPC = zeros(Ns,Nx,Ny);
VERTINT_COND = zeros(Ns,Nx,Ny);

% Read in variables
for i = 1:Ns
    C = CCN_NAMES{i};
    G = GCCN_NAMES{i};
    
    h5_dset = sprintf('/%s/%s/pcprr',C,G);
    fprintf('Reading: %s\n', h5_dset);
    % put the different longitudes in the columns and the different
    % latitudes in the rows (ie, save the transpose of the 2D field from
    % the hdf5 file dataset.
    TEMP = hdf5read(h5_fin, h5_dset);
    PCPRR(i,:,:) = TEMP';
    
    h5_dset = sprintf('/%s/%s/ctop_tempc',C,G);
    fprintf('Reading: %s\n', h5_dset);
    TEMP = hdf5read(h5_fin, h5_dset);
    CLOUDTOP_TEMPC(i,:,:) = TEMP';
    
    h5_dset = sprintf('/%s/%s/vint_cond',C,G);
    fprintf('Reading: %s\n', h5_dset);
    TEMP = hdf5read(h5_fin, h5_dset);
    VERTINT_COND(i,:,:) = TEMP';
end

% Plot the sample data - 4 panel plots each one showing one variable with
% all CCN, SST combinations.

% color map and contour levels for PCPRR
CmapPCPRR = colormap('jet');
CmapPCPRR(1,:) = [ 1 1 1 ];  % change the min value to white
ClevsPCPRR = (0:0.05:2);
%ClevsPCPRR = [ (0:0.2:1.0) (1.3:0.3:2.8) 3.2 4.4 5 6 7 8 9 10 ];
%ClevsPCPRR = [ 0 0.001 0.01 0.1 1.0 10 ];
PtitlePCPRR = sprintf('Precipitation Rate (mm/hr)\n%s', StimeStr);

% color map and contour levels for VERTINT_COND
CmapVERTINT_COND = colormap('jet');
CmapVERTINT_COND(1,:) = [ 1 1 1 ];  % change the min value to white
ClevsVERTINT_COND = (0:0.05:2);
%ClevsVERTINT_COND = [ (0:0.2:1.0) (1.3:0.3:2.8) (3.2:0.4:4.8) 5.4 5.8 6.2 ];
%ClevsVERTINT_COND = [ 0 0.001 0.01 0.1 1 10 ];
PtitleVERTINT_COND = sprintf('Vertically Integrated Condensate (mm)\n%s', StimeStr);

% color map and contour levels for CLOUDTOP_TEMPC
CmapCLOUDTOP_TEMPC = flipud(colormap('gray'));  % want light colors (taller clouds) with colder values
%CmapCLOUDTOP_TEMPC(1,:) = [ 1 1 1 ];  % change the min value to white
%ClevsCLOUDTOP_TEMPC = [ (0:1:10) (12:2:20) ];
ClevsCLOUDTOP_TEMPC = (0:0.5:20);
PtitleCLOUDTOP_TEMPC = sprintf('Cloud Top Temperature (degrees C)\n%s', StimeStr);

close; % issuing the colomap command opens a figure


OutFile = sprintf('%s/pcprr.sample.fig', PlotDir);
PlotMultiPanelSample(PCPRR,Lon,Lat,CCN,GCCN, CmapPCPRR, ClevsPCPRR, PtitlePCPRR, OutFile);
OutFile = sprintf('%s/vint_cloud.sample.fig', PlotDir);
PlotMultiPanelSample(VERTINT_COND,Lon,Lat,CCN,GCCN, CmapVERTINT_COND, ClevsVERTINT_COND, PtitleVERTINT_COND, OutFile);

% REVU placed a very cold temperature (-248.31 deg C) in where there is missing data (no clouds?)
% Change these to extra warm to simulate no cloud regions.
CLOUDTOP_TEMPC(CLOUDTOP_TEMPC <= 0) = 30;
OutFile = sprintf('%s/ctop_tempc.sample.fig', PlotDir);
PlotMultiPanelSample(CLOUDTOP_TEMPC,Lon,Lat,CCN,GCCN, CmapCLOUDTOP_TEMPC, ClevsCLOUDTOP_TEMPC, PtitleCLOUDTOP_TEMPC, OutFile);
