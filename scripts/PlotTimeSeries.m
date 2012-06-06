% Script to plot time series
%

%Exps = { 'z.atex250m.100km.ccn0050.sst298'
%'z.atex250m.100km.ccn0050.sst303'
%'z.atex250m.100km.ccn0100.sst298'
%'z.atex250m.100km.ccn0100.sst303'
%'z.atex250m.100km.ccn0200.sst298'
%'z.atex250m.100km.ccn0200.sst303'
%'z.atex250m.100km.ccn0400.sst298'
%'z.atex250m.100km.ccn0400.sst303'
%'z.atex250m.100km.ccn0800.sst298'
%'z.atex250m.100km.ccn0800.sst303'
%'z.atex250m.100km.ccn1200.sst298'
%'z.atex250m.100km.ccn1200.sst303'
%'z.atex250m.100km.ccn1600.sst298'
%'z.atex250m.100km.ccn1600.sst303' };

Exps = { 'z.atex250m.100km.ccn0050.sst298'
'z.atex250m.100km.ccn0050.sst303'
'z.atex250m.100km.ccn0200.sst298'
'z.atex250m.100km.ccn0200.sst303'
'z.atex250m.100km.ccn0800.sst298'
'z.atex250m.100km.ccn0800.sst303'
'z.atex250m.100km.ccn1600.sst298'
'z.atex250m.100km.ccn1600.sst303' };

CCN = [ 50 50 200 200 800 800 1600 1600 ];
SST = [ 298 303 298 303 298 303 298 303 ];

%Set298 = [ 1 3 5 7 ];
%Set303 = [ 2 4 6 8 ];
Set298 = [ 1 7 ];
Set303 = [ 2 8 ];

% Read in the data, place the time series in the columns of arrays
for i = 1:length(Exps)
    h5_file = sprintf('DIAG/%s/TimeSeries.h5', char(Exps(i)));
    fprintf('Reading time series file: %s\n', h5_file);
    
    AllAvgLwp(:,i) = hdf5read(h5_file,'/AvgLWP');
    AllAvgR2C(:,i) = hdf5read(h5_file,'/Avg_RainToCloud');
    AllRavg2Cavg(:,i) = hdf5read(h5_file,'RainAvgToCloudAvg');
end

% creates sets of data to show on a single plot
AvgLwp_298 = AllAvgLwp(:,Set298);
AvgR2C_298 = AllAvgR2C(:,Set298);
Ravg2Cavg_298 = AllRavg2Cavg(:,Set298);
C298 = CCN(Set298);

AvgLwp_303 = AllAvgLwp(:,Set303);
AvgR2C_303 = AllAvgR2C(:,Set303);
Ravg2Cavg_303 = AllRavg2Cavg(:,Set303);
C303 = CCN(Set303);

%Lstyles = { '-k', '--b', ':r', '-.g' };
Lstyles = { '-k', '-.g' };

PlotTseriesSet(AvgLwp_298, C298, 'Time Series: Average LWP, SST: 298K', 'LWP (mm)', Lstyles, 'South', 'DIAG/TsAvgLwp_298K.jpg');
PlotTseriesSet(AvgLwp_303, C303, 'Time Series: Average LWP, SST: 303K', 'LWP (mm)', Lstyles, 'NorthEast', 'DIAG/TsAvgLwp_303K.jpg');

PlotTseriesSet(AvgR2C_298, C298, 'Time Series: Average Rain/Cloud, SST: 298K', 'Ratio: Rain/Cloud', Lstyles, 'NorthWest', 'DIAG/TsAvgR2C_298K.jpg');
PlotTseriesSet(AvgR2C_303, C303, 'Time Series: Average Rain/Cloud, SST: 303K', 'Ratio: Rain/Cloud', Lstyles, 'NorthWest', 'DIAG/TsAvgR2C_303K.jpg');

PlotTseriesSet(Ravg2Cavg_298, C298, 'Time Series: Average Rain/Cloud, SST: 298K', 'Ratio: Rain/Cloud', Lstyles, 'NorthEast', 'DIAG/TsRavg2Cavg_298K.jpg');
PlotTseriesSet(Ravg2Cavg_303, C303, 'Time Series: Average Rain/Cloud, SST: 303K', 'Ratio: Rain/Cloud', Lstyles, 'NorthWest', 'DIAG/TsRavg2Cavg_303K.jpg');

