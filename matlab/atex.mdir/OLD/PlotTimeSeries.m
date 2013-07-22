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

% Read in the data, place the time series in the columns of arrays
for i = 1:length(Exps)
    h5_file = sprintf('DIAG/%s/TimeSeries.h5', char(Exps(i)));
    fprintf('Reading time series file: %s\n', h5_file);
    
    AllAvgLwp(:,i) = hdf5read(h5_file,'/AvgLWP');
    AllAvgR2C(:,i) = hdf5read(h5_file,'/Avg_RainToCloud');
    AllRavg2Cavg(:,i) = hdf5read(h5_file,'RainAvgToCloudAvg');
    AllAvgCloud(:,i) = hdf5read(h5_file,'/AvgCloud');
    AllAvgRain(:,i) = hdf5read(h5_file,'/AvgRain');
end

% create sets of data to show on a single plot
% Two plots: one for SST = 298, CCN = 50, 1600
%            other for SST = 303, CCN = 50, 1600
Set298 = [ 1 7 ];
Set303 = [ 2 8 ];

AvgLwp_298 = AllAvgLwp(:,Set298);
AvgR2C_298 = AllAvgR2C(:,Set298);
Ravg2Cavg_298 = AllRavg2Cavg(:,Set298);
AvgCloud_298 = AllAvgCloud(:,Set298);
AvgRain_298 = AllAvgRain(:,Set298);
C298 = CCN(Set298);
S298 = SST(Set298);

AvgLwp_303 = AllAvgLwp(:,Set303);
AvgR2C_303 = AllAvgR2C(:,Set303);
Ravg2Cavg_303 = AllRavg2Cavg(:,Set303);
AvgCloud_303 = AllAvgCloud(:,Set303);
AvgRain_303 = AllAvgRain(:,Set303);
C303 = CCN(Set303);
S303 = SST(Set303);

%Lstyles = { '-k', '--b', ':r', '-.g' };
Lstyles = { '-k', '-.g' };

PlotTseriesSet(AvgLwp_298, C298, S298, 'Time Series: Average LWP, SST: 298K', [ 0 0.25 ], 'LWP (mm)', Lstyles, 'NorthEast', 'DIAG/TsAvgLwp_298K.jpg');
PlotTseriesSet(AvgLwp_303, C303, S303, 'Time Series: Average LWP, SST: 303K', [ 0 0.25 ], 'LWP (mm)', Lstyles, 'NorthEast', 'DIAG/TsAvgLwp_303K.jpg');

PlotTseriesSet(AvgR2C_298, C298, S298, 'Time Series: Average Rain/Cloud, SST: 298K', [ 0 2e5 ], 'Ratio: Rain/Cloud', Lstyles, 'NorthEast', 'DIAG/TsAvgR2C_298K.jpg');
PlotTseriesSet(AvgR2C_303, C303, S303, 'Time Series: Average Rain/Cloud, SST: 303K', [ 0 2e5 ], 'Ratio: Rain/Cloud', Lstyles, 'NorthEast', 'DIAG/TsAvgR2C_303K.jpg');

PlotTseriesSet(Ravg2Cavg_298, C298, S298, 'Time Series: Average Rain/Cloud, SST: 298K', [ 0 1 ], 'Ratio: Rain/Cloud', Lstyles, 'NorthEast', 'DIAG/TsRavg2Cavg_298K.jpg');
PlotTseriesSet(Ravg2Cavg_303, C303, S303, 'Time Series: Average Rain/Cloud, SST: 303K', [ 0 1 ], 'Ratio: Rain/Cloud', Lstyles, 'NorthEast', 'DIAG/TsRavg2Cavg_303K.jpg');

PlotTseriesSet(AvgCloud_298, C298, S298, 'Time Series: Average CWP, SST: 298K', [ 0 0.25 ], 'CWP (mm)', Lstyles, 'NorthEast', 'DIAG/TsAvgCloud_298K.jpg');
PlotTseriesSet(AvgCloud_303, C303, S303, 'Time Series: Average CWP, SST: 303K', [ 0 0.25 ], 'CWP (mm)', Lstyles, 'NorthEast', 'DIAG/TsAvgCloud_303K.jpg');

PlotTseriesSet(AvgRain_298, C298, S298, 'Time Series: Average RWP, SST: 298K', [ 0 0.25 ], 'RWP (mm)', Lstyles, 'NorthEast', 'DIAG/TsAvgRain_298K.jpg');
PlotTseriesSet(AvgRain_303, C303, S303, 'Time Series: Average RWP, SST: 303K', [ 0 0.25 ], 'RWP (mm)', Lstyles, 'NorthEast', 'DIAG/TsAvgRain_303K.jpg');


%%% % creates sets of data to show on a single plot
%%% % One plot: combinations of SST = 298,303 and CCN = 50, 1600
%%% Dset = [ 1 7 2 8 ];
%%% 
%%% AvgLwp = AllAvgLwp(:,Dset);
%%% AvgR2C = AllAvgR2C(:,Dset);
%%% Ravg2Cavg = AllRavg2Cavg(:,Dset);
%%% C = CCN(Dset);
%%% S = SST(Dset);
%%% 
%%% Lstyles = { '-k', '--k', '-g', '--g' };
%%% 
%%% PlotTseriesSet(AvgLwp, C, S, 'Time Series: Average LWP', [ 0 0.25 ], 'LWP (mm)', Lstyles, 'NorthEast', 'DIAG/TsAvgLwp.jpg');
%%% 
%%% PlotTseriesSet(AvgR2C, C, S, 'Time Series: Average Rain/Cloud', [ 0 2e5 ], 'Ratio: Rain/Cloud', Lstyles, 'NorthWest', 'DIAG/TsAvgR2C.jpg');
%%% 
%%% PlotTseriesSet(Ravg2Cavg, C, S, 'Time Series: Average Rain/Cloud', [ 0 1 ], 'Ratio: Rain/Cloud', Lstyles, 'NorthWest', 'DIAG/TsRavg2Cavg.jpg');
%%% 
