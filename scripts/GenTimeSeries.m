% Script to create time series plots for two quantities
%   1. domain averaged LWP (VERTINT_COND works since no ice in these sims)
%   2. ratio of vertically integrated rain water to vertically integrated cloud
%      water.
%

clear;

% Read in the PCPRR data
%
% After reading in, the dimensions will be (x,y,t)

Exps = { 'z.atex250m.100km.ccn0050.sst298',
'z.atex250m.100km.ccn0050.sst303',
'z.atex250m.100km.ccn0100.sst298',
'z.atex250m.100km.ccn0100.sst303',
'z.atex250m.100km.ccn0200.sst298',
'z.atex250m.100km.ccn0200.sst303',
'z.atex250m.100km.ccn0400.sst298',
'z.atex250m.100km.ccn0400.sst303',
'z.atex250m.100km.ccn0800.sst298',
'z.atex250m.100km.ccn0800.sst303',
'z.atex250m.100km.ccn1200.sst298',
'z.atex250m.100km.ccn1200.sst303',
'z.atex250m.100km.ccn1600.sst298',
'z.atex250m.100km.ccn1600.sst303' };

% LWP --> VERTINT_COND
% Rain --> VERTINT_RAIN
% Cloud --> VERTINT_CLOUD
for i = 1:size(Exps,1)
  h5_fin_cond = sprintf('REVU/%s/VERTINT_COND.h5',char(Exps(i)));
  h5_fin_rain = sprintf('REVU/%s/VERTINT_RAIN.h5',char(Exps(i)));
  h5_fin_cloud = sprintf('REVU/%s/VERTINT_CLOUD.h5',char(Exps(i)));
  h5_fout = sprintf('DIAG/%s/TimeSeries.h5',char(Exps(i)));
  fprintf('Generating Time Series:\n');
  fprintf('  Input LWP file: %s\n',h5_fin_cond);
  fprintf('  Input rain file: %s\n',h5_fin_rain);
  fprintf('  Input cloud file: %s\n',h5_fin_cloud);
  fprintf('  Output file: %s\n', h5_fout);
  fprintf('\n');

  % Read in data
  LWP   = hdf5read(h5_fin_cond, '/VERTINT_COND');
  RAIN  = hdf5read(h5_fin_rain, '/VERTINT_RAIN');
  CLOUD = hdf5read(h5_fin_cloud, '/VERTINT_CLOUD');

  % LWP time series is simply average of domain at each time step
  [ AVG_LWP ] = GenTseriesAvg2d(LWP);

  % ratio of rain to cloud can be computed two ways:
  %   1. An average of the ratios of each element in the 2D fields
  %   2. The ratio of the domain averages of the 2D fields
  %
  % Do both and save them in the output file

  AVG_R2C = GenTseriesAvg2d(RAIN ./ CLOUD); % average of ratios
  RAVG2CAVG = GenTseriesAvg2d(RAIN) ./ GenTseriesAvg2d(CLOUD); % ratio of averages

  % Save the histograms
  hdf5write(h5_fout, '/AvgLWP', AVG_LWP);
  hdf5write(h5_fout, '/Avg_RainToCloud', AVG_R2C, 'WriteMode', 'append');
  hdf5write(h5_fout, '/RainAvgToCloudAvg', RAVG2CAVG, 'WriteMode', 'append');
end
