% Script to generate sets of plots from SW data
%
%   1) Time series of top of model domain average values
%   2) Time, domain average of top of model - bar graph

clear;

SST = [ 293 298 303 ];
LegLocs = { 'NorthWest' 'NorthEast' 'NorthEast' };

%CCN = [ 50 100 200 400 800 1200 1600 ]; % too cluttered
%Lcolors = { 'k', 'm', 'b', 'c', 'g', 'y', 'r' };
CCN = [ 100 400 800 1600 ];
Lcolors = { 'b', 'g', 'y', 'r' };

CNTL_SST = 298;
CNTL_CCN = 50;

Ns = length(SST);
Nc = length(CCN);

YlimTS = [ -200 200 ];
YlabTS = 'Radiative Flux (W/m^2)';

Times = (0:5/60:36); % 36 hrs in 5 min increments

% Create the output directory if it doesn't exist
OutDir = 'PLOTS';
if (exist(OutDir,'dir') == 0)
  mkdir(OutDir);
end

% read in the control profile
Hfile = sprintf('DIAGS/swup_tdavg_ATEX_C%04d_S%03d.h5',CNTL_CCN,CNTL_SST);
fprintf('Reading control HDF5 file: %s\n\n', Hfile);
CNTL_SWUP_DOMAVG = hdf5read(Hfile, '/swup');
Z = size(CNTL_SWUP_DOMAVG,1) - 1;  % top level is all zeros (boundary)
CNTL_SWUP = squeeze(CNTL_SWUP_DOMAVG(Z,:));

% Each set: all CCN levels for a single given SST
for i = 1:Ns
  fprintf('Generating plot for SST = %d\n', SST(i));
  Ptitle = sprintf('SW Forcing Difference: Top of Model, SST = %d', SST(i));

  % collect the latent heat data
  for j = 1:Nc
    Hfile = sprintf('DIAGS/swup_tdavg_ATEX_C%04d_S%03d.h5',CCN(j),SST(i));
    fprintf('  Reading HDF5 file: %s\n', Hfile);
    SWUP_DOMAVG = hdf5read(Hfile, '/swup');
    SWUP(j,:) = squeeze(SWUP_DOMAVG(Z,:)); % time series

    % subtract off the control series
    SWUP_DIFF(j,:) = SWUP(j,:) - CNTL_SWUP;

    % generate the legend text
    LegText(j) = { sprintf('CCN: %d/cc', CCN(j)) };
  end
  fprintf('\n');

  % plot it
  OutFile = sprintf('%s/swup_S%03d.jpg', OutDir, SST(i));
  fprintf('  Saving plot in file: %s\n\n', OutFile);
  PlotTseriesSet2( Times, SWUP_DIFF, Ptitle, YlimTS, YlabTS, Lcolors, LegText, LegLocs(i), OutFile )
end
