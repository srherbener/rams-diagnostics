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

YlimTS = [ -150 150 ];
% For no diff plots: YlimTS = [ 0 300 ];
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
%  Below: i is the index into SST values
%         j is the index into CCN values
%
% For the time, domain averages place these numbers in a matrix with different
% SST values on the rows and different CCN values on the columns. This will facilitate
% drawing a bar graph with the different CCN values for a given SST grouped together.
k = 0;
for i = 1:Ns
  fprintf('Generating plot for SST = %d\n', SST(i));
  Ptitle = sprintf('Upward SW Flux Difference, Top of Model, SST: %dK', SST(i));

  % collect the latent heat data
  for j = 1:Nc
    Hfile = sprintf('DIAGS/swup_tdavg_ATEX_C%04d_S%03d.h5',CCN(j),SST(i));
    fprintf('  Reading HDF5 file: %s\n', Hfile);
    SWUP_DOMAVG = hdf5read(Hfile, '/swup');
    SWUP(j,:) = squeeze(SWUP_DOMAVG(Z,:)); % time series using avarages on level Z

    % subtract off the control series
    SWUP_DIFF(j,:) = SWUP(j,:) - CNTL_SWUP;

    % generate the legend text
    LegText(j) = { sprintf('CCN: %d/cc', CCN(j)) };
  end
  fprintf('\n');

  % Compute the time average of the data in SWUP
  %   SWUP: rows are CCN value, columns are time so we just want to take the
  %   average of each row.
  SWUP_TDAVG(i,:) = mean(SWUP,2);
  Xtext(i) = { sprintf('SST: %dK', SST(i)) };

  % plot it
  OutFile = sprintf('%s/swup_S%03d.jpg', OutDir, SST(i));
  fprintf('  Saving plot in file: %s\n\n', OutFile);
  PlotTseriesSet( Times, SWUP_DIFF, Ptitle, YlimTS, YlabTS, Lcolors, LegText, LegLocs(i), OutFile );
end

% Plot a bar graph showing the relationship of the time,domain averaged top-of-model SW
% Note that LegText is already set to what we need from the time series plots above
fprintf('Generating bar chart for time domain averaged SW\n');
Fig = figure;

bar(SWUP_TDAVG);
set(gca, 'FontSize', 20);
title('Upward SW Flux: Top of Model');
set(gca, 'XtickLabel', Xtext);
ylabel(YlabTS);
ylim([ 0 150 ]);
legend(LegText);
legend boxoff;

Pfile = 'PLOTS/swup_tdavg.jpg';
fprintf('  Saving bar chart in: %s\n', Pfile);
saveas(Fig,Pfile);
close(Fig);
fprintf('\n');
