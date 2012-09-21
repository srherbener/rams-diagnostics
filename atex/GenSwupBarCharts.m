function [ ] = GenSwupBarCharts(ConfigFile)
% GenSwupBarCharts generate bar charts showing impact of SST and CCN concentration of TOA SW upward flux
%

Config = ReadConfig(ConfigFile);

ControlCase = Config.ControlCase;
Pname   = Config.Pexp.Ename;
Tdir = Config.TsavgDir;
Pdir = Config.PlotDir;

SST = [ 293 298 303 ];
CCN = [ 100 400 800 1600 ];

Ns = length(SST);
Nc = length(CCN);

% Create the output directory if it doesn't exist
if (exist(Pdir,'dir') == 0)
  mkdir(Pdir);
end

% Each set: all CCN levels for a single given SST
%  Below: i is the index into SST values
%         j is the index into CCN values
%
% For the time, domain averages place these numbers in a matrix with different
% SST values on the rows and different CCN values on the columns. This will facilitate
% drawing a bar graph with the different CCN values for a given SST grouped together.
k = 0;
for i = 1:Ns
  % collect the latent heat data
  for j = 1:Nc
    Hfile = sprintf('%s/hda_swup_z.atex.ccn%04d.sst%d.gcn10m5.h5', Tdir, CCN(j),SST(i));
    fprintf('  Reading HDF5 file: %s\n', Hfile);
    SWUP_DOMAVG = hdf5read(Hfile, 'hda_swup');

    % top level is all zeros so use top - 1 level
    % SWUP_DOMAVG is (x,y,z,t)
    SWUP(j,:) = squeeze(SWUP_DOMAVG(:,:,end-1,:));

    % generate the legend text
    if (i == 1)
      LegText{j} = sprintf('CCN: %d/cc', CCN(j));
    end
  end
  fprintf('\n');

LegText

  % Compute the time average of the data in SWUP
  %   SWUP: rows are CCN value, columns are time so we just want to take the
  %   average of each row.
  SWUP_TDAVG(i,:) = mean(SWUP,2);
  Xtext{i} = sprintf('SST: %dK', SST(i));
end

% Plot a bar graph showing the relationship of the time,domain averaged top-of-model SW
% Note that LegText is already set to what we need from the time series plots above
fprintf('Generating bar chart for time domain averaged SW\n');
Fig = figure;

bar(SWUP_TDAVG);
set(gca, 'FontSize', 20);
title('Upward SW Flux: Top of Model');
set(gca, 'XtickLabel', Xtext);
ylabel('Radiative Flux (W/m^2)');
ylim([ 0 150 ]);
legend(LegText);
legend boxoff;

Pfile = sprintf('%s/swup_tdavg.jpg', Pdir);
fprintf('  Saving bar chart in: %s\n', Pfile);
saveas(Fig,Pfile);
close(Fig);
fprintf('\n');

end
