function [ ] = PlotMaxWindAll(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

Adir = Config.AzavgDir;

Pdir = Config.PlotDir;
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

Cases = {
 'SST_CONST'
 'SST_OBS'
 'SST_OBS_UBMIN'
 };

LegText = {
 'SST CONST'
 'SST OBS'
 'SST OBS UBMIN'
 };

% read in the Vt data
Nc = length(Cases);
for icase = 1:Nc
  Hfile = sprintf('../%s/%s/speed_t_TSD_3GRIDS.h5', Cases{icase}, Adir);
  Hdset = 'speed_t';
  fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);

  % S is (r,z,t)
  S = squeeze(hdf5read(Hfile, Hdset));
  if (icase == 1)
    [ Nr Nz Nt ] = size(S);
    TIMES = 1:Nt;
    SPEED_T = zeros([ Nc Nr Nz Nt ]);
  end
  SPEED_T(icase,:,:,:) = S;
end

% ST is organized as (c,r,z,t) so take the max along the middle two dimensions
% Generate the time series of the maximum of the near surface Vt. Use Z = 2
% since Z = 1 is below surface in RAMS.
% Also, skip the first radial entry. There is a miscalculation in the
% speed_t diagnostic due to TS Debby travelling over the Verde Islands.
R1 = 2;
R2 = Nr;
ST_MAX = squeeze(max(squeeze(SPEED_T(:,R1:R2,2,:)), [] ,2));

% NHC Best Track (every six hours) data
% time step 1 from the simulation is where the NHC data starts
% NHC wind is in mph, multiply by 0.447 to convert to m/s
NHC_WIND  = [ 35   35   35   35   35   40   45   50   50   50   50   50   50   ] * 0.447;
NHC_TIMES = (1:6:73);

% plot
OutFile = sprintf('%s/TsDebbyWindAll.jpg', Pdir);

FigWind = figure;
set(gca, 'FontSize', 18);
SimST = plot(TIMES, ST_MAX, 'linewi', 2);
hold on;
NhcST = plot(NHC_TIMES, NHC_WIND, '+k', 'linewi', 2);
title('TS Debby Maximum Tangential Wind Speed');
xlabel('Time');
set(gca,'xtick', (13:24:61));
set(gca,'xticklabel', { 'Aug22:0Z', 'Aug23:0Z', 'Aug24:0Z' });
ylabel('Wind Speed (m/s)');
ylim([ 0 25 ]);
legend(SimST, LegText, 'Location', 'SouthEast');

fprintf('Writing: %s\n', OutFile);
saveas(FigWind, OutFile);
close(FigWind);
