function [ ] = PlotSlpAll(ConfigFile)

[ Config ] = ReadConfig(ConfigFile);

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

% read in the sea level pressure
Nc = length(Cases);
for icase = 1:Nc
  Hfile = sprintf('../%s/AzAveragedData/sea_press_TSD_3GRIDS.h5', Cases{icase});
  Hdset = 'sea_press';
  fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);

  % S is (r,t)
  S = squeeze(hdf5read(Hfile, Hdset));
  if (icase == 1)
    [ Nr Nt ] = size(S);
    TIMES = 1:Nt;
    SPRESS = zeros([ Nc Nr Nt ]);
  end
  SPRESS(icase,:,:) = S;
end

% SPRESS is organized as (c,r,t) so take the minimum along the first dimension
% generate the time series of the minimum SLP of the horizontal domain
R1 = 1;
R2 = Nr;
SLP_MIN = squeeze(min(SPRESS(:,R1:R2,:),[],2));

% NHC Best Track (every six hours) data
% time step 1 from the simulation is where the NHC data starts
NHC_SLP = [ 1007 1007 1007 1007 1005 1003 1002 1001 1001 1000  999 1000 1000 ];
NHC_TIMES = (1:6:73);


% plot
OutFile = sprintf('%s/TsDebbyPressAll.jpg', Pdir);

FigPress = figure;
set(gca, 'FontSize', 18);
SimSLP = plot(TIMES, SLP_MIN, 'linewi', 2);
hold on;
NhcSLP = plot(NHC_TIMES, NHC_SLP, '+b', 'linewi', 2);
title('TS Debby Minumum SLP');
xlabel('Time');
set(gca,'xtick', (13:24:61));
set(gca,'xticklabel', { 'Aug22:0Z', 'Aug23:0Z', 'Aug24:0Z' });
ylabel('Pressure (mb)');
ylim([ 990 1010 ]);
legend(SimSLP, LegText, 'Location', 'NorthEast');

fprintf('Writing: %s\n', OutFile);
saveas(FigPress, OutFile);
close(FigPress);
