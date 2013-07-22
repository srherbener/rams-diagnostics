function [] = MakeKevtRmwFig(ConfigFile)
% MakeKevtRmwFig - merge plots showing KE-Vt and RMW

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

PlotDir = Config.PlotDir;

% grab the storm phases and storm regions figures from the plot dir
FileList{1} = sprintf('%s/KeVt.fig', PlotDir);
FileList{2} = sprintf('%s/TS_rmw.fig', PlotDir);

OutFile = sprintf('%s/KevtRmw.fig', PlotDir);

% Merge the plots stack in one column
[ Fig, Splots ] = figmerge(FileList, [ 2 1 ]);

% Get the handles of the subplots (Splots from figmerge call
% will contain stale subplot handles).
SKevt = subplot(2,1,1);
set(SKevt, 'OuterPosition', [ 0.05 0.55 0.9 0.4 ]);
legend(SKevt, 'show', 'Location', 'NorthWest');
legend boxoff;

SRmw  = subplot(2,1,2);
set(SRmw,  'OuterPosition', [ 0.05 0.05 0.9 0.4 ]);
legend(SRmw, 'show', 'Location', 'NorthEast');
legend boxoff;

% Resize the entire figure
Fpos = get(gcf, 'Position');
Fpos(3) = 250;
Fpos(4) = 500;
set(gcf, 'Position', Fpos);

saveas(Fig, OutFile);
close(Fig);

end
