function [] = MakeSamplingFigure(ConfigFile)
% MakeSamplingFigure - merge plots showing RI, SS time slices and EW, CO, SS radial regions

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

PlotDir = Config.PlotDir;
FigDir = Config.FigDir;

% make sure output directory exists
if (exist(FigDir, 'dir') ~= 7)
  mkdir(FigDir);
end

% grab the storm phases and storm regions figures from the plot dir
FileList{1} = sprintf('%s/StormPhases.fig', PlotDir);
FileList{2} = sprintf('%s/StormRegions.fig', PlotDir);

OutFile = sprintf('%s/SamplingFig.jpg', FigDir);

[ Fig, Splots ] = figmerge(FileList, [ 1 2 ]);
% reverse gray map
colormap(flipud(colormap('gray')));

% Get the handles of the subplots, note Splots from figmerge call
% will contain stale subplot handles.
SPhases = subplot(1,2,1);
SRegion = subplot(1,2,2);

% Create a dummy axes that covers the whole plotting region.
% This can be used to place text, titles, etc.
DAxes = axes('Position', [ 0 0 1 1 ], 'Visible', 'off');
text(0.5, 0.95, 'Temporal and Spatial Averaging', 'HorizontalAlignment', 'center', 'FontSize', 24);

% Adjust aspect ratio of subplots
set(SPhases, 'PlotBoxAspectRatio', [ 1 1 1 ]);
set(SRegion, 'PlotBoxAspectRatio', [ 1 1 1 ]);

saveas(Fig, OutFile);
close(Fig);

end
