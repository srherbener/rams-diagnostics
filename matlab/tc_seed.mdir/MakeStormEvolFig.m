function [] = MakeStormEvolFig(ConfigFile)
% MakeStormEvolFig - merge plots showing Ke-Vt and RMW

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

PlotDir = Config.PlotDir;
FigDir = Config.FigDir;

% make sure output directory exists
if (exist(FigDir, 'dir') ~= 7)
  mkdir(FigDir);
end

% grab the Ke-Vt and RMW plots
FileList{1} = sprintf('%s/KeVt.fig', PlotDir);
FileList{2} = sprintf('%s/TS_rmw.fig', PlotDir);

OutFile = sprintf('%s/StormEvolution.jpg', FigDir);

% merge in the figures
[ Fig, Splots ] = figmerge(FileList, [ 1 2 ]);

% Get the handles of the subplots, note Splots from figmerge call
% will contain stale subplot handles.
clear Splots;
KVplot = subplot(1,2,1);
RMWplot = subplot(1,2,2);

% Create a dummy axis that covers the whole plotting region.
% This can be used to place text, titles, etc.
DAxis = axes('Position', [ 0 0 1 1 ], 'Visible', 'off');
text(0.5, 0.95, 'Storm Evolution', 'HorizontalAlignment', 'center', 'FontSize', 24);

% Adjust aspect ratio of subplots
set(KVplot, 'PlotBoxAspectRatio', [ 1 1 1 ]);
set(RMWplot, 'PlotBoxAspectRatio', [ 1 1 1 ]);

% Adjust the position of the subplots
% The phases plot needs to move a tad to the left
% The region plot needs to move a tad to the right
% These are both controlled by the first entry, 'x' (position -> [ x y w h ])
KV_pos = get(KVplot, 'Position');
RMW_pos = get(RMWplot, 'Position');
KV_pos(1) = KV_pos(1) - 0.02;  % move left
RMW_pos(1) = RMW_pos(1) + 0.03;  % move right
set(KVplot, 'Position', KV_pos);
set(RMWplot, 'Position', RMW_pos);

saveas(Fig, OutFile);
close(Fig);

%%% % Use TileFigs to create an image with two plots side by side
%%% TempFile = sprintf('%s.%s', tempname, 'jpg');
%%% TileFigs(FileList, [ 1 2 ], TempFile);
%%% Img = imread(TempFile);
%%% 
%%% Fig = figure;
%%% 
%%% imshow(Img);
%%% 
%%% saveas(Fig, OutFile);
%%% close(Fig);
%%% 
%%% delete(TempFile);

end
