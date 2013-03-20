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

% Get the handles of the subplots, note Splots from figmerge call
% will contain stale subplot handles.
SKevt = subplot(2,1,1);
SRmw  = subplot(2,1,2);

% Create a dummy axes that covers the whole plotting region.
% This can be used to place text, titles, etc.
Daxes = axes('Position', [ 0 0 1 1 ], 'Visible', 'off');

% Adjust aspect ratio axes within the subplots
set(SKevt, 'DataAspectRatio', [  10 1 1 ]);
set(SRmw,  'DataAspectRatio', [ 1.6 1 1 ]);


% Adjust the position of the subplots
% The phases plot needs to move a tad to the left
% The region plot needs to move a tad to the right
% These are both controlled by the first entry, 'x' (position -> [ x y w h ])
%SK_pos = get(SKevt, 'Position');
%SR_pos = get(SRmw,  'Position');
%SK_pos(1) = SK_pos(1) - 0.02;  % move left
%SR_pos(1) = SR_pos(1) + 0.03;  % move right
%set(SKevt, 'Position', SK_pos);
%set(SRmw,  'Position', SR_pos);

saveas(Fig, OutFile);
close(Fig);

end
