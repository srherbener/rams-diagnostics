function [ ] = GenLegendPlots(ConfigFile)
% GenLegendPlots generate plots that show only a legend (for use with multi-panel plots)

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
UndefVal = Config.UndefVal;

Pdir = Config.PlotDir;

% make sure output directory exists
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

Fsize = 20;
Lwidth = 2;

% list of specs for each legend
% format for each legend spec:
%   { PlotSet Ncols Wscale Hscale OutFile   }
%
%   Ncols - number of columns for arranging legend entries
%   Wscale, Hscale - factors for cropping jpeg image
%      Wscale is multiplied by the image width
%      Hscale is multiplied by the image height
LegendList = {
  { 'CO_UP_S293_G10M5_ALL' 3 0.95 0.10 'Legend_CO_CCN.jpg' }
  { 'CO_UP_C0050_G10M5'    3 0.50 0.10 'Legend_CO_SST.jpg' }
  { 'CG_UP_C50_C1600_S298' 2 0.85 0.10 'Legend_CG_CCN.jpg' }

  };
Nleg = length(LegendList);

for ileg = 1:Nleg
  PSname  = LegendList{ileg}{1};
  Ncols   = LegendList{ileg}{2};
  Wscale  = LegendList{ileg}{3};
  Hscale  = LegendList{ileg}{4};
  OutFile = LegendList{ileg}{5};

  PSnum = find(strcmp(PSname, { Config.PlotSets(:).Name }));
  if (isempty(PSnum))
    fprintf('WARNING: Could not find a match for plot set: %s, legend number: %s\n', PSname, ileg); 
    fprintf('         Skipping this legend plot.\n');
    fprintf('\n');
    continue;
  end
  OutFile = sprintf('%s/%s', Pdir, OutFile);

  fprintf('Generating legend plot:\n');
  fprintf('  Plot Set: %s --> %d\n', PSname, PSnum);

  % create a plot with legend only:
  %   1) plot dummy curves
  %   2) turn off visibility of the plot
  %   3) plot legend
  PlotSet = Config.PlotSets(PSnum);
  clear LegText;
  clear LineColors;
  clear LineStyles;
  clear Gscales;
  Ncases = PlotSet.Ncases;
  for icase = 1:Ncases
    LegText{icase}    = PlotSet.Cases(icase).Legend;
    LineColors{icase} = PlotSet.Cases(icase).Lcolor;
    LineStyles{icase} = PlotSet.Cases(icase).Lstyle;
    Gscales{icase}    = PlotSet.Cases(icase).Lgscale;
  end

  % dummy data - set up to use Plot2dSet
  X = repmat( [ 1:10 ], [ Ncases 1 ]);
  Y = ones([ Ncases length(X) ] ) .* X;

  % dummy plot
  Fig = figure;

  Ptitle = '';
  Pmarkers = '';
  Xlabel = '';
  Ylabel = '';
  LegLoc = 'none';
  AxisProps(1).Name = 'Position'; % position axes across entire figure extent
  AxisProps(1).Val = [ 0 0 1 1 ];
  AddMeas = 'none'; 
  Plot2dSet(X, Y, Ptitle, Pmarkers, Xlabel, Ylabel, LineColors, LineStyles, Gscales, LegText, LegLoc, AxisProps, AddMeas, Fig);

  % turn off visibility of plot
  set(gca, 'Visible', 'off');
  Lines = findobj(gca, 'Type', 'line');
  for i = 1:length(Lines)
    set(Lines(i), 'Visible', 'off');
  end

  % add in legend, position in lower left corner
  legend(LegText, 'FontSize', Fsize, 'Orientation', 'horizontal', 'Location', 'SouthWest');

  % crop the jpeg image
  OrigPaperPos = get(gcf, 'PaperPosition');
  PaperPos(1) = OrigPaperPos(1);
  PaperPos(2) = OrigPaperPos(2);
  PaperPos(3) = OrigPaperPos(3) * Wscale;
  PaperPos(4) = OrigPaperPos(4) * Hscale;
  set(gcf, 'PaperPosition', PaperPos);

  fprintf('  Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  set(gcf, 'PaperPosition', OrigPaperPos); % restore paper position
  close(Fig);
  fprintf('\n');
end
