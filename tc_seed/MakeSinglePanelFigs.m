function [] = MakeSinglePanelFigs(ConfigFile)
% MakeSinglePanelFigs - covert matlab fig to JPEG

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

PlotDir = Config.PlotDir;
FigDir = Config.FigDir;

% make sure output directory exists
if (exist(FigDir, 'dir') ~= 7)
  mkdir(FigDir);
end

% list of drawings to translate
FigList = { 'InitVortex'
            'SampleCcnProf'
            'SpinUpVtSpl' };

for i = 1:length(FigList)
  InFile  = sprintf('%s/%s.fig', PlotDir, FigList{i});
  OutFile = sprintf('%s/%s.jpg', FigDir,  FigList{i});

  fprintf('MATLAB figure file: %s\n', InFile);

  OrigFig = openfig(InFile);
  OrigPlots = findall(0, 'type', 'axes');

  Fig = figure;
  for i = 1:length(OrigPlots)
    copyobj(OrigPlots(i), Fig);
  end
  saveas(Fig, OutFile);
  close(Fig);

  close(OrigFig);

  fprintf('Writing file: %s\n', OutFile);
  fprintf('\n');
end

end
