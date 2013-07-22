function [ ] = MakeBlankPlot(ConfigFile)
%MakeBlankPlot function to create a blank plot - useful for tiling plots to make a figure

Config = ReadConfig(ConfigFile);

Pdir = Config.PlotDir;
% make sure output directory exists
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

OutFile = sprintf('%s/Blank.jpg', Pdir);
fprintf('Writing: %s\n', OutFile);

Fig = figure;
saveas(Fig, OutFile);
close(Fig);

end
