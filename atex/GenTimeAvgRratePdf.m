function [ ] = GenTimeAvgRratePdf(ConfigFile)
% GenTimeAvgRratePdf generate time averages of the rain rate (pcprr) histogram data

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);

Tdir = Config.TsavgDir;
OutDir = 'DIAGS';

RrateHistVar = 'hist_pcprr';

% Create the output directory if it doesn't exist
if (exist(OutDir,'dir') == 0)
  mkdir(OutDir);
end

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  RrateHistFile = sprintf('%s/%s_%s.h5', Tdir, RrateHistVar, Case);
  OutFile = sprintf('%s/tavg_%s_%s.h5', OutDir, RrateHistVar, Case);

  fprintf('***************************************************************\n');
  fprintf('Generating time averaged rain rate PDF:\n');
  fprintf('  Case: %s\n', Case);
  fprintf('  Input rain rain histogram file: %s\n', RrateHistFile);
  fprintf('  Output file: %s\n', OutFile);
  fprintf('\n');

  % The histogram is organized as (x,y,z,t) with x being the bins, y and z being dummy
  % coordinates (size = 1), and t the time steps. Use squeeze to compress into a 2D array
  % (x,t) that has bins and time steps.
  RRHIST = squeeze(hdf5read(RrateHistFile, RrateHistVar));
  BINS = hdf5read(RrateHistFile, 'x_coords');

  % The number of points selected for the histogram at each time step can vary so
  % to get a mean PDF, first convert each histogram to a PDF, then take the averages
  % of the PDF entries.
  % It's possible that a histogram has all zero counts. In this case, the sum will
  % be zero, which can be used to flag this case since there are no negative counts.
  % Exclude the all zeros histogram from the final averaging since it will throw
  % off keeping the sum of the averaged PDF equalling one.
  SUMS = sum(RRHIST,1);
  [ Nb, Nt ] = size(RRHIST);
  ipdf = 0;
  for i = 1:Nt
    if (SUMS(i) ~= 0)
      ipdf = ipdf + 1;
      PDF(:,ipdf) = RRHIST(:,i) / SUMS(i);
    end
  end
  AVGPDF = mean(PDF,2);

fprintf('DEBUG: Nt: %d, ipdf: %d\n', Nt, ipdf);

  hdf5write(OutFile, '/AvgPdf', AVGPDF);
  hdf5write(OutFile, '/Bins', BINS, 'WriteMode', 'append');

end
