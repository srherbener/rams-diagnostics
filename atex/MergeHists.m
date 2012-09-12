% Script to merge together histograms from GenRrateHist.m
%

clear;

% Read in the light, medium and heavy histograms for each experiment and
% merge them together into one dataset.

%Exps = [ 'z.atex250m.100km.ccn0050.sst298' ];
Exps = [ 'z.atex250m.100km.ccn0050.sst298'
'z.atex250m.100km.ccn0050.sst303'
'z.atex250m.100km.ccn0100.sst298'
'z.atex250m.100km.ccn0100.sst303'
'z.atex250m.100km.ccn0200.sst298'
'z.atex250m.100km.ccn0200.sst303'
'z.atex250m.100km.ccn0400.sst298'
'z.atex250m.100km.ccn0400.sst303'
'z.atex250m.100km.ccn0800.sst298'
'z.atex250m.100km.ccn0800.sst303'
'z.atex250m.100km.ccn1200.sst298'
'z.atex250m.100km.ccn1200.sst303'
'z.atex250m.100km.ccn1600.sst298'
'z.atex250m.100km.ccn1600.sst303' ];

% make these agree with order of the experments
CCN = [ 50 50 100 100 200 200 400 400 800 800 1200 1200 1600 1600 ]';

SST = [ 298 303 298 303 298 303 298 303 298 303 298 303 298 303 ]';


for i = 1:size(Exps,1)
  h5_fin = sprintf('DIAG/%s/PCPRR.hist.h5',Exps(i,:));
  fprintf('Mergeing Histograms:\n');
  fprintf('  Input file: %s\n',h5_fin);
  fprintf('\n');

  % Read in histograms
  Hists = hdf5read(h5_fin, '/Hists');
  Bins = hdf5read(h5_fin, '/Bins');
  Npts = hdf5read(h5_fin, '/Npts');
  
  % preallocate the output arrays - this assumes all histograms read out of
  % the hdf5 files are the same size.
  if (i == 1)
      Nexps = size(Exps,1);
      [ Nbins, Nt ] = size(Hists);
      H = zeros(Nbins,Nexps);
      B = zeros(Nbins,Nexps);
      N = zeros(1,Nexps);
  end
  
  % Compute the accumulated, across time, histogram; and scale by the
  % number of time points (Ncols) so that the resulting histogram is scaled
  % by the total points (Nx * Ny * Nt) of the original data.
  [ Nrows, Ncols ] = size(Hists);
  Hacc = sum(Hists,2);  % sum up the rows --> accumulate each bin across
                        % the time steps
  Npts = Npts * Ncols;
  
  % Record each experiments accumulated histrogram in an array where the
  % rows are the bins and the columns are the experiments.
  H(:,i) = Hacc;
  N(i) = Npts;
  B(:,i) = Bins;
end

% Save the histograms
h5_fout = sprintf('DIAG/PCPRR.hist.allexps.h5');
fprintf('Saving merged histograms in: %s\n', h5_fout);

hdf5write(h5_fout, '/Hists', H);
hdf5write(h5_fout, '/Bins', B, 'WriteMode', 'append');
hdf5write(h5_fout, '/Npts', N, 'WriteMode', 'append');
hdf5write(h5_fout, '/CcnConc', CCN, 'WriteMode', 'append');
hdf5write(h5_fout, '/Sst', SST, 'WriteMode', 'append');

