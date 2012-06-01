% Script to create rainfall intensity histograms
%

clear;

% Read in the PCPRR data
%
% After reading in, the dimensions will be (x,y,t)

Exps = [ 'z.atex250m.100km.ccn0050.sst298',
'z.atex250m.100km.ccn0050.sst303',
'z.atex250m.100km.ccn0100.sst298',
'z.atex250m.100km.ccn0100.sst303',
'z.atex250m.100km.ccn0200.sst298',
'z.atex250m.100km.ccn0200.sst303',
'z.atex250m.100km.ccn0400.sst298',
'z.atex250m.100km.ccn0400.sst303',
'z.atex250m.100km.ccn0800.sst298',
'z.atex250m.100km.ccn0800.sst303',
'z.atex250m.100km.ccn1200.sst298',
'z.atex250m.100km.ccn1200.sst303',
'z.atex250m.100km.ccn1600.sst298',
'z.atex250m.100km.ccn1600.sst303' ];

Bins_l = (0.01:0.01:0.1);
Bins_m = (0.1:0.1:1);
Bins_h = (1:1:20);

% Rain rate is in the REVU var PCPRR
for i = 1:size(Exps,1)
  h5_fin = sprintf('REVU/%s/PCPRR.h5',Exps(i,:));
  h5_fout = sprintf('DIAG/%s/PCPRR.hist.h5',Exps(i,:));
  fprintf('Generating Histograms:\n');
  fprintf('  Input file: %s\n',h5_fin);
  fprintf('  Output file: %s\n', h5_fout);
  fprintf('\n');

  % Read in rain rate data
  PCPRR = hdf5read(h5_fin, '/PCPRR');

  % Record number of points in horzontal domain
  [ Nx, Ny, Nt ] = size(PCPRR); 
  Npts = Nx * Ny;

  % Generate histograms for PCPRR
  Hists_l = GenHist2d(PCPRR, Bins_l, min(Bins_l), max(Bins_l));
  Hists_m = GenHist2d(PCPRR, Bins_m, min(Bins_m), max(Bins_m));
  Hists_h = GenHist2d(PCPRR, Bins_h, min(Bins_h), max(Bins_h));

  % Save the histograms
  hdf5write(h5_fout, '/Hists_l', Hists_l);
  hdf5write(h5_fout, '/Bins_l', Bins_l, 'WriteMode', 'append');

  hdf5write(h5_fout, '/Hists_m', Hists_m, 'WriteMode', 'append');
  hdf5write(h5_fout, '/Bins_m', Bins_m, 'WriteMode', 'append');

  hdf5write(h5_fout, '/Hists_h', Hists_h, 'WriteMode', 'append');
  hdf5write(h5_fout, '/Bins_h', Bins_h, 'WriteMode', 'append');

  hdf5write(h5_fout, '/Npts', Npts, 'WriteMode', 'append');
end
