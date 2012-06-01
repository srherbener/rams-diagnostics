% Script to create rainfall intensity histograms
%

clear;

% Read in the PCPRR data
%
% After reading in, the dimensions will be (x,y,t)

Exps = [ 'z.atex250m.100km.ccn0050.sst298' ],
%%% Exps = [ 'z.atex250m.100km.ccn0050.sst298',
%%% 'z.atex250m.100km.ccn0050.sst303',
%%% 'z.atex250m.100km.ccn0100.sst298',
%%% 'z.atex250m.100km.ccn0100.sst303',
%%% 'z.atex250m.100km.ccn0200.sst298',
%%% 'z.atex250m.100km.ccn0200.sst303',
%%% 'z.atex250m.100km.ccn0400.sst298',
%%% 'z.atex250m.100km.ccn0400.sst303',
%%% 'z.atex250m.100km.ccn0800.sst298',
%%% 'z.atex250m.100km.ccn0800.sst303',
%%% 'z.atex250m.100km.ccn1200.sst298',
%%% 'z.atex250m.100km.ccn1200.sst303',
%%% 'z.atex250m.100km.ccn1600.sst298',
%%% 'z.atex250m.100km.ccn1600.sst303' ];

Clevs_l = ( 0.0:0.03:0.3 );
Clevs_m = ( 0.0:0.008:0.08 );
Clevs_h = ( 0.0:0.001:0.01 );

Times = ( 0:432 );  % shift down one so that the first time step is time zero

% Rain rate is in the REVU var PCPRR
for i = 1:size(Exps,1)
  h5_fin = sprintf('DIAG/%s/PCPRR.hist.h5',Exps(i,:));
  Pfile_l = sprintf('PLOTS/%s/PCPRR.hist.light.jpg',Exps(i,:));
  Pfile_m = sprintf('PLOTS/%s/PCPRR.hist.medium.jpg',Exps(i,:));
  Pfile_h = sprintf('PLOTS/%s/PCPRR.hist.heavy.jpg',Exps(i,:));
  fprintf('Generating Histogram Plots:\n');
  fprintf('  Input file: %s\n',h5_fin);
  fprintf('  Output files:\n');
  fprintf('    %s\n', Pfile_l);
  fprintf('    %s\n', Pfile_m);
  fprintf('    %s\n', Pfile_h);
  fprintf('\n');

  Ptitle_l = { 'Fractional Area - light precip', Exps(i,:) };
  Ptitle_m = { 'Fractional Area - medium precip', Exps(i,:) };
  Ptitle_h = { 'Fractional Area - heavy precip', Exps(i,:) };

  % read in the histogram data
  Hist_l = hdf5read(h5_fin, '/Hists_l');
  Hist_m = hdf5read(h5_fin, '/Hists_m');
  Hist_h = hdf5read(h5_fin, '/Hists_h');

  Bins_l = hdf5read(h5_fin, '/Bins_l');
  Bins_m = hdf5read(h5_fin, '/Bins_m');
  Bins_h = hdf5read(h5_fin, '/Bins_h');

  % generate the plots
  PlotHist2d(Hist_l, Times, Bins_l, Clevs_l, Ptitle_l, Pfile_l);
  PlotHist2d(Hist_m, Times, Bins_m, Clevs_m, Ptitle_m, Pfile_m);
  PlotHist2d(Hist_h, Times, Bins_h, Clevs_h, Ptitle_h, Pfile_h);
end
