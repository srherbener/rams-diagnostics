% Script to create rainfall intensity PDF
%

clear;

% Read in the PCPRR data
%
% After reading in, the dimensions will be (x,y,t)

h5file = '/Users/steveherbener/Downloads/PCPRR.h5';
PCPRR = hdf5read(h5file, '/PCPRR');

% Generate PDF for PCPRR
Bins = [0.001, 0.01, 0.1, 1, 10, 100 ];
RrHists = GenHist2d(PCPRR, Bins, min(Bins), max(Bins));

% Save the histograms
hdf5write('PCPRR.hist.h5', '/Hists', RrHists);
hdf5write('PCPRR.hist.h5', '/Bins', Bins, 'WriteMode', 'append');