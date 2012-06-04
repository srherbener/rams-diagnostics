% Script to plot selected histograms on a single figure.

clear;

% Read in the histgram data

h5_fin = 'DIAG/PCPRR.hist.allexps.h5';

Hists = hdf5read(h5_fin, '/Hists');
Bins = hdf5read(h5_fin, '/Bins');
CCN = hdf5read(h5_fin, '/CcnConc');
SST = hdf5read(h5_fin, '/Sst');
Npts = hdf5read(h5_fin, '/Npts');

% Pick out two sets of data for plotting:
%
%    SST = 298, CCN = 50, 200, 800, 1600
%    SST = 303, CCN = 50, 200, 800, 1600
%
% The data are organized as:
%    Experiment number     SST     CCN
%
%           1              298      50
%           2              303      50
%           3              298     100
%           4              303     100
%           5              298     200
%           6              303     200
%           7              298     400
%           8              303     400
%           9              298     800
%          10              303     800
%          11              298    1200
%          12              303    1200
%          13              298    1600
%          14              303    1600
% 
% So the first set consists of Epx numbers: 1 5  9 13
% and the second set:                       2 6 10 14

Set1 = [ 1 5  9 13 ];
Set2 = [ 2 6 10 14 ];

H1 = Hists(:,Set1);
B1 = Bins(:,Set1);
C1 = CCN(Set1);
S1 = SST(Set1);
N1 = Npts(Set1);

H2 = Hists(:,Set2);
B2 = Bins(:,Set2);
C2 = CCN(Set2);
S2 = SST(Set2);
N2 = Npts(Set2);

Lstyles = { '-k', '--b', ':r', '-.g' };

PlotHistSet(H1, B1, C1, S1, N1, Lstyles, 'DIAG/PCPRR.hist1.jpg');
PlotHistSet(H2, B2, C2, S2, N2, Lstyles, 'DIAG/PCPRR.hist2.jpg');



