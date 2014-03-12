function [ ] = GenCfCdBarPlots(ConfigFile)
% GenCfCdBarPlots generate cloud fraction, max cloud depth bar plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
Ddir = Config.DiagDir;
Pdir = Config.PlotDir;

% Variables are all linear arrays
%   ith entry has corresponding SST(i), CCN(i), GCCN(i) values associated with it

InFile = sprintf('%s/cf_cd_bar_graph_data.h5', Ddir);
CF_AVG   = squeeze(hdf5read(InFile, 'cf_avg'));
COT_AVG  = squeeze(hdf5read(InFile, 'cot_avg'));
COT_L0P01_AVG  = squeeze(hdf5read(InFile, 'cot_avg_lwp_0p01'));
COT_L0P10_AVG  = squeeze(hdf5read(InFile, 'cot_avg_lwp_0p10'));
COT_L1P00_AVG  = squeeze(hdf5read(InFile, 'cot_avg_lwp_1p00'));
COT_LB1_AVG  = squeeze(hdf5read(InFile, 'cot_avg_lwp_b1'));
COT_LB2_AVG  = squeeze(hdf5read(InFile, 'cot_avg_lwp_b2'));
COT_LB3_AVG  = squeeze(hdf5read(InFile, 'cot_avg_lwp_b3'));
COT_LB4_AVG  = squeeze(hdf5read(InFile, 'cot_avg_lwp_b4'));
CD_AVG   = squeeze(hdf5read(InFile, 'cd_avg'));
ACC_PCP  = squeeze(hdf5read(InFile, 'acc_pcp'))/1000; % meters

CCN  = squeeze(hdf5read(InFile, 'ccn'));
SST  = squeeze(hdf5read(InFile, 'sst'));
GCCN = squeeze(hdf5read(InFile, 'gccn'));

% Read in histograms and convert to single values using ReduceHists
PrFile = sprintf('%s/hist_data.h5', Ddir);
PR_HIST = squeeze(hdf5read(PrFile, 'pcprr_hist'));
PR_BINS = squeeze(hdf5read(PrFile, 'pcprr_bins'));
PR_AVG = ReduceHists(PR_HIST, 1, PR_BINS, 'wtmean');

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First look at the GCCN = 1e-5 (essentially off, especially since the mean radius was small at 1um)
% Make vertical bar plots where the SST values are color coded bars, and the different CCN values go along the x-axis

Select = SST == 293 & GCCN == 1e-5; % Sel has 1's where we want to select data out of arrays, zeros elsewhere
Nsel = sum(Select);

X = 1:Nsel; % using numbers allows for even spacing on the bar plot

CFA = zeros([ Nsel 3 ]);
COTA = zeros([ Nsel 3 ]);
COT_L0P01_A = zeros([ Nsel 3 ]);
COT_L0P10_A = zeros([ Nsel 3 ]);
COT_L1P00_A = zeros([ Nsel 3 ]);
COT_LB1_A = zeros([ Nsel 3 ]);
COT_LB2_A = zeros([ Nsel 3 ]);
COT_LB3_A = zeros([ Nsel 3 ]);
COT_LB4_A = zeros([ Nsel 3 ]);
CDA = zeros([ Nsel 3 ]);
ACP = zeros([ Nsel 3 ]);
PRA = zeros([ Nsel 3 ]);

CFA(:,1) = CF_AVG(Select);
COTA(:,1) = COT_AVG(Select);
COT_L0P01_A(:,1) = COT_L0P01_AVG(Select);
COT_L0P10_A(:,1) = COT_L0P10_AVG(Select);
COT_L1P00_A(:,1) = COT_L1P00_AVG(Select);
COT_LB1_A(:,1) = COT_LB1_AVG(Select);
COT_LB2_A(:,1) = COT_LB2_AVG(Select);
COT_LB3_A(:,1) = COT_LB3_AVG(Select);
COT_LB4_A(:,1) = COT_LB4_AVG(Select);
CDA(:,1) = CD_AVG(Select);
ACP(:,1) = ACC_PCP(Select);
PRA(:,1) = PR_AVG(Select);

Select = SST == 298 & GCCN == 1e-5;
CFA(:,2) = CF_AVG(Select);
COTA(:,2) = COT_AVG(Select);
COT_L0P01_A(:,2) = COT_L0P01_AVG(Select);
COT_L0P10_A(:,2) = COT_L0P10_AVG(Select);
COT_L1P00_A(:,2) = COT_L1P00_AVG(Select);
COT_LB1_A(:,2) = COT_LB1_AVG(Select);
COT_LB2_A(:,2) = COT_LB2_AVG(Select);
COT_LB3_A(:,2) = COT_LB3_AVG(Select);
COT_LB4_A(:,2) = COT_LB4_AVG(Select);
CDA(:,2) = CD_AVG(Select);
ACP(:,2) = ACC_PCP(Select);
PRA(:,2) = PR_AVG(Select);

Select = SST == 303 & GCCN == 1e-5;
CFA(:,3) = CF_AVG(Select);
COTA(:,3) = COT_AVG(Select);
COT_L0P01_A(:,3) = COT_L0P01_AVG(Select);
COT_L0P10_A(:,3) = COT_L0P10_AVG(Select);
COT_L1P00_A(:,3) = COT_L1P00_AVG(Select);
COT_LB1_A(:,3) = COT_LB1_AVG(Select);
COT_LB2_A(:,3) = COT_LB2_AVG(Select);
COT_LB3_A(:,3) = COT_LB3_AVG(Select);
COT_LB4_A(:,3) = COT_LB4_AVG(Select);
CDA(:,3) = CD_AVG(Select);
ACP(:,3) = ACC_PCP(Select);
PRA(:,3) = PR_AVG(Select);

% Got x and y data for the bar plot. Now set up the plotting.
iaxis = 1;

% Large font
AxisProps(iaxis).Name = 'FontSize';
AxisProps(iaxis).Val  = 20;
iaxis = iaxis + 1;

% Put the actual CCN values on the tick marks for the x-axis.
X_CCN = CCN(Select);
for i = 1:Nsel
  Xlabels{i} = sprintf('%d', X_CCN(i));
end
AxisProps(iaxis).Name = 'XTickLabel';
AxisProps(iaxis).Val = Xlabels;
iaxis = iaxis + 1;

% Misc plotting

% bar colors - triplet is [ r g b ]
% want black cyan magenta
Bcolors = { [ 0 0 0 ] [ 0 1 1 ] [ 1 0 1 ] };

LegText = { '293 K', '298 K', '303 K' };

Xlabel = 'N_a (#/cc)';
CFlabel = 'Cloud Fraction';
COTlabel = 'Cloud Optical Thickness';
COTlabel_lpw_0p01 = 'Cloud Optical Thickness';
COTlabel_lpw_0p10 = 'Cloud Optical Thickness';
COTlabel_lpw_1p00 = 'Cloud Optical Thickness';
COTlabel_lpw_b1 = 'Cloud Optical Thickness';
COTlabel_lpw_b2 = 'Cloud Optical Thickness';
COTlabel_lpw_b3 = 'Cloud Optical Thickness';
COTlabel_lpw_b4 = 'Cloud Optical Thickness';
CDlabel = 'Max Cloud Depth (m)';
APlabel = 'Accum Precip (m)';
PRlabel = 'Avg Precip Rate (mm/h)';

CFAfile = sprintf('%s/cf_avg_bars.jpg', Pdir);
COTAfile = sprintf('%s/cot_avg_bars.jpg', Pdir);
COT_L0P01_Afile = sprintf('%s/cot_lwp_0p01_avg_bars.jpg', Pdir);
COT_L0P10_Afile = sprintf('%s/cot_lwp_0p10_avg_bars.jpg', Pdir);
COT_L1P00_Afile = sprintf('%s/cot_lwp_1p00_avg_bars.jpg', Pdir);
COT_LB1_Afile = sprintf('%s/cot_lwp_b1_avg_bars.jpg', Pdir);
COT_LB2_Afile = sprintf('%s/cot_lwp_b2_avg_bars.jpg', Pdir);
COT_LB3_Afile = sprintf('%s/cot_lwp_b3_avg_bars.jpg', Pdir);
COT_LB4_Afile = sprintf('%s/cot_lwp_b4_avg_bars.jpg', Pdir);
CDAfile = sprintf('%s/cd_avg_bars.jpg', Pdir);
ACPfile = sprintf('%s/acc_pcp_bars.jpg', Pdir);
PRAfile = sprintf('%s/pcprr_avg_bars.jpg', Pdir);

% Make six plots: CD start end avg, CF start end avg
% PlotBarSet( X, Y, Ptitle, Pmarkers, Xlabel, Ylabel, Bcolors, LegText, LegLoc, AxisProps, OutFile )

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 1.5 ];

PlotBarSet(X, CFA,  '',         { 'a' }, Xlabel, CFlabel,  Bcolors, LegText, 'NorthWest', AxisProps, CFAfile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 5 ];
PlotBarSet(X, COTA, 'All Points',         { 'a' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COTAfile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 8 ];
PlotBarSet(X, COT_L0P01_A, 'LWP > 0.01 mm', { 'b' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_L0P01_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 25 ];
PlotBarSet(X, COT_L0P10_A, 'LWP > 0.10 mm', { 'c' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_L0P10_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 100 ];
PlotBarSet(X, COT_L1P00_A, 'LWP > 1.00 mm', { 'd' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_L1P00_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 0.4 ];
PlotBarSet(X, COT_LB1_A, 'LWP < 0.01 mm', { 'b' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_LB1_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 3 ];
PlotBarSet(X, COT_LB2_A, '0.01 mm <= LWP < 0.10 mm', { 'c' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_LB2_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 15 ];
PlotBarSet(X, COT_LB3_A, '0.10 mm <= LWP < 1.00 mm', { 'd' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_LB3_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 100 ];
PlotBarSet(X, COT_LB4_A, 'LWP >= 1.00 mm', { 'e' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_LB4_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 6000 ];

PlotBarSet(X, CDA, '',         { 'c' }, Xlabel, CDlabel, Bcolors, LegText, 'NorthWest', AxisProps, CDAfile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 400 ];

PlotBarSet(X, ACP, '',         { 'b' }, Xlabel, APlabel, Bcolors, LegText, 'NorthWest', AxisProps, ACPfile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 0.3 ];

PlotBarSet(X, PRA, '',         { 'd' }, Xlabel, PRlabel, Bcolors, LegText, 'NorthWest', AxisProps, PRAfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next look at the GCCN impact (mean radius, 3um)
%
% Have SST = 298 and CCN = 50, 400, 1600; and one more SST = 303, CCN = 400

Select = SST == 298 & GCCN == 1e-5 & (CCN == 50 | CCN == 400 | CCN == 1600);
Nsel = sum(Select);

X = 1:Nsel; % using numbers allows for even spacing on the bar plot

CFA = zeros([ Nsel 4 ]);
COTA = zeros([ Nsel 4 ]);
COT_L0P01_A = zeros([ Nsel 4 ]);
COT_L0P10_A = zeros([ Nsel 4 ]);
COT_L1P00_A = zeros([ Nsel 4 ]);
COT_LB1_A = zeros([ Nsel 4 ]);
COT_LB2_A = zeros([ Nsel 4 ]);
COT_LB3_A = zeros([ Nsel 4 ]);
COT_LB4_A = zeros([ Nsel 4 ]);
CDA = zeros([ Nsel 4 ]);
ACP = zeros([ Nsel 4 ]);
PRA = zeros([ Nsel 4 ]);

CFA(:,1) = CF_AVG(Select);
COTA(:,1) = COT_AVG(Select);
COT_L0P01_A(:,1) = COT_L0P01_AVG(Select);
COT_L0P10_A(:,1) = COT_L0P10_AVG(Select);
COT_L1P00_A(:,1) = COT_L1P00_AVG(Select);
COT_LB1_A(:,1) = COT_LB1_AVG(Select);
COT_LB2_A(:,1) = COT_LB2_AVG(Select);
COT_LB3_A(:,1) = COT_LB3_AVG(Select);
COT_LB4_A(:,1) = COT_LB4_AVG(Select);
CDA(:,1) = CD_AVG(Select);
ACP(:,1) = ACC_PCP(Select);
PRA(:,1) = PR_AVG(Select);

Select = SST == 298 & GCCN == 1e-4;
CFA(:,2) = CF_AVG(Select);
COTA(:,2) = COT_AVG(Select);
COT_L0P01_A(:,2) = COT_L0P01_AVG(Select);
COT_L0P10_A(:,2) = COT_L0P10_AVG(Select);
COT_L1P00_A(:,2) = COT_L1P00_AVG(Select);
COT_LB1_A(:,2) = COT_LB1_AVG(Select);
COT_LB2_A(:,2) = COT_LB2_AVG(Select);
COT_LB3_A(:,2) = COT_LB3_AVG(Select);
COT_LB4_A(:,2) = COT_LB4_AVG(Select);
CDA(:,2) = CD_AVG(Select);
ACP(:,2) = ACC_PCP(Select);
PRA(:,2) = PR_AVG(Select);

Select = SST == 298 & GCCN == 1e-2;
CFA(:,3) = CF_AVG(Select);
COTA(:,3) = COT_AVG(Select);
COT_L0P01_A(:,3) = COT_L0P01_AVG(Select);
COT_L0P10_A(:,3) = COT_L0P10_AVG(Select);
COT_L1P00_A(:,3) = COT_L1P00_AVG(Select);
COT_LB1_A(:,3) = COT_LB1_AVG(Select);
COT_LB2_A(:,3) = COT_LB2_AVG(Select);
COT_LB3_A(:,3) = COT_LB3_AVG(Select);
COT_LB4_A(:,3) = COT_LB4_AVG(Select);
CDA(:,3) = CD_AVG(Select);
ACP(:,3) = ACC_PCP(Select);
PRA(:,3) = PR_AVG(Select);

Select = SST == 298 & GCCN == 1;
CFA(:,4) = CF_AVG(Select);
COTA(:,4) = COT_AVG(Select);
COT_L0P01_A(:,4) = COT_L0P01_AVG(Select);
COT_L0P10_A(:,4) = COT_L0P10_AVG(Select);
COT_L1P00_A(:,4) = COT_L1P00_AVG(Select);
COT_LB1_A(:,4) = COT_LB1_AVG(Select);
COT_LB2_A(:,4) = COT_LB2_AVG(Select);
COT_LB3_A(:,4) = COT_LB3_AVG(Select);
COT_LB4_A(:,4) = COT_LB4_AVG(Select);
CDA(:,4) = CD_AVG(Select);
ACP(:,4) = ACC_PCP(Select);
PRA(:,4) = PR_AVG(Select);

% Got x and y data for the bar plot. Now set up the plotting.
iaxis = 1;

% Large font
AxisProps(iaxis).Name = 'FontSize';
AxisProps(iaxis).Val  = 20;
iaxis = iaxis + 1;

% Put the actual CCN values on the tick marks for the x-axis.
AxisProps(iaxis).Name = 'XTickLabel';
AxisProps(iaxis).Val = { '50' '400' '1600' };
iaxis = iaxis + 1;

% Misc plotting

% bar colors - triplet is [ r g b ]
% want black cyan blue magenta
Bcolors = { [ 0 0 0 ] [ 0 1 1 ] [ 0 0 1 ] [ 1 0 1 ] };

LegText = { '10^-^5/cc', '10^-^4/cc', '10^-^2/cc' '1/cc' };

Xlabel = 'N_a (#/cc)';
CFlabel = 'Cloud Fraction';
COTlabel = 'Cloud Optical Thickness';
CDlabel = 'Max Cloud Depth (m)';
APlabel = 'Accum Precip (m)';
PRlabel = 'Avg Precip Rate (mm/h)';

CFAfile = sprintf('%s/cf_avg_bars_gccn.jpg', Pdir);
COTAfile = sprintf('%s/cot_avg_bars_gccn.jpg', Pdir);
COT_L0P01_Afile = sprintf('%s/cot_lwp_0p01_avg_bars_gccn.jpg', Pdir);
COT_L0P10_Afile = sprintf('%s/cot_lwp_0p10_avg_bars_gccn.jpg', Pdir);
COT_L1P00_Afile = sprintf('%s/cot_lwp_1p00_avg_bars_gccn.jpg', Pdir);
COT_LB1_Afile = sprintf('%s/cot_lwp_b1_avg_bars_gccn.jpg', Pdir);
COT_LB2_Afile = sprintf('%s/cot_lwp_b2_avg_bars_gccn.jpg', Pdir);
COT_LB3_Afile = sprintf('%s/cot_lwp_b3_avg_bars_gccn.jpg', Pdir);
COT_LB4_Afile = sprintf('%s/cot_lwp_b4_avg_bars_gccn.jpg', Pdir);
CDAfile = sprintf('%s/cd_avg_bars_gccn.jpg', Pdir);
ACPfile = sprintf('%s/acc_pcp_bars_gccn.jpg', Pdir);
PRAfile = sprintf('%s/pcprr_avg_bars_gccn.jpg', Pdir);

% Make four plots: CD start and end, CF start and end
% PlotBarSet( X, Y, Ptitle, Pmarkers, Xlabel, Ylabel, Bcolors, LegText, LegLoc, AxisProps, OutFile )

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 1.5 ];

PlotBarSet(X, CFA,  '',         { 'a' }, Xlabel, CFlabel,  Bcolors, LegText, 'NorthWest', AxisProps, CFAfile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 5 ];
PlotBarSet(X, COTA, 'All points', { 'a' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COTAfile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 8 ];
PlotBarSet(X, COT_L0P01_A, 'LWP > 0.01 mm', { 'b' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_L0P01_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 25 ];
PlotBarSet(X, COT_L0P10_A, 'LWP > 0.10 mm', { 'c' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_L0P10_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 100 ];
PlotBarSet(X, COT_L1P00_A, 'LWP > 1.00 mm', { 'd' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_L1P00_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 0.2 ];
PlotBarSet(X, COT_LB1_A, 'LWP < 0.01 mm', { 'b' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_LB1_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 3 ];
PlotBarSet(X, COT_LB2_A, '0.01 mm <= LWP < 0.10 mm', { 'c' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_LB2_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 10 ];
PlotBarSet(X, COT_LB3_A, '0.10 mm <= LWP < 1.00 mm', { 'd' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_LB3_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 100 ];
PlotBarSet(X, COT_LB4_A, 'LWP >= 1.00 mm', { 'e' }, Xlabel, COTlabel, Bcolors, LegText, 'NorthWest', AxisProps, COT_LB4_Afile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 6000 ];

PlotBarSet(X, CDA, '',         { 'c' }, Xlabel, CDlabel, Bcolors, LegText, 'NorthWest', AxisProps, CDAfile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 400 ];

PlotBarSet(X, ACP, '',         { 'b' }, Xlabel, APlabel, Bcolors, LegText, 'NorthWest', AxisProps, ACPfile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 0.3 ];

PlotBarSet(X, PRA, '',         { 'd' }, Xlabel, PRlabel, Bcolors, LegText, 'NorthWest', AxisProps, PRAfile);

end
