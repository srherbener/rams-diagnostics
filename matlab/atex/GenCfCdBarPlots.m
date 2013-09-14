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
CF_START = squeeze(hdf5read(InFile, 'cf_start'));
CF_END   = squeeze(hdf5read(InFile, 'cf_end'));
CF_AVG   = squeeze(hdf5read(InFile, 'cf_avg'));
CD_START = squeeze(hdf5read(InFile, 'cd_start'));
CD_END   = squeeze(hdf5read(InFile, 'cd_end'));
CD_AVG   = squeeze(hdf5read(InFile, 'cd_avg'));

CCN  = squeeze(hdf5read(InFile, 'ccn'));
SST  = squeeze(hdf5read(InFile, 'sst'));
GCCN = squeeze(hdf5read(InFile, 'gccn'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First look at the GCCN = 1e-5 (essentially off, especially since the mean radius was small at 1um)
% Make vertical bar plots where the SST values are color coded bars, and the different CCN values go along the x-axis

Select = SST == 293 & GCCN == 1e-5; % Sel has 1's where we want to select data out of arrays, zeros elsewhere
Nsel = sum(Select);

X = 1:Nsel; % using numbers allows for even spacing on the bar plot

CFS = zeros([ Nsel 3 ]);
CDS = zeros([ Nsel 3 ]);
CFE = zeros([ Nsel 3 ]);
CDE = zeros([ Nsel 3 ]);
CFA = zeros([ Nsel 3 ]);
CDA = zeros([ Nsel 3 ]);

CFS(:,1) = CF_START(Select);
CDS(:,1) = CD_START(Select);
CFE(:,1) = CF_END(Select);
CDE(:,1) = CD_END(Select);
CFA(:,1) = CF_AVG(Select);
CDA(:,1) = CD_AVG(Select);

Select = SST == 298 & GCCN == 1e-5;
CFS(:,2) = CF_START(Select);
CDS(:,2) = CD_START(Select);
CFE(:,2) = CF_END(Select);
CDE(:,2) = CD_END(Select);
CFA(:,2) = CF_AVG(Select);
CDA(:,2) = CD_AVG(Select);

Select = SST == 303 & GCCN == 1e-5;
CFS(:,3) = CF_START(Select);
CDS(:,3) = CD_START(Select);
CFE(:,3) = CF_END(Select);
CDE(:,3) = CD_END(Select);
CFA(:,3) = CF_AVG(Select);
CDA(:,3) = CD_AVG(Select);

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

Xlabel = 'N_c (#/cc)';
CFlabel = 'Cloud Fraction';
CDlabel = 'Max Cloud Depth (m)';

CFSfile = sprintf('%s/cf_start_bars.jpg', Pdir);
CFEfile = sprintf('%s/cf_end_bars.jpg', Pdir);
CDSfile = sprintf('%s/cd_start_bars.jpg', Pdir);
CDEfile = sprintf('%s/cd_end_bars.jpg', Pdir);
CFAfile = sprintf('%s/cf_avg_bars.jpg', Pdir);
CDAfile = sprintf('%s/cd_avg_bars.jpg', Pdir);

% Make six plots: CD start end avg, CF start end avg
% PlotBarSet( X, Y, Ptitle, Pmarkers, Xlabel, Ylabel, Bcolors, LegText, LegLoc, AxisProps, OutFile )

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 1.5 ];

PlotBarSet(X, CFS, 't = 12 h', { 'a' }, Xlabel, CFlabel, Bcolors, LegText, 'NorthWest', AxisProps, CFSfile);
PlotBarSet(X, CFE, 't = 36 h', { 'b' }, Xlabel, CFlabel, Bcolors, LegText, 'NorthWest', AxisProps, CFEfile);
PlotBarSet(X, CFA, '',         { 'a' }, Xlabel, CFlabel, Bcolors, LegText, 'NorthWest', AxisProps, CFAfile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 6000 ];

PlotBarSet(X, CDS, 't = 12 h', { 'a' }, Xlabel, CDlabel, Bcolors, LegText, 'NorthWest', AxisProps, CDSfile);
PlotBarSet(X, CDE, 't = 36 h', { 'b' }, Xlabel, CDlabel, Bcolors, LegText, 'NorthWest', AxisProps, CDEfile);
PlotBarSet(X, CDA, '',         { 'b' }, Xlabel, CDlabel, Bcolors, LegText, 'NorthWest', AxisProps, CDAfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next look at the GCCN impact (mean radius, 3um)
%
% Have SST = 298 and CCN = 50, 400, 1600; and one more SST = 303, CCN = 400

Select = SST == 298 & GCCN == 1e-5 & (CCN == 50 | CCN == 400 | CCN == 1600);
Nsel = sum(Select);

X = 1:Nsel; % using numbers allows for even spacing on the bar plot

CFS = zeros([ Nsel 4 ]);
CDS = zeros([ Nsel 4 ]);
CFE = zeros([ Nsel 4 ]);
CDE = zeros([ Nsel 4 ]);
CFA = zeros([ Nsel 4 ]);
CDA = zeros([ Nsel 4 ]);

CFS(:,1) = CF_START(Select);
CDS(:,1) = CD_START(Select);
CFE(:,1) = CF_END(Select);
CDE(:,1) = CD_END(Select);
CFA(:,1) = CF_AVG(Select);
CDA(:,1) = CD_AVG(Select);

Select = SST == 298 & GCCN == 1e-4;
CFS(:,2) = CF_START(Select);
CDS(:,2) = CD_START(Select);
CFE(:,2) = CF_END(Select);
CDE(:,2) = CD_END(Select);
CFA(:,2) = CF_AVG(Select);
CDA(:,2) = CD_AVG(Select);

Select = SST == 298 & GCCN == 1e-2;
CFS(:,3) = CF_START(Select);
CDS(:,3) = CD_START(Select);
CFE(:,3) = CF_END(Select);
CDE(:,3) = CD_END(Select);
CFA(:,3) = CF_AVG(Select);
CDA(:,3) = CD_AVG(Select);

Select = SST == 298 & GCCN == 1;
CFS(:,4) = CF_START(Select);
CDS(:,4) = CD_START(Select);
CFE(:,4) = CF_END(Select);
CDE(:,4) = CD_END(Select);
CFA(:,4) = CF_AVG(Select);
CDA(:,4) = CD_AVG(Select);

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

Xlabel = 'N_c (#/cc)';
CFlabel = 'Cloud Fraction';
CDlabel = 'Max Cloud Depth (m)';

CFSfile = sprintf('%s/cf_start_bars_gccn.jpg', Pdir);
CFEfile = sprintf('%s/cf_end_bars_gccn.jpg', Pdir);
CDSfile = sprintf('%s/cd_start_bars_gccn.jpg', Pdir);
CDEfile = sprintf('%s/cd_end_bars_gccn.jpg', Pdir);
CFAfile = sprintf('%s/cf_avg_bars_gccn.jpg', Pdir);
CDAfile = sprintf('%s/cd_avg_bars_gccn.jpg', Pdir);

% Make four plots: CD start and end, CF start and end
% PlotBarSet( X, Y, Ptitle, Pmarkers, Xlabel, Ylabel, Bcolors, LegText, LegLoc, AxisProps, OutFile )

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 1.5 ];

PlotBarSet(X, CFS, 't = 12 h', { 'a' }, Xlabel, CFlabel, Bcolors, LegText, 'NorthWest', AxisProps, CFSfile);
PlotBarSet(X, CFE, 't = 36 h', { 'b' }, Xlabel, CFlabel, Bcolors, LegText, 'NorthWest', AxisProps, CFEfile);
PlotBarSet(X, CFA, '',         { 'a' }, Xlabel, CFlabel, Bcolors, LegText, 'NorthWest', AxisProps, CFAfile);

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 6000 ];

PlotBarSet(X, CDS, 't = 12 h', { 'a' }, Xlabel, CDlabel, Bcolors, LegText, 'NorthWest', AxisProps, CDSfile);
PlotBarSet(X, CDE, 't = 36 h', { 'b' }, Xlabel, CDlabel, Bcolors, LegText, 'NorthWest', AxisProps, CDEfile);
PlotBarSet(X, CDA, '',         { 'b' }, Xlabel, CDlabel, Bcolors, LegText, 'NorthWest', AxisProps, CDAfile);

end
