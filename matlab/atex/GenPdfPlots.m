function [ ] = GenPdfPlots(ConfigFile)
% GenPdfPlots generate pdf plots (pcprr so far)

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
Ddir = Config.DiagDir;
Pdir = Config.PlotDir;

% Variables are all linear arrays
%   ith entry has corresponding SST(i), CCN(i), GCCN(i) values associated with it

InFile = sprintf('%s/pdf_data.h5', Ddir);
PCPRR_PDF  = squeeze(hdf5read(InFile, 'pcprr_pdf'));
PCPRR_BINS = squeeze(hdf5read(InFile, 'pcprr_bins'));

CCN  = squeeze(hdf5read(InFile, 'ccn'));
SST  = squeeze(hdf5read(InFile, 'sst'));
GCCN = squeeze(hdf5read(InFile, 'gccn'));

Nb = length(PCPRR_BINS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First look at the GCCN = 1e-5 (essentially off, especially since the mean radius was small at 1um)
% Make vertical bar plots where the SST values are color coded bars, and the different CCN values go along the x-axis

Select = SST == 293 & GCCN == 1e-5; % Sel has 1's where we want to select data out of arrays, zeros elsewhere
Nsel = sum(Select);
PDF_PR293 = PCPRR_PDF(:,Select);

Select = SST == 298 & GCCN == 1e-5;
PDF_PR298 = PCPRR_PDF(:,Select);

Select = SST == 303 & GCCN == 1e-5;
PDF_PR303 = PCPRR_PDF(:,Select);

% Got x and y data for the bar plot. Now set up the plotting.
iaxis = 1;

% Large font
AxisProps(iaxis).Name = 'FontSize';
AxisProps(iaxis).Val  = 20;
iaxis = iaxis + 1;

% log scale on x-axis
AxisProps(iaxis).Name = 'Xscale';
AxisProps(iaxis).Val  = 'log';
iaxis = iaxis + 1;

% Put the actual CCN values on the tick marks for the x-axis.
X_CCN = CCN(Select);
for i = 1:Nsel
  LegText{i} = sprintf('%d/cc', X_CCN(i));
end

Lcolors = {
    'k'
    'g'
    'b'
    'c'
    'r'
    'y'
    'm'
    };

Lstyles = {
    '-'
    '-'
    '-'
    '-'
    '-'
    '-'
    '-'
    };

Gscales = zeros([ Nsel 1 ]);

% Misc plotting

Xlabel = 'Precip Rate (mm/h)';
Ylabel = '';

PR293file = sprintf('%s/pcprr_pdf_S293.jpg', Pdir);
PR298file = sprintf('%s/pcprr_pdf_S298.jpg', Pdir);
PR303file = sprintf('%s/pcprr_pdf_S303.jpg', Pdir);

% Make six plots: CD start end avg, CF start end avg
% Plot2dSet( X, Y, Ptitle, Pmarkers, Xlabel, Ylabel, Lcolors, Lstyles,
%            Gscales, LegText, LegLoc, AxisProps, AddMeas, OutFile )

AxisProps(iaxis).Name = 'Xlim';
AxisProps(iaxis).Val  = [ 0.001 100 ];
iaxis = iaxis + 1;

AxisProps(iaxis).Name = 'Ylim';
AxisProps(iaxis).Val  = [ 0 0.035 ];


BINS = repmat(PCPRR_BINS', [ Nsel 1 ]);

Plot2dSet(BINS, PDF_PR293', 'SST: 293K', { 'c' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', PR293file);
Plot2dSet(BINS, PDF_PR298', 'SST: 298K', { 'f' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', PR298file);
Plot2dSet(BINS, PDF_PR303', 'SST: 303K', { 'i' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', PR303file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next look at the GCCN impact (mean radius, 3um)
%
% Have SST = 298 and CCN = 50, 400, 1600; and one more SST = 303, CCN = 400

% Select = SST == 298 & GCCN == 1e-5 & (CCN == 50 | CCN == 400 | CCN == 1600);
% Nsel = sum(Select);
% 
% X = 1:Nsel; % using numbers allows for even spacing on the bar plot
% 
% CFS = zeros([ Nsel 4 ]);
% CDS = zeros([ Nsel 4 ]);
% CFE = zeros([ Nsel 4 ]);
% CDE = zeros([ Nsel 4 ]);
% CFA = zeros([ Nsel 4 ]);
% CDA = zeros([ Nsel 4 ]);
% ACP = zeros([ Nsel 4 ]);
% 
% CFS(:,1) = CF_START(Select);
% CDS(:,1) = CD_START(Select);
% CFE(:,1) = CF_END(Select);
% CDE(:,1) = CD_END(Select);
% CFA(:,1) = CF_AVG(Select);
% CDA(:,1) = CD_AVG(Select);
% ACP(:,1) = ACC_PCP(Select);
% 
% Select = SST == 298 & GCCN == 1e-4;
% CFS(:,2) = CF_START(Select);
% CDS(:,2) = CD_START(Select);
% CFE(:,2) = CF_END(Select);
% CDE(:,2) = CD_END(Select);
% CFA(:,2) = CF_AVG(Select);
% CDA(:,2) = CD_AVG(Select);
% ACP(:,2) = ACC_PCP(Select);
% 
% Select = SST == 298 & GCCN == 1e-2;
% CFS(:,3) = CF_START(Select);
% CDS(:,3) = CD_START(Select);
% CFE(:,3) = CF_END(Select);
% CDE(:,3) = CD_END(Select);
% CFA(:,3) = CF_AVG(Select);
% CDA(:,3) = CD_AVG(Select);
% ACP(:,3) = ACC_PCP(Select);
% 
% Select = SST == 298 & GCCN == 1;
% CFS(:,4) = CF_START(Select);
% CDS(:,4) = CD_START(Select);
% CFE(:,4) = CF_END(Select);
% CDE(:,4) = CD_END(Select);
% CFA(:,4) = CF_AVG(Select);
% CDA(:,4) = CD_AVG(Select);
% ACP(:,4) = ACC_PCP(Select);
% 
% % Got x and y data for the bar plot. Now set up the plotting.
% iaxis = 1;
% 
% % Large font
% AxisProps(iaxis).Name = 'FontSize';
% AxisProps(iaxis).Val  = 20;
% iaxis = iaxis + 1;
% 
% % Put the actual CCN values on the tick marks for the x-axis.
% AxisProps(iaxis).Name = 'XTickLabel';
% AxisProps(iaxis).Val = { '50' '400' '1600' };
% iaxis = iaxis + 1;
% 
% % Misc plotting
% 
% % bar colors - triplet is [ r g b ]
% % want black cyan blue magenta
% Bcolors = { [ 0 0 0 ] [ 0 1 1 ] [ 0 0 1 ] [ 1 0 1 ] };
% 
% LegText = { '10^-^5/cc', '10^-^4/cc', '10^-^2/cc' '1/cc' };
% 
% Xlabel = 'N_c (#/cc)';
% CFlabel = 'Cloud Fraction';
% CDlabel = 'Max Cloud Depth (m)';
% APlabel = 'Accum Precip (m)';
% 
% CFSfile = sprintf('%s/cf_start_bars_gccn.jpg', Pdir);
% CFEfile = sprintf('%s/cf_end_bars_gccn.jpg', Pdir);
% CDSfile = sprintf('%s/cd_start_bars_gccn.jpg', Pdir);
% CDEfile = sprintf('%s/cd_end_bars_gccn.jpg', Pdir);
% CFAfile = sprintf('%s/cf_avg_bars_gccn.jpg', Pdir);
% CDAfile = sprintf('%s/cd_avg_bars_gccn.jpg', Pdir);
% ACPfile = sprintf('%s/acc_pcp_bars_gccn.jpg', Pdir);
% 
% % Make four plots: CD start and end, CF start and end
% % PlotBarSet( X, Y, Ptitle, Pmarkers, Xlabel, Ylabel, Bcolors, LegText, LegLoc, AxisProps, OutFile )
% 
% AxisProps(iaxis).Name = 'Ylim';
% AxisProps(iaxis).Val  = [ 0 1.5 ];
% 
% PlotBarSet(X, CFS, 't = 12 h', { 'a' }, Xlabel, CFlabel, Bcolors, LegText, 'NorthWest', AxisProps, CFSfile);
% PlotBarSet(X, CFE, 't = 36 h', { 'b' }, Xlabel, CFlabel, Bcolors, LegText, 'NorthWest', AxisProps, CFEfile);
% PlotBarSet(X, CFA, '',         { 'a' }, Xlabel, CFlabel, Bcolors, LegText, 'NorthWest', AxisProps, CFAfile);
% 
% AxisProps(iaxis).Name = 'Ylim';
% AxisProps(iaxis).Val  = [ 0 6000 ];
% 
% PlotBarSet(X, CDS, 't = 12 h', { 'a' }, Xlabel, CDlabel, Bcolors, LegText, 'NorthWest', AxisProps, CDSfile);
% PlotBarSet(X, CDE, 't = 36 h', { 'b' }, Xlabel, CDlabel, Bcolors, LegText, 'NorthWest', AxisProps, CDEfile);
% PlotBarSet(X, CDA, '',         { 'b' }, Xlabel, CDlabel, Bcolors, LegText, 'NorthWest', AxisProps, CDAfile);
% 
% AxisProps(iaxis).Name = 'Ylim';
% AxisProps(iaxis).Val  = [ 0 400 ];
% 
% PlotBarSet(X, ACP, '',         { 'a' }, Xlabel, APlabel, Bcolors, LegText, 'NorthWest', AxisProps, ACPfile);

end
