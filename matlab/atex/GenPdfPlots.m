function [ ] = GenPdfPlots(ConfigFile)
% GenPdfPlots generate pdf plots (pcprr so far)

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
Ddir = Config.DiagDir;
Pdir = Config.PlotDir;

% Variables are all linear arrays
%   ith entry has corresponding SST(i), CCN(i), GCCN(i) values associated with it

InFile = sprintf('%s/hist_data.h5', Ddir);
PR_HIST  = squeeze(hdf5read(InFile, 'pcprr_hist'));
PR_BINS  = squeeze(hdf5read(InFile, 'pcprr_bins'));
ALB_HIST = squeeze(hdf5read(InFile, 'albedo_hist'));
ALB_BINS = squeeze(hdf5read(InFile, 'albedo_bins'));

CCN  = squeeze(hdf5read(InFile, 'ccn'));
SST  = squeeze(hdf5read(InFile, 'sst'));
GCCN = squeeze(hdf5read(InFile, 'gccn'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First look at the GCCN = 1e-5 (essentially off, especially since the mean radius was small at 1um)
% Make vertical bar plots where the SST values are color coded bars, and the different CCN values go along the x-axis

Select = SST == 293 & GCCN == 1e-5; % Sel has 1's where we want to select data out of arrays, zeros elsewhere
Nsel = sum(Select);
PDF_PR293 = PR_HIST(:,Select);
PDF_ALB293 = ALB_HIST(:,Select);

Select = SST == 298 & GCCN == 1e-5;
PDF_PR298 = PR_HIST(:,Select);
PDF_ALB298 = ALB_HIST(:,Select);

Select = SST == 303 & GCCN == 1e-5;
PDF_PR303 = PR_HIST(:,Select);
PDF_ALB303 = ALB_HIST(:,Select);

for i = 1:Nsel
  PDF_PR293(:,i) = PDF_PR293(:,i) ./ sum(PDF_PR293(:,i));
  PDF_PR298(:,i) = PDF_PR298(:,i) ./ sum(PDF_PR298(:,i));
  PDF_PR303(:,i) = PDF_PR303(:,i) ./ sum(PDF_PR303(:,i));

  PDF_ALB293(:,i) = PDF_ALB293(:,i) ./ sum(PDF_ALB293(:,i));
  PDF_ALB298(:,i) = PDF_ALB298(:,i) ./ sum(PDF_ALB298(:,i));
  PDF_ALB303(:,i) = PDF_ALB303(:,i) ./ sum(PDF_ALB303(:,i));
end

% Legend
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

% Plot2dSet( X, Y, Ptitle, Pmarkers, Xlabel, Ylabel, Lcolors, Lstyles,
%            Gscales, LegText, LegLoc, AxisProps, AddMeas, OutFile )

Xlabel = 'Precip Rate (mm/h)';
Ylabel = '';

PR293file = sprintf('%s/pcprr_pdf_S293.jpg', Pdir);
PR298file = sprintf('%s/pcprr_pdf_S298.jpg', Pdir);
PR303file = sprintf('%s/pcprr_pdf_S303.jpg', Pdir);

AxisProps(1).Name = 'FontSize';
AxisProps(1).Val  = 35;

Xvals = [ 0.001 0.01 0.1 1 10 100 ];
AxisProps(2).Name = 'Xtick';
AxisProps(2).Val = Xvals;

AxisProps(3).Name = 'Xscale';
AxisProps(3).Val  = 'log';

AxisProps(4).Name = 'Xlim';
AxisProps(4).Val  = [ 0.001 100 ];

AxisProps(5).Name = 'Ylim';
AxisProps(5).Val  = [ 0 0.035 ];


BINS = repmat(PR_BINS', [ Nsel 1 ]);

Plot2dSet(BINS, PDF_PR293', 'SST: 293K', { 'b' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', PR293file);
Plot2dSet(BINS, PDF_PR298', 'SST: 298K', { 'e' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', PR298file);
Plot2dSet(BINS, PDF_PR303', 'SST: 303K', { 'h' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', PR303file);



% Albedo plots
Xlabel = 'Albedo';
Ylabel = '';

ALB293file = sprintf('%s/albedo_pdf_S293.jpg', Pdir);
ALB298file = sprintf('%s/albedo_pdf_S298.jpg', Pdir);
ALB303file = sprintf('%s/albedo_pdf_S303.jpg', Pdir);

clear AxisProps;

AxisProps(1).Name = 'FontSize';
AxisProps(1).Val  = 35;

AxisProps(2).Name = 'Xlim';
AxisProps(2).Val  = [ 0 1 ];

AxisProps(3).Name = 'Ylim';
AxisProps(3).Val  = [ 0 0.2 ];


BINS = repmat(ALB_BINS', [ Nsel 1 ]);

Plot2dSet(BINS, PDF_ALB293', 'SST: 293K', { 'c' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', ALB293file);
Plot2dSet(BINS, PDF_ALB298', 'SST: 298K', { 'f' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', ALB298file);
Plot2dSet(BINS, PDF_ALB303', 'SST: 303K', { 'i' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', ALB303file);

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
