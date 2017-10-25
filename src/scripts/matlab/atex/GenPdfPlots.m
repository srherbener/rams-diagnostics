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
COT_HIST = squeeze(hdf5read(InFile, 'copt_thick_hist'));
COT_BINS = squeeze(hdf5read(InFile, 'copt_thick_bins'));

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
PDF_COT293 = COT_HIST(:,Select);

Select = SST == 298 & GCCN == 1e-5;
PDF_PR298 = PR_HIST(:,Select);
PDF_ALB298 = ALB_HIST(:,Select);
PDF_COT298 = COT_HIST(:,Select);

Select = SST == 303 & GCCN == 1e-5;
PDF_PR303 = PR_HIST(:,Select);
PDF_ALB303 = ALB_HIST(:,Select);
PDF_COT303 = COT_HIST(:,Select);

for i = 1:Nsel
  PDF_PR293(:,i) = PDF_PR293(:,i) ./ sum(PDF_PR293(:,i));
  PDF_PR298(:,i) = PDF_PR298(:,i) ./ sum(PDF_PR298(:,i));
  PDF_PR303(:,i) = PDF_PR303(:,i) ./ sum(PDF_PR303(:,i));

  PDF_ALB293(:,i) = PDF_ALB293(:,i) ./ sum(PDF_ALB293(:,i));
  PDF_ALB298(:,i) = PDF_ALB298(:,i) ./ sum(PDF_ALB298(:,i));
  PDF_ALB303(:,i) = PDF_ALB303(:,i) ./ sum(PDF_ALB303(:,i));
  
  PDF_COT293(:,i) = PDF_COT293(:,i) ./ sum(PDF_COT293(:,i));
  PDF_COT298(:,i) = PDF_COT298(:,i) ./ sum(PDF_COT298(:,i));
  PDF_COT303(:,i) = PDF_COT303(:,i) ./ sum(PDF_COT303(:,i));
end

% Legend
X_CCN = CCN(Select);
for i = 1:Nsel
  LegText{i} = sprintf('%d/cc', X_CCN(i));
end

Lcolors = {
    'black'
    'green'
    'blue'
    'cyan'
    'red'
    'yellow'
    'magenta'
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

Xlabel = 'Precip Rate (mm h^-^1)';
Ylabel = '';

PR293file = sprintf('%s/pcprr_pdf_S293.jpg', Pdir);
PR298file = sprintf('%s/pcprr_pdf_S298.jpg', Pdir);
PR303file = sprintf('%s/pcprr_pdf_S303.jpg', Pdir);

AxisProps(1).Name = 'FontSize';
AxisProps(1).Val  = 35;

Xvals = [ 0.001 0.01 0.1 1 10 100 1000 ];
AxisProps(2).Name = 'Xtick';
AxisProps(2).Val = Xvals;

AxisProps(3).Name = 'Xscale';
AxisProps(3).Val  = 'log';

AxisProps(4).Name = 'Xlim';
AxisProps(4).Val  = [ 0.001 100 ];

AxisProps(5).Name = 'Ylim';
AxisProps(5).Val  = [ 0 0.1 ];

AxisProps(6).Name = 'Yscale';
AxisProps(6).Val  = 'log';

Yvals = [ 1e-5 1e-3 1e-1 ];
AxisProps(7).Name = 'Ytick';
AxisProps(7).Val = Yvals;


BINS = repmat(PR_BINS', [ Nsel 1 ]);

Plot2dSet(BINS, PDF_PR293', '293K', { 'g' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', PR293file);
Plot2dSet(BINS, PDF_PR298', '298K', { 'h' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', PR298file);
Plot2dSet(BINS, PDF_PR303', '303K', { 'i' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
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
AxisProps(3).Val  = [ 0 0.5 ];

AxisProps(4).Name = 'Yscale';
AxisProps(4).Val  = 'log';

Yvals = [ 1e-5 1e-3 1e-1 ];
AxisProps(5).Name = 'Ytick';
AxisProps(5).Val = Yvals;


BINS = repmat(ALB_BINS', [ Nsel 1 ]);

Plot2dSet(BINS, PDF_ALB293', '293K', { 'c' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', ALB293file);
Plot2dSet(BINS, PDF_ALB298', '298K', { 'f' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', ALB298file);
Plot2dSet(BINS, PDF_ALB303', '303K', { 'i' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', ALB303file);
     

% Cloud optical thickness plots
Xlabel = 'Cloud OT';
Ylabel = '';

COT293file = sprintf('%s/copt_thick_pdf_S293.jpg', Pdir);
COT298file = sprintf('%s/copt_thick_pdf_S298.jpg', Pdir);
COT303file = sprintf('%s/copt_thick_pdf_S303.jpg', Pdir);

clear AxisProps;

AxisProps(1).Name = 'FontSize';
AxisProps(1).Val  = 35;

AxisProps(2).Name = 'Xlim';
AxisProps(2).Val  = [ 0 300 ];

AxisProps(3).Name = 'Ylim';
AxisProps(3).Val  = [ 0 0.8 ];

AxisProps(4).Name = 'Yscale';
AxisProps(4).Val  = 'log';

Yvals = [ 1e-5 1e-3 1e-1 ];
AxisProps(5).Name = 'Ytick';
AxisProps(5).Val = Yvals;


BINS = repmat(COT_BINS', [ Nsel 1 ]);

Plot2dSet(BINS, PDF_COT293', '293K', { 'd' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', COT293file);
Plot2dSet(BINS, PDF_COT298', '298K', { 'e' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', COT298file);
Plot2dSet(BINS, PDF_COT303', '303K', { 'f' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', COT303file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next look at the GCCN impact (mean radius, 3um)
%
% Have SST = 298 and CCN = 50, 400, 1600; and one more SST = 303, CCN = 400

Select = SST == 298 & CCN == 50;
Nsel = sum(Select);
PDF_PR50 = PR_HIST(:,Select);
PDF_ALB50 = ALB_HIST(:,Select);
PDF_COT50 = COT_HIST(:,Select);

Select = SST == 298 & CCN == 400;
PDF_PR400 = PR_HIST(:,Select);
PDF_ALB400 = ALB_HIST(:,Select);
PDF_COT400 = COT_HIST(:,Select);

Select = SST == 298 & CCN == 1600;
PDF_PR1600 = PR_HIST(:,Select);
PDF_ALB1600 = ALB_HIST(:,Select);
PDF_COT1600 = COT_HIST(:,Select);

for i = 1:Nsel
  PDF_PR50(:,i) = PDF_PR50(:,i) ./ sum(PDF_PR50(:,i));
  PDF_PR400(:,i) = PDF_PR400(:,i) ./ sum(PDF_PR400(:,i));
  PDF_PR1600(:,i) = PDF_PR1600(:,i) ./ sum(PDF_PR1600(:,i));

  PDF_ALB50(:,i) = PDF_ALB50(:,i) ./ sum(PDF_ALB50(:,i));
  PDF_ALB400(:,i) = PDF_ALB400(:,i) ./ sum(PDF_ALB400(:,i));
  PDF_ALB1600(:,i) = PDF_ALB1600(:,i) ./ sum(PDF_ALB1600(:,i));
  
  PDF_COT50(:,i) = PDF_COT50(:,i) ./ sum(PDF_COT50(:,i));
  PDF_COT400(:,i) = PDF_COT400(:,i) ./ sum(PDF_COT400(:,i));
  PDF_COT1600(:,i) = PDF_COT1600(:,i) ./ sum(PDF_COT1600(:,i));
end

% Rearrange the columns - PDF's have columns corresponding to:
%    GCCN = 1e-5 1 1e-2 1e-4
% and want
%    GCCN = 1e-5 1e-4 1e-2 1

PDF_PR50 = PDF_PR50(:,[ 1 4 3 2 ]);
PDF_PR400 = PDF_PR400(:,[ 1 4 3 2 ]);
PDF_PR1600 = PDF_PR1600(:,[ 1 4 3 2 ]);

PDF_ALB50 = PDF_ALB50(:,[ 1 4 3 2 ]);
PDF_ALB400 = PDF_ALB400(:,[ 1 4 3 2 ]);
PDF_ALB1600 = PDF_ALB1600(:,[ 1 4 3 2 ]);

PDF_COT50 = PDF_COT50(:,[ 1 4 3 2 ]);
PDF_COT400 = PDF_COT400(:,[ 1 4 3 2 ]);
PDF_COT1600 = PDF_COT1600(:,[ 1 4 3 2 ]);

% Legend

LegText = {
    '10^-^5/cc'
    '10^-^4/cc'
    '10^-^2/cc'
    '1/cc'
    };

Lcolors = {
    'black'
    'cyan'
    'blue'
    'magenta'
    };

Lstyles = {
    '-'
    '-'
    '-'
    '-'
    };

Gscales = zeros([ Nsel 1 ]);

% Plot2dSet( X, Y, Ptitle, Pmarkers, Xlabel, Ylabel, Lcolors, Lstyles,
%            Gscales, LegText, LegLoc, AxisProps, AddMeas, OutFile )

Xlabel = 'Precip Rate (mm h^-^1)';
Ylabel = '';

PR50file = sprintf('%s/pcprr_pdf_C50_S298.jpg', Pdir);
PR400file = sprintf('%s/pcprr_pdf_C400_S298.jpg', Pdir);
PR1600file = sprintf('%s/pcprr_pdf_C1600_S298.jpg', Pdir);

AxisProps(1).Name = 'FontSize';
AxisProps(1).Val  = 35;

Xvals = [ 0.001 0.01 0.1 1 10 100 1000 ];
AxisProps(2).Name = 'Xtick';
AxisProps(2).Val = Xvals;

AxisProps(3).Name = 'Xscale';
AxisProps(3).Val  = 'log';

AxisProps(4).Name = 'Xlim';
AxisProps(4).Val  = [ 0.001 100 ];

AxisProps(5).Name = 'Ylim';
AxisProps(5).Val  = [ 0 0.1 ];

AxisProps(6).Name = 'Yscale';
AxisProps(6).Val  = 'log';

Yvals = [ 1e-5 1e-3 1e-1 ];
AxisProps(7).Name = 'Ytick';
AxisProps(7).Val = Yvals;


BINS = repmat(PR_BINS', [ Nsel 1 ]);

Plot2dSet(BINS, PDF_PR50', '298K 50/cc', { 'g' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', PR50file);
Plot2dSet(BINS, PDF_PR400', '298K 400/cc', { 'h' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', PR400file);
Plot2dSet(BINS, PDF_PR1600', '298K 1600/cc', { 'i' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', PR1600file);



% Albedo plots
Xlabel = 'Albedo';
Ylabel = '';

ALB50file = sprintf('%s/albedo_pdf_C50_S298.jpg', Pdir);
ALB400file = sprintf('%s/albedo_pdf_C400_S298.jpg', Pdir);
ALB1600file = sprintf('%s/albedo_pdf_C1600_S298.jpg', Pdir);

clear AxisProps;

AxisProps(1).Name = 'FontSize';
AxisProps(1).Val  = 35;

AxisProps(2).Name = 'Xlim';
AxisProps(2).Val  = [ 0 1 ];

AxisProps(3).Name = 'Ylim';
AxisProps(3).Val  = [ 0 0.5 ];

AxisProps(4).Name = 'Yscale';
AxisProps(4).Val  = 'log';

Yvals = [ 1e-5 1e-3 1e-1 ];
AxisProps(5).Name = 'Ytick';
AxisProps(5).Val = Yvals;


BINS = repmat(ALB_BINS', [ Nsel 1 ]);

Plot2dSet(BINS, PDF_ALB50', '298K 50/cc', { 'c' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', ALB50file);
Plot2dSet(BINS, PDF_ALB400', '298K 400/cc', { 'f' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', ALB400file);
Plot2dSet(BINS, PDF_ALB1600', '298K 1600/cc', { 'i' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', ALB1600file);
     

% Cloud optical thickness plots
Xlabel = 'Cloud OT';
Ylabel = '';

COT50file = sprintf('%s/copt_thick_pdf_C50_S298.jpg', Pdir);
COT400file = sprintf('%s/copt_thick_pdf_C400_S298.jpg', Pdir);
COT1600file = sprintf('%s/copt_thick_pdf_C1600_S298.jpg', Pdir);

clear AxisProps;

AxisProps(1).Name = 'FontSize';
AxisProps(1).Val  = 35;

AxisProps(2).Name = 'Xlim';
AxisProps(2).Val  = [ 0 300 ];

AxisProps(3).Name = 'Ylim';
AxisProps(3).Val  = [ 0 0.8 ];

AxisProps(4).Name = 'Yscale';
AxisProps(4).Val  = 'log';

Yvals = [ 1e-5 1e-3 1e-1 ];
AxisProps(5).Name = 'Ytick';
AxisProps(5).Val = Yvals;


BINS = repmat(COT_BINS', [ Nsel 1 ]);

Plot2dSet(BINS, PDF_COT50', '298K 50/cc', { 'd' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', COT50file);
Plot2dSet(BINS, PDF_COT400', '298K 400/cc', { 'e' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', COT400file);
Plot2dSet(BINS, PDF_COT1600', '298K 1600/cc', { 'f' }, Xlabel, Ylabel, Lcolors, Lstyles, ...
         Gscales, LegText, 'NorthEast', AxisProps, 'none', COT1600file);


end
