function [ ] = GenHistFiles(ConfigFile)
% GenHistData generate data indexed by SST and CCN concentration for histograms

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Hdir = './HDF5';

    PrName = 'pcprr';
    PrFprefix = 'pcprr';
    CtName = 'cloud_opt_thick';
    CtFprefix = 'cloud_opt_thick';
    LwpName = 'vertint_cond';
    LwpFprefix = 'vint_cond';

    PR_BINS = 10.^(-3:0.02:2); % pcprr bins, even logarithmic spacing from 0.001 to 100.0
    ALB_BINS = 0:0.01:1; % albedo bins
    CT_BINS = 0:300;  % cloud optical thickness bins

    Nprb = length(PR_BINS);
    Nab  = length(ALB_BINS);
    Nctb  = length(CT_BINS);

    g = 0.85; % for albedo calc, scattering asymmetry param (cloud droplets interacting with visible light)

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;
        PrFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, PrFprefix, Case);
        CtFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, CtFprefix, Case);
        LwpFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, LwpFprefix, Case);

        fprintf('***************************************************************\n');
        fprintf('Generating hist data:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input precip rate file: %s\n', PrFile);
        fprintf('    Var name: %s\n', PrName);
        fprintf('  Input cloud optical thickness file: %s\n', CtFile);
        fprintf('    Var name: %s\n', CtName);
        fprintf('  LWP file: %s\n', LwpFile);
        fprintf('    Var name: %s\n', LwpName);
        fprintf('\n');

        % Use nctoolbox so that we can walk through the data one time step at a time.
        % Data in input files are organized as: (t, y, x) representing (time, lat, lon)
        PR_DS = ncgeodataset(PrFile);
        PR_VAR = PR_DS.geovariable(PrName);
        T_VAR = PR_DS.geovariable('t_coords');
        T = T_VAR.data(:);

        Size = PR_VAR.size;
        Nt = Size(1);
        Ny = Size(2);
        Nx = Size(3);

        Ntrim = (Nx-2) * (Ny-2); % used for reshaping 2D field after trimming off the borders

        CT_DS = ncgeodataset(CtFile);
        CT_VAR = CT_DS.geovariable(CtName);
 
        LWP_DS = ncgeodataset(LwpFile);
        LWP_VAR = LWP_DS.geovariable(LwpName);

        % Output will be of the form (Nb, Nt) where Nb is the number of bins in the histogram
        PR_HIST        = zeros([ Nprb Nt ]);
        PR_LWP_B1_HIST = zeros([ Nprb Nt ]);
        PR_LWP_B2_HIST = zeros([ Nprb Nt ]);
        PR_LWP_B3_HIST = zeros([ Nprb Nt ]);
        PR_LWP_B4_HIST = zeros([ Nprb Nt ]);

        CT_HIST        = zeros([ Nctb Nt ]);
        CT_LWP_B1_HIST = zeros([ Nctb Nt ]);
        CT_LWP_B2_HIST = zeros([ Nctb Nt ]);
        CT_LWP_B3_HIST = zeros([ Nctb Nt ]);
        CT_LWP_B4_HIST = zeros([ Nctb Nt ]);

        ALB_HIST        = zeros([ Nab Nt ]);
        ALB_LWP_B1_HIST = zeros([ Nab Nt ]);
        ALB_LWP_B2_HIST = zeros([ Nab Nt ]);
        ALB_LWP_B3_HIST = zeros([ Nab Nt ]);
        ALB_LWP_B4_HIST = zeros([ Nab Nt ]);

        % Walk through each time step forming the histogram counts for that time step. These
        % counts can then be combined in different ways during a post processing step.
        for i = 1:Nt
          % read in the 2D vars and reshape them into vectors
          % RAMS uses the borders for boundaries so trim these off
          PR  = squeeze(PR_VAR.data(i,:,:));
          CT  = squeeze(CT_VAR.data(i,:,:));
          LWP = squeeze(LWP_VAR.data(i,:,:));

          PR = reshape(PR(2:end-1,2:end-1), [ 1  Ntrim ]);
          CT = reshape(CT(2:end-1,2:end-1), [ 1  Ntrim ]);
          LWP = reshape(LWP(2:end-1,2:end-1), [ 1  Ntrim ]);

          % Divide PR and CT into 4 bins according to LWP (mm)
          %    b1:         LWP < 0.01
          %    b2: 0.01 <= LWP < 0.10
          %    b3: 0.10 <= LWP < 1.00
          %    b4: 1.00 <= LWP
          PR_LWP_B1 = PR(LWP < 0.01);
          PR_LWP_B2 = PR(LWP >= 0.01 & LWP < 0.10);
          PR_LWP_B3 = PR(LWP >= 0.10 & LWP < 1.00);
          PR_LWP_B4 = PR(LWP >= 1.00);

          CT_LWP_B1 = CT(LWP < 0.01);
          CT_LWP_B2 = CT(LWP >= 0.01 & LWP < 0.10);
          CT_LWP_B3 = CT(LWP >= 0.10 & LWP < 1.00);
          CT_LWP_B4 = CT(LWP >= 1.00);

          % form the histograms
          PR_HIST(:,i)        = histc(PR, PR_BINS);
          PR_LWP_B1_HIST(:,i) = histc(PR_LWP_B1, PR_BINS);
          PR_LWP_B2_HIST(:,i) = histc(PR_LWP_B2, PR_BINS);
          PR_LWP_B3_HIST(:,i) = histc(PR_LWP_B3, PR_BINS);
          PR_LWP_B4_HIST(:,i) = histc(PR_LWP_B4, PR_BINS);

          CT_HIST(:,i)        = histc(CT, CT_BINS);
          CT_LWP_B1_HIST(:,i) = histc(CT_LWP_B1, CT_BINS);
          CT_LWP_B2_HIST(:,i) = histc(CT_LWP_B2, CT_BINS);
          CT_LWP_B3_HIST(:,i) = histc(CT_LWP_B3, CT_BINS);
          CT_LWP_B4_HIST(:,i) = histc(CT_LWP_B4, CT_BINS);

          % calculate albdeo and form the histograms
          ALB        = ((1 - g) .* CT) ./ (((1 - g) .* CT) + 2);
          ALB_LWP_B1 = ((1 - g) .* CT_LWP_B1) ./ (((1 - g) .* CT_LWP_B1) + 2);
          ALB_LWP_B2 = ((1 - g) .* CT_LWP_B2) ./ (((1 - g) .* CT_LWP_B2) + 2);
          ALB_LWP_B3 = ((1 - g) .* CT_LWP_B3) ./ (((1 - g) .* CT_LWP_B3) + 2);
          ALB_LWP_B4 = ((1 - g) .* CT_LWP_B4) ./ (((1 - g) .* CT_LWP_B4) + 2);
          
          ALB_HIST(:,i)        = histc(ALB, ALB_BINS);
          ALB_LWP_B1_HIST(:,i) = histc(ALB_LWP_B1, ALB_BINS);
          ALB_LWP_B2_HIST(:,i) = histc(ALB_LWP_B2, ALB_BINS);
          ALB_LWP_B3_HIST(:,i) = histc(ALB_LWP_B3, ALB_BINS);
          ALB_LWP_B4_HIST(:,i) = histc(ALB_LWP_B4, ALB_BINS);
        end
 
     
        % output
        OutFile = sprintf('%s/hist_data_%s.h5', Ddir, Case);
        fprintf('Writing: %s\n', OutFile);

        hdf5write(OutFile, '/pcprr_bins',        PR_BINS);
        hdf5write(OutFile, '/pcprr_hist',        PR_HIST,        'WriteMode', 'append');
        hdf5write(OutFile, '/pcprr_lwp_b1_hist', PR_LWP_B1_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/pcprr_lwp_b2_hist', PR_LWP_B2_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/pcprr_lwp_b3_hist', PR_LWP_B3_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/pcprr_lwp_b4_hist', PR_LWP_B4_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, '/cot_bins',        CT_BINS,        'WriteMode', 'append');
        hdf5write(OutFile, '/cot_hist',        CT_HIST,        'WriteMode', 'append');
        hdf5write(OutFile, '/cot_lwp_b1_hist', CT_LWP_B1_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/cot_lwp_b2_hist', CT_LWP_B2_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/cot_lwp_b3_hist', CT_LWP_B3_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/cot_lwp_b4_hist', CT_LWP_B4_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, '/albedo_bins',        ALB_BINS,        'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_hist',        ALB_HIST,        'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_lwp_b1_hist', ALB_LWP_B1_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_lwp_b2_hist', ALB_LWP_B2_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_lwp_b3_hist', ALB_LWP_B3_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_lwp_b4_hist', ALB_LWP_B4_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
        fprintf('\n');
    end
end
