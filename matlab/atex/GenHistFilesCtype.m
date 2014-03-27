function [ ] = GenHistFilesCtype(ConfigFile)
% GenHistDataCtype generate histograms for vars based on cloud types
%
% Two cloud types:
%   stratiform
%   cumuliform

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
    CdName = 'cloud_depth';
    CdFprefix = 'cloud_depth';
    CfName = 'cloud_frac';
    CfFprefix = 'cloud_frac';

    PR_BINS = 10.^(-3:0.02:2); % pcprr bins, even logarithmic spacing from 0.001 to 100.0
    ALB_BINS = 0:0.01:1; % albedo bins
    CT_BINS = 0:300;  % cloud optical thickness bins
    CF_BINS = 0:0.01:1;  % cloud fraction bins

    Nprb = length(PR_BINS);
    Nab  = length(ALB_BINS);
    Nctb  = length(CT_BINS);
    Ncfb  = length(CF_BINS);

    g = 0.85; % for albedo calc, scattering asymmetry param (cloud droplets interacting with visible light)

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;
        PrFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, PrFprefix, Case);
        CtFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, CtFprefix, Case);
        CdFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, CdFprefix, Case);
        CfFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, CfFprefix, Case);
        LwpFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, LwpFprefix, Case);

        fprintf('***************************************************************\n');
        fprintf('Generating hist data:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input precip rate file: %s\n', PrFile);
        fprintf('    Var name: %s\n', PrName);
        fprintf('  Input cloud optical thickness file: %s\n', CtFile);
        fprintf('    Var name: %s\n', CtName);
        fprintf('  Input cloud depth file: %s\n', CdFile);
        fprintf('    Var name: %s\n', CdName);
        fprintf('  Input cloud fraction file: %s\n', CfFile);
        fprintf('    Var name: %s\n', CfName);
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
 
        CD_DS = ncgeodataset(CdFile);
        CD_VAR = CD_DS.geovariable(CdName);
 
        CF_DS = ncgeodataset(CfFile);
        CF_VAR = CF_DS.geovariable(CfName);
 
        LWP_DS = ncgeodataset(LwpFile);
        LWP_VAR = LWP_DS.geovariable(LwpName);

        % Output will be of the form (Nb, Nt) where Nb is the number of bins in the histogram
        PR_HIST        = zeros([ Nprb Nt ]);
        PR_STRAT_HIST = zeros([ Nprb Nt ]);
        PR_CUMUL_HIST = zeros([ Nprb Nt ]);

        CT_HIST        = zeros([ Nctb Nt ]);
        CT_STRAT_HIST = zeros([ Nctb Nt ]);
        CT_CUMUL_HIST = zeros([ Nctb Nt ]);

        CF_HIST        = zeros([ Ncfb Nt ]);
        CF_STRAT_HIST = zeros([ Ncfb Nt ]);
        CF_CUMUL_HIST = zeros([ Ncfb Nt ]);

        ALB_HIST        = zeros([ Nab Nt ]);
        ALB_STRAT_HIST = zeros([ Nab Nt ]);
        ALB_CUMUL_HIST = zeros([ Nab Nt ]);

        % Walk through each time step forming the histogram counts for that time step. These
        % counts can then be combined in different ways during a post processing step.
        for i = 1:Nt
          % read in the 2D vars and reshape them into vectors
          % RAMS uses the borders for boundaries so trim these off
          PR  = squeeze(PR_VAR.data(i,:,:));
          CT  = squeeze(CT_VAR.data(i,:,:));
          CD  = squeeze(CD_VAR.data(i,:,:));
          CF  = squeeze(CF_VAR.data(i,:,:));
          LWP = squeeze(LWP_VAR.data(i,:,:));

          PR = reshape(PR(2:end-1,2:end-1), [ 1  Ntrim ]);
          CT = reshape(CT(2:end-1,2:end-1), [ 1  Ntrim ]);
          CD = reshape(CD(2:end-1,2:end-1), [ 1  Ntrim ]);
          CF = reshape(CF(2:end-1,2:end-1), [ 1  Ntrim ]);
          LWP = reshape(LWP(2:end-1,2:end-1), [ 1  Ntrim ]);

          % Divide PR and CT into 2 bins according to LWP (mm) and cdepth (m)
          %    stratiform: 0.01 <= LWP <= 0.20  AND 50 <= CD <= 1500
          %    cumuliform: LWP > 0.2  OR CD > 1500
          S_SELECT = (LWP >= 0.01 & LWP <= 0.2 & CD >= 50 & CD <= 1500);
          C_SELECT = (LWP > 0.2 | CD > 1500);

          PR_STRAT = PR(S_SELECT);
          PR_CUMUL = PR(C_SELECT);

          CT_STRAT = CT(S_SELECT);
          CT_CUMUL = CT(C_SELECT);

          CF_STRAT = CF(S_SELECT);
          CF_CUMUL = CF(C_SELECT);

          % form the histograms
          PR_HIST(:,i)        = histc(PR, PR_BINS);
          PR_STRAT_HIST(:,i) = histc(PR_STRAT, PR_BINS);
          PR_CUMUL_HIST(:,i) = histc(PR_CUMUL, PR_BINS);

          CT_HIST(:,i)        = histc(CT, CT_BINS);
          CT_STRAT_HIST(:,i) = histc(CT_STRAT, CT_BINS);
          CT_CUMUL_HIST(:,i) = histc(CT_CUMUL, CT_BINS);

          CF_HIST(:,i)        = histc(CF, CF_BINS);
          CF_STRAT_HIST(:,i) = histc(CF_STRAT, CF_BINS);
          CF_CUMUL_HIST(:,i) = histc(CF_CUMUL, CF_BINS);

          % calculate albdeo and form the histograms
          ALB        = ((1 - g) .* CT) ./ (((1 - g) .* CT) + 2);
          ALB_STRAT = ((1 - g) .* CT_STRAT) ./ (((1 - g) .* CT_STRAT) + 2);
          ALB_CUMUL = ((1 - g) .* CT_CUMUL) ./ (((1 - g) .* CT_CUMUL) + 2);
          
          ALB_HIST(:,i)        = histc(ALB, ALB_BINS);
          ALB_STRAT_HIST(:,i) = histc(ALB_STRAT, ALB_BINS);
          ALB_CUMUL_HIST(:,i) = histc(ALB_CUMUL, ALB_BINS);
        end
 
     
        % output
        OutFile = sprintf('%s/hist_data_ctype_%s.h5', Ddir, Case);
        fprintf('Writing: %s\n', OutFile);

        hdf5write(OutFile, '/pcprr_bins',        PR_BINS);
        hdf5write(OutFile, '/pcprr_hist',        PR_HIST,        'WriteMode', 'append');
        hdf5write(OutFile, '/pcprr_strat_hist', PR_STRAT_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/pcprr_cumul_hist', PR_CUMUL_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, '/cot_bins',        CT_BINS,        'WriteMode', 'append');
        hdf5write(OutFile, '/cot_hist',        CT_HIST,        'WriteMode', 'append');
        hdf5write(OutFile, '/cot_strat_hist', CT_STRAT_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/cot_cumul_hist', CT_CUMUL_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, '/cf_bins',        CF_BINS,        'WriteMode', 'append');
        hdf5write(OutFile, '/cf_hist',        CF_HIST,        'WriteMode', 'append');
        hdf5write(OutFile, '/cf_strat_hist', CF_STRAT_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/cf_cumul_hist', CF_CUMUL_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, '/albedo_bins',        ALB_BINS,        'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_hist',        ALB_HIST,        'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_strat_hist', ALB_STRAT_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_cumul_hist', ALB_CUMUL_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
        fprintf('\n');
    end
end
