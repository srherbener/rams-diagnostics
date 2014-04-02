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
    CmName = 'cloud_mask';     % 1 for cloudy column, 0 for clear column
    CmFprefix = 'cloud_mask';

    PR_BINS = 10.^(-5:0.02:2);    % pcprr bins, even log spacing from 1e-5 to 100.0
    ALB_BINS = 0:0.01:1;          % albedo bins
    CT_BINS = 10.^(-4:0.02:2.6);  % cloud optical thickness bins, even log spacing from 1e-4 to ~400
    CD_BINS = 0:100:4000;         % cloud depth bins
    LWP_BINS = 10.^(-4:0.02:1.4); % lwp bins, even log spacing from 1e-4 to ~25

    Nprb = length(PR_BINS);
    Nab  = length(ALB_BINS);
    Nctb = length(CT_BINS);
    Ncd  = length(CD_BINS);
    Nlwp = length(LWP_BINS);

    g = 0.85; % for albedo calc, scattering asymmetry param (cloud droplets interacting with visible light)

    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;
        PrFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, PrFprefix, Case);
        CtFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, CtFprefix, Case);
        CdFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, CdFprefix, Case);
        CmFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, CmFprefix, Case);
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
        fprintf('  Input cloud mask file: %s\n', CmFile);
        fprintf('    Var name: %s\n', CmName);
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
 
        CM_DS = ncgeodataset(CmFile);
        CM_VAR = CM_DS.geovariable(CmName);
 
        LWP_DS = ncgeodataset(LwpFile);
        LWP_VAR = LWP_DS.geovariable(LwpName);

        % Output will be of the form (Nb, Nt) where Nb is the number of bins in the histogram
        PR_HIST       = zeros([ Nprb Nt ]);
        PR_STRAT_HIST = zeros([ Nprb Nt ]);
        PR_SCMIX_HIST = zeros([ Nprb Nt ]);
        PR_CUMUL_HIST = zeros([ Nprb Nt ]);

        CT_HIST       = zeros([ Nctb Nt ]);
        CT_STRAT_HIST = zeros([ Nctb Nt ]);
        CT_SCMIX_HIST = zeros([ Nctb Nt ]);
        CT_CUMUL_HIST = zeros([ Nctb Nt ]);

        CD_HIST       = zeros([ Ncd Nt ]);
        CD_STRAT_HIST = zeros([ Ncd Nt ]);
        CD_SCMIX_HIST = zeros([ Ncd Nt ]);
        CD_CUMUL_HIST = zeros([ Ncd Nt ]);

        LWP_HIST       = zeros([ Nlwp Nt ]);
        LWP_STRAT_HIST = zeros([ Nlwp Nt ]);
        LWP_SCMIX_HIST = zeros([ Nlwp Nt ]);
        LWP_CUMUL_HIST = zeros([ Nlwp Nt ]);

        ALB_HIST       = zeros([ Nab Nt ]);
        ALB_STRAT_HIST = zeros([ Nab Nt ]);
        ALB_SCMIX_HIST = zeros([ Nab Nt ]);
        ALB_CUMUL_HIST = zeros([ Nab Nt ]);

        % Walk through each time step forming the histogram counts for that time step. These
        % counts can then be combined in different ways during a post processing step.
        for i = 1:Nt
          % read in the 2D vars and reshape them into vectors
          % RAMS uses the borders for boundaries so trim these off
          PR  = squeeze(PR_VAR.data(i,:,:));
          CT  = squeeze(CT_VAR.data(i,:,:));
          CD  = squeeze(CD_VAR.data(i,:,:));
          CM  = squeeze(CM_VAR.data(i,:,:));
          LWP = squeeze(LWP_VAR.data(i,:,:));

          PR = reshape(PR(2:end-1,2:end-1), [ 1  Ntrim ]);
          CT = reshape(CT(2:end-1,2:end-1), [ 1  Ntrim ]);
          CD = reshape(CD(2:end-1,2:end-1), [ 1  Ntrim ]);
          CM = reshape(CM(2:end-1,2:end-1), [ 1  Ntrim ]);
          LWP = reshape(LWP(2:end-1,2:end-1), [ 1  Ntrim ]);

          % Form the selection based on precip rate
          %   Exclude clear columns (CM == 0) for lightly precipitating stratiform category
          %   Select only on precip rate for other categories (assume you have clouds)
          %
          %    stratiform with light or no precip:  PR < 0.001 and CM == 1
          %    stratiform: 0.001 <= PR <= 1.0
          %    mix:        1.0 < PR < 10.0            
          %    cumuliform: PR >= 10.0
          %
          LWP1 = 0.01;
          LWP2 = 0.10;
          CD1 = 50;
          CD2 = 1500;
          PR1 = 0.001;
          PR2 = 1.0;
          PR3 = 10.0;

          SNP_SELECT = (PR < PR1 & CM > 0.5);
          S_SELECT = (PR >= PR1 & PR <= PR2);
          M_SELECT = (PR >  PR2 & PR <  PR3);
          C_SELECT = (PR >  PR3);

          PR_STRNP = PR(SNP_SELECT);
          PR_STRAT = PR(S_SELECT);
          PR_SCMIX = PR(M_SELECT);
          PR_CUMUL = PR(C_SELECT);

          CT_STRNP = CT(SNP_SELECT);
          CT_STRAT = CT(S_SELECT);
          CT_SCMIX = CT(M_SELECT);
          CT_CUMUL = CT(C_SELECT);

          CD_STRNP = CD(SNP_SELECT);
          CD_STRAT = CD(S_SELECT);
          CD_SCMIX = CD(M_SELECT);
          CD_CUMUL = CD(C_SELECT);

          LWP_STRNP = LWP(SNP_SELECT);
          LWP_STRAT = LWP(S_SELECT);
          LWP_SCMIX = LWP(M_SELECT);
          LWP_CUMUL = LWP(C_SELECT);

          % form the histograms
          PR_HIST(:,i)       = histc(PR, PR_BINS);
          PR_STRNP_HIST(:,i) = histc(PR_STRNP, PR_BINS);
          PR_STRAT_HIST(:,i) = histc(PR_STRAT, PR_BINS);
          PR_SCMIX_HIST(:,i) = histc(PR_SCMIX, PR_BINS);
          PR_CUMUL_HIST(:,i) = histc(PR_CUMUL, PR_BINS);

          CT_HIST(:,i)       = histc(CT, CT_BINS);
          CT_STRNP_HIST(:,i) = histc(CT_STRNP, CT_BINS);
          CT_STRAT_HIST(:,i) = histc(CT_STRAT, CT_BINS);
          CT_SCMIX_HIST(:,i) = histc(CT_SCMIX, CT_BINS);
          CT_CUMUL_HIST(:,i) = histc(CT_CUMUL, CT_BINS);

          CD_HIST(:,i)       = histc(CD, CD_BINS);
          CD_STRNP_HIST(:,i) = histc(CD_STRNP, CD_BINS);
          CD_STRAT_HIST(:,i) = histc(CD_STRAT, CD_BINS);
          CD_SCMIX_HIST(:,i) = histc(CD_SCMIX, CD_BINS);
          CD_CUMUL_HIST(:,i) = histc(CD_CUMUL, CD_BINS);

          LWP_HIST(:,i)       = histc(LWP, LWP_BINS);
          LWP_STRNP_HIST(:,i) = histc(LWP_STRNP, LWP_BINS);
          LWP_STRAT_HIST(:,i) = histc(LWP_STRAT, LWP_BINS);
          LWP_SCMIX_HIST(:,i) = histc(LWP_SCMIX, LWP_BINS);
          LWP_CUMUL_HIST(:,i) = histc(LWP_CUMUL, LWP_BINS);

          % calculate albdeo and form the histograms
          ALB       = ((1 - g) .* CT      ) ./ (((1 - g) .* CT      ) + 2);
          ALB_STRNP = ((1 - g) .* CT_STRNP) ./ (((1 - g) .* CT_STRNP) + 2);
          ALB_STRAT = ((1 - g) .* CT_STRAT) ./ (((1 - g) .* CT_STRAT) + 2);
          ALB_SCMIX = ((1 - g) .* CT_SCMIX) ./ (((1 - g) .* CT_SCMIX) + 2);
          ALB_CUMUL = ((1 - g) .* CT_CUMUL) ./ (((1 - g) .* CT_CUMUL) + 2);
          
          ALB_HIST(:,i)       = histc(ALB, ALB_BINS);
          ALB_STRNP_HIST(:,i) = histc(ALB_STRNP, ALB_BINS);
          ALB_STRAT_HIST(:,i) = histc(ALB_STRAT, ALB_BINS);
          ALB_SCMIX_HIST(:,i) = histc(ALB_SCMIX, ALB_BINS);
          ALB_CUMUL_HIST(:,i) = histc(ALB_CUMUL, ALB_BINS);
        end
 
     
        % output
        OutFile = sprintf('%s/hist_data_ctype_%s.h5', Ddir, Case);
        fprintf('Writing: %s\n', OutFile);

        hdf5write(OutFile, '/pcprr_bins',       PR_BINS);
        hdf5write(OutFile, '/pcprr_hist',       PR_HIST,       'WriteMode', 'append');
        hdf5write(OutFile, '/pcprr_strnp_hist', PR_STRNP_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/pcprr_strat_hist', PR_STRAT_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/pcprr_scmix_hist', PR_SCMIX_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/pcprr_cumul_hist', PR_CUMUL_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, '/cot_bins',       CT_BINS,       'WriteMode', 'append');
        hdf5write(OutFile, '/cot_hist',       CT_HIST,       'WriteMode', 'append');
        hdf5write(OutFile, '/cot_strnp_hist', CT_STRNP_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/cot_strat_hist', CT_STRAT_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/cot_scmix_hist', CT_SCMIX_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/cot_cumul_hist', CT_CUMUL_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, '/cd_bins',       CD_BINS,       'WriteMode', 'append');
        hdf5write(OutFile, '/cd_hist',       CD_HIST,       'WriteMode', 'append');
        hdf5write(OutFile, '/cd_strnp_hist', CD_STRNP_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/cd_strat_hist', CD_STRAT_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/cd_scmix_hist', CD_SCMIX_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/cd_cumul_hist', CD_CUMUL_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, '/albedo_bins',       ALB_BINS,       'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_hist',       ALB_HIST,       'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_strnp_hist', ALB_STRNP_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_strat_hist', ALB_STRAT_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_scmix_hist', ALB_SCMIX_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/albedo_cumul_hist', ALB_CUMUL_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, '/lwp_bins',       LWP_BINS,       'WriteMode', 'append');
        hdf5write(OutFile, '/lwp_hist',       LWP_HIST,       'WriteMode', 'append');
        hdf5write(OutFile, '/lwp_strnp_hist', LWP_STRNP_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/lwp_strat_hist', LWP_STRAT_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/lwp_scmix_hist', LWP_SCMIX_HIST, 'WriteMode', 'append');
        hdf5write(OutFile, '/lwp_cumul_hist', LWP_CUMUL_HIST, 'WriteMode', 'append');

        hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
        fprintf('\n');
    end
end
