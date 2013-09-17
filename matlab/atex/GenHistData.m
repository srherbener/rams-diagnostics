function [ ] = GenHistData(ConfigFile)
% GenHistData generate data indexed by SST and CCN concentration for histograms

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Tdir = Config.TsavgDir;
    Hdir = './HDF5';

    PrName = 'pcprr';
    PrFprefix = 'pcprr';
    CtName = 'cloud_opt_thick';
    CtFprefix = 'cloud_opt_thick';

    PR_BINS = 10.^(-3:0.02:2); % pcprr bins, even logarithmic spacing from 0.001 to 100.0
    ALB_BINS = 0:0.01:1; % albedo bins
    CT_BINS = 0:300;  % cloud optical thickness bins

    Tstart = 12;
    Tend = 36;

    g = 0.85; % for albedo calc, scattering asymmetry param (cloud droplets interacting with visible light)

    i = 0;
    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;
        PrFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, PrFprefix, Case);
        CtFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, CtFprefix, Case);

        % grab the CCN concentration and SST value from the case name
        CCNstr = regexprep(regexp(Case, 'ccn[0-9]+', 'match'), 'ccn', '');
        SSTstr = regexprep(regexp(Case, 'sst[0-9]+', 'match'), 'sst', '');
        GCCNstr = regexprep(regexp(Case, 'gcn10m[0-9]', 'match'), 'gcn', '');
        GCCNstr = regexprep(GCCNstr, '10m', '1.0e-');

        CCNval = str2num(CCNstr{1});   % regexprep creates a cell array
        SSTval = str2num(SSTstr{1});   % sscanf does not parse leading zeros (eg, '0400') as expected
        GCCNval = str2num(GCCNstr{1});   % sscanf does not parse leading zeros (eg, '0400') as expected

        fprintf('***************************************************************\n');
        fprintf('Generating pdf data:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input precip rate file: %s\n', PrFile);
        fprintf('    Var name: %s\n', PrName);
        fprintf('  Input cloud optical thickness file: %s\n', CtFile);
        fprintf('    Var name: %s\n', PrName);
        fprintf('  Start Time: %d\n', Tstart);
        fprintf('  End Time: %d\n', Tend);
        fprintf('  CCN: %d\n', CCNval);
        fprintf('  SST: %d\n', SSTval);
        fprintf('  GCN: %f\n', GCCNval);
        fprintf('\n');
 
        % Precip rate will be organized as (x,y,t)
        PR = squeeze(hdf5read(PrFile, PrName));
        CT = squeeze(hdf5read(CtFile, CtName));
    
        % Grab time coordinates
        T = hdf5read(PrFile, 't_coords')/3600; % hrs
 
        % Select out the given time step
        T1 = find(T >= Tstart, 1, 'first');
        T2 = find(T <= Tend, 1, 'last');

        % Record data
        i = i + 1;

        % calculate PDF of precip rate over the temporal range, T1 to T2
        PR_SELECT = PR(:,:,T1:T2);
        PR_HIST(:,i) = histc(PR_SELECT(:), PR_BINS);

        % calculate PDFs of cloud thickness and albedo over the temporal range, T1 to T2
        CT_SELECT = CT(:,:,T1:T2);
        CT_HIST(:,i) = histc(CT_SELECT(:), CT_BINS);

        ALB = ((1 - g) .* CT_SELECT) ./ (((1 - g) .* CT_SELECT) + 2);
        ALB_HIST(:,i) = histc(ALB(:), ALB_BINS);

        CCN(i) = CCNval;
        SST(i) = SSTval;
        GCCN(i) = GCCNval;
    end

    % output --> Use REVU format, 4D var, *_coords
    X = 1:i;
    Y = 1;
    Z = 1;
    T = 1;
     
    OutFile = sprintf('%s/hist_data.h5', Ddir);
    fprintf('Writing: %s\n', OutFile);
    hdf5write(OutFile, '/pcprr_hist' ,     PR_HIST);
    hdf5write(OutFile, '/pcprr_bins',      PR_BINS,  'WriteMode', 'append');
    hdf5write(OutFile, '/copt_thick_hist', CT_HIST,  'WriteMode', 'append');
    hdf5write(OutFile, '/copt_thick_bins', CT_BINS,  'WriteMode', 'append');
    hdf5write(OutFile, '/albedo_hist',     ALB_HIST, 'WriteMode', 'append');
    hdf5write(OutFile, '/albedo_bins',     ALB_BINS, 'WriteMode', 'append');
    hdf5write(OutFile, '/ccn',             CCN,      'WriteMode', 'append');
    hdf5write(OutFile, '/sst',             SST,      'WriteMode', 'append');
    hdf5write(OutFile, '/gccn',            GCCN,     'WriteMode', 'append');

    hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
    hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
    hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
    hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
    fprintf('\n');
end
