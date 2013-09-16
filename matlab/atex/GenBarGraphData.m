function [ ] = GenBarGraphData(ConfigFile)
% GenBarGraphData generate data indexed by SST and CCN concentration for bar graphs

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Tdir = Config.TsavgDir;
    Hdir = './HDF5';

    CfName = 'hda_cloud_mask';
    CfFprefix = 'hda_cloud_frac';
    CdName = 'max_cloud_depth';
    CdFprefix = 'max_cdepth';
    TpName = 'totpcp';
    TpFprefix = 'totpcp';

    Tstart = 12;
    Tend = 36;

    i = 0;
    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;
        CfFile = sprintf('%s/%s_%s.h5', Tdir, CfFprefix, Case);
        CdFile = sprintf('%s/%s_%s.h5', Tdir, CdFprefix, Case);
        TpFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, TpFprefix, Case);

        % grab the CCN concentration and SST value from the case name
        CCNstr = regexprep(regexp(Case, 'ccn[0-9]+', 'match'), 'ccn', '');
        SSTstr = regexprep(regexp(Case, 'sst[0-9]+', 'match'), 'sst', '');
        GCCNstr = regexprep(regexp(Case, 'gcn10m[0-9]', 'match'), 'gcn', '');
        GCCNstr = regexprep(GCCNstr, '10m', '1.0e-');

        CCNval = str2num(CCNstr{1});   % regexprep creates a cell array
        SSTval = str2num(SSTstr{1});   % sscanf does not parse leading zeros (eg, '0400') as expected
        GCCNval = str2num(GCCNstr{1});   % sscanf does not parse leading zeros (eg, '0400') as expected

        fprintf('***************************************************************\n');
        fprintf('Generating start end times for CF and CD:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input cloud fraction file: %s\n', CfFile);
        fprintf('    Var name: %s\n', CfName);
        fprintf('  Input cloud depth file: %s\n', CdFile);
        fprintf('    Var name: %s\n', CdName);
        fprintf('  Input total precip file: %s\n', TpFile);
        fprintf('    Var name: %s\n', TpName);
        fprintf('  Start Time: %d\n', Tstart);
        fprintf('  End Time: %d\n', Tend);
        fprintf('  CCN: %d\n', CCNval);
        fprintf('  SST: %d\n', SSTval);
        fprintf('  GCN: %f\n', GCCNval);
        fprintf('\n');
 
        % Cloud frac and cloud depth will be organized as (t)
        % Total precip will be organized as (x,y,t)
        CF = squeeze(hdf5read(CfFile, CfName));
        CD = squeeze(hdf5read(CdFile, CdName));
        TP = squeeze(hdf5read(TpFile, TpName));
    
        % Grab time coordinates
        T = hdf5read(CfFile, 't_coords')/3600; % hrs
 
        % Select out the given time step
        T1 = find(T >= Tstart, 1, 'first');
        T2 = find(T <= Tend, 1, 'last');

        % Record data
        i = i + 1;

        CF_Start(i) = CF(T1);
        CD_Start(i) = CD(T1);

        CF_End(i)   = CF(T2);
        CD_End(i)   = CD(T2);

        CF_Avg(i) = mean(CF(T1:T2));
        CD_Avg(i) = mean(CD(T1:T2));

        % calculate total accumulated precip by subtracting the total domain
        % precip at T1 from that at T2
        TP1 = squeeze(TP(:,:,T1));
        TP2 = squeeze(TP(:,:,T2));
        ACC_PCP(i) = sum(TP2(:)) - sum(TP1(:));

        CCN(i) = CCNval;
        SST(i) = SSTval;
        GCCN(i) = GCCNval;
    end

    % output --> Use REVU format, 4D var, *_coords
    X = 1:i;
    Y = 1;
    Z = 1;
    T = 1;
     
    OutFile = sprintf('%s/cf_cd_bar_graph_data.h5', Ddir);
    fprintf('Writing: %s\n', OutFile);
    hdf5write(OutFile, '/cf_start', CF_Start);
    hdf5write(OutFile, '/cf_end',   CF_End,   'WriteMode', 'append');
    hdf5write(OutFile, '/cf_avg',   CF_Avg,   'WriteMode', 'append');
    hdf5write(OutFile, '/cd_start', CD_Start, 'WriteMode', 'append');
    hdf5write(OutFile, '/cd_end',   CD_End,   'WriteMode', 'append');
    hdf5write(OutFile, '/cd_avg',   CD_Avg,   'WriteMode', 'append');
    hdf5write(OutFile, '/acc_pcp',  ACC_PCP,  'WriteMode', 'append');
    hdf5write(OutFile, '/ccn',      CCN,      'WriteMode', 'append');
    hdf5write(OutFile, '/sst',      SST,      'WriteMode', 'append');
    hdf5write(OutFile, '/gccn',     GCCN,     'WriteMode', 'append');

    hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
    hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
    hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
    hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
    fprintf('\n');
end
