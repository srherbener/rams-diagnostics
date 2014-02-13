function [ ] = GenMomentData(ConfigFile)
% GenMomentData generate w moment and w flux data

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;

    FprefixList = {
      'w_M3'
      'w_speed_flux'
      'w_theta_flux'
      'w_thetav_flux'
      'w_vapor_flux'
      };

    VnameList = {
      'w-w-w'
      'w-speed'
      'w-theta'
      'w-theta_v'
      'w-vapor'
      };

    Tstart = 12;
    Tend = 36;

    i = 0;
%    for icase = 1:length(Config.Cases)
    for icase = 1:1
        Case = Config.Cases(icase).Cname;

      for ivar = 1:length(FprefixList)
        Fname = sprintf('%s/%s_%s.h5', Ddir, FprefixList{ivar}, Case);
        Vname = VnameList{ivar};

        fprintf('***************************************************************\n');
        fprintf('Generating start end times for CF and CD:\n');
        fprintf('  Case: %s\n', Case);
        fprintf('  Input file: %s\n', Fname);
        fprintf('    Var name: %s\n', Vname);
        fprintf('\n');

        TERMS = squeeze(hdf5read(Fname, Vname));
        NPTS = squeeze(hdf5read(Fname, 'num_points'));
        Z = squeeze(hdf5read(Fname, 'z_coords'));
        T = squeeze(hdf5read(Fname, 't_coords')) / 3600;   % hr

        % select and sum up across the time range
        T1 = find(T >= Tstart, 1, 'first');
        T2 = find(T <= Tend, 1, 'last');

        MOMENTS = GenMoments(TERMS, NPTS, T1, T2);

whos
      end
    end
  end


%%%         CfFile = sprintf('%s/%s_%s.h5', Tdir, CfFprefix, Case);
%%%         CdFile = sprintf('%s/%s_%s.h5', Tdir, CdFprefix, Case);
%%%         CotFile = sprintf('%s/%s_%s.h5', Tdir, CotFprefix, Case);
%%%         CotFile_lwp_0p01 = sprintf('%s/%s_%s.h5', Tdir, CotFprefix_lwp_0p01, Case);
%%%         CotFile_lwp_0p10 = sprintf('%s/%s_%s.h5', Tdir, CotFprefix_lwp_0p10, Case);
%%%         CotFile_lwp_1p00 = sprintf('%s/%s_%s.h5', Tdir, CotFprefix_lwp_1p00, Case);
%%%         CotFile_lwp_b1 = sprintf('%s/%s_%s.h5', Tdir, CotFprefix_lwp_b1, Case);
%%%         CotFile_lwp_b2 = sprintf('%s/%s_%s.h5', Tdir, CotFprefix_lwp_b2, Case);
%%%         CotFile_lwp_b3 = sprintf('%s/%s_%s.h5', Tdir, CotFprefix_lwp_b3, Case);
%%%         CotFile_lwp_b4 = sprintf('%s/%s_%s.h5', Tdir, CotFprefix_lwp_b4, Case);
%%%         TpFile = sprintf('%s/%s-%s-AS-1999-02-10-040000-g1.h5', Hdir, TpFprefix, Case);
%%% 
%%%         % grab the CCN concentration and SST value from the case name
%%%         CCNstr = regexprep(regexp(Case, 'ccn[0-9]+', 'match'), 'ccn', '');
%%%         SSTstr = regexprep(regexp(Case, 'sst[0-9]+', 'match'), 'sst', '');
%%%         GCCNstr = regexprep(regexp(Case, 'gcn10m[0-9]', 'match'), 'gcn', '');
%%%         GCCNstr = regexprep(GCCNstr, '10m', '1.0e-');
%%% 
%%%         CCNval = str2num(CCNstr{1});   % regexprep creates a cell array
%%%         SSTval = str2num(SSTstr{1});   % sscanf does not parse leading zeros (eg, '0400') as expected
%%%         GCCNval = str2num(GCCNstr{1});   % sscanf does not parse leading zeros (eg, '0400') as expected
%%% 
%%%         fprintf('  Input cloud depth file: %s\n', CdFile);
%%%         fprintf('    Var name: %s\n', CdName);
%%%         fprintf('  Input cloud optical thickness file: %s\n', CotFile);
%%%         fprintf('  Input cloud optical thickness file: %s\n', CotFile_lwp_0p01);
%%%         fprintf('  Input cloud optical thickness file: %s\n', CotFile_lwp_0p10);
%%%         fprintf('  Input cloud optical thickness file: %s\n', CotFile_lwp_1p00);
%%%         fprintf('  Input cloud optical thickness file: %s\n', CotFile_lwp_b1);
%%%         fprintf('  Input cloud optical thickness file: %s\n', CotFile_lwp_b2);
%%%         fprintf('  Input cloud optical thickness file: %s\n', CotFile_lwp_b3);
%%%         fprintf('  Input cloud optical thickness file: %s\n', CotFile_lwp_b4);
%%%         fprintf('    Var name: %s\n', CotName);
%%%         fprintf('  Input total precip file: %s\n', TpFile);
%%%         fprintf('    Var name: %s\n', TpName);
%%%         fprintf('  Start Time: %d\n', Tstart);
%%%         fprintf('  End Time: %d\n', Tend);
%%%         fprintf('  CCN: %d\n', CCNval);
%%%         fprintf('  SST: %d\n', SSTval);
%%%         fprintf('  GCN: %f\n', GCCNval);
%%%         fprintf('\n');
%%%  
%%%         % Cloud frac and cloud depth will be organized as (t)
%%%         % Total precip will be organized as (x,y,t)
%%%         CF  = squeeze(hdf5read(CfFile, CfName));
%%%         CD  = squeeze(hdf5read(CdFile, CdName));
%%%         COT = squeeze(hdf5read(CotFile, CotName));
%%%         COT_L0P01 = squeeze(hdf5read(CotFile_lwp_0p01, CotName));
%%%         COT_L0P10 = squeeze(hdf5read(CotFile_lwp_0p10, CotName));
%%%         COT_L1P00 = squeeze(hdf5read(CotFile_lwp_1p00, CotName));
%%%         COT_B1 = squeeze(hdf5read(CotFile_lwp_b1, CotName));
%%%         COT_B2 = squeeze(hdf5read(CotFile_lwp_b2, CotName));
%%%         COT_B3 = squeeze(hdf5read(CotFile_lwp_b3, CotName));
%%%         COT_B4 = squeeze(hdf5read(CotFile_lwp_b4, CotName));
%%%         TP  = squeeze(hdf5read(TpFile, TpName));
%%%     
%%%         % Grab time coordinates
%%%         T = hdf5read(CfFile, 't_coords')/3600; % hrs
%%%  
%%%         % Select out the given time step
%%%         T1 = find(T >= Tstart, 1, 'first');
%%%         T2 = find(T <= Tend, 1, 'last');
%%% 
%%%         % Record data
%%%         i = i + 1;
%%% 
%%%         [ CF_Avg(i) CF_Npts(i) ] = CountsToAvg(CF, T1, T2);
%%%         CD_Avg(i) = mean(CD(T1:T2));
%%% 
%%%         [ COT_Avg(i) COT_Npts(i) ] = CountsToAvg(COT, T1, T2);
%%% 
%%%         [ COT_L0P01_Avg(i) COT_L0P01_Npts(i) ] = CountsToAvg(COT_L0P01, T1, T2);
%%%         [ COT_L0P10_Avg(i) COT_L0P10_Npts(i) ] = CountsToAvg(COT_L0P10, T1, T2);
%%%         [ COT_L1P00_Avg(i) COT_L1P00_Npts(i) ] = CountsToAvg(COT_L1P00, T1, T2);
%%% 
%%%         [ COT_B1_Avg(i) COT_B1_Npts(i) ] = CountsToAvg(COT_B1, T1, T2);
%%%         [ COT_B2_Avg(i) COT_B2_Npts(i) ] = CountsToAvg(COT_B2, T1, T2);
%%%         [ COT_B3_Avg(i) COT_B3_Npts(i) ] = CountsToAvg(COT_B3, T1, T2);
%%%         [ COT_B4_Avg(i) COT_B4_Npts(i) ] = CountsToAvg(COT_B4, T1, T2);
%%% 
%%%         % calculate total accumulated precip by subtracting the total domain
%%%         % precip at T1 from that at T2
%%%         TP1 = squeeze(TP(:,:,T1));
%%%         TP2 = squeeze(TP(:,:,T2));
%%%         ACC_PCP(i) = sum(TP2(:)) - sum(TP1(:));
%%% 
%%%         CCN(i) = CCNval;
%%%         SST(i) = SSTval;
%%%         GCCN(i) = GCCNval;
%%%     end
%%% 
%%%     % output --> Use REVU format, 4D var, *_coords
%%%     X = 1:i;
%%%     Y = 1;
%%%     Z = 1;
%%%     T = 1;
%%%      
%%%     OutFile = sprintf('%s/cf_cd_bar_graph_data.h5', Ddir);
%%%     fprintf('Writing: %s\n', OutFile);
%%%     hdf5write(OutFile, '/cf_avg',   CF_Avg);
%%%     hdf5write(OutFile, '/cd_avg',   CD_Avg,   'WriteMode', 'append');
%%%     hdf5write(OutFile, '/cot_avg',  COT_Avg,  'WriteMode', 'append');
%%%     hdf5write(OutFile, '/cot_avg_lwp_0p01',  COT_L0P01_Avg,  'WriteMode', 'append');
%%%     hdf5write(OutFile, '/cot_avg_lwp_0p10',  COT_L0P10_Avg,  'WriteMode', 'append');
%%%     hdf5write(OutFile, '/cot_avg_lwp_1p00',  COT_L1P00_Avg,  'WriteMode', 'append');
%%%     hdf5write(OutFile, '/cot_avg_lwp_b1',  COT_B1_Avg,  'WriteMode', 'append');
%%%     hdf5write(OutFile, '/cot_avg_lwp_b2',  COT_B2_Avg,  'WriteMode', 'append');
%%%     hdf5write(OutFile, '/cot_avg_lwp_b3',  COT_B3_Avg,  'WriteMode', 'append');
%%%     hdf5write(OutFile, '/cot_avg_lwp_b4',  COT_B4_Avg,  'WriteMode', 'append');
%%%     hdf5write(OutFile, '/acc_pcp',  ACC_PCP,  'WriteMode', 'append');
%%%     hdf5write(OutFile, '/ccn',      CCN,      'WriteMode', 'append');
%%%     hdf5write(OutFile, '/sst',      SST,      'WriteMode', 'append');
%%%     hdf5write(OutFile, '/gccn',     GCCN,     'WriteMode', 'append');
%%% 
%%%     hdf5write(OutFile, 'x_coords', X, 'WriteMode', 'append');
%%%     hdf5write(OutFile, 'y_coords', Y, 'WriteMode', 'append');
%%%     hdf5write(OutFile, 'z_coords', Z, 'WriteMode', 'append');
%%%     hdf5write(OutFile, 't_coords', T, 'WriteMode', 'append');
%%%     fprintf('\n');
%%% end
