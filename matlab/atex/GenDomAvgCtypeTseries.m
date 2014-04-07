function [ ] = GenDomAvgCtypeTseries(ConfigFile)
% GenDomAvgCtypeTseries generate time series of domain averages from domain histograms

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;

  FileHeader = 'ATEX PDF data';

  InFprefix = 'hist_data_ctype';
  VarSets = {
    { { 'pcprr'  'pcprr_strnp'  'pcprr_strat'  'pcprr_scmix'  'pcprr_cumul'  } 'pcprr_bins'  'avg_ts_ctype_pcprr'  }
    { { 'cot'    'cot_strnp'    'cot_strat'    'cot_scmix'    'cot_cumul'    } 'cot_bins'    'avg_ts_ctype_cot'    }
    { { 'albedo' 'albedo_strnp' 'albedo_strat' 'albedo_scmix' 'albedo_cumul' } 'albedo_bins' 'avg_ts_ctype_albedo' }
    { { 'cd'     'cd_strnp'     'cd_strat'     'cd_scmix'     'cd_cumul'     } 'cd_bins'     'avg_ts_ctype_cd'     }
    { { 'lwp'    'lwp_strnp'    'lwp_strat'    'lwp_scmix'    'lwp_cumul'    } 'lwp_bins'    'avg_ts_ctype_lwp'    }
    };
  Nset = length(VarSets);

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;
    InFile = sprintf('%s/%s_%s.h5', Ddir, InFprefix, Case);

    fprintf('***************************************************************\n');
    fprintf('Generating domain average data:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('  Input hist file: %s\n', InFile);

    % grab the t coordinates, use these for combinations
    T = squeeze(hdf5read(InFile, 't_coords'));
    Nt = length(T);

    for iset = 1:Nset
      VarList = VarSets{iset}{1};
      BinsName = VarSets{iset}{2};
      OutFprefix = VarSets{iset}{3};
      fprintf('    Variable set:\n');

      % read in the bins and replicate the bins for every time point
      % histc() was used to create histograms which means that BINS
      % are the edges of the histogram bins. The very last bin however
      % contains a count of the sampled data points that exactly match
      % that bin value. Therefore, use the average value between the
      % edges for calculating a domain average for bins 1 through n-1
      % and use the value (instead of the average) for bin n.
      fprintf('      Bins: %s\n', BinsName);
      fprintf('\n');
      BINS = squeeze(hdf5read(InFile, BinsName));       % column vector
      BINAVG = (BINS(1:end-1) + BINS(2:end)) .* 0.5;
      BINAVG = [ BINAVG; BINS(end) ];
      BINS_MAT = repmat(BINAVG, [ 1 Nt ]);

      % Write a header into the file so that a write statement without append mode
      % will be run first which consequently will erase an existing file and replace
      % it with the contents about to be generated here.
      OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);
      hdf5write(OutFile, 'Header', FileHeader); 

      Nvar = length(VarList);
      for ivar = 1:Nvar
        InVarName  = sprintf('%s_hist', VarList{ivar});
        fprintf('      %s\n', InVarName);
        OutAvgName  = sprintf('%s_avg', VarList{ivar});
        OutNptsName = sprintf('%s_npts', VarList{ivar});

        % Record average and number of point for each time step
        % HIST will be organized as (b,t) where b are the bins of the histograms
        HIST = squeeze(hdf5read(InFile, InVarName));
        NPTS = squeeze(sum(HIST, 1));
        SUMS = squeeze(sum((BINS_MAT .* HIST), 1)); 

        AVG = SUMS ./ NPTS;

        % Save as 4D (x,y,z,t) variable
        OutVar = reshape(AVG, [ 1 1 1 Nt ]);
        hdf5write(OutFile, OutAvgName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(NPTS, [ 1 1 1 Nt ]);
        hdf5write(OutFile, OutNptsName, OutVar, 'WriteMode', 'append');
      end
      fprintf('\n');

      % Write out dummy coordinate values to keep ReadSelectXyzt happy
      Xdummy = 1;
      Ydummy = 1;
      Zdummy = 1;

      hdf5write(OutFile, '/x_coords', Xdummy, 'WriteMode', 'append');
      hdf5write(OutFile, '/y_coords', Ydummy, 'WriteMode', 'append');
      hdf5write(OutFile, '/z_coords', Zdummy, 'WriteMode', 'append');
      hdf5write(OutFile, '/t_coords', T,      'WriteMode', 'append');

      fprintf('\n');
    end
    fprintf('\n');
  end
end
