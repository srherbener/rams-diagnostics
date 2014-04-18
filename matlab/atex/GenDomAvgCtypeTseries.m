function [ ] = GenDomAvgCtypeTseries(ConfigFile)
% GenDomAvgCtypeTseries generate time series of domain averages from domain histograms

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;
  Tdir = Config.TsavgDir;

  FileHeader = 'ATEX HDA time series data';

  VarSets = {
    { { 'hda_pcprr'  'hda_pcprr_strnp'  'hda_pcprr_strat'  'hda_pcprr_cumul'  'hda_pcprr_all_cld'  } 'hda_pcprr'           'avg_ts_ctype_pcprr'  }
    { { 'hda_cot'    'hda_cot_strnp'    'hda_cot_strat'    'hda_cot_cumul'    'hda_cot_all_cld'    } 'hda_cloud_opt_thick' 'avg_ts_ctype_cot'    }
    { { 'hda_cdepth' 'hda_cdepth_strnp' 'hda_cdepth_strat' 'hda_cdepth_cumul' 'hda_cdepth_all_cld' } 'hda_cloud_depth'     'avg_ts_ctype_cdepth' }
    { { 'hda_lwp'    'hda_lwp_strnp'    'hda_lwp_strat'    'hda_lwp_cumul'    'hda_lwp_all_cld'    } 'hda_vint_cond'       'avg_ts_ctype_lwp'    }
    };
  Nset = length(VarSets);

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;
    fprintf('***************************************************************\n');
    fprintf('Generating domain average time series:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for iset = 1:Nset
      InFprefixList = VarSets{iset}{1};
      InVarName     = VarSets{iset}{2};
      OutFprefix    = VarSets{iset}{3};

      % Write a header into the file so that a write statement without append mode
      % will be run first which consequently will erase an existing file and replace
      % it with the contents about to be generated here.
      OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);
      fprintf('  Writing: %s\n', OutFile);
      fprintf('\n');
      hdf5write(OutFile, 'Header', FileHeader); 
 
      Nfile = length(InFprefixList);
      for ifile = 1:Nfile
        InFprefix = InFprefixList{ifile};
        InFile = sprintf('%s/%s_%s.h5', Tdir, InFprefix, Case);
        fprintf('    Reading: %s (%s)\n', InFile, InVarName);

        if (ifile == 1)
          % pass the t coordinates on to the output file
          T = squeeze(hdf5read(InFile, 't_coords'));
          Nt = length(T);
        end

        HDATA = squeeze(hdf5read(InFile, InVarName)); % hdata is (2,t): (1,:) --> SUMS, (2,:) --> COUNTS
        SUMS = squeeze(HDATA(1,:));
        NPTS = squeeze(HDATA(2,:));

        AVG = SUMS ./ NPTS;

        % Use the input file prefix (unique names) for the basis of the output variable names
        OutVname = regexprep(InFprefix, '^hda_', '');
        OutAvgName  = sprintf('%s_avg', OutVname);
        OutNptsName = sprintf('%s_npts', OutVname);

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
    end
    fprintf('\n');
  end
end
