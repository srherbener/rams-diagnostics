function [ ] = GenPdfFiles(ConfigFile)
% GenPdfFiles generate pdfs from hist data

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;

  FileHeader = 'ATEX PDF data';

  InFprefix = 'hist_data';
  VarSets = {
    { { 'pcprr'  'pcprr_lwp_b1'  'pcprr_lwp_b2'  'pcprr_lwp_b3'  'pcprr_lwp_b4'  } 'pcprr_bins'  'pdf_pcprr' }
    { { 'cot'    'cot_lwp_b1'    'cot_lwp_b2'    'cot_lwp_b3'    'cot_lwp_b4'    } 'cot_bins'    'pdf_cot' }
    { { 'albedo' 'albedo_lwp_b1' 'albedo_lwp_b2' 'albedo_lwp_b3' 'albedo_lwp_b4' } 'albedo_bins' 'pdf_albedo' }
    { { 'cd' 'cd_lwp_b1' 'cd_lwp_b2' 'cd_lwp_b3' 'cd_lwp_b4' } 'cd_bins' 'pdf_cd' }
    { { 'cf' 'cf_lwp_b1' 'cf_lwp_b2' 'cf_lwp_b3' 'cf_lwp_b4' } 'cf_bins' 'pdf_cf' }
    };
  Nset = length(VarSets);

  TimeSelects = {
    { 12    36    'TALL'   }      
    { 12    13    'TSTART' }
    { 23.5  24.5  'TMID'   }
    { 35    36    'TEND'   }
    };
  Ntsel = length(TimeSelects);

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;
    InFile = sprintf('%s/%s_%s.h5', Ddir, InFprefix, Case);

    fprintf('***************************************************************\n');
    fprintf('Generating pdf data:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('  Input hist file: %s\n', InFile);

    % grab the t coordinates, use these for combinations
    T = squeeze(hdf5read(InFile, 't_coords'))/3600; % hours

    for iset = 1:Nset
      VarList = VarSets{iset}{1};
      BinsName = VarSets{iset}{2};
      OutFprefix = VarSets{iset}{3};
      fprintf('    Variable set:\n');

      % Write a header into the file so that a write statement without append mode
      % will be run first which consequently will erase an existing file and replace
      % it with the contents about to be generated here.
      OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);
      hdf5write(OutFile, 'Header', FileHeader); 

      Nvar = length(VarList);
      for ivar = 1:Nvar
        InVarName  = sprintf('%s_hist', VarList{ivar});
        fprintf('      %s\n', InVarName);
        OutVarName  = sprintf('%s_pdf', VarList{ivar});

        for its = 1:Ntsel
          Tstart = TimeSelects{its}{1};
          Tend   = TimeSelects{its}{2};
          Tname  = TimeSelects{its}{3};

          T1 = find(T >= Tstart, 1, 'first');
          T2 = find(T <= Tend,   1, 'last');

          % sum up bins across selected times, and convert counts to PDF
          % HIST will be organized as (b,t) where b are the bins of the histograms
          HIST = squeeze(hdf5read(InFile, InVarName));
          HIST_SELECT = squeeze(sum(HIST(:,T1:T2), 2)); % reduce to single histogram (vector)
          PDF = HIST_SELECT ./ sum(HIST_SELECT);

          % With the LWP selection it is possible to get no data points selected (all counts
          % equal to zero in histogram). When this happens, get nans in the PDF since the sum
          % is zero. Ie, all the PDF entries are 0/0 --> nan. Leave as nan to signal downstream
          % processes that the PDF is empty.

          OutName = sprintf('%s_%s', OutVarName, Tname);
          hdf5write(OutFile, OutName, PDF, 'WriteMode', 'append');
        end
      end
      fprintf('\n');

      % Save bins as 'x_coords'
      % Write out dummy coordinate values to keep ReadSelectXyzt happy
      fprintf('      Bins: %s\n', BinsName);
      BINS = squeeze(hdf5read(InFile, BinsName));
 
      Ydummy = 1;
      Zdummy = 1;
      Tdummy = 1;

      hdf5write(OutFile, '/x_coords', BINS,   'WriteMode', 'append');
      hdf5write(OutFile, '/y_coords', Ydummy, 'WriteMode', 'append');
      hdf5write(OutFile, '/z_coords', Zdummy, 'WriteMode', 'append');
      hdf5write(OutFile, '/t_coords', Tdummy, 'WriteMode', 'append');

      fprintf('\n');
    end
    fprintf('\n');
  end
end
