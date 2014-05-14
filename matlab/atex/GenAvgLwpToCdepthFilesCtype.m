function [ ] = GenAvgLwpToCdepthFilesCtype(ConfigFile)
% GenAvgLwpToCdepthFilesCtype special case: generate LWP to Cdepth ratio

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Tdir = Config.TsavgDir;
  Ddir = Config.DiagDir;

  FileHeader = 'ATEX PDF data';

  LwpSet    = { 'hda_vint_cond'       { 'lwp'    'lwp_strnp'    'lwp_strat'    'lwp_cumul'    'lwp_all_cld'    'lwp_stall'    } };
  CdepthSet = { 'hda_cloud_depth'     { 'cdepth' 'cdepth_strnp' 'cdepth_strat' 'cdepth_cumul' 'cdepth_all_cld' 'cdepth_stall' } };

  OutFprefix = 'avg_ctype_lwp2cdepth';

  TimeSelects = {
    { 12    36    'TALL'   }      
    { 12    13    'TSTART' }
    { 23.5  24.5  'TMID'   }
    { 35    36    'TEND'   }
    };
  Ntsel = length(TimeSelects);

  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;

    fprintf('***************************************************************\n');
    fprintf('Generating avg data:\n');
    fprintf('  Case: %s\n', Case);

    LwpInVname        = LwpSet{1};
    LwpInFprefixList  = LwpSet{2};

    CdepthInVname        = CdepthSet{1};
    CdepthInFprefixList  = CdepthSet{2};

    fprintf('    Input files:\n');

    % Write a header into the file so that a write statement without append mode
    % will be run first which consequently will erase an existing file and replace
    % it with the contents about to be generated here.
    OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);
    hdf5write(OutFile, 'Header', FileHeader); 

    Nvar = length(LwpInFprefixList);
    for ivar = 1:Nvar
      LwpInFprefix = LwpInFprefixList{ivar};
      LwpInVarName = LwpInVname;

      CdepthInFprefix = CdepthInFprefixList{ivar};
      CdepthInVarName = CdepthInVname;

      OutVarName = regexprep(LwpInFprefix, '^lwp', 'lwp2cdepth');

      LwpInFile = sprintf('%s/hda_%s_%s.h5', Tdir, LwpInFprefix , Case);
      CdepthInFile = sprintf('%s/hda_%s_%s.h5', Tdir, CdepthInFprefix , Case);
      fprintf('      LWP: %s (%s)\n', LwpInFile, LwpInVarName);
      fprintf('      CD:  %s (%s)\n', CdepthInFile, CdepthInVarName);

      % grab the hda data for LWP and Cdepth
      % hda data is: (y,t) where y is size 2
      %   y(1) -> Sum
      %   y(2) -> N 
      %
      % want to end up with:
      %
      %    (Sum(LWP)/Nlwp) / (Sum(CD)/Ncd) = (Sum(LWP)/Sum(CD)) / (Nlpw/Ncd)
      %
      % so just divide the sums and Ns and pass it on through to the rest
      % of the code (that divides sum/n
      %
      LWP_HDA  = squeeze(hdf5read(LwpInFile, LwpInVarName));
      CD_HDA  = squeeze(hdf5read(CdepthInFile, CdepthInVarName));
      T = squeeze(hdf5read(LwpInFile, 't_coords'))/3600; % hours
      Nt = length(T);

      HDA = zeros([ 2 Nt ]);
      SUMS = squeeze(LWP_HDA(1,:)) ./ squeeze(CD_HDA(1,:));
      NPTS = squeeze(LWP_HDA(2,:)) ./ squeeze(CD_HDA(2,:));

      HDA(1,:) = SUMS;
      HDA(2,:) = NPTS;

      OutAvgName  = sprintf('%s_avg', OutVarName);
      OutNptsName = sprintf('%s_npts', OutVarName);

      for its = 1:Ntsel
        Tstart = TimeSelects{its}{1};
        Tend   = TimeSelects{its}{2};
        Tname  = TimeSelects{its}{3};

        T1 = find(T >= Tstart, 1, 'first');
        T2 = find(T <= Tend,   1, 'last');

        % sum up bins across selected times, and convert HDA counts to an average
        [ AVG NPTS ] = CountsToAvg(HDA, T1, T2);

        % With the data selection it is possible to get no data points selected (all counts
        % equal to zero in hda file). When this happens, get nans in Avg since the sum
        % is zero. Ie, the AVG entry is 0/0 --> nan. Leave as nan to signal downstream
        % processes that the AVG is empty.

        OutName = sprintf('%s_%s', OutAvgName, Tname);
        hdf5write(OutFile, OutName, AVG, 'WriteMode', 'append');

        OutName = sprintf('%s_%s', OutNptsName, Tname);
        hdf5write(OutFile, OutName, NPTS, 'WriteMode', 'append');
      end
    end
    fprintf('\n');

    % output coords so that ReadSelectXyzt can handle the output files
    Xdummy = 1;
    Ydummy = 1;
    Zdummy = 1;
    Tdummy = 1;

    fprintf('    Writing: %s\n', OutFile);

    hdf5write(OutFile, '/x_coords', Xdummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', Ydummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Zdummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', Tdummy, 'WriteMode', 'append');

    fprintf('\n');
  end
end
