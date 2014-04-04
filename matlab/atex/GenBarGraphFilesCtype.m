function [ ] = GenBarGraphFilesCtype(ConfigFile)
% GenBarGraphFilesCtype generate data bar graphs from cloud type selected data

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Tdir = Config.TsavgDir;

    VarSets = {
      { % cloud optical thickness
        'avg_ctype_cot'
        {
           'cot_TSTART'
           'cot_TMID'
           'cot_TEND'
           'cot_TALL'
%           'cot_strnp_TSTART'
%           'cot_strnp_TMID'
%           'cot_strnp_TEND'
%           'cot_strnp_TALL'
           'cot_strat_TSTART'
           'cot_strat_TMID'
           'cot_strat_TEND'
           'cot_strat_TALL'
           'cot_scmix_TSTART'
           'cot_scmix_TMID'
           'cot_scmix_TEND'
           'cot_scmix_TALL'
           'cot_cumul_TSTART'
           'cot_cumul_TMID'
           'cot_cumul_TEND'
           'cot_cumul_TALL'
        }
        'bgraph_cot'
      }

      { % cloud depth
        'avg_ctype_cdepth'
        {
           'cdepth_TSTART'
           'cdepth_TMID'
           'cdepth_TEND'
           'cdepth_TALL'
%           'cdepth_strnp_TSTART'
%           'cdepth_strnp_TMID'
%           'cdepth_strnp_TEND'
%           'cdepth_strnp_TALL'
           'cdepth_strat_TSTART'
           'cdepth_strat_TMID'
           'cdepth_strat_TEND'
           'cdepth_strat_TALL'
           'cdepth_scmix_TSTART'
           'cdepth_scmix_TMID'
           'cdepth_scmix_TEND'
           'cdepth_scmix_TALL'
           'cdepth_cumul_TSTART'
           'cdepth_cumul_TMID'
           'cdepth_cumul_TEND'
           'cdepth_cumul_TALL'
        }
        'bgraph_cdepth'
      }

      { % liquid water path
        'avg_ctype_lwp'
        {
           'lwp_TSTART'
           'lwp_TMID'
           'lwp_TEND'
           'lwp_TALL'
%           'lwp_strnp_TSTART'
%           'lwp_strnp_TMID'
%           'lwp_strnp_TEND'
%           'lwp_strnp_TALL'
           'lwp_strat_TSTART'
           'lwp_strat_TMID'
           'lwp_strat_TEND'
           'lwp_strat_TALL'
           'lwp_scmix_TSTART'
           'lwp_scmix_TMID'
           'lwp_scmix_TEND'
           'lwp_scmix_TALL'
           'lwp_cumul_TSTART'
           'lwp_cumul_TMID'
           'lwp_cumul_TEND'
           'lwp_cumul_TALL'
        }
        'bgraph_lwp'
      }

      { % precip rate
        'avg_ctype_pcprr'
        {
           'pcprr_TSTART'
           'pcprr_TMID'
           'pcprr_TEND'
           'pcprr_TALL'
%           'pcprr_strnp_TSTART'
%           'pcprr_strnp_TMID'
%           'pcprr_strnp_TEND'
%           'pcprr_strnp_TALL'
           'pcprr_strat_TSTART'
           'pcprr_strat_TMID'
           'pcprr_strat_TEND'
           'pcprr_strat_TALL'
           'pcprr_scmix_TSTART'
           'pcprr_scmix_TMID'
           'pcprr_scmix_TEND'
           'pcprr_scmix_TALL'
           'pcprr_cumul_TSTART'
           'pcprr_cumul_TMID'
           'pcprr_cumul_TEND'
           'pcprr_cumul_TALL'
        }
        'bgraph_pcprr'
      }

      { % cloud fraction
        'avg_ctype_cfrac'
        {
           'cfrac_TSTART'
           'cfrac_TMID'
           'cfrac_TEND'
           'cfrac_TALL'
           % any selection by cloud type will just
           % create CF = 1.0 for all cases
        }
        'bgraph_cfrac'
      }

      };
    Nvarsets = length(VarSets);

    CcnNames = { 'ccn0050' 'ccn0100' 'ccn0200' 'ccn0400' 'ccn0800' 'ccn1600' };
    CcnVals  = [       50       100       200       400       800      1600  ];
    Nccn = length(CcnNames);

    GccnNames = { 'gcn10m5.1um' 'gcn10m4.3um' 'gcn10m2.3um' 'gcn10m0.3um' }; 
    GccnVals  = [     1e-5          1e-4          1e-2             1      ]; 
    Ngccn = length(GccnNames);

    SstNames = { 'sst293' 'sst298' 'sst303' };
    SstVals  = [     293      298      303  ];
    Nsst = length(SstNames);

    % Organize data from set of files into a 4D array with dimensions
    %   dim 1 -> Var 
    %   dim 2 -> SST
    %   dim 3 -> CCN
    %   dim 4 -> GCCN
    %
    % If a file is missing, put a nan into the output array
    %
    for ivset = 1:length(VarSets)
      InFprefix  = VarSets{ivset}{1};
      VarList    = VarSets{ivset}{2};
      OutFprefix = VarSets{ivset}{3};

      Nvars = length(VarList);
      OutFile = sprintf('%s/%s.h5', Ddir, OutFprefix);

      % Form the HDF5 variable names by inserting "_avg" and "_npts" into
      % the names in VarList
      for ivar = 1:Nvars
        AvgVarList{ivar}  = regexprep(VarList{ivar}, '_T', '_avg_T');
        NptsVarList{ivar} = regexprep(VarList{ivar}, '_T', '_npts_T');
      end

      fprintf('*************************************************************************\n');
      fprintf('Generating bar graph data:\n');
      fprintf('  Variables:\n');
      fprintf('\n');
      for ivar = 1:Nvars
        fprintf('    %s --> %s, %s\n', VarList{ivar}, AvgVarList{ivar}, NptsVarList{ivar});
      end
      fprintf('\n');
      fprintf('  Input files:\n');

      % find all of the combinations of sst, ccn and gccn for each variable
      OutAvgs = nan([ Nvars Nsst Nccn Ngccn ]);  % leave nan in entry if file is missing
      OutNpts = nan([ Nvars Nsst Nccn Ngccn ]);  % leave nan in entry if file is missing
      for isst = 1:Nsst
        for iccn = 1:Nccn
          for igccn = 1:Ngccn
            InFile = sprintf('%s/%s_z.atex.%s.%s.%s.h5', Ddir, InFprefix, CcnNames{iccn}, SstNames{isst}, GccnNames{igccn});
            if (exist(InFile, 'file') == 2)
              fprintf('   %s\n', InFile);

              % use the histogram counts and bin values to calculate an average
              % then place the average into the output array
              for ivar = 1:Nvars
                AVG = squeeze(hdf5read(InFile, AvgVarList{ivar}));
                OutNpts(ivar,isst,iccn,igccn) = squeeze(hdf5read(InFile, NptsVarList{ivar}));

                if (isnan(AVG))
                  OutAvgs(ivar,isst,iccn,igccn) = 0;
                else
                  OutAvgs(ivar,isst,iccn,igccn) = AVG;
                end
              end
            end
          end
        end
      end
      fprintf('\n');

      % dump out the averages plus the coordinate values into the output file
      fprintf('  Writing: %s\n', OutFile);

      hdf5write(OutFile, 'Averages', OutAvgs);
      hdf5write(OutFile, 'Npoints',  OutNpts,  'WriteMode', 'append');

      hdf5write(OutFile, 'VarNames', VarList,   'WriteMode', 'append');
      hdf5write(OutFile, 'SST',      SstVals,  'WriteMode', 'append');
      hdf5write(OutFile, 'CCN',      CcnVals,  'WriteMode', 'append');
      hdf5write(OutFile, 'GCCN',     GccnVals, 'WriteMode', 'append');

      fprintf('\n');
    end
end
