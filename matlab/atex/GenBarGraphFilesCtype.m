function [ ] = GenBarGraphFilesCtype(ConfigFile)
% GenBarGraphFilesCtype generate data bar graphs from cloud type selected data

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Tdir = Config.TsavgDir;
    Hdir = './HDF5';


    VarSets = {
      { % cloud optical thickness
        'pdf_ctype_cot'
        {
          'cot_hist_TSTART'
          'cot_hist_TMID'
          'cot_hist_TEND'
          'cot_hist_TALL'
          'cot_strnp_hist_TSTART'
          'cot_strnp_hist_TMID'
          'cot_strnp_hist_TEND'
          'cot_strnp_hist_TALL'
          'cot_strat_hist_TSTART'
          'cot_strat_hist_TMID'
          'cot_strat_hist_TEND'
          'cot_strat_hist_TALL'
          'cot_scmix_hist_TSTART'
          'cot_scmix_hist_TMID'
          'cot_scmix_hist_TEND'
          'cot_scmix_hist_TALL'
          'cot_cumul_hist_TSTART'
          'cot_cumul_hist_TMID'
          'cot_cumul_hist_TEND'
          'cot_cumul_hist_TALL'
        }
        'x_coords'
        'bgraph_cot'
      }

      { % cloud depth 
        'pdf_ctype_cd'
        {
          'cd_hist_TSTART'
          'cd_hist_TMID'
          'cd_hist_TEND'
          'cd_hist_TALL'
          'cd_strnp_hist_TSTART'
          'cd_strnp_hist_TMID'
          'cd_strnp_hist_TEND'
          'cd_strnp_hist_TALL'
          'cd_strat_hist_TSTART'
          'cd_strat_hist_TMID'
          'cd_strat_hist_TEND'
          'cd_strat_hist_TALL'
          'cd_scmix_hist_TSTART'
          'cd_scmix_hist_TMID'
          'cd_scmix_hist_TEND'
          'cd_scmix_hist_TALL'
          'cd_cumul_hist_TSTART'
          'cd_cumul_hist_TMID'
          'cd_cumul_hist_TEND'
          'cd_cumul_hist_TALL'
        }
        'x_coords'
        'bgraph_cd'
      }

      { % precip rate
        'pdf_ctype_pcprr'
        {
          'pcprr_hist_TSTART'
          'pcprr_hist_TMID'
          'pcprr_hist_TEND'
          'pcprr_hist_TALL'
          'pcprr_strnp_hist_TSTART'
          'pcprr_strnp_hist_TMID'
          'pcprr_strnp_hist_TEND'
          'pcprr_strnp_hist_TALL'
          'pcprr_strat_hist_TSTART'
          'pcprr_strat_hist_TMID'
          'pcprr_strat_hist_TEND'
          'pcprr_strat_hist_TALL'
          'pcprr_scmix_hist_TSTART'
          'pcprr_scmix_hist_TMID'
          'pcprr_scmix_hist_TEND'
          'pcprr_scmix_hist_TALL'
          'pcprr_cumul_hist_TSTART'
          'pcprr_cumul_hist_TMID'
          'pcprr_cumul_hist_TEND'
          'pcprr_cumul_hist_TALL'
        }
        'x_coords'
        'bgraph_pcprr'
      }

      { % albedo
        'pdf_ctype_albedo'
        {
          'albedo_hist_TSTART'
          'albedo_hist_TMID'
          'albedo_hist_TEND'
          'albedo_hist_TALL'
          'albedo_strnp_hist_TSTART'
          'albedo_strnp_hist_TMID'
          'albedo_strnp_hist_TEND'
          'albedo_strnp_hist_TALL'
          'albedo_strat_hist_TSTART'
          'albedo_strat_hist_TMID'
          'albedo_strat_hist_TEND'
          'albedo_strat_hist_TALL'
          'albedo_scmix_hist_TSTART'
          'albedo_scmix_hist_TMID'
          'albedo_scmix_hist_TEND'
          'albedo_scmix_hist_TALL'
          'albedo_cumul_hist_TSTART'
          'albedo_cumul_hist_TMID'
          'albedo_cumul_hist_TEND'
          'albedo_cumul_hist_TALL'
        }
        'x_coords'
        'bgraph_albedo'
      }

      { % lwp
        'pdf_ctype_lwp'
        {
          'lwp_hist_TSTART'
          'lwp_hist_TMID'
          'lwp_hist_TEND'
          'lwp_hist_TALL'
          'lwp_strnp_hist_TSTART'
          'lwp_strnp_hist_TMID'
          'lwp_strnp_hist_TEND'
          'lwp_strnp_hist_TALL'
          'lwp_strat_hist_TSTART'
          'lwp_strat_hist_TMID'
          'lwp_strat_hist_TEND'
          'lwp_strat_hist_TALL'
          'lwp_scmix_hist_TSTART'
          'lwp_scmix_hist_TMID'
          'lwp_scmix_hist_TEND'
          'lwp_scmix_hist_TALL'
          'lwp_cumul_hist_TSTART'
          'lwp_cumul_hist_TMID'
          'lwp_cumul_hist_TEND'
          'lwp_cumul_hist_TALL'
        }
        'x_coords'
        'bgraph_lwp'
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

    % Organize data from set of files into a 3D array with dimensions
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
      BinVname   = VarSets{ivset}{3};
      OutFprefix = VarSets{ivset}{4};

      Nvars = length(VarList);
      OutFile = sprintf('%s/%s.h5', Ddir, OutFprefix);

      fprintf('*************************************************************************\n');
      fprintf('Generating bar graph data:\n');
      fprintf('  Variables:\n');
      fprintf('\n');
      for ivar = 1:Nvars
        fprintf('    %s\n', VarList{ivar});
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
              BINS = squeeze(hdf5read(InFile, BinVname));
              for ivar = 1:Nvars
                Vname = VarList{ivar};
                HIST = squeeze(hdf5read(InFile, Vname));

                N = sum(HIST);
                S = sum(HIST .* BINS);
                if (N == 0)
                  OutAvgs(ivar,isst,iccn,igccn) = 0;
                else
                  OutAvgs(ivar,isst,iccn,igccn) = S / N;
                end

                OutNpts(ivar,isst,iccn,igccn) = N;
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

      hdf5write(OutFile, 'VarNames', VarList,  'WriteMode', 'append');
      hdf5write(OutFile, 'SST',      SstVals,  'WriteMode', 'append');
      hdf5write(OutFile, 'CCN',      CcnVals,  'WriteMode', 'append');
      hdf5write(OutFile, 'GCCN',     GccnVals, 'WriteMode', 'append');

      fprintf('\n');
    end
end
