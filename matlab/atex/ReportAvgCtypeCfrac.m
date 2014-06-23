function [ ] = ReportAvgCtypeCfrac(ConfigFile)
% ReportAvgCtypeCfrac report ctype averages for cloud fraction

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;

  InFprefix = 'avg_ctype_cfrac';

  InVarStrnp = 'cfrac_strnp_avg_TALL';
  InVarStmix = 'cfrac_stmix_avg_TALL';
  InVarStrat = 'cfrac_strat_avg_TALL';

  InVarScmix = 'cfrac_scmix_avg_TALL';

  InVarCumul = 'cfrac_cumul_avg_TALL';

  fprintf('\n');
  fprintf('%40s %10s %10s %10s\n', 'Case', 'ST', 'SC', 'CV');
  fprintf('\n');
  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;
    InFile = sprintf('%s/%s_%s.h5', Ddir, InFprefix, Case);

    STRNP_AVG = hdf5read(InFile, InVarStrnp);
    STMIX_AVG = hdf5read(InFile, InVarStmix);
    STRAT_AVG = hdf5read(InFile, InVarStrat);
    SCMIX_AVG = hdf5read(InFile, InVarScmix);
    CUMUL_AVG = hdf5read(InFile, InVarCumul);

    STALL_AVG = STRNP_AVG + STMIX_AVG + STRAT_AVG;

    ALL_CLOUD = STALL_AVG + SCMIX_AVG + CUMUL_AVG;
    STALL_PERCENT = 100 * (STALL_AVG / ALL_CLOUD);
    SCMIX_PERCENT = 100 * (SCMIX_AVG / ALL_CLOUD);
    CUMUL_PERCENT = 100 * (CUMUL_AVG / ALL_CLOUD);

    fprintf('%40s %10.2f %10.2f %10.2f\n', Case, STALL_PERCENT, SCMIX_PERCENT, CUMUL_PERCENT);
  end
end
