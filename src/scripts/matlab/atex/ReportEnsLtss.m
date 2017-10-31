function [ ] = ReportEnsLtss(ConfigFile)
% ReportEnsLtss report ensemble stats for LTSS

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;

  InFile = sprintf('%s/%s', Ddir, 'ensemble_ts.h5');

  VarList = {
    { 293 'ltss_293_ens_min' 'ltss_293_ens_max' 'ltss_293_ens_avg' }
    { 298 'ltss_298_ens_min' 'ltss_298_ens_max' 'ltss_298_ens_avg' }
    { 303 'ltss_303_ens_min' 'ltss_303_ens_max' 'ltss_303_ens_avg' }
    };

  T1 = 145;
  T2 = 433;

  fprintf('\n');
  fprintf('%15s %15s %15s %15s\n', 'SST', 'Min LTSS', 'Max LTSS', 'Avg LTSS');
  fprintf('\n');
  for ivar = 1:length(VarList)
    SST = VarList{ivar}{1};

    HDATA = squeeze(hdf5read(InFile, VarList{ivar}{2}));
    MIN_LTSS = mean(HDATA(T1:T2));
    HDATA = squeeze(hdf5read(InFile, VarList{ivar}{3}));
    MAX_LTSS = mean(HDATA(T1:T2));
    HDATA = squeeze(hdf5read(InFile, VarList{ivar}{4}));
    AVG_LTSS = mean(HDATA(T1:T2));

    fprintf('%15d %15.2f %15.2f %15.2f\n', SST, MIN_LTSS, MAX_LTSS, AVG_LTSS);
  end
end
