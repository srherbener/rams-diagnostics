function [ ] = GenEnsembleData(ConfigFile)
% GenGenEnsembleData generate ensemble means

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Tdir = Config.TsavgDir;
  Ddir = Config.DiagDir;

  LtssFprefix = 'ltss';
  CfracFprefix = 'hda_ts_cf';

  % ensembles
  %     Name
  %     InVarName
  %     InDir
  %     InFileList
  EnsList = {
    % lower tropospheric static stability, LTSS
    {
    'ltss_293'
    'ltss'
    Tdir
      {
      'ltss_z.atex.ccn0050.sst293.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0100.sst293.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0200.sst293.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0400.sst293.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0800.sst293.gcn10m5.1um.h5'
      'ltss_z.atex.ccn1600.sst293.gcn10m5.1um.h5'
      }
    } 

    {
    'ltss_298'
    'ltss'
    Tdir
      {
      'ltss_z.atex.ccn0050.sst298.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0100.sst298.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0200.sst298.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0400.sst298.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0800.sst298.gcn10m5.1um.h5'
      'ltss_z.atex.ccn1600.sst298.gcn10m5.1um.h5'
      }
    } 

    {
    'ltss_303'
    'ltss'
    Tdir
      {
      'ltss_z.atex.ccn0050.sst303.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0100.sst303.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0200.sst303.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0400.sst303.gcn10m5.1um.h5'
      'ltss_z.atex.ccn0800.sst303.gcn10m5.1um.h5'
      'ltss_z.atex.ccn1600.sst303.gcn10m5.1um.h5'
      }
    } 

    % domain cloud fraction
    {
    'cloud_frac_293'
    'cloud_frac'
    Ddir
      {
      'hda_ts_cf_z.atex.ccn0050.sst293.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0100.sst293.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0200.sst293.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0400.sst293.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0800.sst293.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn1600.sst293.gcn10m5.1um.h5'
      }
    } 

    {
    'cloud_frac_298'
    'cloud_frac'
    Ddir
      {
      'hda_ts_cf_z.atex.ccn0050.sst298.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0100.sst298.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0200.sst298.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0400.sst298.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0800.sst298.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn1600.sst298.gcn10m5.1um.h5'
      }
    } 

    {
    'cloud_frac_303'
    'cloud_frac'
    Ddir
      {
      'hda_ts_cf_z.atex.ccn0050.sst303.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0100.sst303.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0200.sst303.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0400.sst303.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn0800.sst303.gcn10m5.1um.h5'
      'hda_ts_cf_z.atex.ccn1600.sst303.gcn10m5.1um.h5'
      }
    } 

    };
  Nens = length(EnsList);

  % last 12 hours of time series
  Tstart = 145;
  Tend = 433;

  OutFile = sprintf('%s/%s', Ddir, 'ensemble_ts.h5');

  fprintf('***************************************************************\n');
  fprintf('Generating ensemble data:\n');
  fprintf('  Output file: %s\n', OutFile);
  fprintf('\n');

  % write a header into the file so a new file will be created
  hdf5write(OutFile, 'Header', 'ATEX Ensemble Data');

  for iens = 1:Nens
    Name       = EnsList{iens}{1}; 
    InVarName  = EnsList{iens}{2}; 
    InDir      = EnsList{iens}{3}; 
    InFileList = EnsList{iens}{4}; 

    fprintf('  Ensemble: %s\n', Name);
    fprintf('\n');
    % collect ensemble data
    Nfiles = length(InFileList);
    for ifile = 1:Nfiles
      InFile = sprintf('%s/%s', InDir, InFileList{ifile});
      fprintf('  Reading: %s (%s)\n', InFile, InVarName);

      TSERIES = squeeze(hdf5read(InFile, InVarName));
      if (ifile == 1)
        T = squeeze(hdf5read(InFile, 't_coords'));
        Nt = length(T);

        ALL_TS = zeros([ Nfiles Nt ]);
      end
      ALL_TS(ifile, :) = TSERIES;
    end
    fprintf('\n');

    % find min max sum and npoints, and write out to output file
    ENS_MIN = min(ALL_TS, [], 1);
    ENS_MAX = max(ALL_TS, [], 1);
    ENS_AVG = mean(ALL_TS, 1);

    % write out data
    OutName = sprintf('%s_%s', Name, 'ens_min');
    hdf5write(OutFile, OutName, ENS_MIN, 'WriteMode', 'append');

    OutName = sprintf('%s_%s', Name, 'ens_max');
    hdf5write(OutFile, OutName, ENS_MAX, 'WriteMode', 'append');

    OutName = sprintf('%s_%s', Name, 'ens_avg');
    hdf5write(OutFile, OutName, ENS_AVG, 'WriteMode', 'append');
  end

  % output coords so that ReadSelectXyzt can handle the output files
  Xdummy = 1;
  Ydummy = 1;
  Zdummy = 1;

  hdf5write(OutFile, '/x_coords', Xdummy, 'WriteMode', 'append');
  hdf5write(OutFile, '/y_coords', Ydummy, 'WriteMode', 'append');
  hdf5write(OutFile, '/z_coords', Zdummy, 'WriteMode', 'append');
  hdf5write(OutFile, '/t_coords', T,      'WriteMode', 'append');
end
