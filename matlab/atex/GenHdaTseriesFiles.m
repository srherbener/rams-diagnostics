function [ ] = GenHdaTseriesFiles(ConfigFile)
% GenHdaTseriesFiles generate horizontal domain time series

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Tdir = Config.TsavgDir;
  Ddir = Config.DiagDir;

  FileHeader = 'ATEX PDF data';

  VarSets = {
    { 'hda_cloud_frac' 'hda_cloud_mask' 1 2 'hda_ts_cf' 'cloud_frac'  }
    };
  Nset = length(VarSets);


  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;

    for iset = 1:Nset
      InFprefix  = VarSets{iset}{1};
      InVname    = VarSets{iset}{2};
      InSums     = VarSets{iset}{3};
      InNpts     = VarSets{iset}{4};
      OutFprefix = VarSets{iset}{5};
      OutVname   = VarSets{iset}{6};
    
      InFile = sprintf('%s/%s_%s.h5', Tdir, InFprefix, Case);

      fprintf('***************************************************************\n');
      fprintf('Generating hda tseries:\n');
      fprintf('  Case: %s\n', Case);
      fprintf('  Input hda file: %s\n', InFile);
      fprintf('    Variable: %s\n', InVname);
      fprintf('\n');

      % Pull out the sums and counts and calculate the averages
      HDATA = squeeze(hdf5read(InFile, InVname));  % organized (2,t)
      SUMS = squeeze(HDATA(InSums,:));
      NPTS = squeeze(HDATA(InNpts,:));
      HDA = SUMS ./ NPTS;

      % Grab the time coordinates
      T = squeeze(hdf5read(InFile, 't_coords'));
      Nt = length(T);

      % output in GRADS 4D variable format so that GenLinePlot can read it in
      X = 1;
      Y = 1;
      Z = 1;
      
      OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);
      fprintf('  Writing: %s\n', OutFile);

      OutVar = reshape(HDA, [ 1 1 1 Nt ]); % 4D var --> (x,y,z,t)
      hdf5write(OutFile, OutVname, OutVar);

      hdf5write(OutFile, '/x_coords', X, 'WriteMode', 'append');
      hdf5write(OutFile, '/y_coords', Y, 'WriteMode', 'append');
      hdf5write(OutFile, '/z_coords', Z, 'WriteMode', 'append');
      hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');
      
      fprintf('\n');
    end
    fprintf('\n');
  end
end
