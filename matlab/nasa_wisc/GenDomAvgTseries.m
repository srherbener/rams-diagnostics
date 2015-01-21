function [ ] = GenDomAvgTseries(ConfigFile)
% GenDomAvgTseries generate time series of domain averages

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % input_file input_dataset case_name outfile_prefix output_dataset
  VarSets = {
%    { 'HDF5/RCE70_OLD_RI/vint_vapor-a-AS-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE70_OLD_RI' 'avg_precip_water' 'precip_water' }
%    { 'HDF5/MATT/vint_vapor-a-AS-2012-01-01-000000-g1.h5'         'vertint_vapor' 'MATT'         'avg_precip_water' 'precip_water' }
    { 'HDF5/RCE50_RECT/vint_vapor-a-AS-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE50_RECT' 'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE50_OLD_RI/vint_vapor-a-AS-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE50_OLD_RI' 'avg_precip_water' 'precip_water' }
    };
  Nset = length(VarSets);


  fprintf('***************************************************************\n');
  fprintf('Generating domain average time series:\n');

  for iset = 1:Nset
    InFile     = VarSets{iset}{1};
    InVname    = VarSets{iset}{2};
    Case       = VarSets{iset}{3};
    OutFprefix = VarSets{iset}{4};
    OutVname   = VarSets{iset}{5};

    OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);

    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Get entire field from input file
    % Do horizontal domain average
    %   if 2D field then result is a single number
    %   if 3D field then result is a vertical profile
    fprintf('    Reading: %s (%s)\n', InFile, InVname);
 
    HDATA = squeeze(hdf5read(InFile, InVname));
    Z     = squeeze(hdf5read(InFile, 'z_coords'));
    T     = squeeze(hdf5read(InFile, 't_coords'));
    Nz = length(Z);
    Nt = length(T);


    % Data is either (x,y,t) or (x,y,z,t) so calculate means on the first two dimensions and
    % write out the result
    AVG = squeeze(mean(HDATA,1));
    AVG = squeeze(mean(AVG,1));

    % Z will be the appropriate size
    OutVar = reshape(AVG, [ 1 1 Nz Nt ]);

    % Write out dummy coordinate values to keep ReadSelectXyzt happy
    Xdummy = 1;
    Ydummy = 1;

    fprintf('    Writing: %s (%s)\n', OutFile, OutVname);
    fprintf('\n');

    hdf5write(OutFile, OutVname, OutVar); 
    hdf5write(OutFile, '/x_coords', Xdummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', Ydummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Z,      'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T,      'WriteMode', 'append');
  end
end
