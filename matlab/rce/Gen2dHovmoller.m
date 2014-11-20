function [ ] = Gen2dHovmoller(ConfigFile)
% Gen2dHovmoller generate hovmoller data from time series of 2D fields

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
    { 'HDF5/RCE70_OLD_RI/vint_vapor-a-AS-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE70_OLD_RI' 'midy' 1 'hov_precip_water' 'precip_water' }
    { 'HDF5/MATT/vint_vapor-a-AS-2012-01-01-000000-g1.h5'         'vertint_vapor' 'MATT'         'midy' 1 'hov_precip_water' 'precip_water' }
    { 'HDF5/RCE50_NEW_RI/vint_vapor-a-AS-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE50_NEW_RI' 'midy' 1 'hov_precip_water' 'precip_water' }
    { 'HDF5/RCE50_OLD_RI/vint_vapor-a-AS-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE50_OLD_RI' 'midy' 1 'hov_precip_water' 'precip_water' }
    };
  Nset = length(VarSets);


  fprintf('***************************************************************\n');
  fprintf('Generating hovmoller data:\n');

  for iset = 1:Nset
    InFile     = VarSets{iset}{1};
    InVname    = VarSets{iset}{2};
    Case       = VarSets{iset}{3};
    Slice      = VarSets{iset}{4};
    CoordIncKm = VarSets{iset}{5};
    OutFprefix = VarSets{iset}{6};
    OutVname   = VarSets{iset}{7};

    OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);

    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Get entire field from input file
    % If a 2D field, then data will be organized as (x,y,t)
    % Skip over if not a 2D field
    fprintf('    Reading: %s (%s)\n', InFile, InVname);
 
    HDATA = squeeze(hdf5read(InFile, InVname));

    if (ndims(HDATA) ~= 3)
      fprintf('      WARNING: %s is not a 2D field, skipping this variable\n', InVname);
      continue;
    end

    % Data is (x,y,t)
    T     = squeeze(hdf5read(InFile, 't_coords'));
    [ Nx Ny Nt ] = size(HDATA);

    % Grab the slice given in the Slice spec
    % x and y are lon and lat, but for RCE would rather have these in km
    %   use CoordIncKm to form x or y coordinate values start at zero
    % Put output var into 4D organization: (x,y,z,t)
    if (strcmp(Slice, 'midx'))
      Mid = round(Nx/2);
      HV = squeeze(HDATA(Mid,:,:));  % (y,t)
      X = 1;
      Y = [ 0:Ny-1 ] .* CoordIncKm;
      OutVar = reshape(HV, [ 1 Ny 1 Nt ]);
    elseif (strcmp(Slice, 'midy'))
      Mid = round(Ny/2);
      HV = squeeze(HDATA(:,Mid,:));  % (x,t)
      X = [ 0:Nx-1 ] .* CoordIncKm;
      Y = 1;
      OutVar = reshape(HV, [ Nx 1 1 Nt ]);
    end
    Z = 1;


    fprintf('    Writing: %s (%s)\n', OutFile, OutVname);
    fprintf('\n');

    hdf5write(OutFile, OutVname, OutVar); 
    hdf5write(OutFile, '/x_coords', X, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', Y, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Z, 'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');
  end
end
