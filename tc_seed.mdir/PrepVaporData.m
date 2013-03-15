% PrepVaporData
%
% Script to create a set of files that Vapor can ingest
%

clear;

Xmin = -43; % degrees lon
Xmax = -36;

Ymin = 12; % degrees lat
Ymax = 18;

% If Zmin changed from zero, then change the assignment
% of ZDATA(1) to zero in the code below.
Zmin = 0;  % km
Zmax = 18; 

Tmin = 24; % hr
Tmax = 120;

Vdir = 'VAPOR';
if (exist(Vdir, 'dir') ~= 7)
  mikdir(Vdir)
end

InFiles = {
  'cloud-TCS_SD_C2000-AS-1998-08-22-120000-g3.h5'
  'ccn_conc-TCS_SD_C2000-AS-1998-08-22-120000-g3.h5'
  'ccn_rain_mass-TCS_SD_C2000-AS-1998-08-22-120000-g3.h5'
  'cloud-TCS_SD_C2000-AS-1998-08-22-120000-g3.h5'
  };

InVars = {
  'cloud'
  'ccn_concen'
  'ccn_rain_mass'
  'ELEVATION'
  };

VapVars = {
  'cloud'
  'aerosol'
  'aero_rain'
  'ELEVATION'
  };



fprintf('Preparing input files:\n');
% prepare the 3D fields
for i = 1:length(InFiles)
  clear HDATA;

  Hfile = InFiles{i};
  Hdset = InVars{i};

  Vdset = VapVars{i};

  fprintf(' Reading file: %s, Dataset: %s\n', Hfile, Hdset);
  fprintf('\n');
  if (strcmp(Hdset, 'ELEVATION'))
    fprintf ('    Doing ELEVATION, skip read of main variable\n');
  else
    HDATA = hdf5read(Hfile, Hdset);
  end
  X = hdf5read(Hfile, 'x_coords');      % degrees lon
  Y = hdf5read(Hfile, 'y_coords');      % degrees lat
  Z = hdf5read(Hfile, 'z_coords')/1000; % km
  T = hdf5read(Hfile, 't_coords')/3600; % hr

  % pull out the selected data
  X1 = find(X >= Xmin, 1, 'first');
  X2 = find(X <= Xmax, 1, 'last');
  Y1 = find(Y >= Ymin, 1, 'first');
  Y2 = find(Y <= Ymax, 1, 'last');
  Z1 = find(Z >= Zmin, 1, 'first');
  Z2 = find(Z <= Zmax, 1, 'last');
  T1 = find(T >= Tmin, 1, 'first');
  T2 = find(T <= Tmax, 1, 'last');

  XDATA = X(X1:X2);
  YDATA = Y(Y1:Y2);
  ZDATA = Z(Z1:Z2);
  TDATA = T(T1:T2);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fix up coordinate values
  % Change X and Y data to km, round to integers since used integer spacing in RAMS
  % Force ZDATA(1) to zero - change this if Zmin changes
  Xkm = LonToKm(XDATA, 15); % use 15N as reference
  Ykm = LatToKm(YDATA);

  XDATA = round(Xkm - Xkm(1)); % start at zero
  YDATA = round(Ykm - Ykm(1));

  ZDATA(1) = 0; % change this if Zmin changes

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % if we are doing ELEVATION, write out the coordinate data as well
  % as construction VDATA.
  if (strcmp(Hdset, 'ELEVATION'))
    fprintf('Constructing ELEVATION\n');
    Nx = length(XDATA);
    Ny = length(YDATA);
    Nz = length(ZDATA);
    Nt = length(TDATA);
  
    % Reshape ZDATA into the (x,y,z,t) organization, then
    % replicate the ZDATA column into all the x, y, t locations.
    Zvals = reshape(ZDATA, [ 1 1 Nz 1 ]);
    VDATA = repmat(Zvals, [ Nx Ny 1 Nt ]);

    % Write out the x,y and z values into coordinate files
    OutFile = sprintf('%s/XCOORDS', Vdir);
    fprintf('  Creating coordinate file: %s\n', OutFile);
    WriteCoordData(XDATA, OutFile)

    OutFile = sprintf('%s/YCOORDS', Vdir);
    fprintf('  Creating coordinate file: %s\n', OutFile);
    WriteCoordData(YDATA, OutFile)

    OutFile = sprintf('%s/ZCOORDS', Vdir);
    fprintf('  Creating coordinate file: %s\n', OutFile);
    WriteCoordData(ZDATA, OutFile)
  else
    VDATA = squeeze(HDATA(X1:X2,Y1:Y2,Z1:Z2,T1:T2));
  end

  OutFile = sprintf('%s/%s_vap.h5', Vdir, Vdset);
  fprintf('  Writing file: %s, Dataset: %s\n', OutFile, Vdset);
  hdf5write(OutFile, Vdset, VDATA);
  hdf5write(OutFile, 'x_coords', XDATA, 'WriteMode', 'append');
  hdf5write(OutFile, 'y_coords', YDATA, 'WriteMode', 'append');
  hdf5write(OutFile, 'z_coords', ZDATA, 'WriteMode', 'append');
  hdf5write(OutFile, 't_coords', TDATA, 'WriteMode', 'append');

  % attach the dimensions
  fprintf('Attaching dimensions: %s\n', OutFile);
  file_id = H5F.open(OutFile, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
  var_dset_id = H5D.open(file_id, Vdset, 'H5P_DEFAULT');
  x_dset_id = H5D.open(file_id, 'x_coords', 'H5P_DEFAULT');
  y_dset_id = H5D.open(file_id, 'y_coords', 'H5P_DEFAULT');
  z_dset_id = H5D.open(file_id, 'z_coords', 'H5P_DEFAULT');
  t_dset_id = H5D.open(file_id, 't_coords', 'H5P_DEFAULT');
  
  % set x,y,z,t vars as dimensions
  H5DS.set_scale(x_dset_id, 'x');
  H5DS.set_scale(y_dset_id, 'y');
  H5DS.set_scale(z_dset_id, 'z');
  H5DS.set_scale(t_dset_id, 't');
  
  % attach dimensions to
  % the dimension order gets reversed (third argument to
  % the H5DS.attach_scale routine) since the dimensions 
  % get reversed in the HDF5 file.
  H5DS.attach_scale(var_dset_id, x_dset_id, 3);
  H5DS.attach_scale(var_dset_id, y_dset_id, 2);
  H5DS.attach_scale(var_dset_id, z_dset_id, 1);
  H5DS.attach_scale(var_dset_id, t_dset_id, 0);
  
  % clean up
  H5D.close(x_dset_id);
  H5D.close(y_dset_id);
  H5D.close(z_dset_id);
  H5D.close(t_dset_id);
  H5D.close(var_dset_id);
  H5F.close(file_id);
  fprintf('\n');
end

