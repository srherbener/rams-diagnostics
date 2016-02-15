% PrepVaporData
%
% Script to create a set of files that Vapor can ingest
% Created by Steve Herbener
%
% Modified by Leah Grant to work with vapor 2.2.2 -- August 2013
%   essentially what this script does now is to create an ELEVATION variable
%   in its own .h5 file so that VAPOR can correctly plot RAMS output 
%   on a stretched vertical grid.
%
%   ELEVATION_vap.h5 is created with the same x, y, and z dimensions as in
%   the REVU-generated H5 output file. It is a 3-D array of the heights of
%   the model grid levels.
% 
%   The script also creates ASCII files called XCOORDS, YCOORDS, and
%   ZCOORDS, which are the grid points in km (instead of lat/lon) to feed
%   to VAPOR. Calls the function WriteCoordData (These files aren't actually
%   needed any longer, but I left here just in case)
%
%   The four files (ELEVATION_vap.h5, XCOORDS, YCOORDS, ZCOORDS) are
%   created in the directory specified in the variable Vdir below.
%
%   These files must be created before running the build_VAPOR_data script!
%
%   ** User modifications in code below: change dx, dy grid spacings, path
%   (relative to current working directory in Matlab) where relevant VAPOR
%   files are created, path to REVU H5 output, and REVU H5 file name **

clear;


dx = 2.00; dy = 2.00; % grid spacing in km

%Vdir = 'VAPOR_LC1'; % directory
Vdir = 'VAPOR_C2000';

if (exist(Vdir, 'dir') ~= 7)
  mkdir(Vdir)
end

% path to REVU H5 output file
%InFilesPath = 'LC1/HDF5/';
InFilesPath = '../TCS_GN_C2000/HDF5/';
InFiles = {
  'ccn_conc-TCS_GN_C2000-AS-1998-08-22-120000-g3.h5'
  'cloud-TCS_GN_C2000-AS-1998-08-22-120000-g3.h5'
  };

InVars = {
  'ELEVATION'
  };

VapVars = {
  'ELEVATION'
  };

% %%% shouldn't need to modify anything below here %%%
% ----------------------------------------------------

fprintf('Preparing input files:\n');
% prepare the 3D fields
for i = 1:length(InVars) % This should just be ELEVATION

  Hfile = [InFilesPath,InFiles{1}];
  Hdset = InVars{i};

  Vdset = VapVars{i};

  fprintf(' Reading file: %s, Dataset: %s\n', Hfile, Hdset);
  fprintf('\n');
 
  % read coordinates
  XDATA = h5read(Hfile, '/x_coords');  % degrees lon
  YDATA = h5read(Hfile, '/y_coords');  % degrees lat
  ZDATA = h5read(Hfile, '/z_coords');  % m
  
  % re-set X, Y to be in km based on dx and dy
  Nx=length(XDATA); Ny=length(YDATA); Nz=length(ZDATA); % # points
  Xkm = single( [0:Nx-1]*dx ); % create x, y, z arrays in km
  Ykm = single( [0:Ny-1]*dy ); 
  Zkm = single( ZDATA/1000. );


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % if we are doing ELEVATION, write out the coordinate data as well
  % as construction VDATA.
  if (strcmp(Hdset, 'ELEVATION'))
    fprintf('Constructing ELEVATION\n');
  
    % Reshape ZDATA into the (x,y,z) organization, then
    % replicate the ZDATA column into all the x, y locations.
    Zvals = reshape(Zkm, [ 1 1 Nz ]);
    VDATA = repmat(Zvals, [ Nx Ny 1 ]);

    % Write out the x,y and z values into coordinate files
    OutFile = sprintf('%s/XCOORDS', Vdir);
    fprintf('  Creating coordinate file: %s\n', OutFile);
    %WriteCoordData(XDATA, OutFile)
    WriteCoordData(Xkm, OutFile) % write out km in Coords file so 
      % VAPOR will use km instead of lat/lon

    OutFile = sprintf('%s/YCOORDS', Vdir);
    fprintf('  Creating coordinate file: %s\n', OutFile);
    %WriteCoordData(YDATA, OutFile)
    WriteCoordData(Ykm, OutFile)

    OutFile = sprintf('%s/ZCOORDS', Vdir);
    fprintf('  Creating coordinate file: %s\n', OutFile);
    %WriteCoordData(ZDATA, OutFile)
    WriteCoordData(Zkm, OutFile)
  end

  OutFile = sprintf('%s/%s_vap.h5', Vdir, Vdset);
  fprintf('  Writing file: %s, Dataset: %s\n', OutFile, Vdset);
  if (exist(OutFile,'file')==2)
      delete(OutFile)
  end
  h5create(OutFile, ['/',Vdset], size(VDATA),'datatype','single');
  h5create(OutFile,'/x_coords',Nx,'datatype','single');
  h5create(OutFile,'/y_coords',Ny,'datatype','single');
  h5create(OutFile,'/z_coords',Nz,'datatype','single');
  % ELEVATION data -- in km, since dimensions are specified in km
  h5write(OutFile, ['/',Vdset], VDATA);
  % dimensions in orig. REVU file
  h5write(OutFile, '/x_coords', XDATA);
  h5write(OutFile, '/y_coords', YDATA);
  h5write(OutFile, '/z_coords', ZDATA);  

  % attach the dimensions
  fprintf('Attaching dimensions: %s\n', OutFile);
  file_id = H5F.open(OutFile, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
  var_dset_id = H5D.open(file_id, Vdset, 'H5P_DEFAULT');
  x_dset_id = H5D.open(file_id, 'x_coords', 'H5P_DEFAULT');
  y_dset_id = H5D.open(file_id, 'y_coords', 'H5P_DEFAULT');
  z_dset_id = H5D.open(file_id, 'z_coords', 'H5P_DEFAULT');
  
  % set x,y,z vars as dimensions
  H5DS.set_scale(x_dset_id, 'x');
  H5DS.set_scale(y_dset_id, 'y');
  H5DS.set_scale(z_dset_id, 'z');
  
  % attach dimensions to
  % the dimension order gets reversed (third argument to
  % the H5DS.attach_scale routine) since the dimensions 
  % get reversed in the HDF5 file.
%   H5DS.attach_scale(var_dset_id, x_dset_id, 3); % Steve's code
%   H5DS.attach_scale(var_dset_id, y_dset_id, 2);
%   H5DS.attach_scale(var_dset_id, z_dset_id, 1);
%   H5DS.attach_scale(var_dset_id, t_dset_id, 0);
  H5DS.attach_scale(var_dset_id, x_dset_id, 2); % no longer need time
  H5DS.attach_scale(var_dset_id, y_dset_id, 1);
  H5DS.attach_scale(var_dset_id, z_dset_id, 0);
  
  % clean up
  H5D.close(x_dset_id);
  H5D.close(y_dset_id);
  H5D.close(z_dset_id);
  H5D.close(var_dset_id);
  H5F.close(file_id);
  fprintf('\n');
end



% write dimensions in km to h5 files: ELEVATION and original REVU file
for i=1:length(InFiles)+1
  
  if i>length(InFiles) % ELEVATION file
      OutFile = sprintf('%s/%s_vap.h5', Vdir, Vdset);
  else % REVU H5 file
      OutFile = [InFilesPath,InFiles{i}];
  end
  
  % get the variable names
  clear info nVars Varname
  info = h5info(OutFile);
  nVars = length(info.Datasets);
  for n=1:nVars
      Varname{n}=info.Datasets(n).Name;
  end
  
  fprintf('  Writing dimensions in km to file: %s\n', OutFile);
  if (~any(strcmp(Varname,'xkm'))) % these datasets haven't been created yet
    h5create(OutFile, '/xkm', Nx, 'datatype','single');
    h5create(OutFile, '/ykm', Ny, 'datatype','single');
    h5create(OutFile, '/zkm', Nz, 'datatype','single');
  end
  h5write(OutFile, '/xkm', Xkm);
  h5write(OutFile, '/ykm', Ykm);
  h5write(OutFile, '/zkm', Zkm);

  % attach the dimensions
  fprintf('Attaching dimensions: %s\n', OutFile);
  file_id = H5F.open(OutFile, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
  x_dset_id = H5D.open(file_id, 'x_coords', 'H5P_DEFAULT');
  y_dset_id = H5D.open(file_id, 'y_coords', 'H5P_DEFAULT');
  z_dset_id = H5D.open(file_id, 'z_coords', 'H5P_DEFAULT');
  xkm_dset_id = H5D.open(file_id, 'xkm', 'H5P_DEFAULT');
  ykm_dset_id = H5D.open(file_id, 'ykm', 'H5P_DEFAULT');
  zkm_dset_id = H5D.open(file_id, 'zkm', 'H5P_DEFAULT');
  
  H5DS.attach_scale(xkm_dset_id, x_dset_id, 0);
  H5DS.attach_scale(ykm_dset_id, y_dset_id, 0);
  H5DS.attach_scale(zkm_dset_id, z_dset_id, 0);

  % clean up
  H5D.close(x_dset_id);
  H5D.close(y_dset_id);
  H5D.close(z_dset_id);
  % testing...
  H5D.close(xkm_dset_id);
  H5D.close(ykm_dset_id);
  H5D.close(zkm_dset_id);
  H5F.close(file_id);
  fprintf('\n');

end
