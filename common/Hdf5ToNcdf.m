function [ ] = Hdf5ToNcdf (Fname, Vname, OutFile)
% Hdf5ToNcdf read in a variable from an hdf5 file and write it out in a netcdf file

  fprintf('Reading HDF5 file: %s, Dataset: %s\n', Fname, Vname);
  Var = hdf5read(Fname, Vname); 
  X = hdf5read(Fname, 'x_coords');
  Y = hdf5read(Fname, 'y_coords');
  Z = hdf5read(Fname, 'z_coords');
  T = hdf5read(Fname, 't_coords');
  fprintf('\n');

  % Get the sizes of all the dimensions form Var
  [ Nx Ny Nz Nt ] = size(Var);

  fprintf('Writing netcdf data file: %s\n', OutFile);
  % remove an existing file
  if (exist(OutFile, 'file') == 2)
    delete(OutFile);
  end
  nccreate(OutFile, Vname, 'Dimensions', { 'x_coords' Nx 'y_coords' Ny 'z_coords' Nz 't_coords' Nt }, 'DeflateLevel', 1);
  nccreate(OutFile, 'x_coords', 'Dimensions', { 'x_coords' Nx });
  nccreate(OutFile, 'y_coords', 'Dimensions', { 'y_coords' Ny });
  nccreate(OutFile, 'z_coords', 'Dimensions', { 'z_coords' Nz });
  nccreate(OutFile, 't_coords', 'Dimensions', { 't_coords' Nt });

ncdisp(OutFile);

%  ncwrite(OutFile, Vname, Var);
  ncwrite(OutFile, 'x_coords', X);
  ncwrite(OutFile, 'y_coords', Y);
  ncwrite(OutFile, 'z_coords', Z);
  ncwrite(OutFile, 't_coords', T);

end
