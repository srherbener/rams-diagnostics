%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CreateDimsensionsXyzt()
%
% This routine will create 4D dimensions in the HDF5 file.
%
%
function [] = CreateDimensionsXyzt(File, X, Y, Z, T, Xname, Yname, Zname, Tname)

  Nx = length(X);
  Ny = length(Y);
  Nz = length(Z);
  Nt = length(T);

  % write out the coordinate values
  h5create(File, Xname, Nx, 'DataType', 'single');
  h5write(File,  Xname, X);

  h5create(File, Yname, Ny, 'DataType', 'single');
  h5write(File,  Yname, Y);

  h5create(File, Zname, Nz, 'DataType', 'single');
  h5write(File,  Zname, Z);

  h5create(File, Tname, Nt, 'DataType', 'single');
  h5write(File,  Tname, T);

  % Use low level HDF5 routiens to mark these as dimensions
  file_id = H5F.open(File, 'H5F_ACC_RDWR', 'H5P_DEFAULT');

  x_id = H5D.open(file_id, Xname, 'H5P_DEFAULT');
  y_id = H5D.open(file_id, Yname, 'H5P_DEFAULT');
  z_id = H5D.open(file_id, Zname, 'H5P_DEFAULT');
  t_id = H5D.open(file_id, Tname, 'H5P_DEFAULT');

  % set x,y,z,t as dimensions
  H5DS.set_scale(x_id, 'x');
  H5DS.set_scale(y_id, 'y');
  H5DS.set_scale(z_id, 'z');
  H5DS.set_scale(t_id, 't');
  
  % close dimension variables, and file
  H5D.close(x_id);
  H5D.close(y_id);
  H5D.close(z_id);
  H5D.close(t_id);
  H5F.close(file_id);
end
