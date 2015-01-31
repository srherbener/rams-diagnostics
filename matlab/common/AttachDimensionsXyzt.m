%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AttachDimensionsXyzt()
%
% This routine will attach the x,y,z,t dimensions to the
% given file and dataset. It is assumed that the dataset
% is organized as (x,y,z,t).
%
function [] = AttachDimensionsXyzt(File, Dataset, Xname, Yname, Zname, Tname);

  file_id = H5F.open(File, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
  dset_id = H5D.open(file_id, Dataset, 'H5P_DEFAULT'); 

  x_id = H5D.open(file_id, Xname, 'H5P_DEFAULT');
  y_id = H5D.open(file_id, Yname, 'H5P_DEFAULT');
  z_id = H5D.open(file_id, Zname, 'H5P_DEFAULT');
  t_id = H5D.open(file_id, Tname, 'H5P_DEFAULT');

  % Dimensions go in reverse order since the data gets stored that
  % way. MATLAB is column-major, whereas C is used to write the HDF5
  % file making the file storage row-major.
  %
  H5DS.attach_scale(dset_id, x_id, 3);
  H5DS.attach_scale(dset_id, y_id, 2);
  H5DS.attach_scale(dset_id, z_id, 1);
  H5DS.attach_scale(dset_id, t_id, 0);

  H5D.close(x_id);
  H5D.close(y_id);
  H5D.close(z_id);
  H5D.close(t_id);
  H5D.close(dset_id);
  H5F.close(file_id);
end
