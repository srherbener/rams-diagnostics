%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AttachDimensionsXyzt()
%
% This routine will attach the dimensions specified in DimOrder
% to the given file and dataset.
%
function [] = AttachDimensionsXyzt(File, Dataset, DimOrder, Xname, Yname, Zname, Tname);

  DsInfo = h5info(File, Dataset);
  Ndims = length(DsInfo.Dataspace.Size);

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
  Ndims = length(DimOrder);
  for i = Ndims:-1:1
    switch(DimOrder{i})
      case 'x'
        H5DS.attach_scale(dset_id, x_id, Ndims-i);
      case 'y'
        H5DS.attach_scale(dset_id, y_id, Ndims-i);
      case 'z'
        H5DS.attach_scale(dset_id, z_id, Ndims-i);
      case 't'
        H5DS.attach_scale(dset_id, t_id, Ndims-i);
    end
  end

  H5D.close(x_id);
  H5D.close(y_id);
  H5D.close(z_id);
  H5D.close(t_id);
  H5D.close(dset_id);
  H5F.close(file_id);

end
