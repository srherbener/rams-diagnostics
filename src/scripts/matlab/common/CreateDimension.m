%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CreateDimsension()
%
% This routine will create a dimension dataset in the HDF5 file.
%
% Dtype is one of 'x', 'y', 'z', 'p', 't'
%
%
function [] = CreateDimension(File, Dval, Dname, Dtype)

  Dlen = length(Dval);

  % write out the coordinate values
  h5create(File, Dname, Dlen, 'DataType', 'single');
  h5write(File,  Dname, single(Dval));

  % Use low level HDF5 routiens to mark these as dimensions
  file_id = H5F.open(File, 'H5F_ACC_RDWR', 'H5P_DEFAULT');

  d_id = H5D.open(file_id, Dname, 'H5P_DEFAULT');

  % set x,y,z,t as dimensions
  H5DS.set_scale(d_id, Dtype);
  
  % close dimension variables, and file
  H5D.close(d_id);
  H5F.close(file_id);
end
