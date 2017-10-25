%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WriteStringAttribute()
%
% This routine will create and attach a string attribute
% to the given file and dataset.
%
function [] = WriteStringAttribute(File, Dataset, AttrName, AttrString)

  file_id = H5F.open(File, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
  dset_id = H5D.open(file_id, Dataset, 'H5P_DEFAULT'); 

  % Going to use C style strings which are null terminated, so append
  % the null character on the end of AttrString
  String = sprintf('%s%c', AttrString, char(0));
  Slen = length(String);

  atype_id = H5T.copy('H5T_C_S1');
  H5T.set_size(atype_id, Slen);

  % Use a scalar memory space so that a single string is created
  % (as opposed to an array of characters).
  mspc_id = H5S.create('H5S_SCALAR');

  attr_id = H5A.create(dset_id, AttrName, atype_id, mspc_id , 'H5P_DEFAULT');
  H5A.write(attr_id, atype_id, String);

  H5S.close(mspc_id);
  H5A.close(attr_id);
  H5T.close(atype_id);
  H5D.close(dset_id);
  H5F.close(file_id);
end
