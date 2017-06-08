%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CreateDimsensionsXyzt()
%
% This routine will create 4D dimensions in the HDF5 file.
%
%
function [] = CreateDimensionsXyzt(File, X, Y, Z, T, Xname, Yname, Zname, Tname)

  CreateDimension(File, X, Xname, 'x')
  CreateDimension(File, Y, Yname, 'y')
  CreateDimension(File, Z, Zname, 'z')
  CreateDimension(File, T, Tname, 't')

end
