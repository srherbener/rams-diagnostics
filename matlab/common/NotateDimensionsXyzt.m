%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NotateDimensionsXyzt()
%
% This routine will notate the x,y,z,t dimensions in the
% given file according to the COARDS convention. This will allow
% GRADS to read in this file with its sdfopen command.
%
function [] = NotateDimensionsXyzt(File, Xname, Yname, Zname, Tname);

  WriteStringAttribute(File, Xname, 'axis', 'x');
  WriteStringAttribute(File, Xname, 'units', 'degrees_east');
  WriteStringAttribute(File, Xname, 'long_name', 'longitude');
  WriteStringAttribute(File, Xname, 'DimNames', 'x');
  WriteStringAttribute(File, Xname, 'ArrayOrg', 'row major');

  WriteStringAttribute(File, Yname, 'axis', 'y');
  WriteStringAttribute(File, Yname, 'units', 'degrees_north');
  WriteStringAttribute(File, Yname, 'long_name', 'latitude');
  WriteStringAttribute(File, Yname, 'DimNames', 'y');
  WriteStringAttribute(File, Yname, 'ArrayOrg', 'row major');

  WriteStringAttribute(File, Zname, 'axis', 'z');
  WriteStringAttribute(File, Zname, 'units', 'meters');
  WriteStringAttribute(File, Zname, 'long_name', 'sigma-z');
  WriteStringAttribute(File, Zname, 'DimNames', 'z');
  WriteStringAttribute(File, Zname, 'ArrayOrg', 'row major');

  WriteStringAttribute(File, Tname, 'axis', 't');
  WriteStringAttribute(File, Tname, 'units', 'seconds since 2006-08-20 12:00:00 00:00');
  WriteStringAttribute(File, Tname, 'long_name', 'simulation time');
  WriteStringAttribute(File, Tname, 'DimNames', 't');
  WriteStringAttribute(File, Tname, 'ArrayOrg', 'row major');

end
