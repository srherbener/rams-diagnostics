%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NotateVariableXyzt()
%
% This routine will notate the x,y,z,t dimensions in the
% given file according to the COARDS convention. This will allow
% GRADS to read in this file with its sdfopen command.
%
function [] = NotateDimensionsXyzt(File, Vname, Units, LongName, DimNames);

  WriteStringAttribute(File, Vname, 'units', Units);
  WriteStringAttribute(File, Vname, 'long_name', LongName);
  WriteStringAttribute(File, Vname, 'DimNames', DimNames);
  WriteStringAttribute(File, Vname, 'ArrayOrg', 'row major');
end
