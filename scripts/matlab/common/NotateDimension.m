%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NotateDimension()
%
% This routine will notate a dimensions in the given file
% according to the COARDS convention. This will allow
% GRADS to read in this file with its sdfopen command.
%
% Dtype is one of 'x', 'y', 'z', 'p', 't'.
%
function [] = NotateDimension(File, Dname, Dtype);

  switch Dtype
    case 'x'
      WriteStringAttribute(File, Dname, 'axis', 'x');
      WriteStringAttribute(File, Dname, 'units', 'degrees_east');
      WriteStringAttribute(File, Dname, 'long_name', 'longitude');
      WriteStringAttribute(File, Dname, 'DimNames', 'x');
      WriteStringAttribute(File, Dname, 'ArrayOrg', 'row major');

    case 'y'
      WriteStringAttribute(File, Dname, 'axis', 'y');
      WriteStringAttribute(File, Dname, 'units', 'degrees_north');
      WriteStringAttribute(File, Dname, 'long_name', 'latitude');
      WriteStringAttribute(File, Dname, 'DimNames', 'y');
      WriteStringAttribute(File, Dname, 'ArrayOrg', 'row major');

    case 'z'
      WriteStringAttribute(File, Dname, 'axis', 'z');
      WriteStringAttribute(File, Dname, 'units', 'meters');
      WriteStringAttribute(File, Dname, 'long_name', 'sigma-z');
      WriteStringAttribute(File, Dname, 'DimNames', 'z');
      WriteStringAttribute(File, Dname, 'ArrayOrg', 'row major');

    case 'p'
      WriteStringAttribute(File, Dname, 'axis', 'p');
      WriteStringAttribute(File, Dname, 'units', 'millibars');
      WriteStringAttribute(File, Dname, 'long_name', 'pressure');
      WriteStringAttribute(File, Dname, 'DimNames', 'p');
      WriteStringAttribute(File, Dname, 'ArrayOrg', 'row major');

    case 't'
      WriteStringAttribute(File, Dname, 'axis', 't');
      WriteStringAttribute(File, Dname, 'units', 'seconds since 2006-08-20 12:00:00 00:00');
      %WriteStringAttribute(File, Dname, 'units', 'seconds since 2006-08-21 18:00:00 00:00');
      WriteStringAttribute(File, Dname, 'long_name', 'simulation time');
      WriteStringAttribute(File, Dname, 'DimNames', 't');
      WriteStringAttribute(File, Dname, 'ArrayOrg', 'row major');
  end

end
