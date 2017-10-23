%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NotateDimensionsXyzt()
%
% This routine will notate the x,y,z,t dimensions in the
% given file according to the COARDS convention. This will allow
% GRADS to read in this file with its sdfopen command.
%
function [] = NotateDimensionsXyzt(File, Xname, Yname, Zname, Tname);

  NotateDimension(File, Xname, 'x')
  NotateDimension(File, Yname, 'y')
  NotateDimension(File, Zname, 'z')
  NotateDimension(File, Tname, 't')

end
