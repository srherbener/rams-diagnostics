function [ OutData ] = ReadSelectHdf5(Hfile, Vname, SelString)
% ReadSelectHdf5 read variable out of HDF5 file, and select subset based on ranges
%
% Vname needs a leading '/' character.
%
% This function will apply selection based on ranges given in SelString. The contents
% of SelString should be what you would type at the command line for the
% parenthetical section of an array variable. Eg, if you have a 3D variable, and
% want to select all of the first dimension, and index 2 of the 2nd dimension and indices
% 3 through 7 of the third dimension, then at the command line you would type:
%
%    Var(:,2,3:7)
%
% and the equivalent SelString for this function would be:
%
%    '(:,2,3:7)'
%
% Note that if you want to select the entire array, set SelString to a blank
% string ('').
%

  Hdata = squeeze(h5read(Hfile, Vname));
  SelCmd = sprintf('squeeze(Hdata%s)', SelString);
  OutData = eval(SelCmd);
end
