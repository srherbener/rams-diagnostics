function [ Avar, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar ( Hfile, Vname )
%ReadAzavgVar read variable from an HDF5 file containing output from azavg
%
% This function will read in the variable plus its associated coordinate
% data from the given HDF5 file. It is assumed that Hfile contains a
% dataset named '/<Vname>' and that this dataset is organized (after
% reading) as (x,y,z,t).
%    x -> radius
%    y -> dummy dimension (to make GRADS happy)
%    z -> height
%    t -> time
%
% The data from Vname will be copied to the output array Avar and reduced
% to (r,z,t).

% Read in the variable and get rid of the dummy dimension
Hvar = hdf5read(Hfile, Vname);
Avar = squeeze(Hvar);

% Read in the coordinate values
Rcoords = hdf5read(Hfile, '/x_coords');
Zcoords = hdf5read(Hfile, '/z_coords');
Tcoords = hdf5read(Hfile, '/t_coords');


end
