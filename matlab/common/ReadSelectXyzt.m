function [ OutData ] = ReadSelectXyzt(Hfile, Vname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax)
% ReadSelectXyzt read 4D var (x,y,z,t) out of HDF5 file, and select subset based on ranges
%
% This function will apply selection based on ranges given by Xmin, Xmax, ... on the input
% data array (Data) and return the selected subset. X, Y, Z, T contain the coordinate values
% that will be compared against the specified ranges.
%
% The special coordinate names (x_coords, y_coords, z_coords, t_coords) and handled as
% one dimensional arrays.

  Hdata = hdf5read(Hfile, Vname);
  switch Vname
    case 'x_coords'
      I1 = find(Hdata >= Xmin, 1, 'first');
      I2 = find(Hdata <= Xmax, 1, 'last');

      OutData = squeeze(Hdata(I1:I2));
    case 'y_coords'
      I1 = find(Hdata >= Ymin, 1, 'first');
      I2 = find(Hdata <= Ymax, 1, 'last');

      OutData = squeeze(Hdata(I1:I2));
    case 'z_coords'
      I1 = find(Hdata >= Zmin, 1, 'first');
      I2 = find(Hdata <= Zmax, 1, 'last');

      OutData = squeeze(Hdata(I1:I2));
    case 't_coords'
      I1 = find(Hdata >= Tmin, 1, 'first');
      I2 = find(Hdata <= Tmax, 1, 'last');

      OutData = squeeze(Hdata(I1:I2));
    otherwise
      X = hdf5read(Hfile, 'x_coords');
      Y = hdf5read(Hfile, 'y_coords');
      Z = hdf5read(Hfile, 'z_coords');
      T = hdf5read(Hfile, 't_coords');

      X1 = find(X >= Xmin, 1, 'first');
      X2 = find(X <= Xmax, 1, 'last');
      Y1 = find(Y >= Ymin, 1, 'first');
      Y2 = find(Y <= Ymax, 1, 'last');
      Z1 = find(Z >= Zmin, 1, 'first');
      Z2 = find(Z <= Zmax, 1, 'last');
      T1 = find(T >= Tmin, 1, 'first');
      T2 = find(T <= Tmax, 1, 'last');

      OutData = squeeze(Hdata(X1:X2, Y1:Y2, Z1:Z2, T1:T2));
  end
end
