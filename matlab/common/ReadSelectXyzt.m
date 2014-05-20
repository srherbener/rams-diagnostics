function [ OutData ] = ReadSelectXyzt(Hfile, Vname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax)
% ReadSelectXyzt read 4D var (x,y,z,t) out of HDF5 file, and select subset based on ranges
%
% This function will apply selection based on ranges given by Xmin, Xmax, ... on the input
% data array (Data) and return the selected subset. X, Y, Z, T contain the coordinate values
% that will be compared against the specified ranges.
%
% The special coordinate names (x_coords, y_coords, z_coords, t_coords) and handled as
% one dimensional arrays.
%
% Can't guarantee that variables will always be 4D. MATLAB won't allow an extra dimension of
% size 1 to be tagged onto the end. Eg., if you have x,y,z from one time step then you can't
% form the variable into a 4D (x,y,z,t) with t being a dimension with size one.
%
% To handle this case, always squeeze the input variable and require that the user set the
% select range for the unused variables to nans.

  Hdata = squeeze(hdf5read(Hfile, Vname));
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

      % Form a string with the command for selecting the data, then eval that string
      SelString = '';
      if (~isnan(Xmin) & ~isnan(Xmax))
        X1 = find(X >= Xmin, 1, 'first');
        X2 = find(X <= Xmax, 1, 'last');
        if (strcmp(SelString, ''))
          SelString = sprintf('%d:%d', X1, X2);
        else
          SelString = sprintf('%s, %d:%d', SelString, X1, X2);
        end
      end
      if (~isnan(Ymin) & ~isnan(Ymax))
        Y1 = find(Y >= Ymin, 1, 'first');
        Y2 = find(Y <= Ymax, 1, 'last');
        if (strcmp(SelString, ''))
          SelString = sprintf('%d:%d', Y1, Y2);
        else
          SelString = sprintf('%s, %d:%d', SelString, Y1, Y2);
        end
      end
      if (~isnan(Zmin) & ~isnan(Zmax))
        Z1 = find(Z >= Zmin, 1, 'first');
        Z2 = find(Z <= Zmax, 1, 'last');
        if (strcmp(SelString, ''))
          SelString = sprintf('%d:%d', Z1, Z2);
        else
          SelString = sprintf('%s, %d:%d', SelString, Z1, Z2);
        end
      end
      if (~isnan(Tmin) & ~isnan(Tmax))
        T1 = find(T >= Tmin, 1, 'first');
        T2 = find(T <= Tmax, 1, 'last');
        if (strcmp(SelString, ''))
          SelString = sprintf('%d:%d', T1, T2);
        else
          SelString = sprintf('%s, %d:%d', SelString, T1, T2);
        end
      end

      SelCmd = sprintf('squeeze(Hdata(%s))', SelString);
      OutData = eval(SelCmd);
  end
end
