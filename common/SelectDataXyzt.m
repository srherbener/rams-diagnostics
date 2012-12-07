function [ SelectedData ] = SelectDataXyzt(Data, X, Y, Z, T, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax)
% SelectDataXyzt select subset out of 4D var (x,y,z,t) based on ranges
%
% This function will apply selection based on ranges given by Xmin, Xmax, ... on the input
% data array (Data) and return the selected subset. X, Y, Z, T contain the coordinate values
% that will be compared against the specified ranges.

  X1 = find(X >= Xmin, 1, 'first');
  X2 = find(X <= Xmax, 1, 'last');
  Y1 = find(Y >= Ymin, 1, 'first');
  Y2 = find(Y <= Ymax, 1, 'last');
  Z1 = find(Z >= Zmin, 1, 'first');
  Z2 = find(Z <= Zmax, 1, 'last');
  T1 = find(T >= Tmin, 1, 'first');
  T2 = find(T <= Tmax, 1, 'last');

  SelectedData = squeeze(Data(X1:X2, Y1:Y2, Z1:Z2, T1:T2));
end
