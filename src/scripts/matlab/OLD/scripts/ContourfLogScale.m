function [ ] = ContourfLogScale( Fig, X, Y, Z, Clevs )
%ContourfLogScale create contourf plot with log scale contours
%   This function will convert the data given in Z and Clevs to log values
%   and then plot as a filled contour plot.

figure(Fig);
hold on;

Zlog = log10(Z);

% First plot out original values
contourf(X,Y,Z,Clevs);


end

