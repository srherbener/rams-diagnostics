function [ ] = Plot2dMap( Fig, X, Y, Z, Clevs, Cbounds, Xlab, Ylab, Ptitle, SelectData )
%Plot2dMap Create a countour plot of the 2D data given in vector Z.
%   This function will extract the 2D data given in selected columns of
%   D, convert the column data to 2D data and creat contour plots of that
%   data as subplots in a single figure (specified by Fig).
%
%   Vector X contains the values for the x-axis, and vector Y contains the
%   values for the y-axis.
%
%   It is assumned that the data in Z is organized as:
%     first Ncols (length of X) items --> row = 1, columns 1 through Ncols
%     next Ncols items                --> row = 2, columns 1 through Ncols
%     etc.
%
%   SelectData contains [ x1 x2 y1 y2 ] which are used to select a
%   rectangular section out of the entire map.

Fsize = 20;

figure(Fig);

Nrows = length(Y);
Ncols = length(X);
Map = Create2dMap(Z,Nrows,Ncols);

% Trim out the rectangular region defined by SelectData
x1 = SelectData(1);
x2 = SelectData(2);
y1 = SelectData(3);
y2 = SelectData(4);
XP = X(x1:x2);
YP = Y(y1:y2);
MapP = Map(y1:y2,x1:x2);

% Find the largest absolute value of the entries in Dmap for setting the
% colormap axis. Want to center this about zero so that blue represents
% negative values and red represents positive values.

contourf(XP,YP,MapP,Clevs);
set(gca, 'FontSize', Fsize);
shading flat;
title(Ptitle);
xlabel(Xlab);
ylabel(Ylab);
    
caxis(Cbounds);
colormap(redblue);
cbar = colorbar;
set(cbar, 'FontSize', Fsize);
    
end

