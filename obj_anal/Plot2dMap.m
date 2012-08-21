function [ ] = Plot2dMap( Fig, X, Y, Z, Xlab, Ylab, Ptitle )
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

figure(Fig);

Nrows = length(Y);
Ncols = length(X);

% Find the largest absolute value of the entries in Dmap for setting the
% colormap axis. Want to center this about zero so that blue represents
% negative values and red represents positive values.
Clim = max(abs(Z)) * 1.1;  % need the range to be slightly
                            % larger than the actual data
Cbounds = [ -Clim Clim ];

contourf(X,Y,Create2dMap(Z,Nrows,Ncols));
shading flat;
title(Ptitle);
xlabel(Xlab);
ylabel(Ylab);
    
caxis(Cbounds);
colormap(redblue);
colorbar;
    
end

