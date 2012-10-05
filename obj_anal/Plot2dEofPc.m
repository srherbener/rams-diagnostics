function [ ] = Plot2dEofPc(EOF, PC, X, Y, Clevs, Cbounds, E_Title, E_Xlabel, E_Ylabel, P_Title, P_Xlabel, P_Ylabel, SelectData, OutFile)
%Plot2dEofPc plot out a specified EOF and PC of a set of 2D observations
%   This function will plot out the Nth EOF and PC of set of 2D
%   observations. The EOFs and PCs are in the columns of their respective
%   input arrays. X and Y are the coordinate values for the 2D field.
%
%   SelectData contains [ x1 x2 y1 y2 ] which are used to select a rectangular
%   region out of the entire 2D map.

Fig = figure;
subplot(2,1,1);
Plot2dMap(Fig,X,Y,EOF, Clevs, Cbounds, E_Xlabel, E_Ylabel, E_Title, SelectData);
subplot(2,1,2);
plot(PC);
title(P_Title);
xlabel(P_Xlabel);
ylabel(P_Ylabel);

saveas(Fig, OutFile);
close(Fig);

end
