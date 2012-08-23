function [ ] = Plot2dEofPc( EOF, PC, N, X, Y, ObsDescrip, Exp, E_xlabel, E_ylabel, P_ylabel, OutFile )
%Plot2dEofPc plot out a specified EOF and PC of a set of 2D observations
%   This function will plot out the Nth EOF and PC of set of 2D
%   observations. The EOFs and PCs are in the columns of their respective
%   input arrays. X and Y are the coordinate values for the 2D field.

E_title = sprintf('%s EOF%d: %s', ObsDescrip, N, Exp);
P_title = sprintf('%s PC%d: %s', ObsDescrip, N, Exp);

Fig = figure;
subplot(2,1,1);
Plot2dMap(Fig,X,Y,EOF(:,N)',E_xlabel, E_ylabel, E_title);
subplot(2,1,2);
plot(PC(:,N)');
title(P_title);
xlabel('Timestep');
ylabel(P_ylabel);

saveas(Fig, OutFile);
close(Fig);

end
