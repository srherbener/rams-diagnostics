function [ ] = GenPlotEigSpectrum( EV, Nstar, N, ObsDescrip, Exp, OutFile )
%PlotEigSpectrum generate and plot out the eigenvalue spectrum that resulted from an EOF analysis
%  This function will plot out the first N eigenvalues plus error bars to show the
%  statistical significance of EOFs and PCs. Nstar is the estimate of the number
%  of independent samples.
%
%  The plot title gets built out of the two input strings:
%     ObsDescrip --> description of the observation data in the EOF analysis
%     Exp --> name of experiment

% Generate the eigen spectrum
[ Lambda, VarExpl, Err ] = EigenSpectrum(EV, Nstar);

% Fix up the title
ES_title = sprintf('First %d Eigenvalues of %s: %s', N, ObsDescrip, Exp);

Fig = figure;
errorbar(VarExpl(1:N), Err(1:N)); % just show the first N
title(ES_title);
xlabel('Eigenvalue number');
ylabel('% Variance Explained');

saveas(Fig, OutFile);
close(Fig);

end
