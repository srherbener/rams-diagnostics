function [ ] = PlotEigSpectrum( VarExpl, Err, N, Ptitle, OutFile )
%PlotEigSpectrum plot out the eigenvalue spectrum that resulted from an EOF analysis
%  This function will plot out the first N eigenvalues plus error bars to show the
%  statistical significance of EOFs and PCs. Nstar is the estimate of the number
%  of independent samples.
%

Fig = figure;
errorbar(VarExpl(1:N), Err(1:N)); % just show the first N
title(Ptitle);
xlabel('Eigenvalue number');
ylabel('% Variance Explained');

saveas(Fig, OutFile);
close(Fig);

end
