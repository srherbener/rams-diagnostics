function [ ] = PlotPowerSpectrum( Fig, Pspect, RnSpect, Freqs, F95, Ptitle )
%PlotPowerSpectrum Plot the power spectrum with red noise fit
%   This function will plot the power spectrum Pspect with the red noise
%   fit, RnSpect, overlaid including the 95% significance line.

figure(Fig);
semilogy(Freqs,Pspect);
line(Freqs,RnSpect,'Color','r');
line(Freqs,F95*RnSpect,'Color', 'r', 'LineStyle', '--');
title(Ptitle);
xlabel('Frequency (radians/timestep)');
ylabel('log(Phi(omega))');
legend('Spectrum', 'Red Noise Fit', 'Red Noise Fit: 95% Significance');

end

