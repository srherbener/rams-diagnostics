function [ Rspect ] = RedNoiseSpectrum( Freqs, Efold )
%RedNoiseSpectrum Create red noise frequency spectrum
%   This function will create a frequency spectrum based on red noise given
%   the e-folding time and the desired frequencies (omega).

% The power spectrum for red noise is:
%
%  F(omega) = (2*Efold) / (1 + Efold^2*omega^2)
%

Rspect = (2*Efold) ./ (1 + (Efold^2 * Freqs.^2));

% Set the area under the curve to one
Rspect = SetAreaToOne(Freqs, Rspect);

end

