function [ Fx, Freqs ] = GenSpectrum( X, Window )
%GenSpectrum Spectral analysis of input time series
%   This function will perform a spectral analysis on the input data
%   series, X. It will return a complex vector containing the spectrum
%   along with the associated frequencies.
%
%   Window controls what kind of a window is applied to the data series.
%     0 --> Boxcar window
%     1 --> Hann window

LenX = length(X);

% Apply the window to the entire series
if (Window == 0)
    Xwin = X;  % use series as is (boxcar window)
elseif (Window == 1)
    Xwin = HannWindow(X);
else
    fprintf('ERROR: GenSpectrum: Must use "0" or "1" for Window argument\n');
    Fx = nan(1,LenX/2);
    Freqs = nan(1,LenX/2);
    return;
end

% Run fft to get the A and B values (Fx will be complex with real part A,
% and imaginary part B). Note that Fx will have half the length of X due to
% the Nyquist restriction.
Fx = RunFft(Xwin);
LenFx = length(Fx);

% Calculate the associated frequencies using: Freq = 2*pi*k/N. Note that the
% first entry in Mag is the mean value (k = 0, F = 0), the second entry is
% wavenumber 1, etc. Note too that the proper wavenumber/frequency is
% calculated from using the full length of the original series.
Freqs = zeros(1,LenFx);
for i = 1:LenFx
    k = i - 1;  % wavenumber
    Freqs(i) = (2*pi*k)/LenX;
end


end

