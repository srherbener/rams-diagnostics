function [ Fx, Fy, CO, Q, Freqs ] = GenCrossSpectrum( X, Y, Window )
%GenCrossSpectrum Spectral analysis of input time series
%   This function will perform a corss spectral analysis on the input data
%   series, X and Y. It will return a complex vector containing the
%   individual spectra for X and Y (Fx and Fy respectively) plus the
%   cospectrum (CO) and quadrature spectrum (Q), along with the associated
%   frequencies.
%
%   Window controls what kind of a window is applied to the data series.
%     0 --> Boxcar window
%     1 --> Hann window

LenX = length(X);
LenY = length(Y);

if (LenX ~= LenY)
    fprintf('ERROR: GenCrossSpectrum: X and Y must have the same length\n');
    Fx = nan(1,LenX/2);
    Fy = nan(1,LenY/2);
    CO = nan(1,LenX/2);
    Q = nan(1,LenX/2);
    Freqs = nan(1,LenX/2);
    return;
end

% Apply the window to the entire series
if (Window == 0)
    Xwin = X;  % use series as is (boxcar window)
    Ywin = Y;
elseif (Window == 1)
    Xwin = HannWindow(X);
    Ywin = HannWindow(Y);
else
    fprintf('ERROR: GenCrossSpectrum: Must use "0" or "1" for Window argument\n');
    Fx = nan(1,LenX/2);
    Fy = nan(1,LenY/2);
    CO = nan(1,LenX/2);
    Q = nan(1,LenX/2);
    Freqs = nan(1,LenX/2);
    return;
end

% Run fft to get the A and B values (Fx will be complex with real part A,
% and imaginary part B) for Fx and Fy. Note that Fx and Fy will have half
% the length of X and Y due to the Nyquist restriction.
Fx = RunFft(Xwin);
Fy = RunFft(Ywin);
LenF = length(Fx);

% Generate the cross spectral quantities: CO and Q
Fxy = Fx .* conj(Fy);
CO = real(Fxy);
Q = imag(Fxy);

% Calculate the associated frequencies using: Freq = 2*pi*k/N. Note that the
% first entry in Mag is the mean value (k = 0, F = 0), the second entry is
% wavenumber 1, etc. Note too that the proper wavenumber/frequency is
% calculated from using the full length of the original series.
Freqs = zeros(1,LenF);
for i = 1:LenF
    k = i - 1;  % wavenumber
    Freqs(i) = (2*pi*k)/LenX;
end


end

