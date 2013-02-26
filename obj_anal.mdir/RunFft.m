function [ Fx ] = RunFft( X )
%RunFft Run FFT on data series
%   This function will run fft on X and return the first
%   half of the result (from 0 to the nyquist frequency).
%

% Base this on the built in FFT function, fft. fft returns the A's and B's
% for all wave numbers 1 to N. The Nyquist limit is at wave number N/2 so
% remove the second half of the results from fft. (These are just a repeat
% of the first half since they came from the aliased frequencies resulting
% from the N/2 frequencies beyond the Nyquist limit.)

HalfN = length(X) / 2; % Nyquist limit (length of spectrum we want to keep)
FullFft = fft(X);

Fx = FullFft(1:HalfN);

end

