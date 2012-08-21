function [ FS ] = LowPassButter( S, Wc, Ord )
%LowPassButter Apply low pass butterworth filter to input series
%   This function will apply a low pass butterworth fitler to the input
%   series, S, using Wc for the cutoff frequency (expressed in radians per
%   sample) and using the filter order given in Ord.

% Get the coefficients from the built in function: butter. The cutoff
% frequency given to butter, needs to be normalized to pi radians per
% sample. Ie, a number between 0 and 1 where 1 represent pi radians per
% sample.
[ b, a ] = butter(Ord, Wc/pi, 'low');

% Apply the filter to the input series
FS = filtfilt(b,a,S);


end

