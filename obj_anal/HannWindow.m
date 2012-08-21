function [ Xout ] = HannWindow( X )
%HannWindow Apply Hann window onto a data series
%   This function will apply a Hann window onto the input data series X.
%   One instance of the Hann window will be applied over the entire length
%   of X.
%
%   The formula for the window is:
%
%     w(t) = 0.5 * (1 - cos(2*pi*t/T)
%
%   where t is the index into X and T is the length of x.

N = length(X);
t = (1:N);

Hann = 0.5 * (1 - cos((2*pi*t)./N));

Xout = X .* Hann;

end

