function [ Fx, Fy, CO, Q, Freqs ] = GenSplitCrossSpectrum( X, Y, Size, Incr, Window )
%GenSplitSpectrum Generates a "split" spectrum of the input data series, X
%   This function will take the input data X and split it up into a set of
%   smaller data sets (of length "Size"), compute a spectrum for each of
%   these subsets and return the entire set of specta in the matrix Fx
%   (each of Fx is one of the resulting spectra).
%
%   "Incr" says how far apart to space the subsets.
%
%   Example: if length(X) is 1000, Size is 100 and Incr is 50, then 19
%   subsets will be produced:
%     X(1:100)
%     X(51:150)
%     X(101:200)
%     ...
%     X(901:1000)
%
%   "Window" selects a window function:
%     0 --> Boxcar (no window)
%     1 --> Hann Window

if (length(X) ~= length(Y))
    fprintf('ERROR: GenSplitCrossSpectrum: length of X and Y must be equal\n');
    Fx = nan(1,length(X)/2);
    Fy = Fx; % nans
    CO = Fx;
    Q = Fx;
    Freqs = Fx;
    return;
end

N = length(X);
End = N - Size + 1;
NumSpectra = floor((End - 1) / Incr) + 1;
HalfSize = floor(Size/2);

% load up the rows of the output array with the individual spectra so that
% an average spectrum can be calculate with mean (which does column means).
% NOTE it is up to the caller to do the averaging since we need the
% individual spectra for cross spectral analysis.
Fx = zeros(NumSpectra,HalfSize);
Fy = zeros(NumSpectra,HalfSize);
CO = zeros(NumSpectra,HalfSize);
Q  = zeros(NumSpectra,HalfSize);
ir = 1;
for i = 1:Incr:End
    xstart = i;
    xend = xstart + Size - 1;
    Xsubset = X(xstart:xend);
    Ysubset = Y(xstart:xend);
    [ TFx, TFy, TCO, TQ, Freqs ] = GenCrossSpectrum(Xsubset, Ysubset, Window);
    Fx(ir,:) = [ TFx ];
    Fy(ir,:) = [ TFy ];
    CO(ir,:) = [ TCO ];
    Q(ir,:)  = [ TQ ];
    
    ir = ir + 1;
end

end

