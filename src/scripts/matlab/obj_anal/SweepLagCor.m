function [ LagCors, LagRcoeffs, LagRyints ] = SweepLagCor( X, Y, Lags )
%SweepLagCor Run a series of lag correlations sweeping the lag value
%   This function will perform a series of lag correlations using the lag
%   values in the vector Lags, and return the correlation coefficients in
%   the output vector LagCors.

% Want to run through a series of lag correlations between X and Y. Figure
% out how to get the longest samples (of X and Y) so that the samples are
% the same size for all lag correlations, and that size is a maximum.
% Along with this figure out where the starting point is in the independent
% variable (X).

LenX = length(X);
LenY = length(Y);

if (LenX == LenY)

    MaxLag = max(Lags);
    MinLag = min(Lags);
    NumLags = length(Lags);
    
    if (MinLag >= 0)
        % All lag values are positive we always want the dependent
        % variable array to start at 1, and the independent variable array
        % to start at 1 + Lag.
        % Since Start is 1, the Size is just the length of the data arrays
        % minus MaxLag (MaxLag is assumed > MinLag, so MaxLag will be
        % positive.)
        Start = 1;
    else
        % There are negative lag values, so for the greatest negative lag,
        % we want to start at index 1 of the independent variable array which
        % means we want Start + MinLag = 1. (Start is the first index value
        % in the dependent variable array.)
        Start = 1 - MinLag;
    end

    if ((MinLag * MaxLag) > 0)
        % All lags are of the same sign. Then the size of the sample needs
        % to be the length of the data minus the greater of the absolute
        % values of min and max lags.
        Size = length(X) - max(abs(MinLag),abs(MaxLag));
    else
        % min and max lags have opposite signs. Then the size of the sample
        % needs to be the length of the data minus the sum of the absolute
        % values of min and max lags
        Size = length(X) - (abs(MinLag) + abs(MaxLag));
    end

    LagCors = zeros(1,NumLags);
    LagRcoeffs = zeros(1,NumLags);
    LagRyints = zeros(1,NumLags);
    for i = 1:NumLags
        [ LagRcoeffs(i), LagRyints(i), LagCors(i) ] = LagRegFit(X,Y,Lags(i),Start,Size);
    end
else
    fprintf('ERROR: SweepLagCor: lengths of X and Y must be equal\n');
    LagRcoeffs = nan(1,length(X));
    LagRyints = LagRcoeffs;
    LagCors = LagRcoeffs;
end

end

