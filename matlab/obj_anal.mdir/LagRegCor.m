function [ RegCoeff, LagCoeff  ] = LagRegCor( X, Y, Lag, Start, Size )
% LagRegCor calculate lag correlation and lag regression coefficient
% on a time series
%   This function will perform regression/correlation calculations on X and
%   Y where the corresponding X and Y pairs can be offest from each other
%   in time. The args Lag, Start, and Size are used to specify this offset.
%
%   X - independent variable
%
%   Y - dependent variable
%
%   Lag - specifies the offset in the X (independent var) series will have
%         with respect to the the Y (dependent var) series. Lag of -1 means
%         that X will preceed Y by 1 time step; Lag of +3 means that X will
%         follow Y by 3 time steps, etc.
%
%   Start - specifies where the starting point (time step) is located in
%           the Y time series. The X time series starting point will be
%           offset from the location 'Start'.
%
%   Size - specifies how many oncsecutive time points from X and Y to use
%          in the regression/correlation calculation.
%
%   Lag, Start, Size are used as arguments so the caller can do sweeps on
%   the Lag values. If you want to compare say lags going from 1 to 20, you
%   want Size to be the same each time. One could vary the size in an
%   attempt to maximize Size for each lag, but this results in all time
%   points being used for a lag of 0, n-1 time points for a lag of 1, n-2
%   time points for a lag of 2, n-lag time points in general. Using a
%   different nubmer of time points for different lags messes up the
%   attempt to compare the resulting lag correlations.

% set indices into Data for pulling out the two samples
ixstart = Start + Lag;
ixend = ixstart + (Size-1);

iystart = Start;
iyend = iystart + (Size-1);

if ((ixstart < 1) || (ixend > length(X)))
    fprintf ('ERROR: AutoCor: indices into X are not within the bounds of X\n');
    fprintf ('ERROR:   ixstart: %d\n', ixstart);
    fprintf ('ERROR:   ixend: %d\n', ixend);
    fprintf ('ERROR:   bounds of X: 1 %d\n', length(X));
    LagCoeff = NaN;
elseif ((iystart < 1) || (iyend > length(Y)))
    fprintf ('ERROR: AutoCor: indices into Y are not within the bounds of Y\n');
    fprintf ('ERROR:   iystart: %d\n', iystart);
    fprintf ('ERROR:   iyend: %d\n', iyend);
    fprintf ('ERROR:   bounds of Y: 1 %d\n', length(Y));
    LagCoeff = NaN;
else
    % RegCoeff is a dummy variable
    [ RegCoeff, LagCoeff ] = RegCorCoeffs(X(ixstart:ixend),Y(iystart:iyend));
end

end

