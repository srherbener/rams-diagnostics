function [ Tstd ] = RedNoise( R1, Length )
%RedNoise Generate red noise time series
%   This function will generate a time series with Length number of points
%   using the lag-1 autocorrelation coefficient given in R1. The formula
%   used is:
%
%     X(t+1) = R1*x(t) + sqrt(1-R1^2)*e(t)
%       where e(t) is a random value selected from a standard normal
%       distribution (mean = 0, stddev = 1). randn() will supply
%       these random values.
%
%   The returned time series will be standardized.

T = zeros(1,Length); % allocate space outside loop so that run time
                     % doesn't get bogged down resizing the array
                     % over and over again in the loop

T(1) = 0;
B = sqrt(1-R1^2);
for i = 2:Length
   T(i) = (R1*T(i-1)) + (B*randn(1,1));
end

% standardize the time series
Tstd = StandardizeData(T);

end

