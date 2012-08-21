function [ Cov ] = CovNan( X, Y )
%CovNan calculate covariance between X and Y vectors
%   This function returns the covariance between X and Y (which are assumed
%   to be row vectors). The formula is:
%      Cov = (1/(N-1)) * sum ((Xi - MeanX)*(Yi - MeanY))
%
%   NaNs are thrown out

if (length(X) == length(Y))
    Use = ~isnan(X) & ~isnan(Y); % use only the elements where both
                                 % X and Y contain non NaN values               
    Xuse = X(Use);
    Yuse = Y(Use);
    N = length(Xuse); % X and Y are equal size

    if (N > 1)
        MeanX = MeanNan(Xuse);
        MeanY = MeanNan(Yuse);
        Cov = (1/(N-1)) * sum((Xuse-MeanX).*(Yuse-MeanY));
    else
        fprintf('ERROR: CovNan: X and Y vectors must have more than one common pairs of non NaN elements\n');
        fprintf('                     Setting covariance to NaN\n');
        Cov = NaN;
    end
else
    fprintf('ERROR: CovNan: X and Y vectors must be the same length\n');
    fprintf('                     Setting covariance to NaN\n');
    Cov = NaN;
end


end

