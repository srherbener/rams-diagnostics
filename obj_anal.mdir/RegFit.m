function [ RegCoeff, RegYint, CorCoeff ] = RegFit( X, Y )
%RegFit calculate regression line and correlation coefficient
%   This routine will return the regression line (coefficient and
%   y-intercept) and the correlation coefficient from a least squares
%   linear fit. Also, this routine will throw out NaNs.
%
%   Formula for correlation coefficient:
%     r = cov(X,Y) / (std(X)*std(Y))
%
%   Formula for regression coefficient:
%     a1 = r * (std(Y)/std(X))
%
%   Formula for regression y-intercept:
%     a0 = mean(Y) - a1*mean(X)

if (length(X) == length(Y))
    Use = ~isnan(X) & ~isnan(Y); % use only the elements where both
                                 % X and Y contain non NaN values               
    Xuse = X(Use);
    Yuse = Y(Use);
    
    StdDevX = sqrt(VarNan(Xuse));
    StdDevY = sqrt(VarNan(Yuse));
    CovXY = CovNan(Xuse, Yuse);
    
    CorCoeff = CovXY / (StdDevX * StdDevY);
    RegCoeff = CorCoeff * StdDevY / StdDevX;
    RegYint = MeanNan(Y) - RegCoeff*MeanNan(X);
else
    fprintf('ERROR: RegCorCoeffs: X and Y vectors must have the same length\n');
    fprintf('                     Setting results to NaN\n');
    RegCoeff = NaN;
    CorCoeff = NaN;
    RegYint = NaN;
end


end

