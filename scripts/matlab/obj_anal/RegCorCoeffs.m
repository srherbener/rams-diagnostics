function [ RegCoeff, CorCoeff ] = RegCorCoeffs( X, Y )
%RegCorCoeffs calculate regression and correlation coefficients
%   Throws out NaNs
%   Formula for correlation coefficient:
%     r = cov(X,Y) / (std(X)*std(Y))
%   Formula for regression coefficient:
%     a1 = r * (std(Y)/std(X))

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
else
    fprintf('ERROR: RegCorCoeffs: X and Y vectors must have the same length\n');
    fprintf('                     Setting RegCoeff and CorCoeff to NaN\n');
    RegCoeff = NaN;
    CorCoeff = NaN;
end


end

