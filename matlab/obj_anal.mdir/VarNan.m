function [ Variance ] = VarNan( Data )
%VarNan calculate variance (2nd moment) or row vector
%   throw out NaNs

MeanVal = MeanNan(Data);
MeanOfSquares = MeanNan(Data .* Data);
N = length(Data(~isnan(Data)));

if (N > 1)
    Variance = (N / (N-1)) * (MeanOfSquares - (MeanVal*MeanVal));
else
    fprintf('WARNING: VarNan: Need more than one non NaN value in Data\n');
    fprintf('                     Setting variance to NaN\n');
    Variance = NaN;
end

end

