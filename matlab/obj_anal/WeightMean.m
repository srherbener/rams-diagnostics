function [ Wmean, N ] = WeightMean( Data, Weights )
%WeigthMean compute weigthed mean of row vectors
%   Compute weighted mean, throw out NaNs from
%   calculations.
%
%   mean = sum(weight*data) / sum(weight)

if (length(Data) == length(Weights))
    Use = ~isnan(Data) & ~isnan(Weights); % use only the elements where both
                                          % Data and Weights contain non NaN values               
    DataUse = Data(Use);
    WeightsUse = Weights(Use);
    N = length(DataUse);
    
    if (N ~= 0)
        Wmean = sum(DataUse .* WeightsUse) / sum(WeightsUse);
    else
        fprintf('ERROR: WeightMean: Data and Weights vectors must have at least one common pair of non NaN elements\n');
        fprintf('                   Setting weighted mean to NaN\n');
        Wmean = NaN;
    end
else
    fprintf('ERROR: WeightMean: Nubmer of entries in Data and Weights do not match:\n');
    fprintf('                   Number of entries in Data: %d\n', NumD);
    fprintf('                   Number of entries in Weights: %d\n', NumW);
    fprintf('                   Setting weighted mean to NaN\n');
    Wmean = NaN;
end

end

