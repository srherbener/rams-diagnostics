function [ Mean ] = MeanNan( Data )
%MeanNan compute mean value of row vector
%   Throw out NaN's

Use = ~isnan(Data);
SumData = sum(Data(Use));
N = length(Data(Use));

if (N ~= 0)
    Mean = SumData / N;
else
    fprintf('WARNING: MeanNan: All data are NaN, setting mean to NaN\n');
    Mean = NaN;
end

end

