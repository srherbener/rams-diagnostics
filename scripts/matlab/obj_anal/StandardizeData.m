function [ StdData ] = StandardizeData( Data )
%StandardizeData standardize data points
%   standaized data will have mean of zero and std dev of one
%   formula:
%     StdX(i) = (Data(i) - mean(Data)) / stddev(Data)
%
%   throw out NaNs

DataUse = Data(~isnan(Data)); % strip out NaN entries

MeanD = MeanNan(DataUse);
StdDevD = sqrt(VarNan(DataUse));

StdData = (Data - MeanD) / StdDevD; % return original data length

end

