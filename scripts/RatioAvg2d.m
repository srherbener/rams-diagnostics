function [ AvgOfRatios, RatioOfAvgs ] = RatioAvg2d ( VarN, VarD )
%RatioAvg2d calculate domain average of the ratio of two 2D fields
%   This routine will calculate the domain average value of the 2D fields given
%   by VarN and VarD. VarN is the numerator, and VarD is the denominator
%   of the ratios.
%
%   There are two ways to do this calculation so both methods are used and
%   their results passed back to the caller.
%     1. take the average of the elementwise ratios in the 2D fields
%     2. take the ratio of the domain averages of the 2D fields

% Average of ratios:
%
% Want to avoid division by zero, so use find() to eliminate using locations
% that contain zeros in the denominator array.
%
% Locate the non zero elements in the denominator and use this to form a vector
% of ratios of all the paired up elements of VarN and VarD where the VarD value
% is non zero. Then take the mean of that vector.
NonZ = find(VarD ~= 0);
AvgOfRatios = mean(VarN(NonZ) ./ VarD(NonZ));

% Ratio of averages:
%
% Calculate the mean of each 2D field first. mean() calculated averages along only
% one dimension per call, so nest two calls together: one for the rows and the other
% for the columns.
MeanN = mean(mean(VarN,1),2);
MeanD = mean(mean(VarD,1),2);

if (MeanD == 0)
  RatioOfAvgs = NaN;  % use NaN to let caller know that we tried to divide by zero
else
  RatioOfAvgs = MeanN / MeanD;
end

end

