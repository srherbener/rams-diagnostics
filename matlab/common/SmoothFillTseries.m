function [ OutTs ] = SmoothFillTseries( InTs, TsLen, FilterLen )
%SmoothFillTseries Smooth time series and fill to specified length
%   This function will take the input times series and:
%     1. smooth it using a running mean of length FilterLen
%     2. fill in at the end of the time series with nans so that the time
%        series length will match that given in TsLen
%

% smooth with a running mean of length 'Flen'
OutTs = filtfilt(ones(1,FilterLen)/FilterLen,1,double(InTs));
   
InTsLen = length(InTs);
if (InTsLen ~= TsLen)
    % fill with nans at the end
    OutTs(InTsLen+1:TsLen) = nan;
end

end

