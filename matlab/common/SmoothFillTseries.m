function [ OutTs ] = SmoothFillTseries( InTs, TsLen, FilterLen )
%SmoothFillTseries Smooth time series and fill to specified length
%   This function will take the input times series and:
%     1. smooth it using a running mean of length FilterLen
%     2. fill in at the end of the time series with nans so that the time
%        series length will match that given in TsLen
%

% InTs might have nans in it signifying missing data. To handle this case, do the
% filtering with the nan entries removed, then insert the nans back into the output
% filtered series in their respective places. This should work reasonably well given that
% we are doing a running mean filter. This method may not be a good one for other kinds
% of filters.

InTsLen = length(InTs);

Use = ~isnan(InTs);
InTsNoNan = InTs(Use);

% smooth with a running mean of length 'Flen'
% run filter only on the input series with the nans removed
% stick the nans back into the output series in the same locations that they appeared in the input series
% run the filter in both directions in order to create a zero-phase result
Bcoeff = ones(1,FilterLen)./FilterLen;
Acoeff = 1;
TempTs = filter(Bcoeff,Acoeff,double(InTsNoNan));
TempTs2 = filter(Bcoeff,Acoeff,TempTs(end:-1:1));
TempTs = TempTs2(end:-1:1);

OutTs = nan([ 1 InTsLen ]);
OutTs(Use) = TempTs;
   
if (InTsLen ~= TsLen)
    % fill with nans at the end
    OutTs(InTsLen+1:TsLen) = nan;
end

end

