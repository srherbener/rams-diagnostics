function [ OutVector ] = CreateCyclicPattern( NumPts, Select, Period )
%PullOutCyclicData pull out periodic entries from Data vector
%   This function will create a row vector with a repeating pattern of
%   1's and 0's. The vector will be NumPts long, with a repeating pattern
%   specified by Select and Period. The length of the pattern will be
%   Period; and 1's appear wherever the remainder after dividing the vector
%   index by Period is one of the values in Select.
%
%   Example: Say you have monthly data points over multiple years starting
%   with January. Then the points in Data will look like:
%       index     Corresponding Month
%         1             Jan
%         2             Feb
%         3             Mar
%         ...           ...
%        11             Nov
%        12             Dec
%        13             Jan
%        14             Feb
%        ...            ...
%        24             Dec
%        25             Jan
%        ...            ...
%
%   Now say you want to pull out Dec,Jan,Feb data. Then set:
%      Period = 12
%      Select = [0,1,2]
%
%   Then when index is 12, 24, 36, ... the remainder after dividing the
%   index by 12 will be 0 thus selecting all the Dec points. Same type of
%   thing happens for Jan when the remainder is 1, Feb when the remainder
%   is 2.

for i = 1:NumPts
    if (ismember(rem(i,Period),Select))
        OutVector(i) = 1;
    else
        OutVector(i) = 0;
    end
end


end

