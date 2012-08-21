function [ OutVector ] = CreateCyclicIndices( NumPts, Select, Period )
%CreateCyclicIndices create vector indices that point to periodic locations
%
%   This function will create a row vector with a list of index values
%   that point to a repeating pattern of locations from which to select
%   data. The repeating pattern is specified by Select and Period. Period
%   is how many indices to go before repeating, and Select shows which
%   indices to pick for the output vector.
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

j = 0; % index for the output array

for i = 1:NumPts
    if (ismember(rem(i,Period),Select))
        j = j + 1;
        OutVector(j) = i;
    end
end


end

