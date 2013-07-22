function [ OutData ] = PullOutCyclicData( Data, Select, Period )
%PullOutCyclicData pull out periodic entries from Data vector
%   This fuction will select data from the input vector Data and place the
%   selected data in the output vector OutData. The vector Select specifies
%   how to select the data. The function walks through Data in order, finds
%   the remainder of an integer divide of the index used on Data, and
%   selects the data if that index is a member of Select.
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

[ Dummy, NumPts ] = size(Data); % assume row vector
j = 0; % index for output array

for i = 1:NumPts
    if (ismember(rem(i,Period),Select))
        j = j + 1;
        OutData(j) = Data(i);
    end
end


end

