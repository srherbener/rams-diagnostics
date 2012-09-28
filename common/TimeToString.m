function [ Dstring, Tstring ] = TimeToString(BtYr, BtMo, BtDay, BtHr, BtMin, BtSec, OffsetHrs)
% TimeToString function to translate numbers representing year/month/day/hr/min/sec into a day and time string
%
% This function will take numbers representing a base time (year, month, day, hour, minute, second)
% and an additional number that is an offset in hours from the base time and return strings
% containing the corresponding day and time.

% BaseTime will be in units of days.
BaseTime = datenum(BtYr, BtMo, BtDay, BtHr, BtMin, BtSec);

OffsetDays = OffsetHrs / 24;

Dstring = datestr(BaseTime + OffsetDays, 'dd');
Tstring = datestr(BaseTime + OffsetDays, 'HH:MM');

end
