function [ Dstring, Tstring ] = TimeToString(BtYr, BtMo, BtDay, BtHr, BtMin, BtSec, OffsetHrs)
% TimeToString function to translate numbers representing year/month/day/hr/min/sec into a date string
%
% This function will take numbers representing a base time (year, month, day, hour, minute, second)
% and an additional number that is an offset in hours from the base time and return a string
% containing the corresponding date/time.

% BaseTime will be in units of days.
BaseTime = datenum(BtYr, BtMo, BtDay, BtHr, BtMin, BtSec);

OffsetDays = OffsetHrs / 24;

Dstring = datestr(BaseTime + OffsetDays, 'mm/dd/yy');
Tstring = datestr(BaseTime + OffsetDays, 'HH:MM');

end
