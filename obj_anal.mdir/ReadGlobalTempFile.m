function [ Temps, NumRows, NumCols, Weights ] = ReadGlobalTempFile( TempsFile )
%ReadGlobalTempFile reads global temp file and returns data plus weights
%
% Input file is 612 x 2592
%   Rows:
%     1 - Jan, 1958
%     2 - Feb, 1958
%     ...
%     612 - Dec, 2008
%
%   Columns:
%     1 - 87.5S, 2.5E
%     2 - 87.5S, 7.5E
%     ...
%    72 - 87.5S, 357.5E
%    
%    73 - 82.5S, 2.5E
%    74 - 82.5S, 7.5E
%     ...
%   144 - 82.5S, 357.5E
%
%   ...
%   2592 - 87.5N, 357.5E
%
%   First 72 are for 87.5S, going from 2.5E thru 357.5E
%   Next 72 are for 82.5S, "
%   Last 72 are for 87.5N, "
%
%   Increments for lat and lon are both 5 degrees.

% data file is space delimited, very convenient for importdata()
Temps = importdata(TempsFile);

[ NumRows, NumCols ] = size(Temps);

% Walk through one row, weights are cos(latitude)
for i = 1:NumCols
    % want i = 1..72 to return zero from floor() so subtract 1 from i
    % to "push i back one integer"
    LatNum = floor((i-1)/72);
    LatAngle = -87.5 + (LatNum * 5.0);
    %Need to convert degrees to radians
    LatAngle = (LatAngle * pi) / 180.0;
    Weights(i) = cos(LatAngle);
end

end

