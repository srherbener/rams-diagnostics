function [ Eindex, NumPts ] = ReadEnsoFile( EnsoFile )
%ReadEnsoFile read the ENSO index file and load into a vector
%
% Input file format:
%   Line 1 - Jan, 1958
%   Line 2 - Feb, 1958
%   ...
%   Line 612 - Dec, 2008
%
% Each line has one number which is the ENSO index for that month

% importdata() will create a column vector since the lines are treated
% as rows in a matrix and there is one index per line --> transpose
% Eindex after running importdata.

Eindex = importdata(EnsoFile);
Eindex = Eindex';
NumPts = length(Eindex); 

end

