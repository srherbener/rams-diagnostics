function [ Tseries ] = GenTseriesAvg2d( Var )
%GenTseriesAvg2d create time series of domain average of 2D fields
%   This routine will create a time series of the domain average value of
%   the 2D fields given by Var. It is assumed that Var has 3 dimensions
%   organized as (x,y,t).  The x,y portions are the 2D fields, and t is
%   the "series" dimension. The resulting time series is passed back
%   through the variable Tseries.
%

% This turns out to be quite simple to calculate. The matlab mean() function
% will calculate averages along a specified dimension in a multidimensional
% array. mean() can be used to first generate means along the columns of each
% 2D field, then to average these columns into a single number for each time
% step.
%
% squeeze() is used to elimnate the degenerate dimension (size = 1).

Tseries = squeeze(mean(mean(Var,1),2));

end

