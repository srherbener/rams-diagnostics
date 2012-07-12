function [ AvgProfile ] = DomainTimeAvgProfile ( Var )
%DomainTimeAvgProfile calculate domain time average profile of 3D field
%  This function will transform the given time series of a 3D field into
%  a single profile (in z) which was calculated by taking the horizontal
%  domain average at each level at each time step and then averaging across
%  time these resulting profiles.
%
%  The input variable is assumed to be of the form: (x,y,z,t)
%
%  The output variable will be of the form: (z)
%

% All we need to do is average across the x, y, and t dimensions of
% the input variable. This will reduce the result to a 1D vector
% in dimension z.
AvgProfile = squeeze(mean(mean(mean(Var,1),2),4))'; % transpose this so you return a row vector

end

