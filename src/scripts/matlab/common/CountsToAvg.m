function [ Avg, Npts ] = CountsToAvg( HdaData, T1, T2 )
% CountsToAvg will convert the (sum,count) data from the tsavg 'hda' format into
%             the average value (sum/count)
%

% HdaData comes from hda being run on either 2D or 3D field data.
%
%  2D: HdaData
%      (1,t) - sums for each time step
%      (2,t) - counts for each time step
%
%      where t is time
%
%  3: HdaData
%     (1,z,t) - sums
%     (2,z,t) - counts
%
%     where z is height, t is time

  if (ndims(HdaData) == 3)
    Sum = squeeze(sum(HdaData(1,:,T1:T2) ,3));
    Npts = squeeze(sum(HdaData(2,:,T1:T2) ,3));
  elseif (ndims(HdaData) == 2)
    Sum = squeeze(sum(HdaData(1,T1:T2) ,2));
    Npts = squeeze(sum(HdaData(2,T1:T2) ,2));
  else
    % this will produce a nan in Avg
    fprintf('CountsToAvg: ERROR: can only handle 2D or 3D arrays: setting Avg and Npts to zero\n');
    Sum = 0;
    Npts = 0;
  end

  Avg = Sum ./ Npts;

end
