function [ TsAvgOfRatios, TsRatioOfAvgs ] = GenRatioTsAvg2d( VarN, VarD )
%GenRatioTsAvg2d create time series of domain average of the ratio of two 2D fields
%   This routine will create a time series of the domain average value of
%   the 2D fields given by VarN and VarD. It is assumed that VarN and VarD 
%   have 3 dimensions organized as (x,y,t).  The x,y portions are the 2D fields,
%   and t is the "series" dimension. VarN is the numerator, and VarD is the
%   denominator of the ratios.
%
%   There are two ways to do this calculation so both methods are used and
%   their results passed back to the caller.
%     1. take the average of the elementwise ratios in the 2D fields
%     2. take the ratio of the domain averages of the 2D fields

% Want to avoid division by zero, so use find() to eliminate using locations
% that contain zeros in the denominator array.

Nt = size(VarN,3);  % number of time steps
TsAvgOfRatios = zeros(1,Nt);
TsRatioOfAvgs = zeros(1,Nt);

% find returns the linear indexes, so we need to look at the 2D fields at each
% time step separately.
for i = 1:Nt
  % pick out the 2D fields at this time step
  N = squeeze(VarN(:,:,i));
  D = squeeze(VarD(:,:,i));

  % Average of ratios:
  %
  % Locate the non zero elements in the denominator and use this to
  % form a vector of ratios of all the paired up elements of N and
  % D where the D value is non zero. Then take the mean of that vector.
  NonZ = find(D ~= 0);
  TsAvgOfRatios(i) = mean(N(NonZ) ./ D(NonZ));

  % Ratio of averages:
  %
  MeanN = mean(mean(N,1),2);
  MeanD = mean(mean(D,1),2);

  if (MeanD == 0)
    TsRatioOfAvgs(i) = NaN;  % use NaN to let caller know that we tried to divide by zero
  else
    TsRatioOfAvgs(i) = MeanN / MeanD;
  end
end

end

