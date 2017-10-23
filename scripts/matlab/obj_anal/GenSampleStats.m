function [ SampleMeans, SampleStdDevs ] = GenSampleStats( X, NumSamples, SampleSize )
%GenSampleStats Generate sample means and sample standard deviations
%   This function treats X as the population. A sample of SampleSize
%   consecutive elements are taken from X and the mean and standard
%   deviation of that sample are calculated. This is repeated NumSamples
%   times and the resulting mean/stddev pairs are returned in SampleMeans
%   and SampleStdDevs.

N = length(X);
SsizeM1 = SampleSize - 1;
RandiMax = N - SsizeM1; % randi(RandiMax) is being used in the loop
                        % below to select the starting point in X
                        % for a sample. The sample needs to be
                        % SampleSize consecutive points which means
                        % that the indices of X (for the sample) need
                        % to go from the randomly selected starting
                        % value plus (SampleSize - 1).

SampleMeans = zeros(1,NumSamples);
SampleStdDevs = zeros(1,NumSamples);
for i = 1:NumSamples
    Istart = randi(RandiMax);
    Iend = Istart + SsizeM1;
    
    S = X(Istart:Iend);
    SampleMeans(i) = MeanNan(S);
    SampleStdDevs(i) = sqrt(VarNan(S));
end


end

