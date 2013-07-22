function [ Zscore ] = FisherZscore( Rho, R, N )
% FisherZscore generate zscore resulting from the Fisher transformation
%
% This function will generate the z score resulting from the Fisher r to z
% transoformation of the given correlation values. This can be used to test
% hypotheses about the magnitude of the population correlation (argument Rho)
% given the correlation from a linear regression (argument R) that has
% been performed on a sample from the population. The argument N is the
% number of independent samples used in the linear regression.
% 
% The Fisher transform is done using:
%
%  Fz = 1/2 * ln((1+R)/(1-R))
%
% Fz is normally distributed with:
%    mean = 1/2 * ln((1+Rho)/(1-Rho))
%    standard deviation = 1 / sqrt(N-3)
%
% This function returns the z-score obtained from Fz and the generated
% mean and standard deviation as shown above. If N is < 3, a nan is returned
% to signal that the z score could not be produced.
%
% Since Fz is a symmetric distribution and since we are interested in testing
% the magnitude of correlations (which can be negative or positive), this
% function will perform the transformation on the absolute values of Rho and
% R which will always yield the z score on the right hand side of the distribution
% (positive z scores).
%
% Note that R can be an array filled with correlation coefficients from multiple
% linear regressions.

R = abs(R);
Rho = abs(Rho);

Fz = 0.5 * log((1 + R) ./ (1 - R));

MuFz = 0.5 * log((1 + Rho) / (1 - Rho));
if (N > 3)
  SigmaFz = 1 / sqrt(N - 3);
else
  SigmaFz = nan;
end

% generate the z score in the typical fashion using the
% sample mean and std deviation (MuFz and SigmaFz).
Zscore = (Fz - MuFz) / SigmaFz;

end
