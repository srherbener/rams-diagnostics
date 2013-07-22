function [ Lambdas, VarExpl, Errors ] = EigenSpectrum( E, Nstar )
%EigenSpectrum For the eigenvalue spectrum
%   This function will form the eigenvalue spectrum from the input
%   eigenvalues matrix, E, and output the error ranges in the
%   output vector Errors. The specific eigenvalues are output in the
%   vector Lambdas. The error range is given by:
%
%      Error range = %variance explained * sqrt(2/Nstar)

Lambdas = diag(E);
VarExpl = (Lambdas / sum(Lambdas)) * 100; % %variance explained
Errors = VarExpl * sqrt(2/Nstar);

end

