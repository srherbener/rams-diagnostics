function [ Tstat, DegOfFree ] = DiffOfMeansTstat( X1, X2 )
%DiffOfMeansTstat Calculate Student's t-statistic for difference of two
%means
%   This function takes the two data sets given by X1, and X2 and generates
%   the Student's t-statistic for the difference of the means of these data
%   sets. NaN's are thrown out.

X1bar = MeanNan(X1);
X1var = VarNan(X1);
N1 = length(X1(~isnan(X1))); %MeanNan and VarNan throw out nans

X2bar = MeanNan(X2);
X2var = VarNan(X2);
N2 = length(X2(~isnan(X2)));

DegOfFree = N1 + N2 - 2;

if (DegOfFree > 0);
    JointVar = ((N1*X1var) + (N2*X2var)) / DegOfFree;
    Tstat = (X1bar - X2bar) / sqrt(JointVar * (1/N1 + 1/N2));
else
    fprintf ('WARNING: DiffOfMeansTstat: Degrees of Freedom must be > 0: %d\n', DegOfFree);
    fprintf ('         Setting Tstat to NaN\n');
    Tstat = NaN;
end
end

