function [ EOF, PC, Evals ] = EofNan( A, ObsInCols)
%EofNan Performs EOF analysis and handles missing data
%   This function will perform EOF analysis of a data set given in matrix A
%   and return the EOFs and PCs. It will handle the case where A contains
%   missing data (NaNs).
%
%   Since the built in svd command in matlab does not handle NaNs, need to
%   compute the covariance matrix (A*A' or A'*A) with MatMultNan which will
%   handle (and eliminate) NaNs. Then eigenanalye the covariance matrix to
%   get the EOFs and PCs.
%
%   The flag ObsInCols should be set to 1 if the columns represent
%   observed values (as opposed to the columns representing the instances
%   of sampling). If the rows represent observed values (data), then set
%   ObsInCols to 0. For example if the rows are time and the columns space,
%   set ObsInCols since the spatial data (at one time sample) are
%   in the columns.
%

% The eigenvectors of the sample covariance matrix give you the EOFs,
% then the PCs are obtained using: PC = A * EOF
% The eigenvectors of the data covariance give you the PCs,
% then the EOFs are obtained using: EOF = A' * PC

% Figure out which order to multiply A and A', ie pick the one that yields
% the smaller covariance matrix.
% 
% A*A' yields a square matrix with size "number of rows in A"
% A'*A yields a square matrix with size "number of columns in A"
%
% For effeiciency, want to do the analysis on the smaller covairance
% matrix. Keep track of whether or not the transpose was first in the
% matrix multiply so we can sort out which covariance matrix we just
% calculated.
[ Nrows, Ncols ] = size(A);
if (Nrows > Ncols)
    TransposeFirst = true;
    [ C, Counts ] = MatMultNan(A',A);
    % Only keep the covariance value if there were at least half of the
    % element-wise pairs in the dot products used for that dot product.
    % (The MatMultNan routine will toss out pairs that have NaNs during the
    % dot product operations.)
    KeepThreshold = round(length(A(:,1)) / 2);
    
else
    TransposeFirst = false;
    [ C, Counts ] = MatMultNan(A,A');
    KeepThreshold = round(length(A(1,:)) / 2);
end
% Need to divide by the counts to form the covariance
Counts(Counts < KeepThreshold) = 0;
C = C ./ Counts;


if (ObsInCols == 1)
    % rows are samples and columns are data --> the sample covariance
    % matrix is based on A'*A, the data covariance matrix is based on A*A'.
    %
    % NormCols is a function that emulates the newer MATLAB function normc.
    if (TransposeFirst)
        % C is the sample covariance matrix
        % Using Ncols-2 since eigs wants the number of eigenvalues you ask
        % for to be less than the size of the matrix.
        [ EOF, Evals ] = eigs(C,Ncols-2);
        [ RawPC, Counts ] = MatMultNan(A,EOF);
        PC = NormCols(RawPC);
    else
        % C is the data covariance matrix
        [ PC, Evals ] = eigs(C,Nrows-2);
        [ RawEOF, Counts ] = MatMultNan(A',PC);
        EOF = NormCols(RawEOF);
    end
else
    % rows are data and columns are samples --> the sample covariance is
    % based on A*A' and the data covariance is based on A'*A.
    if (TransposeFirst)
        % C is the data covariance matrix
        [ PC, Evals ] = eigs(C,Ncols-2);
        [ RawEOF, Counts ] = MatMultNan(A,PC);
        EOF = NormCols(RawEOF);
    else
        % C is the sample covariance matrix
        [ EOF, Evals ] = eigs(C,Nrows-2);
        [ RawPC, Counts ] = MatMultNan(A',EOF);
        PC = NormCols(RawPC);
    end
end

end

