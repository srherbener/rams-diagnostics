function [ Anormc ] = NormCols( A )
%NormCols Normalize columns of an array.
%   This function will normalize (divide by length) the columns of
%   matrix A. This function is intended to be used in older versions of
%   MATLAB that do not have the normc built-in function.

% sum(A) will sum up the columns of matrix A and places these in a row
% vector, which is ideal for this task. Create a row vector with the
% Euclidean lengths of the columns, then use repmat to expand this into a
% matrix that matches the dimensions of A. Then use element-wise divide to
% finish the normalization task.

% row vector with one value per column of A, where that value is the sqaure
% root of the sum of the column elments.
LenCols = sqrt(sum(A.*A));

% repmat will take the matrix (LenCols is a 1 x NumberOfColumsInA matrix)
% and "tile" it to create another matrix. The second argument is the number
% of rows to create in the tiling, and the third argument is the number of
% columns in the tiling. size(A,1) returns the number of rows in A, so the
% following will replicate the LenCol vector and stack the copies
% vertically to get a matrix that has the corresponding column length value
% for the input matrix A.
LenMat = repmat(LenCols,size(A,1),1);

Anormc = A ./ LenMat;

end

