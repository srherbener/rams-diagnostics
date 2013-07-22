function [ ProdMat, Counts ] = MatMultNan( A, B)
%MatMultNan Matrix multiply with missing data (NaNs)
%   This function will perform a matrix multiply of A * B while handling
%   possible missing data (NaNs) in either of A or B.
%
%   In the output, ProdMat, the NaNs will be replaced with zeros. The
%   Counts array will show, for each element in ProdMat, how many
%   element-wise pairs of the two vectors used to compute that entry were
%   used. A count of zero would mean that no element-wise pairs were found
%   where both values in the pair were not NaNs. Count of 1 means 1 non-NaN
%   pair was found, count of 2 - found 2 non-NaN pairs, etc. It is
%   intended that the Count array can be used by the caller to sort out
%   where "good" data in ProdMat exist. It can also be used to create
%   covariances where you would do an element-wise division of Counts (or
%   Counts - 1).

[ NrA, NcA ] = size(A);
[ NrB, NcB ] = size(B); 

if (NcA == NrB)
    % The length of the rows of A (which is the number of columns in A)
    % match the length of the columns of B (number of rows in B),
    % so we are good to go. The resulting matrix will have dimensions:
    % NrA X NcB.
    %
    % Thanks to Rob Seigel for this algorithm.
    %  1) Do a matrix multiply of matrices that are shaped the same as A
    %  and B, but contain 1's where non-NaN data exist and 0's where NaNs
    %  exist. Doing this will sum up how many element-wise pairs where both
    %  values are non-NaN exist for each element in the output matrix (ie,
    %  all of the combinations of dot products between rows of A and
    %  columns of B). This yeilds the Counts output matrix.
    %
    %  2) Replace all of the NaNs in A and B with zeros and then perform
    %  the matrix multiply. This will exclude the element-wise pairs that
    %  had NaN values from the dot products due to the zeros that are now
    %  in those places.
    Counts = double(~isnan(A)) * double(~isnan(B));
    TempA = A;
    TempA(isnan(A)) = 0;
    TempB = B;
    TempB(isnan(B)) = 0;
    ProdMat = TempA * TempB;
else
    fprintf('ERROR: MatMultNan: the inner dimensions of A * B do not match\n');
    ProdMat = nan(1,1);
end


end

