function [ NoMeans ] = RemoveMeanNan( A, ByColumn )
%RemoveMeanNan remove mean values either by column or by row
%   This function will subtract out mean values from entries in a matrix.
%   It will handle NaNs in the data by tossing them out when computing
%   means, but will leave the original entries in A alone.
%
%   If ByColumn is a 0 then row means are computed and subtracted from
%   every entry in that row. If ByColumn is a 1, then columns means are
%   computed and subtracted from every entry in that column.

[ Nrows, Ncols ] = size(A);
NoMeans = zeros(Nrows, Ncols);

% Walk through A either column by column or row by row subtracting out the
% corresponding mean as you go. i indexes rows, j indexes columns
if (ByColumn == 1)
    % Column means
    for j = 1:Ncols
        Mean = MeanNan(A(:,j));
        for i = 1:Nrows
            NoMeans(i,j) = A(i,j) - Mean;
        end
    end
else
    % Remove row means
    for i = 1:Nrows
        Mean = MeanNan(A(i,:));
        for j = 1:Ncols
            NoMeans(i,j) = A(i,j) - Mean;
        end
    end
end


end

