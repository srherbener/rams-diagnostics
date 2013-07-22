function [ Rmean ] = RunningMean( X, Size )
%RunningMean Generate running mean
%   This function will generate a running mean of X using a window size
%   specified by Size. The endpoints of the series are handled by repeating
%   the end data as many times as necessary. This allows the returned
%   running mean to be the same size as the input X.

N = length(X);
HalfS = floor(Size/2);
Rmean = zeros(1,N);
for i = 1:N
    jstart = i - HalfS;
    jend = jstart + Size - 1;
    TempMean = 0;
    for j = jstart:jend
        if (j < 1)
            TempMean = TempMean + X(1);
        elseif (j > N)
            TempMean = TempMean + X(N);
        else
            TempMean = TempMean + X(j);
        end
    end
    Rmean(i) = TempMean / Size;
end


end

