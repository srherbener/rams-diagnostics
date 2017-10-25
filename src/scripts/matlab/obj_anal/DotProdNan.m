function [ Dprod, NumUsed ] = DotProdNan( X, Y, Weights, UseWeights )
%DotProdNan performs dot product on vectors X and Y, throws out NaNs
%   This function will perform a dot product on X and Y.
%   NaNs will be thrown out of the calculation.
%   If UseWeights is non zero, then the dot product will be weighted by the
%   values in the vector Weigths. Note Weights will not be checked for NaNs
%   so the caller is responsible for making Weights be free of NaNs.


if (length(X) == length(Y))

    Use = ~isnan(X) & ~isnan(Y); % Pick out element-wise pairs of X,Y data
                                 % where both elements are not NaN
    NumUsed = sum(Use);
    if (NumUsed > 0)
        if (UseWeights ~= 0)
            if (length(X) == length(Weights))
                Dprod = sum(X(Use) .* Y(Use) .* Weights(Use));
            else
                fprintf('ERROR: DotProdNan: the length of the Weights vector must be equal to the input vectors X and Y\n');
                Dprod = NaN;
                NumUsed = 0;
            end  
        else
            Dprod = sum(X(Use) .* Y(Use));
        end
    else
        % There are no element-wise pairs in X,Y that have both elements
        % not equal to NaN. Set the result to zero. The return value
        % NumUsed will allow the caller to decide what to do with the
        % result when there were no element-wise pairs selected.
        Dprod = 0;
    end
else
    fprintf('ERROR: DotProdNan: the lengths of the input vectors X and Y must be equal\n');
    Dprod = NaN;
    NumUsed = 0;
end

end

