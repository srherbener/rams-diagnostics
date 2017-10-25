function [ Map ] = Create2dMap( X, Nrows, Ncols )
%Create2dMap Takes 2D data stored in a 1D vector and outputs as 2D array
%   This routine takes the 2D data that has been stored in a 1D array (for
%   the purpose of something like EOF analysis) and loads it into a 2D
%   array (for the purpose of something like plotting).
%
%   The dimensions of the 2D data are paseed in through Nrows (number of
%   rows) and Ncols (number of columns). The organization of the data in X
%   is expected to be:
%      First Ncols entries: 1 through Ncols columns for row number 1
%      Next Ncols entries: 1 through Ncols columns for row number 2
%      etc.

% reshape can be used for this. Reshape wants to take elements column-wise
% from its input and place elements column-wise first into its output. This
% is opposite of the organization we have in X. What we want is row-wise
% first action. This can be done by telling reshape to build the transpose
% of the result we want, and then transpose that result before outputting
% it.

Map = reshape(X,Ncols,Nrows)'; % Note the transpose operator on reshape()


end

