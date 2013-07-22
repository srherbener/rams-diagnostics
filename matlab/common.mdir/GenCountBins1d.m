function [ OUT_COUNTS, XL, XU ] = GenCountBins1d(Counts, Xbins, Xmin, Xmax, Xsize)
% GenCountBins1d combine bins from counts (histogram fashion)
%
%   Counts - array holding the count data
%   Xbins - vector holding the bin edge vaules for the x dimension
%   Xsize - integers that say how many adjacent bins in the x-direction
%           that are to be used to combine into the output bins
%
% This function will create bins out of the count data.
% The argument Counts is an array organized as (x,t) where
%
%   x - bins
%   t - time
%
% The arguments Xsize is an integer that says how many adjacent bins
% from the input data are to be combined into a single bin for the output
% data. This can be used to reduce the number of bins if desired. Say there
% are 100 Xbins in the input, then if Xsize is set to 5 there will be
% 20 Xbins in the output. The first output bin contains the sum of the counts
% from input bins 1 through 5; the second output bin contains the sum of the
% counts from input bins 6 through 10; etc.
%
% Output:
%
%  OUT_COUNTS - combined counts in the new bin structure
%  XL - lower boundaries of each Xbin
%  XU - upper boundaries of each Xbin
%

% The data in Counts will be divided up into the bins as follows:
%
%    For each bin, except the last one, the entry in Counts in bin location i
%    represents the count for grid cells selected from the range:
%
%         Bins(i) <= X < Bins(i+1)
%
%    The last bin, the entry represents the count for grid cells
%    selected by:
%
%            X == Bins(end)
%
% Note it's pretty unlikely that the last bin has any non-zero counts since
% X would have had to have been exactly equal to the bin edge value. Because
% of this, there are really effectively length(Bins)-1 bins. (we can just
% add the last bin to the prior bin and git rid of the final bin)

[ Nx, Nt ] = size(Counts);

% Create the bin edge data
TempXL = Xbins(1:end-1);
TempXU = Xbins(2:end);

% In the Counts array, add in the counts from the last bin into the next to
% last bin. This means that we need to add
% the "band" around the ends of array (where the -1's are located in the
% example below) to the entries adjacent to that band.
%
%  Counts = 
%    1  4  7         1 4 7
%    2  5  8   --->  2 5 8
%    3  6  9         2 5 8
%   -1 -1 -1
%
TempC = Counts(1:end-1,:); % extract all but the last bin
TempC(end,:,:) = TempC(end,:,:) + Counts(end,:); % Add the last bin in Counts to last bin in TempC

% Form the output by combining the bins. This depends on the new number of bins
% Nb being an integer. NR and NT are now organized as (x,t).
%
X1 = find(TempXL >= Xmin, 1, 'first');
X2 = find(TempXU <= Xmax, 1, 'last');

% Form the indices to pick out the proper data spaced at the interval defined
% by Xsize.
UseX = X1:Xsize:X2;
NewNx = length(UseX);

OffsetX = Xsize - 1;
XL = TempXL(UseX);
XU = TempXU(UseX+OffsetX);

OUT_COUNTS = zeros(NewNx,Nt);
for i = 1:Xsize
  OffsetX = i - 1;
  OUT_COUNTS = OUT_COUNTS + TempC(UseX+OffsetX,:);
end

end
