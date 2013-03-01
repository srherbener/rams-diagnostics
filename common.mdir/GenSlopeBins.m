function [ NR, NT, BL, BU ] = GenSlopeBins(Counts, Bins, BinGroupSize)
% GenSlopeBins generate bins for slope data from POP diagnostic
%
%   Counts - array holding the count datat from the HDF5 pop file
%   Bins - vector holding the bin edge vaules
%   BinGroupSize - integer that says how many adjacent bins are to
%                  be used to combine into the output bins
%
% This function will create bins out of the count data from the pop diagnostic.
% The argument Counts is an array organized as (x,y,z,t) where
%
%   x - bins
%   y - size = 2 --> (Nt,Nr)
%   z - dummy dimension
%   t - time
%
% The argument BinGroupSize is an integer that says how many adjacent bins
% from the input data are to be combined into a single bin for the output
% data. This can be used to reduce the number of bins if desired. Say there
% are 100 bins in the input, then if BinGroupSize is set to 5 there will be
% 20 bins in the output. The first output bin contains the sum of the counts
% from input bins 1 through 5; the second output bin contains the sum of the
% counts from input bins 6 through 10; etc.
%
% Output:
%
%  NR - total number of grid cells for each bin where it is raining
%  NT - total number of grid cells for each bin
%  BL - lower boundaries of each bin
%  BU - upper boundaries of each bin
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

% Create the bin edge data
TempBL = Bins(1:end-1);
TempBU = Bins(2:end);

% Split up the Counts array into the corresponding NR and NT arrays, and add
% in the counts from the last bin into the next to last bin.
TempNR = squeeze(Counts(1:end-1,2,:,:)); % extract all but the last bin
TempNT = squeeze(Counts(1:end-1,1,:,:));
TempNR(end,:) = TempNR(end,:) + squeeze(Counts(end,2,:,:))'; % add in the last bin to
TempNT(end,:) = TempNT(end,:) + squeeze(Counts(end,2,:,:))'; % then end of NR and NT

% Form the output by combining the bins. This depends on the new number of bins
% Nb being an integer.

% Form the indices to pick out the proper data spaced at the interval defined
% by BinGroupSize.
OrigNb = length(TempBL);
Use = 1:BinGroupSize:OrigNb;

Offset = BinGroupSize - 1;
BL = TempBL(Use);
BU = TempBU(Use+Offset);

NR = TempNR(Use,:);
NT = TempNT(Use,:);
for i = 2:BinGroupSize
  Offset = i - 1;
  NR = NR + TempNR(Use+Offset,:);
  NT = NT + TempNT(Use+Offset,:);
end

end
