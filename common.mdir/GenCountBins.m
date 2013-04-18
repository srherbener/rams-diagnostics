function [ OUT_COUNTS, XL, XU, YL, YU ] = GenCountBins(Counts, Xbins, Ybins, Xmin, Xmax, Xsize, Ymin, Ymax, Ysize)
% GenCountBins combine bins from counts (histogram fashion)
%
%   Counts - array holding the count data
%   Xbins - vector holding the bin edge vaules for the x dimension
%   Ybins - vector holding the bin edge vaules for the y dimension
%   Xsize, Ysize - integers that say how many adjacent bins (in the x- and
%                  y-directions that are to be used to combine into the output bins
%
% This function will create bins out of the count data.
% The argument Counts is an array organized as (x,y,z,t) where
%
%   x - bins
%   y - bins
%   z - levels
%   t - time
%
% The arguments Xsize, Ysize are integers that say how many adjacent bins
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
%  YL - lower boundaries of each Ybin
%  YU - upper boundaries of each Ybin
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

[ Nx, Ny, Nz, Nt ] = size(Counts);

% Create the bin edge data
TempXL = Xbins(1:end-1);
TempXU = Xbins(2:end);
TempYL = Ybins(1:end-1);
TempYU = Ybins(2:end);

% In the Counts array, add in the counts from the last bin into the next to
% last bin (in both the x and y directions). This means that we need to add
% the "band" around the ends of array (where the -1's are located in the
% example below) to the entries adjacent to that band.
%
%  Counts = 
%    1  4  7 -1         1 4 6
%    2  5  8 -1   --->  2 5 7
%    3  6  9 -1         2 5 6
%   -1 -1 -1 -1
%
TempC = squeeze(Counts(1:end-1,1:end-1,:,:)); % extract all but the last bin

TempC(:,end,:) = TempC(:,end,:) + reshape(Counts(1:end-1,end,:,:), [ Nx-1 1 Nt ] );
TempC(end,:,:) = TempC(end,:,:) + reshape(Counts(end,1:end-1,:,:), [ 1 Ny-1 Nt ] );
TempC(end,end,:) = TempC(end,end,:) + reshape(Counts(end,end,:,:), [ 1 1 Nt ]);

% Form the output by combining the bins. This depends on the new number of bins
% Nb being an integer. NR and NT are now organized as (x,y,t).
%
X1 = find(TempXL >= Xmin, 1, 'first');
X2 = find(TempXU <= Xmax, 1, 'last');
Y1 = find(TempYL >= Ymin, 1, 'first');
Y2 = find(TempYU <= Ymax, 1, 'last');

% Form the indices to pick out the proper data spaced at the interval defined
% by (Xsize,Ysize).
UseX = X1:Xsize:X2;
UseY = Y1:Ysize:Y2;
NewNx = length(UseX);
NewNy = length(UseY);

OffsetX = Xsize - 1;
OffsetY = Ysize - 1;
XL = TempXL(UseX);
XU = TempXU(UseX+OffsetX);
YL = TempYL(UseY);
YU = TempYU(UseY+OffsetY);

OUT_COUNTS = zeros(NewNx,NewNy,Nt);
for i = 1:Xsize
  OffsetX = i - 1;
  for j = 1:Ysize
    OffsetY = j - 1;
    OUT_COUNTS = OUT_COUNTS + TempC(UseX+OffsetX,UseY+OffsetY,:);
  end
end

end
