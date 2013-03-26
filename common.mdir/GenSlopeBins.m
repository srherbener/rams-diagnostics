function [ NR, NT, XL, XU, YL, YU ] = GenSlopeBins(Counts, Xbins, Ybins, Xsize, Ysize)
% GenSlopeBins generate bins for slope data from POP diagnostic
%
%   Counts - array holding the count datat from the HDF5 pop file
%   Xbins - vector holding the bin edge vaules for the x dimension
%   Ybins - vector holding the bin edge vaules for the y dimension
%   Xsize, Ysize - integers that say how many adjacent bins (in the x- and
%                  y-directions that are to be used to combine into the output bins
%
% This function will create bins out of the count data from the pop diagnostic.
% The argument Counts is an array organized as (x,y,z,t) where
%
%   x - bins
%   y - bins
%   z - (Nt,Nr)  (note size of z dimension is 2)
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
%  NR - total number of grid cells for each bin where it is raining
%  NT - total number of grid cells for each bin
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

[ Nx, Ny, Nz, Ntime ] = size(Counts);

% Create the bin edge data
TempXL = Xbins(1:end-1);
TempXU = Xbins(2:end);
TempYL = Ybins(1:end-1);
TempYU = Ybins(2:end);

% Split up the Counts array into the corresponding NR and NT arrays, and add
% in the counts from the last bin into the next to last bin (in both the x and y
% directions).
%
% Note that we need to add in a "band" around the ends of array (where the -1's
% are located in the example below).
%
%  Counts = 
%    1  4  7 -1         1 4 6
%    2  5  8 -1   --->  2 5 7
%    3  6  9 -1         2 5 6
%   -1 -1 -1 -1
%
% Note that the -1 in the lower right got added in resulting in the 
% '9' entry getting a -1 added three times.
%
TempNR = squeeze(Counts(1:end-1,1:end-1,2,:)); % extract all but the last bin
TempNT = squeeze(Counts(1:end-1,1:end-1,1,:));

TempNR(:,end,:) = TempNR(:,end,:) + reshape(Counts(1:end-1,end,2,:), [ Nx-1 1 Ntime ] );
TempNR(end,:,:) = TempNR(end,:,:) + reshape(Counts(end,1:end-1,2,:), [ 1 Ny-1 Ntime ] );
TempNR(end,end,:) = TempNR(end,end,:) + reshape(Counts(end,end,2,:), [ 1 1 Ntime ]);

TempNT(:,end,:) = TempNT(:,end,:) + reshape(Counts(1:end-1,end,1,:), [ Nx-1 1 Ntime ] );
TempNT(end,:,:) = TempNT(end,:,:) + reshape(Counts(end,1:end-1,1,:), [ 1 Ny-1 Ntime ] );
TempNT(end,end,:) = TempNT(end,end,:) + reshape(Counts(end,end,1,:), [ 1 1 Ntime ]);

% Form the output by combining the bins. This depends on the new number of bins
% Nb being an integer. NR and NT are now organized as (x,y,t).

% Form the indices to pick out the proper data spaced at the interval defined
% by (Xsize,Ysize).
OrigNx = length(TempXL);
OrigNy = length(TempYL);
UseX = 1:Xsize:OrigNx;
UseY = 1:Ysize:OrigNy;
NewNx = length(UseX);
NewNy = length(UseY);

OffsetX = Xsize - 1;
OffsetY = Ysize - 1;
XL = TempXL(UseX);
XU = TempXU(UseX+OffsetX);
YL = TempYL(UseY);
YU = TempYU(UseY+OffsetY);

NR = zeros(NewNx,NewNy,Ntime);
NT = zeros(NewNx,NewNy,Ntime);
for i = 1:Xsize
  OffsetX = i - 1;
  for j = 1:Ysize
    OffsetY = j - 1;
    NR = NR + TempNR(UseX+OffsetX,UseY+OffsetY,:);
    NT = NT + TempNT(UseX+OffsetX,UseY+OffsetY,:);
  end
end

end
