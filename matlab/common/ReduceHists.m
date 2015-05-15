function [ Hreduced ] = ReduceHists( Hdata, Hdim, Bins, Method, Param )
% ReduceHists reduces histogram bins to a single number
%
% Method: 'wtmean' - weighted mean: estimates mean of original data
%
%         'max'    - max: estimates the mode of the histogram
%
%         'farea'  - "fractional" area of histogram counts. This is like the median
%                    measurement except you can specify (Param argument) what portion
%                    of the distribution is to the left and right of the value. Setting
%                    Param to 0.5 returns the median value (1/2 of the distribution is
%                    to the left and 1/2 to the right). Setting Param to 0.25 gives the
%                    value where 1/4 of the distribution is to the left, and 3/4 to the
%                    right.
%
%         'pct'    - Percentile of bin values corresponding to the non-zero counts,
%                    percentile is specified in Param argument. This is another
%                    way to do what farea is doing. This seems to produce "block-ier"
%                    (discretize looking) results.
%

% Bin values are the edges of the bins where count(i) are the
% values, x, that were: bin(i) <= x < bin(i+1). The count in
% final bin, number n, is for x == bin(n). For the calculations
% below, associate the counts with the average of the bin edges,
% and for count n, just use bin n value.
Nb = length(Bins);
BinVals = (Bins(1:end-1) + Bins(2:end)) .* 0.5;
BinVals(Nb) = Bins(Nb);

% Create an array the is the same size as Hdata and contains the
% bin values going down the dimension with the histogram counts.
% This array will be used by the wtmean and pct90 methods.
%
% The bin value array (BINS) can be formed from an array that has
% the same number of dimensions as Hdata with all dimension sizes
% equal to one except for the dimension with the histogram counts, which has
% its size set to the same size in Hdata. Note that this array will
% have the same total length as the Bins array and by copying in
% the Bins array into this array will automatically place one copy
% of the bin values running along the dimension that corresponds to
% that of the histogram counts in Hdata. Then repmat can be used to
% tile the bin values into the same size of Hdata.

Hsize = size(Hdata);

% dimension sizes for array with all dims sizes equal to one
% except for the dimension with the histograms
InitSize = ones([ 1 length(Hsize) ]);
InitSize(Hdim) = Hsize(Hdim);

% dimesions sized for tiling with repmat where all dim sizes match
% those of Hdata except for the dimension we placed the bin values
% (which has a repeat factor of 1).
RepFactor = Hsize;
RepFactor(Hdim) = 1;

BINS = ones(InitSize);
BINS(:) = BinVals;
BINS = repmat(BINS, RepFactor);

% BINS now is now the same size as Hdata with the bin values
% running down the same dimension as the histogram counts.
%

switch Method
    case 'pct'
        % Copy the BINS array into a new array. Set the locations where
        % the counts (Hdata) has zeros to nans. Then run the prctile()
        % command to get the 90th percentile of the bin value distributions.
        BINS_WITH_COUNTS = BINS;
        BINS_WITH_COUNTS(Hdata == 0) = nan;
        Hreduced = squeeze(prctile(BINS_WITH_COUNTS, [ Param ], Hdim));

    case 'wtmean'
        % Weighted mean - this estmates the mean of the original data
        
        % Want to form the average of the bin values weighted by their
        % corresponding counts in the histogram. Ie,
        %
        %    sum(bin(i) * count(i)) / sum(count(i)
        %
        % To do this efficiently, need a array with the same dimensions as
        % Hdata (histograms) with the bin values going along the dimension
        % containing the counts (Hdim). This array can be used then to do
        % an element by element multiply to form the numerator sum. (Doing
        % the element by element multiply will be way faster than writing a
        % loop to walk through each element of Hdata.)
        %
        % It's possible to get all counts equal to zero since we can filter
        % out input data before making the counts. In this case both SUM_BC
        % and SUM_C will be zero and we will proceed to do a 0/0 operation.
        % This results in a nan which is what we want so we can exclude this
        % case from downstream analysis.
        SUM_BC = squeeze(sum(Hdata .* BINS, Hdim));
        SUM_C = squeeze(sum(Hdata, Hdim));
        Hreduced = SUM_BC ./ SUM_C;

    case 'max'
        % Max value
        %
        % Find the indices of the max values along the bins dimension.
        % Then use that result to create an array that contains the bin
        % value corresponding to the index where the max count
        % existed.
        [ MAXC_V, MAXC_I ] = max(Hdata, [], Hdim);
        Hreduced = squeeze(BinVals(MAXC_I));

    case 'farea'
        % fraction of area
        %
        % Define the area under the histogram curve assuming that each bin
        % has a unit width. Then the total area is the sum of the histogram
        % counts.
        %
        % Form the fractional area that we are looking for as:
        %   Param .* (sum of histogram counts)
        %
        % Form a cumulative sum of the bin counts and find the two bins
        % where the cumulative sum values bracket the fractional area
        % value. Then do a linear interpolation of the corresponding bin values
        % for the output values.

        % Form the cumulative sum array, plus an array matching the size of
        % Hdata with the fractional area value repeated along the
        % histogram dimension.
        CSUM = cumsum(Hdata, Hdim);
        PSUM = Param .* sum(Hdata, Hdim);
        DIFF = CSUM - repmat(PSUM, InitSize);

        % After subtracting PSUM from CSUM, the final negative entry along
        % the histogram dimension will be the left side of the interval where
        % PSUM lies, and the first postivie entry is the right side of the
        % interval where PSUM lies.
        LEFT = double(DIFF < 0);   % convert logical types to double
        RIGHT = double(DIFF >= 0);  % for subsequent operations

        % RIGHT has 1's where DIFF was >= 0. Do a cumulative sum and zero
        % out all the entries greater than 1 which will leave the first
        % positive entry with a 1, all others zero.
        RIGHT = cumsum(RIGHT,Hdim);
        RIGHT(RIGHT > 1) = 0;

        % LEFT side of the interval is a little trickier since cumsum
        % goes in the wrong direction. A cumsum going from right to left
        % can be accomplished by subtracting the cumsum from the sum+1 and
        % zeroing out the result where LEFT was also zero (ie, multiply by LEFT).
        LEFT_CS = cumsum(LEFT,Hdim);
        LEFT_S  = repmat(sum(LEFT,Hdim)+1, InitSize);
        LEFT = (LEFT_S - LEFT_CS) .* LEFT;
        LEFT(LEFT > 1) = 0;

        % Now RIGHT has a 1 on the right side of the interval where
        % Param*Sum lies, and zero's elsewhere. Ditto for LEFT with
        % respect to the left side of the interval.
        B1  = squeeze(sum(LEFT  .* BINS,Hdim));
        B2  = squeeze(sum(RIGHT .* BINS,Hdim));
        CS1 = squeeze(sum(LEFT  .* CSUM,Hdim));
        CS2 = squeeze(sum(RIGHT .* CSUM,Hdim));

        % Return the linear interpolation
        Hreduced = B1 + ((B2-B1) .* ((squeeze(PSUM)-CS1) ./ (CS2-CS1)));

        % It is possible for PSUM_HD1 to fall to the left of the first
        % entry in CSUM or to the right of the last entry of CSUM. If this
        % happens either LEFT or RIGHT will be all zeros. Check for these
        % cases and set the output to NAN to denote that a solution wasn't
        % attainable.
        SUM_LEFT  = squeeze(sum(LEFT,Hdim));
        SUM_RIGHT = squeeze(sum(RIGHT,Hdim));
        Hreduced(SUM_LEFT == 0) = nan;
        Hreduced(SUM_RIGHT == 0) = nan;

    otherwise
        fprintf('WARNING: ReduceHists: Unrecognized reduction method: %s\n', Method);
        fprintf('WARNING:              Setting output to all nans\n');
        Hreduced = nan(size(Hdata));
end


end
