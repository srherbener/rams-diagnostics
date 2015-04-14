function [ Hreduced ] = ReduceHists( Hdata, Hdim, Bins, Method, Param )
% ReduceHists reduces histogram bins to a single number
%
% Method: 'wtmean' - weighted mean
%         'max'    - max
%         'com'    - center of mass
%         'pct'    - Percentile of bin values corresponding
%                    to the non-zero counts, percentile is
%                    specified in Param argument.
%


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
BINS(:) = Bins;
BINS = repmat(BINS, RepFactor);

% BINS now is now the same size as Hdata with the bin values
% running down the same dimension as the histogram counts.
%

switch Method
    case 'pct'
        % Copy the Bins array into a new array. Set the locations where
        % the counts (Hdata) has zeros to nans. Then run the prctile()
        % command to get the 90th percentile of the bin value distributions.
        BINS_WITH_COUNTS = BINS;
        BINS_WITH_COUNTS(Hdata == 0) = nan;
        Hreduced = squeeze(prctile(BINS_WITH_COUNTS, [ Param ], Hdim));

    case 'wtmean'
        % Weighted mean
        
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
        Hreduced = squeeze(Bins(MAXC_I));

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
        PSUM_HD1 = Param .* sum(Hdata, Hdim);
        PSUM = repmat(PSUM_HD1, InitSize);
        DIFF = CSUM - PSUM;

        % After subtracting PSUM from CSUM, the final negative entry along
        % the histogram dimension will be the left side of the interval where
        % PSUM lies, and the first postivie entry is the right side of the
        % interval where PSUM lies.
        NEG = DIFF < 0;
        POS = DIFF >= 0;

        % POS has 1's where DIFF was >= 0. Do a cumulative sum and zero
        % out all the entries greater than 1 which will leave the first
        % positive entry with a 1, all others zero.
        RIGHT = cumsum(POS,Hdim);
        RIGHT(RIGHT > 1) = 0;

        % LEFT side of the interval is a little trickier since cumsum
        % goes in the wrong direction. A cumsum going from right to left
        % can be accomplished by subtracting the cumsum from the sum+1 and
        % zeroing out the result where NEG was also zero (ie, multiply by NEG).
        NEG_CS = cumsum(NEG,Hdim);
        NEG_S  = repmat(sum(NEG,2)+1, InitSize);
        LEFT = (NEG_S-NEG_CS) .* NEG;
        LEFT(LEFT > 1) = 0;

        % Now RIGHT has a 1 on the right side of the interval where
        % Param*Sum lies, and zero's elsewhere. Ditto for LEFT with
        % respect to the left side of the interval.
        B1  = squeeze(sum(LEFT  .* BINS,Hdim));
        B2  = squeeze(sum(RIGHT .* BINS,Hdim));
        CS1 = squeeze(sum(LEFT  .* CSUM,Hdim));
        CS2 = squeeze(sum(RIGHT .* CSUM,Hdim));
       
        % Return the linear interpolation
        Hreduced = B1 + ((B2-B1) .* ((squeeze(PSUM_HD1)-CS1) ./ (CS2-CS1)));

    otherwise
        fprintf('WARNING: ReduceHists: Unrecognized reduction method: %s\n', Method);
        fprintf('WARNING:              Setting output to all nans\n');
        Hreduced = nan(size(Hdata));
end


end
