function [ Hreduced ] = ReduceHists( Hdata, Hdim, Bins, Method )
% ReduceHists reduces histogram bins to a single number
%
% Method: 'wtmean' - weighted mean
%         'max'    - max
%         'com'    - center of mass
%

switch Method
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
        % The bin value array (BINS) can be formed from an array that has
        % the same dimensions as Hdata with all dimension sizes equal to
        % one except for the dimension with the histogram counts, which has
        % its size set to the same size in Hdata. Note that this array will
        % have the same total length as the Bins array and by copying in
        % the Bins array into this array will automatically place one copy
        % of the bin values running along the dimension that corresponds to
        % that of the histogram counts in Hdata. The repmat can be used to
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
        % Create an array that contains ones where there are non zero
        % counts in the corresponding histogram, and nans where all bins
        % in the corresponding historgram contain zero counts. This array
        % can then be multiplied by the the result of the max calculations
        % in order to place nans where ever the histograms of all
        % zero counts exist.
        HAVE_COUNTS = squeeze(sum(Hdata, Hdim));
        HAVE_COUNTS(HAVE_COUNTS ~= 0) = 1.0;
        HAVE_COUNTS(HAVE_COUNTS == 0) = nan; 
 
        % Find the indices of the max values along the bins dimension (number 2)
        % Then use that result to create an array organized as (r,z,t) that
        % contains the bin value corresponding to the index where the max count
        % existed.
        [ MAXC_V, MAXC_I ] = max(Hdata, [], Hdim);
        Hreduced = squeeze(Bins(MAXC_I));
        if (size(HAVE_COUNTS) == size(Hreduced'))
          % catch the case where you have a column vector and a row vector
          Hreduced = Hreduced' .* HAVE_COUNTS;
        else
          Hreduced = Hreduced .* HAVE_COUNTS;
        end
    case 'com'
        % Center of mass
        %
        % Create an array that contains ones where there are non zero
        % counts in the corresponding histogram, and nans where all bins
        % in the corresponding historgram contain zero counts. This array
        % can then be multiplied by the the result of the max calculations
        % in order to place nans where ever the histograms of all
        % zero counts exist.
        HAVE_COUNTS = squeeze(sum(Hdata, Hdim));
        HAVE_COUNTS(HAVE_COUNTS ~= 0) = 1.0;
        HAVE_COUNTS(HAVE_COUNTS == 0) = nan;
        
        % For every histogram (along dimension Hdim) find the index nearest where
        % the area under the curve is half the total area under the curve. Then
        % use these indices to place the corresponding bin values into the result.
        % This will track where the "center of mass" of the histogram distribution
        % lies.
        %
        % Accomplish this by subtracting off 1/2 the total sum from the cumulative
        % sum, taking the absolute value of the result, and locating the indices
        % where the minimum occurs (closest differece to zero). Then use these
        % indices to build a result that contains the corresponding bin
        % values.
        CSUM = cumsum(Hdata, Hdim);
        SUM = sum(Hdata, Hdim);
        
        % Expand SUM to the same shape as CSUM with copies of the sum
        % along the histogram dimension.
        Hsize = size(Hdata);
        RepFactor = ones([ 1 length(Hsize) ]);
        RepFactor(Hdim) = Hsize(Hdim);
        SUM = repmat(SUM, RepFactor);
        
        CSHALF = abs(CSUM - (SUM/2));
        [ HALF_V, HALF_I ] = min(CSHALF, [], Hdim);
        Hreduced = squeeze(Bins(HALF_I));
        if (size(HAVE_COUNTS) == size(Hreduced'))
          % catch the case where you have a column vector and a row vector
          Hreduced = Hreduced' .* HAVE_COUNTS;
        else
          Hreduced = Hreduced .* HAVE_COUNTS;
        end
        
    otherwise
        fprintf('WARNING: ReduceHists: Unrecognized reduction method: %s\n', Method);
        fprintf('WARNING:              Setting output to all nans\n');
        Hreduced = nan(size(Hdata));
end


end
