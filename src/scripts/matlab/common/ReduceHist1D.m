function [ Hreduced ] = ReduceHist1D( Counts, Bins, Method, Param )
% ReduceHist1D reduces histogram to a single number
%
% Counts and Bins are vectors
%
% Method: 'wtmean' - weighted mean (results in same as taking
%                    average of data before turned into a histogram)
%         'max'    - max
%         'farea'  - fractional area treating bin values having unit
%                    width, fractional area is specified in Param
%                    argument (0 < Param < 1)
%         'pct'    - Percentile of bin values corresponding
%                    to the non-zero counts, percentile is
%                    specified in Param argument (0 < Param < 100)
%


switch Method
    case 'wtmean'
        % Weighted mean
        
        % Want to form the average of the bin values weighted by their
        % corresponding counts in the histogram. Ie,
        %
        %    sum(bin(i) * count(i)) / sum(count(i)
        %
        Hreduced = sum(Bins .* Counts) / sum(Counts);

    case 'max'
        % Max value
        %
        % Find the index of the max count value and return the corresponding
        % bin value.
        [ MAXC_V, MAXC_I ] = max(Counts);
        Hreduced = Bins(MAXC_I);

    case 'pct'
        % Copy the Bins array into a new array. Set the locations where
        % the counts are zero to nans. Then run the prctile() command to get
        % the percentile of the nonzero bin values.
        BinsWithCounts = Bins;
        BinsWithCounts(Counts == 0) = nan;
        Hreduced = squeeze(prctile(BinsWithCounts, Param));

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

        CumCounts = cumsum(Counts);
        Pval = Param * sum(Counts);

        I1 = find(CumCounts <  Pval, 1, 'last');
        I2 = find(CumCounts >= Pval, 1, 'first');

        if (~isempty(I1) && ~isempty(I2))
          B1 = Bins(I1);
          B2 = Bins(I2);
          CC1 = CumCounts(I1);
          CC2 = CumCounts(I2);
         
          % Return the linear interpolation
          Hreduced = B1 + ((B2-B1) .* ((Pval-CC1) ./ (CC2-CC1)));
        else
          Hreduced = nan;
        end

    otherwise
        fprintf('WARNING: ReduceHist1D: Unrecognized reduction method: %s\n', Method);
        fprintf('WARNING:               Setting output to NaN\n');
        Hreduced = nan;
end


end
