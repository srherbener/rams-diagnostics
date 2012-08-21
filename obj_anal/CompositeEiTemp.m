function [ CmeanPlus, CmeanMinus ] = CompositeEiTemp( Eindex, Temps )
%CompositeEiTemp calculate mean temperature for +/- 1.0 std dev ENSO index
%
%   This function looks at the ENSO index values in Eindex,
%   selects out data points where ENSO index is < -1.0 and data points where
%   ENSO index is > 1.0, and calculates and returns the mean Temperature for
%   those two bins.
%      CmeanPlus is the composite mean for ENSO index > 1.0
%      CmeanMinus is the composite mena for ENSO index < -1.0

if (length(Eindex) == length(Temps))
    ip = 0; % index for "plus" data (ENSO index > 1.0)
    im = 0; % index for "minus" data (ENSO index < -1.0)
    % Grab the data sets for both "plus" and "minus" ENSO index
    for i = 1:length(Eindex)
        if (Eindex(i) > 1.0)
            ip = ip + 1;
            DataPlus(ip) = Temps(i);
        else
            if (Eindex(i) < -1.0)
                im = im + 1;
                DataMinus(im) = Temps(i);
            end
        end
    end
    
    CmeanPlus = MeanNan(DataPlus);
    CmeanMinus = MeanNan(DataMinus);
else
    fprintf ('ERROR: CompositeEiTemp: Number of elements in Eindex and Temps must match\n');
    fprintf ('                        Number of elements, Eindex: %d\n', NumE);
    fprintf ('                        Number of elements, Temps: %d\n', NumT);
    CmeanPlus = NaN;
    CmeanMinus = NaN;
end


end

