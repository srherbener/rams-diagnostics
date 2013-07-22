function [ CohSqr, CSmag, CSphase ] = GenCoherence ( Fx, Fy, CO, Q )
%GenCoherence Generate the coherence of a cross spectrum
%   This function will take the complex Fourier coefficients in Fx, Fy, the
%   real coefficients in CO and Q and generate the coherence. The inputs
%   are either row vectors or matrices containing multiple spectra; one
%   spectrum per each row in the matrices.

% One way or another averaging must be done to get a usable coherence
% value. This is so because the coherence at each frequency without
% averaging always comes out to be equal to one. Therefore, some form of
% averaging must be done before calculating the coherence numbers.
%
% If there is just one realization of the spectrum (inputs are row vectors)
% then average neighboring frequency values.
%
% If there are more than one realization of the spectrum (inputs are
% arrays) then average each frequency over the multiple realizations at
% just that particular frequency.

% Multiply every element of F by it's complex conjugate. If F was a row
% vector, then you are done. If F was a matrix, then finish by taking the
% mean of the matrix (which will average the columns are result in a row
% vector).

if (length(Fx(:,1)) == 1)
    % Have just one realization of the spectrum --> average neighboring
    % frequencies. Avaerage two neighboring frequencies together (according
    % to Hartmann's notes) where freq(i) and freq(i+1) get averaged. Doing
    % it this way will drop the highest frequecy instead of the lowest and
    % typically the lowest frequency is of more interest.
    
    CxSqr = Fx .* conj(Fx);
    CySqr = Fy .* conj(Fy);
    Fxy = Fx .* conj(Fy);
    
    N = length(Fx);
    Nm1 = N - 1;
    MeanFxy = Fxy(1:Nm1) + Fxy(2:N);
    Px = CxSqr(1:Nm1) + CxSqr(2:N);
    Py = CySqr(1:Nm1) + CySqr(2:N);
    
    % Coherence
    CohSqr = (MeanFxy .* conj(MeanFxy)) ./ (Px .* Py);
    
    % magnitude and phase of the averaged cross spectrum
    CSmag = abs(MeanFxy);
    CSphase = angle(MeanFxy);
else
    % Have multiple realizations of the spectrum --> average over the
    % multiple realizations.
    
    % Calculate Cx^2 and Cy^2 from Fx, Fy. Need to form a mean of squares
    % for the coherence formula. Ie, square the values first, then take
    % their mean. After multiplying by the conjugate values, the columns
    % will contain all of the Ck values for a given frequency. Use the
    % built-it function mean since it does column means to complete the
    % mean of squares calculation.
    CxSqr = Fx .* conj(Fx);
    CySqr = Fy .* conj(Fy);
    CxSqr = mean(CxSqr);
    CySqr = mean(CySqr);
    
    % The CO^2 and Q^2 values however, need to be calculated as the square
    % of means (ie, calculate the means first, then square the means). Keep
    % the mean values so that the magnitude and phase of the cross spectrum
    % can be calculated.
    COmean = mean(CO);
    Qmean = mean(Q);
    
    COsqr = COmean.^2;
    Qsqr = Qmean.^2;
    
    % Coherence squared is now the sum of CO^2 and Q^2 divided by the
    % product of Cx^2 and Cy^2.
    CohSqr = (COsqr + Qsqr) ./ (CxSqr .* CySqr);

    % Use COmean and Qmean to get the magnitude and phase of the cross
    % spectrum
    CSmag = sqrt(COsqr + Qsqr);
    CSphase = atan2(Qmean, COmean);
end


end

