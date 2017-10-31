function [ RnSpect, Efold, R1 ] = GenRnoiseSpectrum( S, DeltaT, Freqs )
%GenRnoiseSpectrum Generate red noise spectrum
%   This function will generate the associated red noise spectrum for
%   series S and for the frequency values in Freqs

% Need the e-folding time associated with the data which is given by:
%
%   Efold = -DeltaT / ln(R1), where R1 is the lag-1 autocorrelation
%

N = length(S);

[ Lag1_a1, Lag1_yint, R1 ] = LagRegFit(S,S,1,1,N-1);
Efold = -DeltaT / log(R1);

[ RnSpect ] = RedNoiseSpectrum(Freqs, Efold);

end

