function [ ] = PrintSigPeaks( CohSqr, Csig, Freqs, Phase, Ptitle )
%PrintSigPeaks Print out info about the significant peaks of a cross spectrum
%   This function will take the Coh^2 spectrum, compare against the
%   significance level given in Csig, and print out the frequency and phase
%   angle for all frequences that rise above the significance level.


fprintf('%s\n', Ptitle);
fprintf('  %15s %15s\n', 'Frequency', 'Phase Shift');
SigFreqs = Freqs(CohSqr > Csig);
SigPhase = Phase(CohSqr > Csig);
for i = 1:length(SigFreqs)
    fprintf('  %15.3f %15.3f\n', SigFreqs(i), SigPhase(i));
end
fprintf('\n');

end

