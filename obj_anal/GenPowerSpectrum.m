function [ Pspect ] = GenPowerSpectrum ( F )
%GenPowerSpectrum Generate the power spectrum from the Fourier coefficients
%   This function will take the complex Fourier coefficients in F and
%   create the associated power spectrum. F can be a row vector or a matrix
%   where each row in F is an individual set of Fourier coefficients. In
%   this case, the power spectra of each row is calculated and the ouptut
%   is set to the average of these power spectra.

% Multiply every element of F by it's complex conjugate. If F was a row
% vector, then you are done. If F was a matrix, then finish by taking the
% mean of the matrix (which will average the columns are result in a row
% vector).

% Want to end up with 1/2 * Ck^2 for the power numbers. Ck^2 is obtained
% from multiplying the elements of F by their complex conjugates.

Pspect = (F .* conj(F)) * 0.5; % Form the 1/2 * Ck^2 values

if (length(F(:,1)) > 1)
    % more than one row in F --> average the power numbers in the columns
    Pspect = mean(Pspect);
end

end

