function [ Moments ] = GenMoments( Terms, NumPoints, T1, T2 )
% GenMoments generate moment values from gen_moments output
%
% Input
%    Terms - 4D array (Nterms, Norder, Nz, Nt) containing summation terms
%    NumPoints - 2D array (Nz, Nt) containing N values corresponding to each term
%    T1, T2 - time dimension indices denoting range of time steps for selecting
%             sum terms for the final momement calculations
%
% Output
%    Moments - 2D array (Nz, Norder)

  % First sum up data across the selected time indicies
  Mterms = squeeze(sum(Terms(:,:,:,T1:T2), 4));
  if (ndims(Mterms) == 2)
    % Mterms had size = [ 1 1 Nz Nt ] before the squeeze(sum(Terms... command above.
    % Squeeze eliminated the first two dimensions. Put them back in so the subsequent
    % code will be happy.
    Mterms = reshape(Mterms, [ 1 1 length(Mterms) ]);
  end
  Npts = squeeze(sum(NumPoints(:,T1:T2), 2));

  % Mterms is (Nterms, Norder, Nz)
  % Npts is (Nz)
  %
  % Output is same size as Mterms. Organize Moments in the same
  % fashion and Mterms.
  [ Nterms Norder Nz ] = size(Mterms);
  Moments = zeros([ Nterms Norder Nz ]);

  % Create a map between combinations of variable numbers and indices
  % in the columns of Mterms. This is done to facilitate looking up
  % terms out of the Mterms array.
  i1 = 0;
  i2 = 0;
  i3 = 0;
  i4 = 0;
  TermMap = containers.Map; % like a PERL hash
  for iv1 = 1:Norder
    i1 = i1 + 1;
    key = sprintf('V%d', iv1);
    TermMap(key) = i1;

    for iv2 = iv1+1:Norder
      i2 = i2 + 1;
      key = sprintf('V%dV%d', iv1, iv2);
      TermMap(key) = i2;

      for iv3 = iv2+1:Norder
        i3 = i3 + 1;
        key = sprintf('V%dV%dV%d', iv1, iv2, iv3);
        TermMap(key) = i3;

        for iv4 = iv3+1:Norder
          i4 = i4 + 1;
          key = sprintf('V%dV%dV%dV%d', iv1, iv2, iv3, iv4);
          TermMap(key) = i4;
        end
      end
    end
  end

  % Walk through the variables again, this time computing the moments.
  % Use the TermMap to help pick out the necessary terms for the moment
  % calculations.
  for iv1 = 1:Norder
    key = sprintf('V%d', iv1);
    Moments(TermMap(key), 1, :) = CalcM1(Mterms, Npts, iv1);
    for iv2 = iv1+1:Norder
      key = sprintf('V%dV%d', iv1, iv2);
      Moments(TermMap(key), 2, :) = CalcM2(Mterms, Npts, TermMap, iv1, iv2);
      for iv3 = iv2+1:Norder
        key = sprintf('V%dV%dV%d', iv1, iv2, iv3);
        Moments(TermMap(key), 3, :) = CalcM3(Mterms, Npts, TermMap, iv1, iv2, iv3);
        for iv4 = iv3+1:Norder
          key = sprintf('V%dV%dV%dV%d', iv1, iv2, iv3, iv4);
          Moments(TermMap(key), 4, :) = CalcM4(Mterms, Npts, TermMap, iv1, iv2, iv3, iv4);
        end
      end
    end
  end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CalcM1()
%
% This function calculates the mean of the variable indicated by V1.

function [ M1 ] = CalcM1( Terms, Npts, V1 )

  S1 = squeeze(Terms(V1, 1, :));
  M1 = S1 ./ Npts;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CalcM2()
%
% This function calculates the co-variance of the variables indicated by V1 and V2.

function [ M2 ] = CalcM2( Terms, Npts, TermMap, V1, V2)

  key = sprintf('V%dV%d', V1, V2);
  V12 = TermMap(key);

  S1 = squeeze(Terms(V1, 1, :));
  S2 = squeeze(Terms(V2, 1, :));
  S12 = squeeze(Terms(V12, 2, :));

  M2 = (S12 - (S1 .* S2 ./ Npts)) ./ Npts;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CalcM3()
%
% This function calculates the 3rd moment of the variables indicated by V1, V2 and V3.

function [ M3 ] = CalcM3( Terms, Npts, TermMap, V1, V2, V3)

  key = sprintf('V%dV%d', V1, V2);
  V12 = TermMap(key);

  key = sprintf('V%dV%d', V1, V3);
  V13 = TermMap(key);

  key = sprintf('V%dV%d', V2, V3);
  V23 = TermMap(key);

  key = sprintf('V%dV%dV%d', V1, V2, V3);
  V123 = TermMap(key);

  S1 = squeeze(Terms(V1, 1, :));
  S2 = squeeze(Terms(V2, 1, :));
  S3 = squeeze(Terms(V3, 1, :));
  S12 = squeeze(Terms(V12, 2, :));
  S13 = squeeze(Terms(V13, 2, :));
  S23 = squeeze(Terms(V23, 2, :));
  S123 = squeeze(Terms(V123, 3, :));

  M3 = (S123 - ((S1.*S23 + S2.*S13 + S3.*S12) ./ Npts) + ((2.*S1.* S2.*S3) ./ (Npts .^ 2))) ./ Npts;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CalcM4()
%
% This function calculates the 4th moment of the variables indicated by V1, V2, V3 and V4.

function [ M4 ] = CalcM4( Terms, Npts, TermMap, V1, V2, V3, V4)

  key = sprintf('V%dV%d', V1, V2);
  V12 = TermMap(key);

  key = sprintf('V%dV%d', V1, V3);
  V13 = TermMap(key);

  key = sprintf('V%dV%d', V1, V4);
  V14 = TermMap(key);

  key = sprintf('V%dV%d', V2, V3);
  V23 = TermMap(key);

  key = sprintf('V%dV%d', V2, V4);
  V24 = TermMap(key);

  key = sprintf('V%dV%d', V3, V4);
  V34 = TermMap(key);

  key = sprintf('V%dV%dV%d', V1, V2, V3);
  V123 = TermMap(key);

  key = sprintf('V%dV%dV%d', V1, V2, V4);
  V124 = TermMap(key);

  key = sprintf('V%dV%dV%d', V1, V3, V4);
  V134 = TermMap(key);

  key = sprintf('V%dV%dV%d', V2, V3, V4);
  V234 = TermMap(key);

  key = sprintf('V%dV%dV%dV%d', V1, V2, V3, V4);
  V1234 = TermMap(key);

  S1 = squeeze(Terms(V1, 1, :));
  S2 = squeeze(Terms(V2, 1, :));
  S3 = squeeze(Terms(V3, 1, :));
  S4 = squeeze(Terms(V4, 1, :));

  S12 = squeeze(Terms(V12, 2, :));
  S13 = squeeze(Terms(V13, 2, :));
  S14 = squeeze(Terms(V14, 2, :));
  S23 = squeeze(Terms(V23, 2, :));
  S24 = squeeze(Terms(V24, 2, :));
  S34 = squeeze(Terms(V34, 2, :));

  S123 = squeeze(Terms(V123, 3, :));
  S124 = squeeze(Terms(V124, 3, :));
  S134 = squeeze(Terms(V134, 3, :));
  S234 = squeeze(Terms(V234, 3, :));

  S1234 = squeeze(Terms(V1234, 4, :));

  M4 = (S1234 ...
       - ((S1.*S234 + S2.*S134 + S3.*S124 + S4.*S123) ./ Npts) ...
       + ((S1.*S2.*S34 + S1.*S3.*S24 + S1.*S4.*S23 + S2.*S3.*S14 + S2.*S4.*S13 + S3.*S4.*S12) ./ (Npts .^ 2)) ...
       - ((3.*S1.* S2.*S3.*S4) ./ (Npts .^ 3)) ...
       ) ./ Npts;

end
