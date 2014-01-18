function [ Hcombined ] = CombineBins( Harray, Hdim, GroupSize )
% CombineBins reduces the length of a histogram by combining adjacent bins
%
% Harray is an array that contains histograms.
%
% Hdim is an integer denoting the dimension (in an array) of which 
% contains the histograms.
%
% GroupSize is an integer denoting how many adjacent bins become a single
% bin in the output histogram.
%

  % Algorithm:
  %   Extract vectors that describe indicies of Harray going in linear order (ind2sub)
  %   Along the dimension that contains the histograms, rewrite the corresponding vector
  %     so that it can be given to accumarray() for the purpose of adding up adjacent
  %     elements along that dimension.
  %   Form the output array using accumarray().

  Hsize = size(Harray);
  Npts = prod(Hsize);
  Nndims = ndims(Harray);

  % Don't allow a GroupSize that is bigger than the selected dimension size.
  if (GroupSize > Hsize(Hdim))
    fprintf('ERROR: CombineBins: GroupSize (%d) must be <= size (%d) of dimension selected by Hdim (%d)\n', GroupSize, Hsize(Hdim), Hdim);
    Hcombined = nan;
    return;
  end

  % For now, support up to 4 dimensions.
  % For more dimensions, add vectors to the ind2sub and Indices = [  ] commands below.
  if (ndims(Harray) > 4)
    fprintf('ERROR: CombineBins: only support up to 4 dimensional histogram arrays.\n');
    Hcombined = nan;
    return;
  end
  [ i j k l ] = ind2sub(Hsize, 1:Npts);
  Indices = [ i; j; k; l ];

  % Change the vector (row) in Indices corresponding to the dimension given in Hdim 
  % so that it will collapse that dimension when accumarra()y is called.
  %
  % Temporarily reduce the dimension size to an even multiple of GroupSize. Then figure out
  % the pattern for selecting adjacent blocks (of size GroupSize) out of the dimension. If the original
  % dimension size didn't divide evenly by GroupSize, then tack on the extras needed to combine the remaining
  % entries into an extra output element.

  % Form the access pattern though the size of the selected dimension.
  % Note that if Hsize(Hdim) / GroupSize does not divide evenly, that we will need to tag on
  % an extra piece at the end of the vector.
  RepLength = floor(Hsize(Hdim) / GroupSize);
  VectRem = mod(Hsize(Hdim), GroupSize);
  NewVector = 1:RepLength;
  NewVector = repmat(NewVector, [ GroupSize 1 ]);
  NewVector = [ NewVector(:)' ones([ 1 VectRem ]).*(RepLength+1) ];

  % Replicate the pattern within the context of where its dimension fits into the total array.
  RowRep = prod([ 1 Hsize(1:Hdim-1) ]);
  ColRep = prod([ 1 Hsize(Hdim+1:end) ]);

  NewVector = repmat(NewVector, [ RowRep ColRep ]);
  NewVector = NewVector(:)';

  % Replace the row in Indices corresponding to Hdim with NewVector
  Indices(Hdim,:) = NewVector;

  % Create output
  Hcombined = accumarray(Indices', Harray(:));

end
