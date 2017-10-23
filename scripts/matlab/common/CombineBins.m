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
  % so that it will collapse that dimension when accumarray() is called. At this point,
  % Indices contains the sequence of index values corresponding to walking through Harray
  % in column-major order (ie, matches the order given by Harray(:)). What we want to do
  % is change the indices in the row of Indices so that when Indices is given to accumarray()
  % adjacent elements along Hdim will be added together. Say Harray is 3x3 and Hdim is 2 and
  % group size is 2. Indices starts out containing:
  %
  %     1 2 3 1 2 3 1 2 3
  %     1 1 1 2 2 2 3 3 3
  %
  % We want Indices to be changed to:
  %
  %     1 2 3 1 2 3 1 2 3
  %     1 1 1 1 1 1 2 2 2
  %
  % which tells accumarray() to add together the first two elements along the second dimension,
  % (second row of Indices) to form the 1st element of the output array (in its second dimension)
  % and to add the 3rd elements of Harray, 2nd dimension to form the 2nd element in the output. If
  % Harray started out as:
  %
  %      1 4 7
  %      2 5 8
  %      3 6 9
  %
  % then using the above example, the output of accumarray() will be: 
  %
  %      5 7
  %      7 8
  %      9 9
  %
  % where the first two columns (1st two elements of the 2nd dimension) of Harray got added together.
  %
  % Temporarily reduce the dimension size to an even multiple of GroupSize. Then figure out
  % the pattern for selecting adjacent blocks (of size GroupSize) out of the dimension. If the original
  % dimension size didn't divide evenly by GroupSize, then tack on the extras needed to combine the remaining
  % entries into an extra output element.

  % Form the access pattern through the size of the selected dimension.
  % Note that if Hsize(Hdim) / GroupSize does not divide evenly, that we will need to tag on
  % an extra piece at the end of the vector.
  % Example:
  %   Let's say Hdim = 2 and the size of dimension 2 is 10.
  %
  %   If group size is 3, then we want NewVector to be: 1 1 1 2 2 2 3 3 3 4
  %   which says to accumarray() to add the first 3 elements together (1 1 1) to form the first
  %   element of the output array, the next 3 elements together (2 2 2) to form the second element
  %   of the output array, then next three together (3 3 3) to form the third element of the output,
  %   and add the last element (4) to form the fourth element of the output.
  %
  %   If group size is 2, then we want NewVector to be: 1 1 2 2 3 3 4 4 5 5
  RepLength = floor(Hsize(Hdim) / GroupSize);
  VectRem = mod(Hsize(Hdim), GroupSize);
  NewVector = 1:RepLength;
  NewVector = repmat(NewVector, [ GroupSize 1 ]);
  NewVector = [ NewVector(:)' ones([ 1 VectRem ]).*(RepLength+1) ];

  % Replicate the pattern within the context of where its dimension fits into the total array.
  % Note that it is important for NewVector to be a row vector during this section. This makes
  % the pattern come out correct when changing the output of repmat() to a 1D vector in the
  % NewVector(:) statement.
  RowRep = prod([ 1 Hsize(1:Hdim-1) ]);
  ColRep = prod([ 1 Hsize(Hdim+1:end) ]);

  NewVector = repmat(NewVector, [ RowRep ColRep ]);
  NewVector = NewVector(:)';

  % Replace the row in Indices corresponding to Hdim with NewVector
  Indices(Hdim,:) = NewVector;

  % Create output
  Hcombined = accumarray(Indices', Harray(:));

end
