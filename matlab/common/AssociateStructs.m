%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AssociateStructs
%
% This function will take names specified in one struct array (Base) and find
% these names in another struct array (List), then record the matching indeces
% in the output structure OutStruct. Base is first copied to OutStruct so that
% the caller can replace the struct passed in by List. This gets around the fact
% that MATLAB uses call-by-value making it so that you cannot directly change
% an argument to this function.
%
% Type is used to enable Base to have multiple associations in it. The
% naming scheme for the structur elements needs to adhere to the following
% convection for this routine to work.
%
% 1) Base is an array of structures with an element named '<Type>name'
%     eg, if Type is 'PS', then Base need an element named 'PSname'
% 2) The value stored in Base(i).<Type>name will be searched for in List
%    where List is an array of structures with and element named 'Name'
%     eg, List(:).Name will be checked to see if it matches Base(i).PSname
% 3) The index in List where the match occurred will be entered into an element
%    in OutStruct (copy of Base) named '<Type>num'. A zero is assigned if
%    no match was made.
%     eg, if List(3).Name matched Base(2).PSname, then OutStruct(2).PSnum will
%     be set to 3.
%
% A copy of the input array with the new field (containing the index of the match)
% is returned to facilitate the replacement of Base in the calling routine.
% Since MATLAB uses pass-by-value it does not work to do the assignment
% here in this routine.
% 
function [ OutStruct ] = AssociateStructs(Base, List, Type, Ltype, Btype)

  % copy Base so the caller can replace Base
  OutStruct = Base;

  for i = 1:length(Base)
    BaseName = sprintf('%sname', Type);
    BaseNum  = sprintf('%snum', Type);

    TestName = Base(i).(BaseName);
    TestArray = { List(:).Name };
    
    Index = find(strcmp(TestName, TestArray));
    if (isempty(Index))
      fprintf('WARNING: Could not find a match for %s "%s" specified in %s number %d\n', Ltype, TestName, Btype, i);
      Index = 0;
    end

    OutStruct(i).(BaseNum) = Index;
  end
end
