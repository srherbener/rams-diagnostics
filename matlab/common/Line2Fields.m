%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Line2Fields
%
% This function will parse the string given in Line using the
% delimiter character given in Delim and separate Line into fields
% which get passed back to the caller in the array Fields.
%
function [ Fields ] = Line2Fields ( Line, Delim )

  Remain = Line;
  Fields{1} = '';

  i = 1;
  while (~isempty(Remain))
    [ Fields{i}, Remain ] = strtok(Remain, Delim);
    i = i + 1;
  end

  % Fields will end up with blank string in its last entry as
  % a result of the above loop. If there is more than one element
  % in Fields, then strip off this last entry since it is an artifact
  if ((strcmp(Fields{end},'')) && (length(Fields) > 1))
    Fields = Fields(1:end-1);
  end
end
