function [ Fields ] = Line2Fields ( Line, Delim )
%Line2Fields separate the string into delimited fields
%
% This function will parse the string given in Line using the
% delimiter character given in Delim and separate Line into fields
% which get passed back to the caller in the array Fields.

Remain = Line;

i = 1;
while (~isempty(Remain))
    [ Fields{i}, Remain ] = strtok(Remain, Delim);
    i = i + 1;
end

end
