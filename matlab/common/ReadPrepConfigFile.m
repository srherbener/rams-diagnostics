function [ Cstmts ] = ReadPrepConfigFile(Cfile)
% ReadPrepConfigFile function to read config file and reformat into a list
% of statements.

  % Read file line by line. For each line:
  %   strip off comments ('#' to end of line)
  %   put words into a list until 'End' is found
  %   when 'End' is found
  %    terminate the current list
  %    start a new list
  istmnt = 1;
  itoken = 1;
  GotLine = 1;
  
  fid = fopen(char(Cfile));
  while (GotLine == 1)

    InLine = fgetl(fid);
    if (ischar(InLine))
      % strip off comments
      InLine = regexprep(InLine, '#.*', '');
      
      % break line into tokens
      [ Tokens ] = Line2Fields(InLine, ' ');

      % Add tokens to Cstmts until 'End' is encountered. When that happens,
      % start a new statment.
      for i = 1:length(Tokens)
        Tkn = Tokens{i};
        if (~isempty(Tkn))
          if (strcmp(Tkn, 'End'))
            istmnt = istmnt + 1;
            itoken = 1;
          else
            Cstmts{istmnt}{itoken} = Tkn;
            itoken = itoken + 1;
          end
        end
      end
    else
      GotLine = 0;
    end
  end
  fclose(fid);
  
end
