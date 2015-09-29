function [] = RunGenFigures(FilePattern)
% RunGenFigures - function to run GenFigures() on the list of files represented by the
%                 pattern in FilePattern.

  % Expand FilePattern into a list of file names.
  LsCmd = sprintf('ls -1 %s', FilePattern);
  [ Rcode OutVal ] = system(LsCmd);

  if (Rcode ~= 0)
    fprintf('ERROR: Given file name pattern (%s) did not match any files.\n', FilePattern);
  else
    % Got a list of files: Arrange file names into a list and run GenFigures
    % on each one.
    ConfigList = strsplit(OutVal);

    for i = 1:length(ConfigList)
      Config = ConfigList{i};
      if (~strcmp(Config, ''))
        GenFigures(Config);
      end
    end
  end
end
