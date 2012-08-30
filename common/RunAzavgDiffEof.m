function [ ] = RunAzavgDiffEof( ConfigFile )
% RunAzavgDiffEof run a set of azavg difference EOF analyses.
%
% This routine will read instructions out of 'ConfigFile' and run the
% corresponding EOF analyses.
%

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

fprintf('DEBUG: AzavgDir: %s\n', Config.AzavgDir);
fprintf('DEBUG: TsavgDir: %s\n', Config.TsavgDir);

for icase = 1:length(Config.Cases)
  fprintf('DEBUG: Case: %s\n', Config.Cases{icase});
  for itdir = 1: length(Config.Tdirs)
    fprintf('DEBUG:   Tdir: %s\n', Config.Tdirs{itdir});
    for ieof = 1: length(Config.AzavgEof)
      fprintf('DEBUG:    Eof: %s %.1f %.1f %.1f %.1f\n', Config.AzavgEof(ieof).Var, Config.AzavgEof(ieof).Rmin, Config.AzavgEof(ieof).Rmax, Config.AzavgEof(ieof).Zmin, Config.AzavgEof(ieof).Zmax);
    end
  end
end

end
