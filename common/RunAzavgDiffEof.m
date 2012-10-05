function [ ] = RunAzavgDiffEof( ConfigFile )
% RunAzavgDiffEof run a set of azavg difference EOF analyses.
%
% This routine will read instructions out of 'ConfigFile' and run the
% corresponding EOF analyses.
%

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

EofDir = Config.EofDir;
ControlCase = Config.ControlCase;

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  for ieof = 1: length(Config.AzavgEof)
    % Don't run if the case matches the control
    if (~strcmp(Case, ControlCase))
      fprintf('***********************************************************************\n');
      fprintf('Running azavg EOF: \n');
      fprintf('  Name: %s\n', Config.AzavgEof(ieof).Name);
      fprintf('  Variable: %s\n', Config.AzavgEof(ieof).Rvar);
      fprintf('  Case: %s\n', Case);
      fprintf('  Control: %s\n', ControlCase);
      fprintf('  Data selection:\n');
      fprintf('    Rmin: %.2f\n', Config.AzavgEof(ieof).Rmin);
      fprintf('    Rmax: %.2f\n', Config.AzavgEof(ieof).Rmax);
      fprintf('    Zmin: %.2f\n', Config.AzavgEof(ieof).Zmin);
      fprintf('    Zmax: %.2f\n', Config.AzavgEof(ieof).Zmax);
      fprintf('    Tmin: %.2f\n', Config.AzavgEof(ieof).Tmin);
      fprintf('    Tmax: %.2f\n', Config.AzavgEof(ieof).Tmax);
      fprintf('  Undefined Value: %.2f\n', Config.UndefVal);
      fprintf('  N-star: %d\n', Config.AzavgEofConfig.Nstar);
      fprintf('\n');

      Name = Config.AzavgEof(ieof).Name;
      Vname = Config.AzavgEof(ieof).Rvar;
      InFile1 = sprintf('%s/%s_%s.h5', Config.AzavgDir, Name, Case);
      InFile2 = sprintf('%s/%s_%s.h5', Config.AzavgDir, Name, ControlCase);
      OutFile = sprintf('%s/%s_%s.h5', EofDir, Name, Case);

      SelectData = [ Config.AzavgEof(ieof).Rmin Config.AzavgEof(ieof).Rmax Config.AzavgEof(ieof).Zmin Config.AzavgEof(ieof).Zmax Config.AzavgEof(ieof).Tmin Config.AzavgEof(ieof).Tmax ];

      AzavgDiffEof(InFile1, InFile2, Vname, OutFile, SelectData, Config.UndefVal, Config.AzavgEofConfig.Nstar);
      fprintf('\n');
    end
  end
end

end
