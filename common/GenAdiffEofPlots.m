function [ ] = GenAdiffEofPlots( ConfigFile )
% GenAdiffEofPlots create a set of azavg difference EOF plots
%
% This routine will read instructions out of 'ConfigFile' and create the
% corresponding plots.
%

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  Pcase = Config.Cases(icase).Pname;
  for ieof = 1: length(Config.AzavgEofPlots)
    % Don't run if the case matches the control
    if (~strcmp(Case, Config.AzavgEofConfig.Control))

      Var = Config.AzavgEofPlots(ieof).Var;
      Vname = Config.AzavgEofPlots(ieof).Name;
      Vunits = Config.AzavgEofPlots(ieof).Units;
      EofNum = Config.AzavgEofPlots(ieof).Num;
      Clim = Config.AzavgEofPlots(ieof).Clim;
      Cinc = Config.AzavgEofPlots(ieof).Cinc;
      Rmin = Config.AzavgEofPlots(ieof).Rmin;
      Rmax = Config.AzavgEofPlots(ieof).Rmax;
      Zmin = Config.AzavgEofPlots(ieof).Zmin;
      Zmax = Config.AzavgEofPlots(ieof).Zmax;

      fprintf('***********************************************************************\n');
      fprintf('Generating azavg EOF plots: \n');
      fprintf('  Variable: %s\n', Var);
      fprintf('  Case: %s\n', Case);
      fprintf('  Contour levels:\n');
      fprintf('    Limit: %.2f\n', Clim);
      fprintf('    Increment: %.2f\n', Cinc);
      fprintf('  Data selection:\n');
      fprintf('    Rmin: %.2f\n', Rmin);
      fprintf('    Rmax: %.2f\n', Rmax);
      fprintf('    Zmin: %.2f\n', Zmin);
      fprintf('    Zmax: %.2f\n', Zmax);
      fprintf('\n');
  
      InFile = sprintf('%s/%s_%s.h5', Config.AzavgEofConfig.Dir, Var, Case);
      EofOutFile = sprintf('%s/EOF_%s_%s.jpg', Config.PlotDir, Var, Case);
      EsOutFile = sprintf('%s/ES_%s_%s.jpg', Config.PlotDir, Var, Case);
  
      SelectData = [ Rmin Rmax Zmin Zmax ];
  
      Clevs = (-Clim:Cinc:Clim);
      Cbounds = [ -Clim Clim ];
  
      AdiffEofPlot(InFile, EofOutFile, EsOutFile, Vname, Vunits, Pcase, EofNum, SelectData, Clevs, Cbounds);
  
      fprintf('\n');
    end
  end
end

end
