function [ ] = GenAdiffEofPlots( ConfigFile )
% GenAdiffEofPlots create a set of azavg difference EOF plots
%
% This routine will read instructions out of 'ConfigFile' and create the
% corresponding plots.
%

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

EofDir = Config.EofDir;
ControlCase = Config.ControlCase;
NumEv = Config.AzavgEofConfig.NumEv;

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  Pcase = Config.Cases(icase).Pname;
  for ieof = 1: length(Config.AzavgEofPlots)
    % Don't run if the case matches the control
    if (~strcmp(Case, ControlCase))

      Fprefix = Config.AzavgEofPlots(ieof).Fprefix;
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
  
      InFile = sprintf('%s/%s_%s.h5', EofDir, Fprefix, Case);
      EofOutFile = sprintf('%s/EOF%d_%s_%s.jpg', Config.PlotDir, EofNum, Fprefix, Case);
      PcOutFile = sprintf('%s/PC%d_%s_%s.jpg', Config.PlotDir, EofNum, Fprefix, Case);
      EsOutFile = sprintf('%s/ES%d_%s_%s.jpg', Config.PlotDir, EofNum, Fprefix, Case);
  
      SelectData = [ Rmin Rmax Zmin Zmax ];
  
      Clevs = (-Clim:Cinc:Clim);
      Cbounds = [ -Clim Clim ];
  
      AdiffEofPlot(InFile, EofOutFile, PcOutFile, EsOutFile, Vname, Vunits, Pcase, EofNum, NumEv, SelectData, Clevs, Cbounds);
  
      fprintf('\n');
    end
  end
end

end
