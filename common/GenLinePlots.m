function [ ] = GenLinePlots(ConfigFile)
% GenLinePlots function generate line plots spec'd in a config file

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;

Pname   = Config.ExpName;

Pdir = Config.PlotDir;

% For smoothing, length of a running mean
Flen = 5;

% make the plots
for iplot = 1:length(Config.LinePlots)
    clear Xall;
    clear Yall;
    clear LegText;
    clear LineSpecs;
    clear AxisProps;

    % check associations
    ixv = Config.LinePlots(iplot).XVnum;
    iyv = Config.LinePlots(iplot).YVnum;
    ids = Config.LinePlots(iplot).DSnum;
    ips = Config.LinePlots(iplot).PSnum;
    if (ixv == 0)
      fprintf('WARNING: skipping LinePlot number %d due to no associated PlotVar\n', iplot)
      continue;
    end
    if (iyv == 0)
      fprintf('WARNING: skipping LinePlot number %d due to no associated PlotVar\n', iplot)
      continue;
    end
    if (ids == 0)
      fprintf('WARNING: skipping LinePlot number %d due to no associated PlotDselect\n', iplot)
      continue;
    end
    if (ips == 0)
      fprintf('WARNING: skipping TwoDimPlot number %d due to no associated PlotSet\n', iplot)
      continue;
    end

    Smooth = Config.LinePlots(iplot).Smooth;
    AddMeas = Config.LinePlots(iplot).AddMeas;
    %Ptitle = sprintf('%s: %s', Pname, Config.LinePlots(iplot).Title);
    Ptitle = sprintf('%s', Config.LinePlots(iplot).Title);
    LegLoc = Config.LinePlots(iplot).LegLoc;
    OutFile = sprintf('%s/%s', Pdir, Config.LinePlots(iplot).OutFile);

    AxisProps(1).Name = 'FontSize';
    AxisProps(1).Val = 20; 

    % X variable, axis specs
    Xvname   = Config.PlotVars(ixv).Var;
    Xlabel   = sprintf('%s (%s)', Config.PlotVars(ixv).Label, Config.PlotVars(ixv).Units);
    Xfprefix = Config.PlotVars(ixv).Fprefix;
    Xscale   = Config.PlotVars(ixv).Scale;
    AxisProps(2).Name = 'Xlim';
    AxisProps(2).Val = [ Config.PlotVars(ixv).Min Config.PlotVars(ixv).Max ]; 

    % Y variable, axis specs
    Yvname   = Config.PlotVars(iyv).Var;
    Ylabel   = sprintf('%s (%s)', Config.PlotVars(iyv).Label, Config.PlotVars(iyv).Units);
    Yfprefix = Config.PlotVars(iyv).Fprefix;
    Yscale   = Config.PlotVars(iyv).Scale;
    AxisProps(3).Name = 'Ylim';
    AxisProps(3).Val = [ Config.PlotVars(iyv).Min Config.PlotVars(iyv).Max ]; 

    % Data selection specs
    Xmin = Config.PlotDselects(ids).Xmin;
    Xmax = Config.PlotDselects(ids).Xmax;
    Ymin = Config.PlotDselects(ids).Ymin;
    Ymax = Config.PlotDselects(ids).Ymax;
    Zmin = Config.PlotDselects(ids).Zmin;
    Zmax = Config.PlotDselects(ids).Zmax;
    Tmin = Config.PlotDselects(ids).Tmin;
    Tmax = Config.PlotDselects(ids).Tmax;
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end

    for icase = 1:Config.PlotSets(ips).Ncases
      Case = Config.PlotSets(ips).Cases(icase).Cname;
      LegText{icase} = Config.PlotSets(ips).Cases(icase).Legend;
      LineSpecs{icase} = Config.PlotSets(ips).Cases(icase).Lspec;

      Xfile = sprintf('%s_%s.h5', Xfprefix, Case);
      fprintf('Reading HDF5 file: %s\n', Xfile);
      Xdata = ReadSelectXyzt(Xfile, Xvname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);
      if ((strcmp(Smooth,'x') == 0) || (strcmp(Smooth,'xy') == 0))
        Xdata = SmoothFillTseries(Xdata, length(Xdata), Flen);
      end
      Xall(icase,:) = Xdata * Xscale;
  
      Yfile = sprintf('%s_%s.h5', Yfprefix, Case);
      fprintf('Reading HDF5 file: %s\n', Yfile);
      Ydata = ReadSelectXyzt(Yfile, Yvname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);
      if ((strcmp(Smooth,'y') == 0) || (strcmp(Smooth,'xy') == 0))
        Ydata = SmoothFillTseries(Ydata, length(Ydata), Flen);
      end
      Yall(icase,:) = Ydata * Yscale;

      % If AddMeas is one of the known strings, calculate the appropriate measurement and
      % append that value to the legend text so that it will show up on the plot
      % associated with the correct curve.
      Ltext = '';
      if (strcmp(AddMeas, 'IntArea'))
        Area = trapz(Xdata*Xscale, Ydata*Yscale);
        Ltext = sprintf('%.2e', Area);
        fprintf('  Area under curve: %s\n', Ltext);
      end
      if (strcmp(AddMeas, 'IntVol'))
        % x values are the radii
        % each radial band has a volume approx. equal to:
        %    area of radial band X average height of curve in that radial band
        %      --> PI * (R(outside)^2 - R(inside)^2) X (H(outside)+H(inside))/2
        % Total volume is sum of volumes of each radial band.
        Rsq = (Xdata*Xscale).^2;
        H = (Ydata*Yscale).*0.5;
        RsqDiff = Rsq(2:end) - Rsq(1:end-1);
        Havg = H(2:end) + H(1:end-1);
        Vol = sum(pi * (RsqDiff .* Havg));
        Ltext = sprintf('%.2e', Vol);
        fprintf('  Volume under rotated radial surface: %s\n', Ltext);
      end

      if (~strcmp(Ltext,''))
        LegText{icase} = sprintf('%s (%s)', LegText{icase}, Ltext);
      end
    end

    fprintf('\n');
    fprintf('Writing plot file: %s\n', OutFile);
    Plot2dSet( Xall, Yall, Ptitle, Xlabel, Ylabel, LineSpecs, LegText, LegLoc, AxisProps, OutFile );
    fprintf('\n');
end

end
