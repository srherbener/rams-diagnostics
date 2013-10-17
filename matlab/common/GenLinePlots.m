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
    clear LineColors;
    clear LineStyles;
    clear LineGscales;
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
      fprintf('WARNING: skipping LinePlot number %d due to no associated PlotSet\n', iplot)
      continue;
    end

    Smooth = Config.LinePlots(iplot).Smooth;
    AddMeas = Config.LinePlots(iplot).AddMeas;
    Ptitle = Config.LinePlots(iplot).Title.Main;
    Pmarkers = Config.LinePlots(iplot).Title.Pmarkers;
    LegLoc = Config.LinePlots(iplot).LegLoc;
    OutFile = sprintf('%s/%s', Pdir, Config.LinePlots(iplot).OutFile);

    i_ap = 1;
    AxisProps(i_ap).Name = 'FontSize';
    AxisProps(i_ap).Val = 35; 
    i_ap = i_ap + 1;

    % X variable, axis specs
    Xvname   = Config.PlotVars(ixv).Var;
    if (strcmp(Config.PlotVars(ixv).Units, ' '))
      Xlabel   = sprintf('%s', Config.PlotVars(ixv).Label);
    else
      Xlabel   = sprintf('%s (%s)', Config.PlotVars(ixv).Label, Config.PlotVars(ixv).Units);
    end
    Xfprefix = Config.PlotVars(ixv).Fprefix;
    Xscale   = Config.PlotVars(ixv).Scale;
    if (Config.PlotVars(ixv).Min > Config.PlotVars(ixv).Max)
      AxisProps(i_ap).Name = 'Xlim';
      AxisProps(i_ap).Val = [ Config.PlotVars(ixv).Max Config.PlotVars(ixv).Min ]; 
      i_ap = i_ap + 1;

      AxisProps(i_ap).Name = 'XDir';
      AxisProps(i_ap).Val = 'Reverse';
      i_ap = i_ap + 1;
    else
      AxisProps(i_ap).Name = 'Xlim';
      AxisProps(i_ap).Val = [ Config.PlotVars(ixv).Min Config.PlotVars(ixv).Max ]; 
      i_ap = i_ap + 1;
    end

    % Y variable, axis specs
    Yvname   = Config.PlotVars(iyv).Var;
    if (strcmp(Config.PlotVars(iyv).Units, ' '))
      Ylabel   = sprintf('%s', Config.PlotVars(iyv).Label);
    else
      Ylabel   = sprintf('%s (%s)', Config.PlotVars(iyv).Label, Config.PlotVars(iyv).Units);
    end
    Yfprefix = Config.PlotVars(iyv).Fprefix;
    Yscale   = Config.PlotVars(iyv).Scale;
    if (Config.PlotVars(iyv).Min > Config.PlotVars(iyv).Max)
      AxisProps(i_ap).Name = 'Ylim';
      AxisProps(i_ap).Val = [ Config.PlotVars(iyv).Max Config.PlotVars(iyv).Min ]; 
      i_ap = i_ap + 1;

      AxisProps(i_ap).Name = 'YDir';
      AxisProps(i_ap).Val = 'Reverse';
      i_ap = i_ap + 1;
    else
      AxisProps(i_ap).Name = 'Ylim';
      AxisProps(i_ap).Val = [ Config.PlotVars(iyv).Min Config.PlotVars(iyv).Max ]; 
      i_ap = i_ap + 1;
    end

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
      LineColors{icase} = Config.PlotSets(ips).Cases(icase).Lcolor;
      LineStyles{icase} = Config.PlotSets(ips).Cases(icase).Lstyle;
      LineGscales(icase) = Config.PlotSets(ips).Cases(icase).Lgscale;

      Xfile = sprintf('%s_%s.h5', Xfprefix, Case);
      fprintf('Reading HDF5 file: %s\n', Xfile);
      Xdata = ReadSelectXyzt(Xfile, Xvname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);
      Xdata(Xdata == UndefVal) = nan;
      if (strcmp(Smooth,'x') || strcmp(Smooth,'xy'))
        Xdata = SmoothFillTseries(Xdata, length(Xdata), Flen);
      end
      Xall(icase,:) = Xdata * Xscale;
  
      Yfile = sprintf('%s_%s.h5', Yfprefix, Case);
      fprintf('Reading HDF5 file: %s\n', Yfile);
      Ydata = ReadSelectXyzt(Yfile, Yvname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);
      Ydata(Ydata == UndefVal) = nan;
      if (strcmp(Smooth,'y') || strcmp(Smooth,'xy'))
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
        if (size(RsqDiff) == size(Havg'))
          % if one vector is row and the other column, make them both the same
          Havg = Havg';
        end
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
    Plot2dSet( Xall, Yall, Ptitle, Pmarkers, Xlabel, Ylabel, LineColors, LineStyles, LineGscales, LegText, LegLoc, AxisProps, AddMeas, OutFile );
    fprintf('\n');
end

end
