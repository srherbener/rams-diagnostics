function [ ] = GenFigures(ConfigFile)
% GenFigures function to generate figures

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;

Pdir = Config.PlotDir;

% For smoothing, length of a running mean
Flen = 5;

% make the plots
for iplot = 1:length(Config.LinePlots)
    clear DataSpecs;
    clear AxesSpecs;
    clear LegendSpec;

    % check associations
    ipd = Config.LinePlots(iplot).PDnum;
    ipa = Config.LinePlots(iplot).PAnum;
    ids = Config.LinePlots(iplot).DSnum;
    ips = Config.LinePlots(iplot).PSnum;
    if (ipd == 0)
      fprintf('WARNING: skipping LinePlot number %d due to no associated PlotData\n', iplot)
      continue;
    end
    if (ipa == 0)
      fprintf('WARNING: skipping LinePlot number %d due to no associated PlotAxes\n', iplot)
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

    XAshow = Config.LinePlots(iplot).XAshow;
    YAshow = Config.LinePlots(iplot).YAshow;
    Smooth = Config.LinePlots(iplot).Smooth;
    AxesSpecs.Title = Config.LinePlots(iplot).Title.Main;
    AxesSpecs.Pmarker = Config.LinePlots(iplot).Title.Pmarkers;
    LegendSpec.Loc = Config.LinePlots(iplot).LegLoc;
    OutFile = sprintf('%s/%s', Pdir, Config.LinePlots(iplot).OutFile);

    i_ap = 1;
    AxesSpecs.Props(i_ap).Name = 'FontSize';
    AxesSpecs.Props(i_ap).Val = 25;
    i_ap = i_ap + 1;

    % X variable, axis specs
    Xvname   = Config.PlotVars(ipd).Var;
    if (XAshow > 0)
      % show x-axis
      if (strcmp(Config.PlotVars(ipd).Units, ' '))
        AxesSpecs.Xlabel   = sprintf('%s', Config.PlotVars(ipd).Label);
      else
        AxesSpecs.Xlabel   = sprintf('%s (%s)', Config.PlotVars(ipd).Label, Config.PlotVars(ipd).Units);
      end
    else
      % "hide" x-axis: no label, erase tick labels, keep tick marks
      AxesSpecs.Xlabel = '';
      AxesSpecs.Props(i_ap).Name = 'XTickLabel';
      AxesSpecs.Props(i_ap).Val  = {};
      i_ap = i_ap + 1;
    end
    Xfprefix = Config.PlotVars(ipd).Fprefix;
    Xscale   = Config.PlotVars(ipd).Scale;
    Xoffset  = Config.PlotVars(ipd).Offset;
    Xamin    = Config.PlotAxes(ixa).Min;
    Xamax    = Config.PlotAxes(ixa).Max;
    Xticks   = Config.PlotAxes(ixa).Ticks;
    if (Xamin > Xamax)
      AxesSpecs.Props(i_ap).Name = 'Xlim';
      AxesSpecs.Props(i_ap).Val = [ Xamax Xamin ]; 
      i_ap = i_ap + 1;

      AxesSpecs.Props(i_ap).Name = 'XDir';
      AxesSpecs.Props(i_ap).Val = 'Reverse';
      i_ap = i_ap + 1;
    else
      AxesSpecs.Props(i_ap).Name = 'Xlim';
      AxesSpecs.Props(i_ap).Val = [ Xamin Xamax ]; 
      i_ap = i_ap + 1;
    end
    % axis scale
    AxesSpecs.Props(i_ap).Name = 'Xscale';
    AxesSpecs.Props(i_ap).Val  = Config.PlotAxes(ixa).Scale;
    i_ap = i_ap + 1;
    % axis tick marks
    if (~isempty(Xticks))
      AxesSpecs.Props(i_ap).Name = 'XTick';
      AxesSpecs.Props(i_ap).Val  = Xticks;
      i_ap = i_ap + 1;
    end

    % Y variable, axis specs
    Yvname   = Config.PlotVars(ipa).Var;
    if (YAshow > 0)
      % show y-axis
      if (strcmp(Config.PlotVars(ipa).Units, ' '))
        Ylabel   = sprintf('%s', Config.PlotVars(ipa).Label);
      else
        Ylabel   = sprintf('%s (%s)', Config.PlotVars(ipa).Label, Config.PlotVars(ipa).Units);
      end
    else
      % "hide" y-axis: no label, erase tick labels, keep tick marks
      Ylabel = '';
      AxesSpecs.Props(i_ap).Name = 'YTickLabel';
      AxesSpecs.Props(i_ap).Val  = {};
      i_ap = i_ap + 1;
    end
    Yfprefix = Config.PlotVars(ipa).Fprefix;
    Yscale   = Config.PlotVars(ipa).Scale;
    Yoffset  = Config.PlotVars(ipa).Offset;
    Yamin    = Config.PlotAxes(iya).Min;
    Yamax    = Config.PlotAxes(iya).Max;
    Yticks   = Config.PlotAxes(iya).Ticks;
    if (Yamin > Yamax)
      AxesSpecs.Props(i_ap).Name = 'Ylim';
      AxesSpecs.Props(i_ap).Val = [ Yamax Yamin ]; 
      i_ap = i_ap + 1;

      AxesSpecs.Props(i_ap).Name = 'YDir';
      AxesSpecs.Props(i_ap).Val = 'Reverse';
      i_ap = i_ap + 1;
    else
      AxesSpecs.Props(i_ap).Name = 'Ylim';
      AxesSpecs.Props(i_ap).Val = [ Yamin Yamax ]; 
      i_ap = i_ap + 1;
    end
    % axis scale
    AxesSpecs.Props(i_ap).Name = 'Yscale';
    AxesSpecs.Props(i_ap).Val  = Config.PlotAxes(iya).Scale;
    i_ap = i_ap + 1;
    % axis tick marks
    if (~isempty(Yticks))
      AxesSpecs.Props(i_ap).Name = 'YTick';
      AxesSpecs.Props(i_ap).Val  = Yticks;
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
      Xzoom = Config.PlotSets(ips).Cases(icase).Xzoom;
      Yzoom = Config.PlotSets(ips).Cases(icase).Yzoom;

      Xfile = sprintf('%s_%s.h5', Xfprefix, Case);
      fprintf('Reading HDF5 file: %s\n', Xfile);
      Xdata = ReadSelectXyzt(Xfile, Xvname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax) .* Xzoom;
      Xdata(Xdata == UndefVal) = nan;
      if (strcmp(Smooth,'x') || strcmp(Smooth,'xy'))
        Xdata = SmoothFillTseries(Xdata, length(Xdata), Flen);
      end
      Xall(icase,:) = (Xdata .* Xscale) + Xoffset;
  
      Yfile = sprintf('%s_%s.h5', Yfprefix, Case);
      fprintf('Reading HDF5 file: %s\n', Yfile);
      Ydata = ReadSelectXyzt(Yfile, Yvname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax) .* Yzoom;
      Ydata(Ydata == UndefVal) = nan;
      if (strcmp(Smooth,'y') || strcmp(Smooth,'xy'))
        Ydata = SmoothFillTseries(Ydata, length(Ydata), Flen);
      end
      Yall(icase,:) = (Ydata .* Yscale) + Yoffset;

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
    Fig = figure;
    set(Fig, 'Visible', 'off');
    Plot2dSet( Xall, Yall, Ptitle, Pmarkers, Xlabel, Ylabel, LineColors, LineStyles, LineGscales, LegText, LegLoc, AxisProps, AddMeas, Fig );
    saveas(Fig, OutFile);
    close(Fig); 
    fprintf('\n');
end

end
