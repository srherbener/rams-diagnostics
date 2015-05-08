function [ Config AxesSpecs ] = GenFigures(ConfigFile)
% GenFigures function to generate figures

  [ Config ] = ReadFigureConfig(ConfigFile);
  UndefVal = Config.UndefVal;

  fprintf('Generating Figures: %s\n', ConfigFile);

  for i_panel = 1:length(Config.FigPanels)
    % clear out plot structs
    clear LegendSpecs;

    Pname  = Config.FigPanels(i_panel).Name;
    Smooth = Config.FigPanels(i_panel).Smooth;

    fprintf('  Panel: %s\n', Pname);

    % Check that plot data sets and axes got associated
    ipa = Config.FigPanels(i_panel).PAnum;
    ips = Config.FigPanels(i_panel).PSnum;
    if (ipa == 0)
      fprintf('WARNING: skipping FigPanel number %d due to no associated PlotAxes\n', i_panel)
      continue;
    end
    if (ips == 0)
      fprintf('WARNING: skipping FigPanel number %d due to no associated PlotSet\n', i_panel)
      continue;
    end

    % Fill in the AxesSpecs structure
    [ AxesSpecs ] = GenAxesSpecs(Config, i_panel, ipa);

    % Fill in the DataSpecs structure
    [ DataSpecs ] = GenDataSpecs(Config, i_panel, ips);

    % Legend specs
    LegSpecs.Loc   = Config.FigPanels(i_panel).LegLoc;
    LegSpecs.Fsize = Config.FigPanels(i_panel).LegFsize;


  end
  fprintf('\n');

%%%  % For smoothing, length of a running mean
%%%  Flen = 5;
%%%  
%%%  % make the plots
%%%  for iplot = 1:length(Config.LinePlots)
%%%  
%%%  
%%%      % axis scale
%%%      AxesSpecs.Props(i_ap).Name = 'Yscale';
%%%      AxesSpecs.Props(i_ap).Val  = Config.PlotAxes(iya).Scale;
%%%      i_ap = i_ap + 1;
%%%      % axis tick marks
%%%      if (~isempty(Yticks))
%%%        AxesSpecs.Props(i_ap).Name = 'YTick';
%%%        AxesSpecs.Props(i_ap).Val  = Yticks;
%%%        i_ap = i_ap + 1;
%%%      end
%%%  
%%%      % Data selection specs
%%%      Xmin = Config.PlotDselects(ids).Xmin;
%%%      Xmax = Config.PlotDselects(ids).Xmax;
%%%      Ymin = Config.PlotDselects(ids).Ymin;
%%%      Ymax = Config.PlotDselects(ids).Ymax;
%%%      Zmin = Config.PlotDselects(ids).Zmin;
%%%      Zmax = Config.PlotDselects(ids).Zmax;
%%%      Tmin = Config.PlotDselects(ids).Tmin;
%%%      Tmax = Config.PlotDselects(ids).Tmax;
%%%      
%%%      % make sure output directory exists
%%%      if (exist(Pdir, 'dir') ~= 7)
%%%          mkdir(Pdir);
%%%      end
%%%  
%%%      for icase = 1:Config.PlotSets(ips).Ncases
%%%        Case = Config.PlotSets(ips).Cases(icase).Cname;
%%%        LegText{icase} = Config.PlotSets(ips).Cases(icase).Legend;
%%%        LineColors{icase} = Config.PlotSets(ips).Cases(icase).Lcolor;
%%%        LineStyles{icase} = Config.PlotSets(ips).Cases(icase).Lstyle;
%%%        LineGscales(icase) = Config.PlotSets(ips).Cases(icase).Lgscale;
%%%        Xzoom = Config.PlotSets(ips).Cases(icase).Xzoom;
%%%        Yzoom = Config.PlotSets(ips).Cases(icase).Yzoom;
%%%  
%%%        Xfile = sprintf('%s_%s.h5', Xfprefix, Case);
%%%        fprintf('Reading HDF5 file: %s\n', Xfile);
%%%        Xdata = ReadSelectXyzt(Xfile, Xvname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax) .* Xzoom;
%%%        Xdata(Xdata == UndefVal) = nan;
%%%        if (strcmp(Smooth,'x') || strcmp(Smooth,'xy'))
%%%          Xdata = SmoothFillTseries(Xdata, length(Xdata), Flen);
%%%        end
%%%        Xall(icase,:) = (Xdata .* Xscale) + Xoffset;
%%%    
%%%        Yfile = sprintf('%s_%s.h5', Yfprefix, Case);
%%%        fprintf('Reading HDF5 file: %s\n', Yfile);
%%%        Ydata = ReadSelectXyzt(Yfile, Yvname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax) .* Yzoom;
%%%        Ydata(Ydata == UndefVal) = nan;
%%%        if (strcmp(Smooth,'y') || strcmp(Smooth,'xy'))
%%%          Ydata = SmoothFillTseries(Ydata, length(Ydata), Flen);
%%%        end
%%%        Yall(icase,:) = (Ydata .* Yscale) + Yoffset;
%%%  
%%%        % If AddMeas is one of the known strings, calculate the appropriate measurement and
%%%        % append that value to the legend text so that it will show up on the plot
%%%        % associated with the correct curve.
%%%        Ltext = '';
%%%        if (strcmp(AddMeas, 'IntArea'))
%%%          Area = trapz(Xdata*Xscale, Ydata*Yscale);
%%%          Ltext = sprintf('%.2e', Area);
%%%          fprintf('  Area under curve: %s\n', Ltext);
%%%        end
%%%        if (strcmp(AddMeas, 'IntVol'))
%%%          % x values are the radii
%%%          % each radial band has a volume approx. equal to:
%%%          %    area of radial band X average height of curve in that radial band
%%%          %      --> PI * (R(outside)^2 - R(inside)^2) X (H(outside)+H(inside))/2
%%%          % Total volume is sum of volumes of each radial band.
%%%          Rsq = (Xdata*Xscale).^2;
%%%          H = (Ydata*Yscale).*0.5;
%%%          RsqDiff = Rsq(2:end) - Rsq(1:end-1);
%%%          Havg = H(2:end) + H(1:end-1);
%%%          if (size(RsqDiff) == size(Havg'))
%%%            % if one vector is row and the other column, make them both the same
%%%            Havg = Havg';
%%%          end
%%%          Vol = sum(pi * (RsqDiff .* Havg));
%%%          Ltext = sprintf('%.2e', Vol);
%%%          fprintf('  Volume under rotated radial surface: %s\n', Ltext);
%%%        end
%%%  
%%%        if (~strcmp(Ltext,''))
%%%          LegText{icase} = sprintf('%s (%s)', LegText{icase}, Ltext);
%%%        end
%%%      end
%%%  
%%%      fprintf('\n');
%%%      fprintf('Writing plot file: %s\n', OutFile);
%%%      Fig = figure;
%%%      set(Fig, 'Visible', 'off');
%%%      Plot2dSet( Xall, Yall, Ptitle, Pmarkers, Xlabel, Ylabel, LineColors, LineStyles, LineGscales, LegText, LegLoc, AxisProps, AddMeas, Fig );
%%%      saveas(Fig, OutFile);
%%%      close(Fig); 
%%%      fprintf('\n');
%%%  end

end
