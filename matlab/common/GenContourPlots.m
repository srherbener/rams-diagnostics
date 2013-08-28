function [ ] = GenContourPlots(ConfigFile)
% GenContourPlots function generate contour plots spec'd in a config file

[ Config ] = ReadConfig(ConfigFile);
Pdir = Config.PlotDir;

% make the plots
for iplot = 1:length(Config.ContourPlots)

    % check associations
    ixv = Config.ContourPlots(iplot).XVnum;
    iyv = Config.ContourPlots(iplot).YVnum;
    izv = Config.ContourPlots(iplot).ZVnum;
    ids = Config.ContourPlots(iplot).DSnum;
    ips = Config.ContourPlots(iplot).PSnum;
    if (ixv == 0)
      fprintf('WARNING: skipping ContourPlot number %d due to no associated PlotVar (x)\n', iplot)
      continue;
    end
    if (iyv == 0)
      fprintf('WARNING: skipping ContourPlot number %d due to no associated PlotVar (y)\n', iplot)
      continue;
    end
    if (izv == 0)
      fprintf('WARNING: skipping ContourPlot number %d due to no associated PlotVar (z)\n', iplot)
      continue;
    end
    if (ids == 0)
      fprintf('WARNING: skipping ContourPlot number %d due to no associated PlotDataSelect\n', iplot)
      continue;
    end
    if (ips == 0)
      fprintf('WARNING: skipping ContourPlot number %d due to no associated PlotSet\n', iplot)
      continue;
    end

    Ptitle      = Config.ContourPlots(iplot).Title.Main;
    Pmarkers    = Config.ContourPlots(iplot).Title.Pmarkers;
    Cfill       = Config.ContourPlots(iplot).Fill;
    Cbar        = Config.ContourPlots(iplot).Cbar;
    Cnlevs      = Config.ContourPlots(iplot).Cnlevs;
    Cmap        = Config.ContourPlots(iplot).Cmap;
    OfilePrefix = sprintf('%s/%s', Pdir, Config.ContourPlots(iplot).OutFprefix);

    Npm = length(Pmarkers);

    AxisProps(1).Name = 'FontSize';
    AxisProps(1).Val = 25; 

    % X variable contains x axis values
    Xvname   = Config.PlotVars(ixv).Var;
    Xfprefix = Config.PlotVars(ixv).Fprefix;
    if (strcmp(Config.PlotVars(ixv).Units, ' '))
      Xlabel   = sprintf('%s', Config.PlotVars(ixv).Label);
    else
      Xlabel   = sprintf('%s (%s)', Config.PlotVars(ixv).Label, Config.PlotVars(ixv).Units);
    end

    % Y variable contains y axix values
    Yvname   = Config.PlotVars(iyv).Var;
    Yfprefix = Config.PlotVars(iyv).Fprefix;
    if (strcmp(Config.PlotVars(iyv).Units, ' '))
      Ylabel   = sprintf('%s', Config.PlotVars(iyv).Label);
    else
      Ylabel   = sprintf('%s (%s)', Config.PlotVars(iyv).Label, Config.PlotVars(iyv).Units);
    end

    % Z variable contains contour data
    % use min and max for caxis range
    Zvname   = Config.PlotVars(izv).Var;
    Zfprefix = Config.PlotVars(izv).Fprefix;
    Crange   = [ Config.PlotVars(izv).Min Config.PlotVars(izv).Max ];

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
      fprintf('************************************************************************\n');
      Case = Config.PlotSets(ips).Cases(icase).Cname;
      OutFile = sprintf('%s_%s.jpg', OfilePrefix, Case);

      Hfile = sprintf('%s_%s.h5', Xfprefix, Case);
      fprintf('Reading HDF5 file: %s, Dataset: %s\n', Hfile, Xvname);
      X = ReadSelectXyzt(Hfile, Xvname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);

      Hfile = sprintf('%s_%s.h5', Yfprefix, Case);
      fprintf('Reading HDF5 file: %s, Dataset: %s\n', Hfile, Yvname);
      Y = ReadSelectXyzt(Hfile, Yvname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);

      Hfile = sprintf('%s_%s.h5', Zfprefix, Case);
      fprintf('Reading HDF5 file: %s, Dataset: %s\n', Hfile, Zvname);
      Z = ReadSelectXyzt(Hfile, Zvname, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);

      fprintf('\n');
      fprintf('Writing plot file: %s\n', OutFile);
      if (Npm < 1) then
        Pmark = '';
      else
        % cycle through the list of panel markers as needed
        ipm = mod(icase,Npm);
        if (ipm == 0)
          ipm = Npm;
        end
        if (isempty(Pmarkers{ipm}))
          Pmark = '';
        else
          Pmark = Pmarkers{ipm};
        end
      end
      PlotContour( X, Y, Z, Ptitle, Pmark, Xlabel, Ylabel, Cfill, Cbar, Cmap, Crange, Cnlevs, AxisProps, OutFile );
      fprintf('\n');
    end
end

end
