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
    ixs = Config.LinePlots(iplot).XSnum;
    iys = Config.LinePlots(iplot).YSnum;
    iss = Config.LinePlots(iplot).SSnum;
    ips = Config.LinePlots(iplot).PSnum;
    if (ixs == 0)
      fprintf('WARNING: skipping LinePlot number %d due to no associated Xspec\n', iplot)
      continue;
    end
    if (iys == 0)
      fprintf('WARNING: skipping LinePlot number %d due to no associated Yspec\n', iplot)
      continue;
    end
    if (iss == 0)
      fprintf('WARNING: skipping LinePlot number %d due to no associated Sspec\n', iplot)
      continue;
    end
    if (ips == 0)
      fprintf('WARNING: skipping TwoDimPlot number %d due to no associated PlotSet\n', iplot)
      continue;
    end

    Smooth = Config.LinePlots(iplot).Smooth;
    Ptitle = sprintf('%s: %s', Pname, Config.LinePlots(iplot).Title);
    LegLoc = Config.LinePlots(iplot).LegLoc;
    OutFile = sprintf('%s/%s', Pdir, Config.LinePlots(iplot).OutFile);

    AxisProps(1).Name = 'FontSize';
    AxisProps(1).Val = 20; 

    % X variable, axis specs
    Xvname   = Config.PlotVars(ixs).Var;
    Xlabel   = sprintf('%s (%s)', Config.PlotVars(ixs).Label, Config.PlotVars(ixs).Units);
    Xfprefix = Config.PlotVars(ixs).Fprefix;
    Xscale   = Config.PlotVars(ixs).Scale;
    AxisProps(2).Name = 'Xlim';
    AxisProps(2).Val = [ Config.PlotVars(ixs).Min Config.PlotVars(ixs).Max ]; 

    % Y variable, axis specs
    Yvname   = Config.PlotVars(iys).Var;
    Ylabel   = sprintf('%s (%s)', Config.PlotVars(iys).Label, Config.PlotVars(iys).Units);
    Yfprefix = Config.PlotVars(iys).Fprefix;
    Yscale   = Config.PlotVars(iys).Scale;
    AxisProps(3).Name = 'Ylim';
    AxisProps(3).Val = [ Config.PlotVars(iys).Min Config.PlotVars(iys).Max ]; 

    % Data selection specs
    Xmin = Config.Sspecs(iss).Xmin;
    Xmax = Config.Sspecs(iss).Xmax;
    Ymin = Config.Sspecs(iss).Ymin;
    Ymax = Config.Sspecs(iss).Ymax;
    Zmin = Config.Sspecs(iss).Zmin;
    Zmax = Config.Sspecs(iss).Zmax;
    Tmin = Config.Sspecs(iss).Tmin;
    Tmax = Config.Sspecs(iss).Tmax;
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end

    for icase = 1:Config.PlotSets(ips).Ncases
      Case = Config.PlotSets(ips).Cases(icase).Cname;
      LegText(icase) = { Config.PlotSets(ips).Cases(icase).Legend };
      LineSpecs(icase) = { Config.PlotSets(ips).Cases(icase).Lspec };

      Xfile = sprintf('%s_%s.h5', Xfprefix, Case);
      fprintf('Reading HDF5 file: %s\n', Xfile);
      Xdata = hdf5read(Xfile, Xvname) * Xscale;
      X = hdf5read(Xfile, '/x_coords');
      Y = hdf5read(Xfile, '/y_coords');
      Z = hdf5read(Xfile, '/z_coords');
      T = hdf5read(Xfile, '/t_coords');
      Xdata = SelectDataXyzt(Xdata, X, Y, Z, T, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);
      if ((strcmp(Smooth,'x') == 0) || (strcmp(Smooth,'xy') == 0))
        Xdata = SmoothFillTseries(Xdata, length(Xdata), Flen);
      end
      Xall(icase,:) = Xdata;
  
      Yfile = sprintf('%s_%s.h5', Yfprefix, Case);
      fprintf('Reading HDF5 file: %s\n', Yfile);
      Ydata = hdf5read(Yfile, Yvname) * Yscale;
      X = hdf5read(Yfile, '/x_coords');
      Y = hdf5read(Yfile, '/y_coords');
      Z = hdf5read(Yfile, '/z_coords');
      T = hdf5read(Yfile, '/t_coords');
      Ydata = SelectDataXyzt(Ydata, X, Y, Z, T, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, Tmin, Tmax);
      if ((strcmp(Smooth,'y') == 0) || (strcmp(Smooth,'xy') == 0))
        Ydata = SmoothFillTseries(Ydata, length(Ydata), Flen);
      end
      Yall(icase,:) = Ydata;
    end

    fprintf('\n');
    fprintf('Writing plot file: %s\n', OutFile);
    Plot2dSet( Xall, Yall, Ptitle, Xlabel, Ylabel, LineSpecs, LegText, LegLoc, AxisProps, OutFile );
    fprintf('\n');
end

end
