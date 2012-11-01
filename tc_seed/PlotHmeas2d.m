function [ ] = PlotHmeas2d(ConfigFile)
% PlotHmeas2d function to plot histogram measured data

[ Config ] = ReadConfig(ConfigFile);

InDir = Config.DiagDir;
OutDir = Config.PlotDir;
ControlCase = Config.ControlCase;

% Make sure output directory exists
if (exist(OutDir, 'dir') ~= 7)
    mkdir(OutDir);
end

% big font
Fsize = 20;

% Get some colormaps
Fig = figure;
cmap_def = colormap;
cmap_rb = redblue;
close(Fig);


% make the plots
% steal the config for the 3d plots, simply ignore the isosurface spec
for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  CaseTitle = regexprep(Case, '_', '-');

  for ihmplot = 1: length(Config.HmeasPlot3d)
    Name    = Config.HmeasPlot3d(ihmplot).Name;
    Vname   = Config.HmeasPlot3d(ihmplot).Var;
    Fprefix = Config.HmeasPlot3d(ihmplot).Fprefix;
    Units   = Config.HmeasPlot3d(ihmplot).Units;
    Descrip = Config.HmeasPlot3d(ihmplot).Descrip;
    Ptype   = Config.HmeasPlot3d(ihmplot).Ptype;
    Vtype   = Config.HmeasPlot3d(ihmplot).Vtype;

    Rmin = Config.HmeasPlot3d(ihmplot).Rmin;
    Rmax = Config.HmeasPlot3d(ihmplot).Rmax;
    Zmin = Config.HmeasPlot3d(ihmplot).Zmin;
    Zmax = Config.HmeasPlot3d(ihmplot).Zmax;
    Tmin = Config.HmeasPlot3d(ihmplot).Tmin;
    Tmax = Config.HmeasPlot3d(ihmplot).Tmax;
    Cmin = Config.HmeasPlot3d(ihmplot).Cmin;
    Cmax = Config.HmeasPlot3d(ihmplot).Cmax;

    % skip this iteration if we are on the control case and
    % we are doing a 'diff' plot
    if (strcmp(Ptype, 'diff') && strcmp(Case, ControlCase))
      continue;
    end
 
    fprintf('***********************************************************************\n');
    fprintf('Plotting Histogram Measurements, time mean data:\n');
    fprintf('  Name: %s\n', Name);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('  Plot type: %s\n', Ptype);
    fprintf('\n');

    InFile = sprintf('%s/hmeas_%s_%s.h5', InDir, Fprefix, Case);
    OutFile = sprintf('%s/hmeas2d_%s_%s.jpg', OutDir, Name, Case);
    
    % Read in the histogram data. HDATA will be organized as (r,z) where
    %    r --> radius
    %    z --> heights
    Hdset = sprintf('/%s_%s_tmean', Vname, Vtype);
    fprintf('Reading file: %s, Dataset: %s\n', InFile, Hdset);
    fprintf('\n');
    R = hdf5read(InFile, '/x_coords') / 1000; % km
    Z = hdf5read(InFile, '/z_coords') / 1000; % km
    HDATA = hdf5read(InFile, Hdset);

    % If doing a 'diff' plot, read in the control data
    if (strcmp(Ptype, 'diff'))
      InFile = sprintf('%s/hmeas_%s_%s.h5', InDir, Fprefix, ControlCase);
      Hdset = sprintf('/%s_%s_tmean', Vname, Vtype);
      fprintf('Reading file: %s, Dataset: %s\n', InFile, Hdset);
      fprintf('\n');
      CNTL_HDATA = hdf5read(InFile, Hdset);
    end
    
    % find the indices that go with the specified limits 
    R1 = find(R >= Rmin, 1, 'first');
    R2 = find(R <= Rmax, 1, 'last');
    Z1 = find(Z >= Zmin, 1, 'first');
    Z2 = find(Z <= Zmax, 1, 'last');

    Rvals = R(R1:R2);
    Zvals = Z(Z1:Z2);
    
    Pdata = squeeze(HDATA(R1:R2,Z1:Z2));
    
    % switch (r,z) to (z,r)
    Pdata = permute(Pdata, [ 2 1 ]);

    % If doing a 'diff' plot, subtract off the control data
    if (strcmp(Ptype, 'diff'))
      fprintf('Subtracting off control case: %s\n', ControlCase);
      fprintf('\n');
      CntlPdata = squeeze(CNTL_HDATA(R1:R2,Z1:Z2));
      CntlPdata = permute(CntlPdata, [ 2 1 ]);
      Pdata = Pdata - CntlPdata;
    end
    
    % Plot 
    Ptitle = sprintf('%s, %s (%s)', CaseTitle, Descrip, Units);
    Xlabel = 'Radius (km)';
    Ylabel = 'Height (km)';

    Fig = figure;
    P = contourf(Rvals,Zvals,Pdata);
    set(gca, 'FontSize', Fsize);
    if (strcmp(Ptype, 'diff'))
      colormap(cmap_rb);
    else
      colormap(cmap_def);
    end
    shading flat;
    cbar = colorbar;
    set(cbar, 'FontSize', Fsize);
    caxis([ Cmin Cmax ]);
    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);
    
    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');
    saveas(Fig, OutFile);
    close(Fig);
    
  end
end

end
