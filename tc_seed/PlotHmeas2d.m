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

% Temperature file 
TempFile = 'tempc';
TempVar = 'tempc';

% Get some colormaps
Fig = figure;
cmap_def = colormap;
cmap_rb = redblue;
close(Fig);

% make the plots
% steal the config for the 3d plots, simply ignore the isosurface spec
% use the Tmin, Tmax to do time averaging between a specified interval
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

    Flevel = Config.HmeasPlot3d(ihmplot).Flevel;

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
    fprintf('  Time averaging bounds: %.2f to %.2f\n', Tmin, Tmax);
    fprintf('\n');

    InFile = sprintf('%s/hmeas_%s_%s.h5', InDir, Fprefix, Case);
    OutFile = sprintf('%s/hmeas2d_%s_%s.jpg', OutDir, Name, Case);
    
    % Read in the histogram data. HDATA will be organized as (r,z) where
    %    r --> radius
    %    z --> heights
    Hdset = sprintf('/%s_%s', Vname, Vtype);
    fprintf('Reading file: %s, Dataset: %s\n', InFile, Hdset);
    fprintf('\n');
    R = hdf5read(InFile, '/x_coords') / 1000; % km
    Z = hdf5read(InFile, '/z_coords') / 1000; % km
    T = hdf5read(InFile, '/t_coords') / 3600; % hr
    HDATA = hdf5read(InFile, Hdset);

    % if adding a freezing level line, generate it from the data in
    % the temperature file
    if (Flevel == 1)
      InFile = sprintf('%s/hmeas_%s_%s.h5', InDir, TempFile, Case);
      Hdset = sprintf('/%s_%s', TempVar, Vtype);
      fprintf('Reading file: %s, Dataset: %s\n', InFile, Hdset);
      fprintf('\n');
      TEMP = hdf5read(InFile, Hdset) / 1000; % km
    end

    % If doing a 'diff' plot, read in the control data
    if (strcmp(Ptype, 'diff'))
      InFile = sprintf('%s/hmeas_%s_%s.h5', InDir, Fprefix, ControlCase);
      Hdset = sprintf('/%s_%s', Vname, Vtype);
      fprintf('Reading file: %s, Dataset: %s\n', InFile, Hdset);
      fprintf('\n');
      CNTL_HDATA = hdf5read(InFile, Hdset);
      T_CNTL = hdf5read(InFile, '/t_coords') / 3600; % hr
    end
    
    % find the indices that go with the specified limits 
    R1 = find(R >= Rmin, 1, 'first');
    R2 = find(R <= Rmax, 1, 'last');
    Z1 = find(Z >= Zmin, 1, 'first');
    Z2 = find(Z <= Zmax, 1, 'last');
    T1 = find(T >= Tmin, 1, 'first');
    T2 = find(T <= Tmax, 1, 'last');
    if (strcmp(Ptype, 'diff'))
      T1_CNTL = find(T_CNTL >= Tmin, 1, 'first');
      T2_CNTL = find(T_CNTL <= Tmax, 1, 'last');
    end

    Rvals = R(R1:R2);
    Zvals = Z(Z1:Z2);

    % Gererate the freezing level from the time average of the temperature.
    if (Flevel == 1)
      TempMean = squeeze(nanmean(TEMP(R1:R2,Z1:Z2,T1:T2),3));

      % TempMean is now (r,z)
      % Calculate a line representing the freezing level
      % Don't use >= 0 in the find command because the wtmean method
      % can result in multiple zeros in a profile
      [ Nr, Nz ] = size(TempMean);
      Flvals = zeros(1,Nr);
      for ir = 1:Nr
        i = find(TempMean(ir,:) > 0, 1, 'last');
        if (length(i) == 0)
          i = 1;
        end
        FLvals(ir) = Zvals(i);
      end
    end
    
    Pdata = squeeze(nanmean(HDATA(R1:R2,Z1:Z2,T1:T2),3));
    
    % switch (r,z) to (z,r)
    Pdata = permute(Pdata, [ 2 1 ]);

    % If doing a 'diff' plot, subtract off the control data
    if (strcmp(Ptype, 'diff'))
      fprintf('Subtracting off control case: %s\n', ControlCase);
      fprintf('\n');
      CntlPdata = squeeze(nanmean(CNTL_HDATA(R1:R2,Z1:Z2,T1_CNTL:T2_CNTL),3));
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

    if (Flevel == 1)
      line(Rvals, FLvals, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2.0);
      text(double(Rvals(end)), double(FLvals(end)), 'FL', 'FontSize', Fsize, ...
           'Color', 'k', 'HorizontalAlignment', 'Right', ...
           'VerticalAlignment', 'Top');
    end
    
    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');
    saveas(Fig, OutFile);
    close(Fig);
    
  end
end

end
