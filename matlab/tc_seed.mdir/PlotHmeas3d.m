function [ ] = PlotHmeas3d(ConfigFile)
% PlotHmeas3d function to plot histogram measured data

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

% make the plots
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

    Isurf = Config.HmeasPlot3d(ihmplot).Isurf;
    
    % skip this iteration if we are on the control case and
    % we are doing a 'diff' plot
    if (strcmp(Ptype, 'diff') && strcmp(Case, ControlCase))
      continue;
    end
 
    fprintf('***********************************************************************\n');
    fprintf('Plotting Histogram Measurements:\n');
    fprintf('  Name: %s\n', Name);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('  Plot type: %s\n', Ptype);
    fprintf('\n');

    InFile = sprintf('%s/hmeas_%s_%s.h5', InDir, Fprefix, Case);
    OutFile = sprintf('%s/hmeas3d_%s_%s.jpg', OutDir, Name, Case);
    
    % Read in the histogram data. HDATA will be organized as (r,z,t) where
    %    r --> radius
    %    z --> heights
    %    t --> time
    Hdset = sprintf('/%s_%s', Vname, Vtype);
    fprintf('Reading file: %s, Dataset: %s\n', InFile, Hdset);
    fprintf('\n');
    R = hdf5read(InFile, '/x_coords');
    Z = hdf5read(InFile, '/z_coords');
    T = hdf5read(InFile, '/t_coords');
    HDATA = hdf5read(InFile, Hdset);

    % If doing a 'diff' plot, read in the control data
    if (strcmp(Ptype, 'diff'))
      InFile = sprintf('%s/hmeas_%s_%s.h5', InDir, Fprefix, ControlCase);
      Hdset = sprintf('/%s_%s', Vname, Vtype);
      fprintf('Reading file: %s, Dataset: %s\n', InFile, Hdset);
      fprintf('\n');
      CNTL_HDATA = hdf5read(InFile, Hdset);
    end
    
    % Change time to hours, radius and height to km
    R = R / 1000;
    T = T / 3600;
    Z = Z / 1000;
    
    % find the indices that go with the specified limits 
    R1 = find(R >= Rmin, 1, 'first');
    R2 = find(R <= Rmax, 1, 'last');
    Z1 = find(Z >= Zmin, 1, 'first');
    Z2 = find(Z <= Zmax, 1, 'last');
    T1 = find(T >= Tmin, 1, 'first');
    T2 = find(T <= Tmax, 1, 'last');

    Tvals = T(T1:T2);
    Rvals = R(R1:R2);
    Zvals = Z(Z1:Z2);
    
    Pdata = double(squeeze(HDATA(R1:R2,Z1:Z2,T1:T2)));
    
    % switch (r,z,t) to (r,t,z)
    Pdata = permute(Pdata, [ 1 3 2 ]);

    % If doing a 'diff' plot, subtract off the control data
    if (strcmp(Ptype, 'diff'))
      fprintf('Subtracting off control case: %s\n', ControlCase);
      fprintf('\n');
      CntlT1 = find(T >= Tmin, 1, 'first');
      CntlT2 = find(T <= Tmax, 1, 'last');
      CntlPdata = double(squeeze(CNTL_HDATA(R1:R2,Z1:Z2,CntlT1:CntlT2)));
      CntlPdata = permute(CntlPdata, [ 1 3 2 ]);
      Pdata = Pdata - CntlPdata;
    end
    
    % Plot 
    Ptitle = sprintf('%s, %s: Surface: %.1f (%s)', CaseTitle, Descrip, Isurf, Units);
    Xlabel = 'Time (hr)';
    Ylabel = 'Radius (km)';
    Zlabel = 'Height (km)';

    [ TG, RG, ZG ] = meshgrid(Tvals, Rvals, Zvals);

    Fig = figure;
    %surf(TG,RG,Pdata);
    P = patch(isosurface(TG,RG,ZG,Pdata,Isurf));
    set(gca, 'FontSize', Fsize);
    isonormals(TG,RG,ZG,Pdata,P);
    set(P,'FaceColor', 'Red', 'EdgeColor', 'None');
    view(3);
    camlight;
    lighting gouraud;
    %shading flat;
    %cbar = colorbar;
    %set(cbar, 'FontSize', Fsize);
    %caxis([ Cmin Cmax ]);
    grid on;
    axis([ Tmin Tmax Rmin Rmax Zmin Zmax Cmin Cmax ]);
    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);
    zlabel(Zlabel);
    
    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');
    saveas(Fig, OutFile);
    close(Fig);
    
  end
end

end
