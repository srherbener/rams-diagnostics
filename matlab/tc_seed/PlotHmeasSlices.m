function [ ] = PlotHmeasSlices(ConfigFile)
% PlotHmeasSlices function to plot slice through the 3D hmeas data

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

for ihmplot = 1: length(Config.HmeasSlicePlots)
  Name    = Config.HmeasSlicePlots(ihmplot).Name;
  Vname   = Config.HmeasSlicePlots(ihmplot).Var;
  Fprefix = Config.HmeasSlicePlots(ihmplot).Fprefix;
  Units   = Config.HmeasSlicePlots(ihmplot).Units;
  Descrip = Config.HmeasSlicePlots(ihmplot).Descrip;
  Ptype   = Config.HmeasSlicePlots(ihmplot).Ptype;
  Vtype   = Config.HmeasSlicePlots(ihmplot).Vtype;
  Stype   = Config.HmeasSlicePlots(ihmplot).Stype;

  S1 = Config.HmeasSlicePlots(ihmplot).S1;
  S2 = Config.HmeasSlicePlots(ihmplot).S2;
  S3 = Config.HmeasSlicePlots(ihmplot).S3;

  Rmin = Config.HmeasSlicePlots(ihmplot).Rmin;
  Rmax = Config.HmeasSlicePlots(ihmplot).Rmax;
  Zmin = Config.HmeasSlicePlots(ihmplot).Zmin;
  Zmax = Config.HmeasSlicePlots(ihmplot).Zmax;
  Tmin = Config.HmeasSlicePlots(ihmplot).Tmin;
  Tmax = Config.HmeasSlicePlots(ihmplot).Tmax;
  Cmin = Config.HmeasSlicePlots(ihmplot).Cmin;
  Cmax = Config.HmeasSlicePlots(ihmplot).Cmax;

  ips = Config.HmeasSlicePlots(ihmplot).PSnum;
  if (ips == 0)
    fprintf('WARNING: skipping ProfPlot number %d due to no associated PlotSet\n', iplot)
  else
    for icase = 1:Config.PlotSets(ips).Ncases
      Case = Config.PlotSets(ips).Cases(icase).Cname;
      CaseLegend = Config.PlotSets(ips).Cases(icase).Legend;

      % skip this iteration if we are on the control case and
      % we are doing a 'diff' plot
      if (strcmp(Ptype, 'diff') && strcmp(Case, ControlCase))
        continue;
      end
     
      fprintf('***********************************************************************\n');
      fprintf('Plotting Histogram Measurements by Slice:\n');
      fprintf('  Name: %s\n', Name);
      fprintf('  Case: %s\n', Case);
      fprintf('  Variable: %s\n', Vname);
      fprintf('  Plot type: %s\n', Ptype);
      fprintf('  Variable type: %s\n', Vtype);
      fprintf('  Slice type: %s\n', Stype);
      fprintf('\n');
    
      InFile = sprintf('%s/hmeas_%s_%s.h5', InDir, Fprefix, Case);
      
      % Read in the histogram data. HDATA will be organized as (r,z,t) where
      %    r --> radius
      %    z --> heights
      %    t --> time
      Hdset = sprintf('/%s_%s', Vname, Vtype);
      fprintf('Reading file: %s, Dataset: %s\n', InFile, Hdset);
      fprintf('\n');
      R = hdf5read(InFile, '/x_coords') / 1000; % change to km
      Z = hdf5read(InFile, '/z_coords') / 1000; % change to km
      T = hdf5read(InFile, '/t_coords') / 3600; % change to hr
      HDATA = hdf5read(InFile, Hdset);
    
      % If doing a 'diff' plot, read in the control data
      if (strcmp(Ptype, 'diff'))
        InFile = sprintf('%s/hmeas_%s_%s.h5', InDir, Fprefix, ControlCase);
        Hdset = sprintf('/%s_%s', Vname, Vtype);
        fprintf('Reading file: %s, Dataset: %s\n', InFile, Hdset);
        fprintf('\n');
        CNTL_HDATA = hdf5read(InFile, Hdset);
        CNTL_T = hdf5read(InFile, '/t_coords') / 3600; % change to hrs
      end
      
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
      
      Pdata = squeeze(HDATA(R1:R2,Z1:Z2,T1:T2));
      
      % If doing a 'diff' plot, subtract off the control data
      if (strcmp(Ptype, 'diff'))
        fprintf('Subtracting off control case: %s\n', ControlCase);
        fprintf('\n');
        CntlT1 = find(CNTL_T >= Tmin, 1, 'first');
        CntlT2 = find(CNTL_T <= Tmax, 1, 'last');
        CntlPdata = squeeze(CNTL_HDATA(R1:R2,Z1:Z2,CntlT1:CntlT2));
        Pdata = Pdata - CntlPdata;
      end

      % Smooth the data
      Pdata = smooth3(Pdata);
      
      % Pdata is (r,z,t) at this point. Cut out the slices now.
      %   Stype    Xaxis    Yaxis     Org of Sdata
      %    R       time     height       (z,t)
      %    T      radius    height       (z,r)
      %    Z      radius     time        (t,r)
      if (strcmp(Stype, 'R'))
        X = Tvals;
        Y = Zvals;
        Xlabel = 'Time (hr)';
        Ylabel = 'Height (km)';
        Sunits = 'km';
        S1ind = find(R >= S1, 1, 'first');
        S2ind = find(R >= S2, 1, 'first');
        S3ind = find(R >= S3, 1, 'first');
        Sdata1 = squeeze(Pdata(S1ind,:,:));
        Sdata2 = squeeze(Pdata(S2ind,:,:));
        Sdata3 = squeeze(Pdata(S3ind,:,:));
      end
      if (strcmp(Stype, 'T'))
        X = Rvals;
        Y = Zvals;
        Xlabel = 'Radius (km)';
        Ylabel = 'Height (km)';
        Sunits = 'hr';
        S1ind = find(T >= S1, 1, 'first');
        S2ind = find(T >= S2, 1, 'first');
        S3ind = find(T >= S3, 1, 'first');
        Sdata1 = squeeze(Pdata(:,:,S1ind))'; % note transpose
        Sdata2 = squeeze(Pdata(:,:,S2ind))';
        Sdata3 = squeeze(Pdata(:,:,S3ind))';
      end
      if (strcmp(Stype, 'Z'))
        X = Rvals;
        Y = Tvals;
        Xlabel = 'Radius (km)';
        Ylabel = 'Time (hr)';
        Sunits = 'km';
        S1ind = find(Z >= S1, 1, 'first');
        S2ind = find(Z >= S2, 1, 'first');
        S3ind = find(Z >= S3, 1, 'first');
        Sdata1 = squeeze(Pdata(:,S1ind,:))'; % note transpose
        Sdata2 = squeeze(Pdata(:,S2ind,:))';
        Sdata3 = squeeze(Pdata(:,S3ind,:))';
      end

      % Get some colormaps
      Fig = figure;
      cmap_def = colormap;
      cmap_rb = redblue;
      close(Fig);

      if (strcmp(Ptype, 'diff'))
        cmap = cmap_rb;
      else
        cmap = cmap_def;
      end
    
      % Plot  - 3 slices
      Fig = figure;
      Ptitle = sprintf('%s (%s): %s, %s = %.1f %s', Descrip, Units, CaseLegend, Stype, S1, Sunits);
      OutFile = sprintf('%s/hmslice_%s_S1_%s.jpg', OutDir, Name, Case);
      contourf(X,Y,Sdata1);
      set(gca, 'FontSize', Fsize);
      shading flat;
      colormap(cmap);
      colorbar;
      caxis([ Cmin Cmax ]);
      title(Ptitle);
      xlabel(Xlabel);
      ylabel(Ylabel);
      fprintf('Writing file: %s\n', OutFile);
      fprintf('\n');
      saveas(Fig, OutFile);
      close(Fig);
    
      Ptitle = sprintf('%s (%s): %s, %s = %.1f %s', Descrip, Units, CaseLegend, Stype, S2, Sunits);
      OutFile = sprintf('%s/hmslice_%s_S2_%s.jpg', OutDir, Name, Case);
      contourf(X,Y,Sdata2);
      set(gca, 'FontSize', Fsize);
      shading flat;
      colormap(cmap);
      colorbar;
      caxis([ Cmin Cmax ]);
      title(Ptitle);
      xlabel(Xlabel);
      ylabel(Ylabel);
      fprintf('Writing file: %s\n', OutFile);
      fprintf('\n');
      saveas(Fig, OutFile);
      close(Fig);
    
      Ptitle = sprintf('%s (%s): %s, %s = %.1f %s', Descrip, Units, CaseLegend, Stype, S3, Sunits);
      OutFile = sprintf('%s/hmslice_%s_S3_%s.jpg', OutDir, Name, Case);
      contourf(X,Y,Sdata3);
      set(gca, 'FontSize', Fsize);
      shading flat;
      colormap(cmap);
      colorbar;
      caxis([ Cmin Cmax ]);
      title(Ptitle);
      xlabel(Xlabel);
      ylabel(Ylabel);
      fprintf('Writing file: %s\n', OutFile);
      fprintf('\n');
      saveas(Fig, OutFile);
      close(Fig);

    end % for cases from plot set
  end % have a valid PlotSet
end % for HmeasSlicePlots

end
