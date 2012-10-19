function [ ] = PlotHmeas(ConfigFile)
% PlotHmeas function to plot total KE vs max azimuthally averaged tangential wind

[ Config ] = ReadConfig(ConfigFile);

InDir = Config.DiagDir;
OutDir = Config.PlotDir;

% Make sure output directory exists
if (exist(OutDir, 'dir') ~= 7)
    mkdir(OutDir);
end

Vtype = 'com'; % Center of mass measurement
Rval = 50;  % radius = 50km
Tval = 72;  % time = 72hrs

Rmin = 0;
Rmax = 250;
Tmin = 24;
Tmax = 144;

Fsize = 20;

% make the plots
for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  for ihmeas = 1: length(Config.Hmeas)
    Name = Config.Hmeas(ihmeas).Name;
    Vname = Config.Hmeas(ihmeas).Rvar;
    
    fprintf('***********************************************************************\n');
    fprintf('Plotting Histogram Measurements:\n');
    fprintf('  Name: %s\n', Name);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    InFile = sprintf('%s/%s_%s.h5', InDir, Name, Case);
    OutFileT = sprintf('%s/tseries_%s_%s.jpg', OutDir, Name, Case);
    OutFileR = sprintf('%s/rseries_%s_%s.jpg', OutDir, Name, Case);
    
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
    
    % Change time to hours, radius and height to km
    R = R / 1000;
    T = T / 3600;
    Z = Z / 1000;
    
    % Set limits according to the quantity we are viewing
    if (strcmp(Vname, 'ccn_conc'))
        clims = [ 0 150 ];
        Zmin = 0;
        Zmax = 2.5;
    end
    if (strcmp(Vname, 'w'))
        clims = [ -0.4 0.4 ];
        Zmin = 0;
        Zmax = 10;
    end
    if (strcmp(Vname, 'lh_vapt'))
        clims = [ -35 35 ];
        Zmin = 0;
        Zmax = 20;
    end
    if (strcmp(Vname, 'lh_frzt'))
        clims = [ -1 1 ];
        Zmin = 0;
        Zmax = 20;
    end
    
    % Strip off the dimension according to Ptype
    Rindex = find(R <= Rval, 1, 'last');
    Tindex = find(T <= Tval, 1, 'last');
    
    R1 = find(R >= Rmin, 1, 'first');
    R2 = find(R <= Rmax, 1, 'last');
    Z1 = find(Z >= Zmin, 1, 'first');
    Z2 = find(Z <= Zmax, 1, 'last');
    T1 = find(T >= Tmin, 1, 'first');
    T2 = find(T <= Tmax, 1, 'last');
    
    Tvals = T(T1:T2);
    Rvals = R(R1:R2);
    Zvals = Z(Z1:Z2);
    
    PdataT = squeeze(HDATA(Rindex,Z1:Z2,T1:T2));
    PdataR = squeeze(HDATA(R1:R2,Z1:Z2,Tindex))';
    
    Ylabel = 'Height (km)';
    
    % Plot time series
    Ptitle = sprintf('Time: %s, R = %dkm: %s', Case, Rval, Vname);
    Ptitle = regexprep(Ptitle, '_', '-');
    Xlabel = 'Time (hr)';
    Fig = figure;
    contourf(Tvals,Zvals,PdataT);
    set(gca, 'FontSize', Fsize);
    if (~ strcmp(Vname, 'ccn_conc'))
        colormap(redblue);
    else
        colormap('default');
    end
    shading flat;
    cbar = colorbar;
    set(cbar, 'FontSize', Fsize);
    caxis(clims);
    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);
    
    saveas(Fig, OutFileT);
    close(Fig);
    
    % Plot radial view
    Ptitle = sprintf('Radius: %s, T = %dhr: %s', Case, Tval, Vname);
    Ptitle = regexprep(Ptitle, '_', '-');
    Xlabel = 'Radius (km)';
    Fig = figure;
    contourf(Rvals,Zvals,PdataR);
    set(gca, 'FontSize', Fsize);
    if (~ strcmp(Vname, 'ccn_conc'))
        colormap(redblue);
    else
        colormap('default');
    end
    shading flat;
    cbar = colorbar;
    set(cbar, 'FontSize', Fsize);
    caxis(clims);
    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);
    
    saveas(Fig, OutFileR);
    close(Fig);

  end
end

end
