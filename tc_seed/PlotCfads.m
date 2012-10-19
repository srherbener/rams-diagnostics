function [ ] = PlotCfads(ConfigFile)
% PlotCfads function to plot CFADS of azimuthally averaged histogram data

[ Config ] = ReadConfig(ConfigFile);

InDir = Config.AzavgDir;
OutDir = Config.PlotDir;

% Make sure output directory exists
if (exist(OutDir, 'dir') ~= 7)
    mkdir(OutDir);
end

Rval = 50;  % radius = 50km
Tval = 72;  % time = 72hrs

% make the plots
for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  for ihmeas = 1: length(Config.Hmeas)
    Fprefix = Config.Hmeas(ihmeas).Fprefix;
    Vname = Config.Hmeas(ihmeas).Rvar;
    
    fprintf('***********************************************************************\n');
    fprintf('Plotting CFAD:\n');
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    InFile = sprintf('%s/%s_%s.h5', InDir, Fprefix, Case);
    OutFile = sprintf('%s/cfad_%s_%s.jpg', OutDir, Fprefix, Case);
    
    % Read in the histogram data. HDATA will be organized as (r,b,z,t) where
    %    r --> radius
    %    b --> bins
    %    z --> heights
    %    t --> time
    Hdset = sprintf('/%s', Vname);
    fprintf('Reading file: %s, Dataset: %s\n', InFile, Hdset);
    fprintf('\n');
    R = hdf5read(InFile, '/x_coords');
    B = hdf5read(InFile, '/y_coords');
    Z = hdf5read(InFile, '/z_coords');
    T = hdf5read(InFile, '/t_coords');
    HDATA = hdf5read(InFile, Hdset);
    
    % Change time to hours, radius and height to km
    R = R / 1000;
    T = T / 3600;
    Z = Z / 1000;
    
    % Set limits according to the quantity we are viewing
    if (strcmp(Vname, 'ccn_conc'))
        Xlabel = 'CCN Concentration (#/cc)';
        clims = [ 0 300 ];
        Zmin = 0;
        Zmax = 2.5;
        Bmin = 0;
        Bmax = 150;
    end
    if (strcmp(Vname, 'w'))
        Xlabel = 'Velocity (m/s)';
        clims = [ 0 150 ];
        Zmin = 0;
        Zmax = 10;
        Bmin = -0.4;
        Bmax = 0.4;
    end
    if (strcmp(Vname, 'lh_vapt'))
        Xlabel = 'Heating Rate (K/hr)';
        clims = [ 0 100 ];
        Zmin = 0;
        Zmax = 5;
        Bmin = -5;
        Bmax = 0;
    end
    if (strcmp(Vname, 'lh_frzt'))
        Xlabel = 'Heating Rate (K/hr)';
        clims = [ 0  300 ];
        Zmin = 0;
        Zmax = 10;
        Bmin = -0.3;
        Bmax = 0.3;
    end
    
    % Strip off the dimension according to Ptype
    Rindex = find(R <= Rval, 1, 'last');
    Tindex = find(T <= Tval, 1, 'last');
    
    B1 = find(B >= Bmin, 1, 'first');
    B2 = find(B <= Bmax, 1, 'last');
    Z1 = find(Z >= Zmin, 1, 'first');
    Z2 = find(Z <= Zmax, 1, 'last');
    
    Bvals = B(B1:B2);
    Zvals = Z(Z1:Z2);
    
    Pdata = squeeze(HDATA(Rindex,B1:B2,Z1:Z2,Tindex))';
    
    Ylabel = 'Height (km)';
    Fsize = 20;
    
    % Plot time series
    Ptitle = sprintf('CFAD: %s, %dhr, %dkm: %s', Case, Tval, Rval, Vname);
    Ptitle = regexprep(Ptitle, '_', '-');
    Fig = figure;
    contourf(Bvals,Zvals,Pdata);
    set(gca, 'FontSize', Fsize);
    shading flat;
    cbar = colorbar;
    set(cbar, 'FontSize', Fsize);
    caxis(clims);
    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);
    
    saveas(Fig, OutFile);
    close(Fig);
  end
end

end
