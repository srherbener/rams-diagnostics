function [ ] = PlotRsProfs(ConfigFile)
% PlotRsProfs generate radial series of vertical profile plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
Pname   = Config.ExpName;

Ddir = Config.DiagDir;
Pdir = Config.PlotDir;

% Grab the colormap
Fig = figure;
DefCmap = colormap;
close(Fig);

Fsize = 20;

% Find and replace underscores in Ptitle, Ylabel with blank spaces
% Steal the ProfTsPlots config for now
for iplot = 1:length(Config.ProfTsPlots)
    clear Profs;
    clear LegText;
    clear Lspecs;

    Fprefix = Config.ProfTsPlots(iplot).Fprefix;
    Var = Config.ProfTsPlots(iplot).Var;

    % config for axes
    Cmin = Config.ProfTsPlots(iplot).Cmin;
    Cmax = Config.ProfTsPlots(iplot).Cmax;
    Zmin = Config.ProfTsPlots(iplot).Zmin;
    Zmax = Config.ProfTsPlots(iplot).Zmax;

    % plot specs
    Pspec = Config.ProfTsPlots(iplot).Pspec;

    Ptitle = sprintf('%s: %s', Pname, Config.ProfTsPlots(iplot).Title);
    Tlabel = Config.ProfTsPlots(iplot).Tlabel;
    Zlabel = Config.ProfTsPlots(iplot).Zlabel;
    OutFileBase = sprintf('%s/%s', Pdir, Config.ProfTsPlots(iplot).OutFileBase);
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end

    ips = Config.ProfTsPlots(iplot).PSnum;
    if (ips == 0)
      fprintf('WARNING: skipping ProfPlot number %d due to no associated PlotSet\n', iplot)
    else
      for icase = 1:Config.PlotSets(ips).Ncases
        Case = Config.PlotSets(ips).Cases(icase).Cname;

        % Var is organized (z,t) in the file (profile time series).
        % Height holds the z values, Time holds sim time in seconds.
        Hfile = sprintf('%s/%s_%s.h5', Ddir, Fprefix, Case);
        Hdset = sprintf('/ProfRs_%s', Var);
        fprintf('  HDF5 file: %s, dataset: %s\n', Hfile, Hdset);
        PROF_TS = squeeze(hdf5read(Hfile, Hdset));
        Z = hdf5read(Hfile, 'Height')/1000; % km
        R = hdf5read(Hfile, 'Radius')/1000; % km

        % Trim off the selected z range
        % Each profile goes into a row of LHV
        Z1 = find(Z >= Zmin, 1, 'first');
        Z2 = find(Z <= Zmax, 1, 'last');
        Zvals = Z(Z1:Z2);
        Pdata = squeeze(PROF_TS(:,Z1:Z2))';

        % Create plot
        PtitleCase = sprintf('%s: %s', Ptitle, regexprep(Case, '_', '-'));
        Fig = figure;

        cmap = DefCmap;
        if (strcmp(Pspec, 'min_white'))
          % change first entry (min value) to white
          cmap(1,:) = [ 1 1 1 ];
        end
        if (strcmp(Pspec, 'max_white'))
          % change last entry (max value) to white
          cmap(end,:) = [ 1 1 1 ];
        end

        CF = contourf(R, Zvals, Pdata);
        if (length(CF(:)) == 0)
          fprintf('WARNING: Constant Zdata - attempting to render a constant value contour plot\n');
          % change the first value so that data is not constant
          % pick a value that will hopefully map to the same color as the rest of the map
          %   matlab tries to split into 20 contour levels by default so pick a value that
          %   will likely not jump an adjacent contour level
          Pdata(1) = Pdata(1) + (Cmax - Cmin)/200;
          contourf(T, Zvals, Pdata);
        end
        set(gca, 'FontSize', Fsize)
        colormap(cmap);
        shading flat;
        cbar = colorbar;
        set(cbar, 'FontSize', Fsize)
        caxis([ Cmin Cmax ]);
        title(PtitleCase);
        xlabel(Tlabel);
        ylabel(Zlabel);

        OutFile = sprintf('%s_%s.jpg', OutFileBase, Case);
        fprintf('Writing plot file: %s\n', OutFile);
        saveas(Fig, OutFile);
        fprintf('\n');

        close(Fig);
      end
    end
end
