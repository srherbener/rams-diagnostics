function [ ] = PlotTsProfs(ConfigFile)
% PlotTsProfs generate time series of vertical profile plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
%
[ Config ] = ReadConfig(ConfigFile);

Ddir = Config.DiagDir;
Pdir = Config.PlotDir;

ControlCase = Config.ControlCase;

% Temperature file 
TempFprefix1 = 'pmeas_tempc';
TempVar = 'tempc';

% Grab the colormap
Fig = figure;
DefCmap = colormap;
RbCmap = redblue;
GrayCmap = colormap('gray');
RevGrayCmap = flipud(GrayCmap);
close(Fig);


Fsize = 20;

% Find and replace underscores in Ptitle, Ylabel with blank spaces
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
    Pspec  = Config.ProfTsPlots(iplot).Pspec;
    Flevel = Config.ProfTsPlots(iplot).Flevel;
    Ptype  = Config.ProfTsPlots(iplot).Ptype;

    Ptitle = sprintf('%s', Config.ProfTsPlots(iplot).Title);
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
        Legend = Config.PlotSets(ips).Cases(icase).Legend;

        % Var is organized (z,t) in the file (profile time series).
        % Height holds the z values, Time holds sim time in seconds.
        Hfile = sprintf('%s/%s_%s.h5', Ddir, Fprefix, Case);
        Hdset = sprintf('/ProfTs_%s', Var);
        fprintf('  HDF5 file: %s, dataset: %s\n', Hfile, Hdset);
        PROF_TS = squeeze(hdf5read(Hfile, Hdset));
        Z = hdf5read(Hfile, 'Height')/1000; % km
        T = hdf5read(Hfile, 'Time')/3600;   % hr

        % If doing a diff plot then subtract off the 'control' case
        if (strcmp(Ptype, 'diff'))
            Hfile = sprintf('%s/%s_%s.h5', Ddir, Fprefix, ControlCase);
            Hdset = sprintf('/ProfTs_%s', Var);
            fprintf('  HDF5 file: %s, dataset: %s\n', Hfile, Hdset);
            PROF_TS_CNTL = squeeze(hdf5read(Hfile, Hdset));
            
            PROF_TS = PROF_TS - PROF_TS_CNTL;
        end

        % Trim off the selected z range
        % Each profile goes into a row of LHV
        Z1 = find(Z >= Zmin, 1, 'first');
        Z2 = find(Z <= Zmax, 1, 'last');
        Zvals = Z(Z1:Z2);
        Pdata = squeeze(PROF_TS(Z1:Z2,:));

        % if adding a freezing level line, generate it from the data in
        % the temperature file
        if (Flevel == 1)
          TempFprefix2 = '';
          if (~isempty(strfind(Fprefix,'AR_RI')))
              TempFprefix2 = 'AR_RI';
          else
              if (~isempty(strfind(Fprefix,'EW_RI')))
                  TempFprefix2 = 'EW_RI';
              end
          end
          Hfile = sprintf('%s/%s_%s_%s.h5', Ddir, TempFprefix1, TempFprefix2, Case);
          Hdset = sprintf('/ProfTs_%s', TempVar);
          fprintf('  HDF5 file: %s, dataset: %s\n', Hfile, Hdset);
          TEMP = squeeze(hdf5read(Hfile, Hdset));
          fprintf('\n');

          % TEMP is now (z,t)
          % Calculate a line representing the freezing level
          % Don't use >= 0 in the find command because the wtmean method
          % can result in multiple zeros in a profile
          [ Nz, Nt ] = size(TEMP);
          FLvals = zeros(1,Nt);
          for it = 1:Nt
            i = find(TEMP(:,it) > 0, 1, 'last');
            if (isempty(i))
              i = 1;
            end
            FLvals(it) = Zvals(i);
          end       
        end

        % Create plot
        PtitleCase = sprintf('%s: %s', Ptitle, Legend);
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
        if (strcmp(Pspec, 'redblue'))
          % change last entry (max value) to white
          cmap = RbCmap;
        end
        if (strcmp(Pspec, 'revgray'))
          % change last entry (max value) to white
          cmap = RevGrayCmap;
        end
        if (strcmp(Pspec, 'gray'))
          % change last entry (max value) to white
          cmap = GrayCmap;
        end

        CF = contourf(T, Zvals, Pdata);
        if (isempty(CF(:)))
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

        if (Flevel == 1)
          line(T, FLvals, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2.0);
          text(double(T(end)), double(FLvals(end)), 'FL', 'FontSize', Fsize, ...
               'Color', 'k', 'HorizontalAlignment', 'Right', ...
               'VerticalAlignment', 'Top');
        end

        OutFile = sprintf('%s_%s.fig', OutFileBase, Case);
        fprintf('Writing plot file: %s\n', OutFile);
        saveas(Fig, OutFile);
        fprintf('\n');

        close(Fig);
      end
    end
end
