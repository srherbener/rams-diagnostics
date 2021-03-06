function [ ] = PlotRsProfs(ConfigFile)
% PlotRsProfs generate radial series of vertical profile plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
%
[ Config ] = ReadConfig(ConfigFile);

Ddir = Config.DiagDir;
Pdir = Config.PlotDir;

ControlCase = Config.ControlCase;

% Temperature file 
TempFprefix1 = 'pmeas_tempc_twp4';
TempVar = 'tempc';

% Grab the colormap
Fig = figure;
DefCmap = colormap;
RbCmap = redblue;
GrayCmap = colormap('gray');
RevGrayCmap = flipud(GrayCmap);
close(Fig);

Ncont = 10; % number of contours

Yticks = [ 2 6 10 14 ];
Yticklabels = { '2' '6' '10' '14' };

% Find and replace underscores in Ptitle, Ylabel with blank spaces
for iplot = 1:length(Config.ProfRsPlots)
    clear Profs;
    clear LegText;
    clear Lspecs;

    Fprefix = Config.ProfRsPlots(iplot).Fprefix;
    Var = Config.ProfRsPlots(iplot).Var;
    Scale = Config.ProfRsPlots(iplot).Scale;

    % config for axes
    Cmin = Config.ProfRsPlots(iplot).Cmin;
    Cmax = Config.ProfRsPlots(iplot).Cmax;
    Zmin = Config.ProfRsPlots(iplot).Zmin;
    Zmax = Config.ProfRsPlots(iplot).Zmax;

    % plot specs
    Pspec  = Config.ProfRsPlots(iplot).Pspec;
    Flevel = Config.ProfRsPlots(iplot).Flevel;
    Ptype  = Config.ProfRsPlots(iplot).Ptype;

    Ptitle = Config.ProfRsPlots(iplot).Title.Main;
    Pmarkers = Config.ProfRsPlots(iplot).Title.Pmarkers;
    Tlabel = Config.ProfRsPlots(iplot).Tlabel;
    Zlabel = Config.ProfRsPlots(iplot).Zlabel;
    OutFileBase = sprintf('%s/%s', Pdir, Config.ProfRsPlots(iplot).OutFileBase);
    
    PanelTitle = ~isempty(Pmarkers);
    if (PanelTitle)
        Fsize = 35;
     else
        Fsize = 25;
    end
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end

    ips = Config.ProfRsPlots(iplot).PSnum;
    if (ips == 0)
      fprintf('WARNING: skipping ProfPlot number %d due to no associated PlotSet\n', iplot)
    else
      for icase = 1:Config.PlotSets(ips).Ncases
        Case = Config.PlotSets(ips).Cases(icase).Cname;
        Legend = Config.PlotSets(ips).Cases(icase).Legend;

        % Var is organized (z,t) in the file (profile time series).
        % Height holds the z values, Time holds sim time in seconds.
        Hfile = sprintf('%s/%s_%s.h5', Ddir, Fprefix, Case);
        Hdset = sprintf('/ProfRs_%s', Var);
        fprintf('  HDF5 file: %s, dataset: %s\n', Hfile, Hdset);
        PROF_RS = squeeze(hdf5read(Hfile, Hdset)) .* Scale;
        Z = hdf5read(Hfile, 'Height')/1000; % km
        R = hdf5read(Hfile, 'Radius')/1000; % km
        
        % If doing a diff plot then subtract off the 'control' case
        if (strcmp(Ptype, 'diff'))
            Hfile = sprintf('%s/%s_%s.h5', Ddir, Fprefix, ControlCase);
            Hdset = sprintf('/ProfRs_%s', Var);
            fprintf('  HDF5 file: %s, dataset: %s\n', Hfile, Hdset);
            PROF_RS_CNTL = squeeze(hdf5read(Hfile, Hdset)) .* Scale;
            
            PROF_RS = PROF_RS - PROF_RS_CNTL;
        end

        % Trim off the selected z range
        % Each profile goes into a row of LHV
        Z1 = find(Z >= Zmin, 1, 'first');
        Z2 = find(Z <= Zmax, 1, 'last');
        Zvals = Z(Z1:Z2);
        Pdata = squeeze(PROF_RS(:,Z1:Z2))';
        
        % if adding a freezing level line, generate it from the data in
        % the temperature file
        if (Flevel == 1)
          TempFprefix2 = '';
          if (~isempty(strfind(Fprefix,'AR_RI')))
              TempFprefix2 = 'AR_RI';
          else
            if (~isempty(strfind(Fprefix,'SC_RI')))
                TempFprefix2 = 'SC_RI';
            else
              if (~isempty(strfind(Fprefix,'RB_RI')))
                  TempFprefix2 = 'RB_RI';
              else
                if (~isempty(strfind(Fprefix,'FF_RI')))
                    TempFprefix2 = 'FF_RI';
                else
                  if (~isempty(strfind(Fprefix,'SC_SS')))
                      TempFprefix2 = 'SC_SS';
                  else
                    if (~isempty(strfind(Fprefix,'RB_SS')))
                        TempFprefix2 = 'RB_SS';
                    else
                      if (~isempty(strfind(Fprefix,'FF_SS')))
                          TempFprefix2 = 'FF_SS';
                      end
                    end
                  end
                end
              end
            end
          end
          Hfile = sprintf('%s/%s_%s_%s.h5', Ddir, TempFprefix1, TempFprefix2, Case);
          Hdset = sprintf('/ProfRs_%s', TempVar);
          fprintf('  HDF5 file: %s, dataset: %s\n', Hfile, Hdset);
          TEMP = squeeze(hdf5read(Hfile, Hdset));
          fprintf('\n');

          % TEMP is now (r,z)
          % Calculate a line representing the freezing level
          % Don't use >= 0 in the find command because the wtmean method
          % can result in multiple zeros in a profile
          [ Nr, Nz ] = size(TEMP);
          FLvals = zeros(1,Nr);
          for ir = 1:Nr
            i = find(TEMP(ir,:) > 0, 1, 'last');
            if (isempty(i))
              i = 1;
            end
            FLvals(ir) = Zvals(i);
          end       
        end

        % Create plot
        if (PanelTitle)
            if (isempty(Ptitle))
              PtitleCase = sprintf('(%s) %s', Pmarkers{icase}, Legend);
            else
              PtitleCase = sprintf('(%s) %s: %s', Pmarkers{icase}, Ptitle, Legend);
            end
        else
            PtitleCase = sprintf('%s: %s', Ptitle, Legend);
        end
        
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

        CF = contourf(R, Zvals, Pdata, Ncont);
        if (isempty(CF(:)))
          fprintf('WARNING: Constant Zdata - attempting to render a constant value contour plot\n');
          % change the first value so that data is not constant
          % pick a value that will hopefully map to the same color as the rest of the map
          %   matlab tries to split into 20 contour levels by default so pick a value that
          %   will likely not jump an adjacent contour level
          Pdata(1) = Pdata(1) + (Cmax - Cmin)/200;
          contourf(T, Zvals, Pdata);
        end
        if (strcmp(Ptype, 'diff'))
            % plot the zero contour
            hold on;
            contour(R, Zvals, Pdata, [ 0 0 ], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 5);
        end
        set(gca, 'FontSize', Fsize);
        set(gca, 'LineWidth', 2);
        set(gca, 'TickLength', [ 0.025 0.025 ]);
        set(gca, 'YTick', Yticks);
        set(gca, 'YTickLabel', Yticklabels);
        colormap(cmap);
        shading flat;
        cbar = colorbar;
%        Clabels = num2str(get(cbar, 'Ytick'));
%        Clabels = regexp(Clabels, '[ ][ ]*', 'split');
%        set(cbar, 'YtickLabel', Clabels);
        caxis([ Cmin Cmax ]);
        set(cbar, 'FontSize', Fsize);
        
        if (PanelTitle)
          % The title is in a box that adjusts to the amount of characters in
          % the title. Ie, it doesn't do any good to do Left/Center/Right
          % alignment. But, the entire box can be moved to the left side of the
          % plot.
          T = title(PtitleCase);
          set(T, 'Units', 'Normalized');
          set(T, 'HorizontalAlignment', 'Left');
          Tpos = get(T, 'Position');
          Tpos(1) = 0; % line up with left edge of plot area
          set(T, 'Position', Tpos);
        else
            title(PtitleCase);
        end
        xlabel(Tlabel);
        ylabel(Zlabel);
        
        if (Flevel == 1)
          line(R, FLvals, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2.0);
          text(double(R(end)), double(FLvals(end)), 'FL', 'FontSize', Fsize, ...
               'Color', 'k', 'HorizontalAlignment', 'Right', ...
               'VerticalAlignment', 'Top');
        end

        if (PanelTitle)
          % Fix up the positioning
          Ppos = get(gca, 'Position'); % position of plot area
          Ppos(1) = Ppos(1) * 0.95;
          Ppos(2) = Ppos(2) * 0.95;
          Ppos(3) = Ppos(3) * 0.85;
          Ppos(4) = Ppos(4) * 0.85;
          set(gca, 'Position', Ppos);
        end

        OutFile = sprintf('%s_%s.jpg', OutFileBase, Case);
        fprintf('Writing plot file: %s\n', OutFile);
        saveas(Fig, OutFile);
        fprintf('\n');

        close(Fig);
      end
    end
end
