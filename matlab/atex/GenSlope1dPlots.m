function [ ] = GenSlope1dPlots(ConfigFile)
% GenSlope1dPlots generate slope (POP) plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);

Pdir = Config.PlotDir;
% make sure output directory exists
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

Fsize = 20;

for iplot = 1:length(Config.Slope1dPlots)
    clear SLOPES;

    InDir = Config.Slope1dPlots(iplot).InDir;
    Svar  = Config.Slope1dPlots(iplot).Svar;
    Cvar  = Config.Slope1dPlots(iplot).Cvar;
    XLvar = Config.Slope1dPlots(iplot).XLvar;
    XUvar = Config.Slope1dPlots(iplot).XUvar;
    Tvar  = Config.Slope1dPlots(iplot).Tvar;
    Nsamp = Config.Slope1dPlots(iplot).Nsamp;
    Tcor  = Config.Slope1dPlots(iplot).Tcor;
    Xmin  = Config.Slope1dPlots(iplot).Xmin;
    Xmax  = Config.Slope1dPlots(iplot).Xmax;
    Tmin  = Config.Slope1dPlots(iplot).Tmin;
    Tmax  = Config.Slope1dPlots(iplot).Tmax;
    Cmin  = Config.Slope1dPlots(iplot).Cmin;
    Cmax  = Config.Slope1dPlots(iplot).Cmax;

    Ptitle = sprintf('%s', Config.Slope1dPlots(iplot).Title);
    Xlabel = Config.Slope1dPlots(iplot).Xlabel;
    Ylabel = Config.Slope1dPlots(iplot).Ylabel;
    OutFile = sprintf('%s/%s', Pdir, Config.Slope1dPlots(iplot).OutFile);

    fprintf('Generating slope plot:\n');
    fprintf('  Input directory: %s\n', InDir);
    fprintf('  Time range for data selection (hr):\n');
    fprintf('    Tmin: %.2f\n', Tmin);
    fprintf('    Tmax: %.2f\n', Tmax);
    fprintf('  X Range:\n');
    fprintf('    Xmin: %.2f\n', Xmin);
    fprintf('    Xmax: %.2f\n', Xmax);
    fprintf('  Parameters for Fisher z-score screening:\n');
    fprintf('    Number of samples: %d\n', Nsamp);
    fprintf('    Test correlation: %.2f\n', Tcor);
    fprintf('\n');
    
    % Read in the files and collect into a 2D array (x,y,t). (x,t) comes
    % from the contents of the files, and (y) comes from set of files. The
    % Legend entry on the PlotSet holds the y coordinate values.
    ips = Config.Slope1dPlots(iplot).PSnum;
    
    if (ips == 0)
      fprintf('WARNING: skipping Slope1dPlot number %d due to no associated PlotSet\n', iplot)
      continue;
    end
    
    Ny = Config.PlotSets(ips).Ncases;
    Ycenters = zeros([1 Ny]);
    for icase = 1:Ny
      Case = Config.PlotSets(ips).Cases(icase).Cname;
      Yval = Config.PlotSets(ips).Cases(icase).Legend; 
      Ycenters(icase) = sscanf(Yval, '%d');
      Hfile = sprintf('%s/%s.h5', InDir, Case);
      fprintf('Reading: %s\n', Hfile);
      
      Slopes = hdf5read(Hfile, Svar);
      Cor    = hdf5read(Hfile, Cvar);
      if (icase == 1)
        % read in bin edge values and time
        XL = hdf5read(Hfile, XLvar);
        XU = hdf5read(Hfile, XUvar);
        T  = hdf5read(Hfile, Tvar)/3600; % hr
        
        % initialize the composite variables
        Nx = length(XL);
        Nt = length(T);
        S = zeros([Nx Ny Nt]);
        C = zeros([Nx Ny Nt]);
      end
      
      S(:,icase,:) = Slopes;
      C(:,icase,:) = Cor;
    end
    fprintf('\n');

    % find the indices for selecting out a time range
    T1 = find(T >= Tmin, 1, 'first');
    T2 = find(T <= Tmax, 1, 'last');
    
    X1 = find(XL >= Xmin, 1, 'first');
    X2 = find(XU <= Xmax, 1, 'last');

    % Select out the slope data based on the time range
    S = S(X1:X2,:,T1:T2);
    C = C(X1:X2,:,T1:T2);

    % Select out the slopes that have >= 95% statistical significance
    % Use the fisher z-score to test the null hypothesis that the
    % population correlation is <= Tcor. We want to reject the null
    % hypothesis since we want to keep the slopes from regression fits
    % that indicate that the population correlation > Tcor. For 2-tailed
    % test at 95% significance, need the z-score >= 1.96. Place nans
    % in S where the z-score is < 1.96 so that these slopes do not get
    % figured into the time average.
    S(FisherZscore(Tcor, C, Nsamp) < 1.96) = nan;

    % Save time average of slopes
    % Put the LTSS in the rows, so that in the contour plot
    % the LTSS will end up being the y-axis
    SLOPES = nanmean(S,3)';

    % Use pcolor to make the plot. This will look much better than using contourf
    % since contourf will try to interpolate between values in the grids when it
    % is trying to build the contour lines. Pcolor will simple pick a color for
    % each box that represents the slope value.
    %
    % pcolor(X,Y,C)
    %
    % Pcolor treats the X and Y values as if they were the edges of individual
    % grid cells. This means that the last row and column in C are not used.
    %
    % In order to get all of the data in SLOPES plotted, need to
    %
    %   Form X using the union of XL and XU. Ditto for Y
    %
    %   Pad SLOPES with an extra column and row at the ends of the dimensions.


    % use the bin edges for the plot, bin centers for the tick marks.
    X = [ XL(X1:X2)' XU(X2) ];
    Xcenters = (XL + XU)' / 2;
    Xcenters = Xcenters(X1:X2);
    
    Yshift = (Ycenters(2) - Ycenters(1)) / 2; % assume even spacing
    Y = Ycenters - Yshift;
    Y(end+1) = Ycenters(end) + Yshift;    

    % add an extra row and column at the ends of the rows and columns
    SLOPES(:,end+1) = nan;
    SLOPES(end+1,:) = nan;

    Fig = figure;
    pcolor(X, Y, SLOPES);
    shading flat;
    set(gca, 'FontSize', Fsize);
    %set(gca, 'XTick', X);
    set(gca, 'YTick', Ycenters);
    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);
    cbar = colorbar;
    set(cbar, 'FontSize', Fsize);
    colormap(redblue);
    caxis([ Cmin Cmax ]);
    set(gca, 'Layer', 'top'); % place the axes (tick marks, box, etc) on top of the plot

    fprintf('  Writing plot file: %s\n', OutFile);
    saveas(Fig, OutFile);
    close(Fig);
    fprintf('\n');
end
