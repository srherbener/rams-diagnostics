function [ ] = GenSlopePlots(ConfigFile)
% GenSlopePlots generate slope (POP) plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
UndefVal = Config.UndefVal;

Ddir = Config.DiagDir;
Pdir = Config.PlotDir;
% make sure output directory exists
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

Fsize = 20;

for iplot = 1:length(Config.SlopePlots)
    clear SLOPES;
    clear Y;

    Hfile = Config.SlopePlots(iplot).InFile;
    Svar  = Config.SlopePlots(iplot).Svar;
    Cvar  = Config.SlopePlots(iplot).Cvar;
    XLvar = Config.SlopePlots(iplot).XLvar;
    XUvar = Config.SlopePlots(iplot).XUvar;
    YLvar = Config.SlopePlots(iplot).YLvar;
    YUvar = Config.SlopePlots(iplot).YUvar;
    Tvar  = Config.SlopePlots(iplot).Tvar;
    Nsamp = Config.SlopePlots(iplot).Nsamp;
    Tcor  = Config.SlopePlots(iplot).Tcor;
    Tmin  = Config.SlopePlots(iplot).Tmin;
    Tmax  = Config.SlopePlots(iplot).Tmax;
    Cmin  = Config.SlopePlots(iplot).Cmin;
    Cmax  = Config.SlopePlots(iplot).Cmax;

    Ptitle = sprintf('%s', Config.SlopePlots(iplot).Title);
    Xlabel = Config.SlopePlots(iplot).Xlabel;
    Ylabel = Config.SlopePlots(iplot).Ylabel;
    OutFile = sprintf('%s/%s', Pdir, Config.SlopePlots(iplot).OutFile);

   fprintf('Generating slope plot:\n');
   fprintf('  Input file: %s\n', Hfile);
   fprintf('  Time range for data selection (hr):\n');
   fprintf('    Tmin: %.2f\n', Tmin);
   fprintf('    Tmax: %.2f\n', Tmax);
   fprintf('  Parameters for Fisher z-score screening:\n');
   fprintf('    Number of samples: %d\n', Nsamp);
   fprintf('    Test correlation: %.2f\n', Tcor);
   fprintf('\n');

    % Slopes and correlation coefficients are organized (x,y,t)
    % If on the first case, also read in lower and upper bin edges, and time
    S = hdf5read(Hfile, Svar);
    C = hdf5read(Hfile, Cvar);

    % read in bin edge values and time
    XL = hdf5read(Hfile, XLvar);
    XU = hdf5read(Hfile, XUvar);
    YL = hdf5read(Hfile, YLvar);
    YU = hdf5read(Hfile, YUvar);
    T  = hdf5read(Hfile, Tvar)/3600; % hr

    % find the indices for selecting out a time range
    T1 = find(T >= Tmin, 1, 'first');
    T2 = find(T <= Tmax, 1, 'last');

    % Select out the slope data based on the time range
    S = S(:,:,T1:T2);
    C = C(:,:,T1:T2);

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


    % use the bin edges as is
    X = [ XL(1:end); XU(end) ]';
    Y = [ YL(1:end); YU(end) ]';

    % add an extra row and column at the ends of the rows and columns
    SLOPES(:,end+1) = nan;
    SLOPES(end+1,:) = nan;

    fprintf('  Writing plot file: %s\n', OutFile);
    Fig = figure;

    pcolor(X, Y, SLOPES);
    shading flat;
    set(gca, 'FontSize', Fsize);
    set(gca, 'XTick', X);
    set(gca, 'YTick', Y);
    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);
    cbar = colorbar;
    set(cbar, 'FontSize', Fsize);
    colormap(redblue);
    caxis([ Cmin Cmax ]);
    set(gca, 'Layer', 'top'); % place the axes (tick marks, box, etc) on top of the plot

    saveas(Fig, OutFile);
    close(Fig);
    fprintf('\n');
end
