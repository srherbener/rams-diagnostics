function [ ] = GenPcolorPlots(ConfigFile)
% GenPcolorPlots generate psuedo color plots (color coded tiles)

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
UndefVal = Config.UndefVal;

Pdir = Config.PlotDir;
% make sure output directory exists
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

Fsize = 25;

for icase = 1:length(Config.Cases)
  Case  = Config.Cases(icase).Cname;
  Pname = Config.Cases(icase).Pname;
  for iplot = 1:length(Config.PcolorPlots)
    InFile  = sprintf('%s_%s.h5', Config.PcolorPlots(iplot).Fprefix, Case);
    Cmin    = Config.PcolorPlots(iplot).Cmin;
    Cmax    = Config.PcolorPlots(iplot).Cmax;
    Cticks  = Config.PcolorPlots(iplot).Cticks;
    Cscale  = Config.PcolorPlots(iplot).Cscale;
    Ptitle  = Config.PcolorPlots(iplot).Title.Main;
    Pmarkers = Config.PcolorPlots(iplot).Title.Pmarkers;
    Xlabel  = Config.PcolorPlots(iplot).Xlabel;
    Ylabel  = Config.PcolorPlots(iplot).Ylabel;
    OutFile = sprintf('%s/%s_%s.jpg', Pdir, Config.PcolorPlots(iplot).OutFile, Case);

    fprintf('********************************************************************\n');
    fprintf('Generating psuedo color plot:\n');
    fprintf('  Input file: %s\n', InFile);
    fprintf('  Plot color range: [ %.2f %.2f ]\n', Cmin, Cmax);
    fprintf('  Plot color scale: %s\n', Cscale);
    fprintf('  Output File: %s\n', OutFile);
    fprintf('\n');

    % Read in the counts and do the data selection (along with combining bins)
    COUNTS = hdf5read(InFile, 'COUNTS');
    XL     = hdf5read(InFile, 'XL');
    XU     = hdf5read(InFile, 'XU');
    YL     = hdf5read(InFile, 'YL');
    YU     = hdf5read(InFile, 'YU');
 
    % pcolor wants doubles
    PDATA = double(COUNTS)'; % transpose is for plotting pursposes

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
    PX = [ XL(1:end); XU(end) ]';
    PY = [ YL(1:end); YU(end) ]';

    % add an extra row and column at the ends of the rows and columns
    PDATA(:,end+1) = nan;
    PDATA(end+1,:) = nan;

    % convert zeros into nans so they don't show up on the plot
    PDATA(PDATA == 0) = nan;

    Crange = [ Cmin Cmax ];

    % convert to logarithmic color scale if requested
    if (strcmp(Cscale, 'log'))
      PDATA = log10(PDATA);
      Crange = log10(Crange);
      Cticks = log10(Cticks);
      for i = 1:length(Cticks)
        CtickLabels{i} = sprintf('10^%d', Cticks(i));
      end
    end

    Fig = figure;

    pcolor(PX, PY, PDATA);
    shading flat;
    set(gca, 'FontSize', Fsize);
    set(gca, 'Layer', 'top'); % place the axes (tick marks, box, etc) on top of the plot
% kludge to get cloud depth working
if (~isempty(regexp(InFile, '_cdepth_')))
  fprintf('WARNING: have cloud depth, adjusting x ticks\n');
  set(gca, 'XTick', [ 1000 3000 ]);
end
    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);
    cbar = colorbar;
    set(cbar, 'FontSize', Fsize);
    caxis(Crange);
    if (strcmp(Cscale, 'log'))
      set(cbar, 'YTick', Cticks);
      set(cbar, 'YTickLabel', CtickLabels);
    end

    % fix the position
    Ppos = get(gca, 'Position'); % position of plot area
    Ppos(1) = Ppos(1) * 1.00;
    Ppos(2) = Ppos(2) * 1.00;
    Ppos(3) = Ppos(3) * 0.90;
    Ppos(4) = Ppos(4) * 1.00;
    set(gca, 'Position', Ppos);
 
    saveas(Fig, OutFile);
    close(Fig);
  end
end

end
