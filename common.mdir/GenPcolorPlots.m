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

Fsize = 30;

for icase = 1:length(Config.Cases)
  Case  = Config.Cases(icase).Cname;
  Pname = Config.Cases(icase).Pname;
  for iplot = 1:length(Config.PcolorPlots)
    InFile  = sprintf('%s_%s.h5', Config.PcolorPlots(iplot).Fprefix, Case);
    Var     = Config.PcolorPlots(iplot).Var;
    Xmin    = Config.PcolorPlots(iplot).Xmin;
    Xmax    = Config.PcolorPlots(iplot).Xmax;
    Xgroup  = Config.PcolorPlots(iplot).Xgroup;
    Ymin    = Config.PcolorPlots(iplot).Ymin;
    Ymax    = Config.PcolorPlots(iplot).Ymax;
    Ygroup  = Config.PcolorPlots(iplot).Ygroup;
    Tmin    = Config.PcolorPlots(iplot).Tmin;
    Tmax    = Config.PcolorPlots(iplot).Tmax;
    Cmin    = Config.PcolorPlots(iplot).Cmin;
    Cmax    = Config.PcolorPlots(iplot).Cmax;
    Title   = Config.PcolorPlots(iplot).Title;
    Xlabel  = Config.PcolorPlots(iplot).Xlabel;
    Ylabel  = Config.PcolorPlots(iplot).Ylabel;
    OutFile = sprintf('%s/%s_%s.jpg', Pdir, Config.PcolorPlots(iplot).OutFile, Case);

    if (strcmp(Title, ' '))
      Title = regexprep(Pname,'_', ' ');
    end

   fprintf('********************************************************************\n');
   fprintf('Generating psuedo color plot:\n');
   fprintf('  Input file: %s\n', InFile);
   fprintf('  Input variable: %s\n', Var);
   fprintf('  Data selection:\n');
   fprintf('    Xrange, Xgroup: [ %.2f %.2f], %d\n', Xmin, Xmax, Xgroup);
   fprintf('    Yrange, Ygroup: [ %.2f %.2f], %d\n', Ymin, Ymax, Ygroup);
   fprintf('    Trange: [ %.2f %.2f]\n', Tmin, Tmax);
   fprintf('  Plot color range: [ %.2f %.2f ]\n', Cmin, Cmax);
   fprintf('  Output File: %s\n', OutFile);
   fprintf('\n');

   % Read in the counts and do the data selection (along with combining bins)
   C = hdf5read(InFile, Var);
   X = hdf5read(InFile, 'x_coords');
   Y = hdf5read(InFile, 'y_coords');
   T = hdf5read(InFile, 't_coords')/3600; % hr

   T1 = find(T >= Tmin, 1, 'first');
   T2 = find(T <= Tmax, 1, 'last');

   [ COUNTS, XL, XU, YL, YU ] = GenCountBins(C(:,:,:,T1:T2), X, Y, Xmin, Xmax, Xgroup, Ymin, Ymax, Ygroup);
   TIMES = T(T1:T2);

   % If the time range was more than a single point, then sum up counts across the time
   % dimension. Time will be the last dimension.
   if ((T2-T1) > 0)
     COUNTS = sum(COUNTS, ndims(COUNTS));
   end

   % Turn the counts into a pdf
   PDATA = double(COUNTS ./ sum(COUNTS(:)))'; % transpose is for plotting pursposes

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

    Fig = figure;

    pcolor(PX, PY, PDATA);
    shading flat;
    set(gca, 'FontSize', Fsize);
    set(gca, 'Layer', 'top'); % place the axes (tick marks, box, etc) on top of the plot
    title(Title);
    xlabel(Xlabel);
    ylabel(Ylabel);
    cbar = colorbar;
    set(cbar, 'FontSize', Fsize);
    caxis([ Cmin Cmax ]);
 
    saveas(Fig, OutFile);
    close(Fig);
  end
end

end
