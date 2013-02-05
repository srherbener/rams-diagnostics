function [ ] = PlotVintCond(ConfigFile)
% PlotVintCond plot hovmoller (time vs radius) of vertically integrated condensate 
%

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

Fsize = 20;

% limit view into plot
Rmin = 0;
Rmax = 150; %km
Tmin = 24;
Tmax = 144;

% color/contour limits
Cmin = 0;
Cmax = 30;

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  CaseTitle = regexprep(Case, '.*_', '');
  fprintf('Plotting VintCond for case: %s\n', Case)
  fprintf('\n');
  
  % read in the vint_cond and coordinate data
  % data is (x,y,z,t) where x is radius, y and z are dummy dimensions
  % VC will be (r,t) after the squeeze(hdf5read())
  Hfile = sprintf('%s/vint_cond_%s.h5', AzavgDir, Case);
  fprintf('Reading HDF5 file: %s\n', Hfile);
  VC = squeeze(hdf5read(Hfile, '/vint_cond'));
  R = squeeze(hdf5read(Hfile, '/x_coords')) / 1000; % in km
  T = squeeze(hdf5read(Hfile, '/t_coords')) / 3600; % in hr

  % trim off radial and vertical ranges
  R1 = find(R >= Rmin, 1, 'first');
  R2 = find(R <= Rmax, 1, 'last');
  T1 = find(T >= Tmin, 1, 'first');
  T2 = find(T <= Tmax, 1, 'last');

  % Build data for the plot
  % Convert undefined values to nan so they are excluded from the plot
  % Take average over time (dimension number 3)
  Rvals = R(R1:R2);
  Tvals = T(T1:T2);
  Pdata = VC(R1:R2,T1:T2);
  Pdata (Pdata  == UndefVal) = nan;

  % plot
  Pfile  = sprintf('%s/vint_cond_%s.fig', PlotDir, Case);
  Ptitle = sprintf('Vertically Integrated Condensate (mm): %s', CaseTitle);
  Xlabel = sprintf('Radius (km)');
  Ylabel = sprintf('Time (hr)');
   
  Fig = figure;
  
  contourf(Rvals, Tvals, Pdata');
  set(gca,'FontSize', Fsize);
  colormap(flipud(colormap('gray')));
  shading flat;
  caxis([ Cmin Cmax ]);
  title(Ptitle);
  xlabel(Xlabel);
  ylabel(Ylabel);
  cbar = colorbar;
  set(cbar, 'FontSize', Fsize);

  line([ 40 40 ], [ Tmin Tmax ], 'Color' , 'r', 'LineStyle', '--', 'LineWidth', 2);
  line([ 70 70 ], [ Tmin Tmax ], 'Color' , 'r', 'LineStyle', '--', 'LineWidth', 2);
  text(20,  Tmin, 'EW', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', 'FontSize', Fsize, 'Color', 'k');
  text(55,  Tmin, 'CO', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', 'FontSize', Fsize, 'Color', 'k');
  text(110, Tmin, 'RB', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', 'FontSize', Fsize, 'Color', 'k');
  
  fprintf('Writing plot file: %s\n', Pfile);
  saveas(Fig, Pfile);
  close(Fig);

end

end
