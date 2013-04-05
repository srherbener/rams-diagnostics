function [ ] = PlotVintCond(ConfigFile)
% PlotVintCond plot hovmoller (time vs radius) of vertically integrated condensate 
%

[ Config ] = ReadConfig(ConfigFile);

UndefVal = Config.UndefVal;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

Fsize = 40;

% limit view into plot
Rmin = 0;
Rmax = 150; %km
Tmin = 24;
Tmax = 144;

% color/contour limits
Cmin = 0;
Cmax = 30;

%PanelMarkers = { 'd' 'e' 'f' 'g' };
PanelMarkers = { 'b' 'e' 'f' 'c' };

% Tick marks
Xticks = [ 50 100 ];
Xticklabels = { '50' '100' };

Yticks = [ 40 80 120 ];
Yticklabels = { '40' '80' '120' };

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  Pcase = Config.Cases(icase).Pname;
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
  Pfile  = sprintf('%s/vint_cond_%s.jpg', PlotDir, Case);
  Ptitle = sprintf('%s) %s', PanelMarkers{icase}, Pcase);
  Xlabel = sprintf('Radius (km)');
  Ylabel = sprintf('Time (hr)');
   
  Fig = figure;
  
  contourf(Rvals, Tvals, Pdata');
  set(gca,'FontSize', Fsize);
  set(gca, 'XTick', Xticks);
  set(gca, 'XTickLabel', Xticklabels);
  set(gca, 'YTick', Yticks);
  set(gca, 'YTickLabel', Yticklabels);
  colormap(flipud(colormap('gray')));
  shading flat;
  caxis([ Cmin Cmax ]);
  
  % The title is in a box that adjusts to the amount of characters in
  % the title. Ie, it doesn't do any good to do Left/Center/Right
  % alignment. But, the entire box can be moved to the left side of the
  % plot.
  T = title(Ptitle);
  set(T, 'Units', 'Normalized');
  set(T, 'HorizontalAlignment', 'Left');
  Tpos = get(T, 'Position');
  Tpos(1) = 0; % line up with left edge of plot area
  set(T, 'Position', Tpos);

  xlabel(Xlabel);
  ylabel(Ylabel);
  cbar = colorbar;
  set(cbar, 'FontSize', Fsize);

  line([ 40 40 ], [ Tmin Tmax ], 'Color' , 'k', 'LineStyle', '--', 'LineWidth', 2);
  line([ 70 70 ], [ Tmin Tmax ], 'Color' , 'k', 'LineStyle', '--', 'LineWidth', 2);
  text(20,  Tmin, 'SC', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', 'FontSize', 30, 'Color', 'k');
  text(55,  Tmin, 'RB', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', 'FontSize', 30, 'Color', 'k');
  text(110, Tmin, 'FF', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', 'FontSize', 30, 'Color', 'k');
  
  % Fix up the positioning
  Ppos = get(gca, 'Position'); % position of plot area
  Ppos(1) = Ppos(1) * 1.00;
  Ppos(2) = Ppos(2) * 1.00;
  Ppos(3) = Ppos(3) * 0.85;
  Ppos(4) = Ppos(4) * 0.85;
  set(gca, 'Position', Ppos);
  
  fprintf('Writing plot file: %s\n', Pfile);
  saveas(Fig, Pfile);
  close(Fig);

end

end
