function [ ] = PlotSlp(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);
  
  Tdir = Config.TsavgDir;

  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Case = 'TSD_DRY_DUST';
  SimLabel = 'DD';
  
  % read in the sea level pressure
  Hfile = sprintf('%s/min_sea_press_%s.h5', Tdir, Case);
  Hdset = '/min_sea_press';
  fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);
  PRESS = squeeze(h5read(Hfile, Hdset));
  TIMES = (squeeze(h5read(Hfile, '/t_coords')) ./ 3600) - 42.0;  % hr
  
  % NHC Best Track (every six hours) data
  % time step 1 from the simulation is where the NHC data starts
  % Aug 22, 2006, 6Z through Aug 24, 2006, 18Z (11 points)
  NHC_PRESS = [ 1007 1007 1005 1003 1002 1001 1001 1000  999 1000 1000 ];
  NHC_TIMES = (0:6:60);
  
  % plot
  OutFile = sprintf('%s/TsDebbyPress_%s.jpg', Pdir, Case);

  Fsize = 22;
  LineW = 2;

  LegendFsize = 15;

  Xrange = [  -6   66 ];
  Yrange = [ 995 1010 ];
  
  FigPress = figure;
  SimPress = plot(TIMES, PRESS, 'LineWidth', LineW);
  set(gca, 'FontSize', Fsize);
  hold on;
  NhcPress = plot(NHC_TIMES, NHC_PRESS, '+k', 'LineWidth', LineW);
  
%  T = title('(b)');
%  set(T, 'Units', 'Normalized');
%  set(T, 'HorizontalAlignment', 'Left');
%  Tpos = get(T, 'Position');
%  Tpos(1) = 0; % line up with left edge of plot area
%  set(T, 'Position', Tpos);
  title('Minimum SLP');
  
  xlabel('Time');
  xlim(Xrange);
  set(gca,'xtick', (6:24:54));
  set(gca,'xticklabel', { 'Aug22:12Z', 'Aug23:12Z', 'Aug24:12Z' });
  ylabel('Pressure (mb)');
  ylim(Yrange);
  legend([ NhcPress SimPress ], 'NHC Best Track', SimLabel, 'Location', 'NorthEast');
  
  fprintf('Writing: %s\n', OutFile);
  saveas(FigPress, OutFile);
  close(FigPress);

end
