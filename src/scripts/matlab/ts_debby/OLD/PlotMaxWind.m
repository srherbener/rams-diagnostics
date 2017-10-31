function [ ] = PlotMaxWind(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);
  
  Tdir = Config.TsavgDir;
  
  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  Case = 'TSD_DRY_DUST';
  SimLabel = 'DD';
  Sim25Label = 'DD (25m)';
  
  % read in the tangential wind
  Hfile = sprintf('%s/max_speed10m_%s.h5', Tdir, Case);
  Hdset = '/max_speed10m';
  fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);
  SPEED = squeeze(h5read(Hfile, Hdset));
  TIMES = (squeeze(h5read(Hfile, '/t_coords')) ./3600) - 42.0; %hr
  
  % read in the 25m data
  Hfile = sprintf('%s/max_speed25m_%s.h5', Tdir, Case);
  Hdset = '/max_speed25m';
  fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);
  SPEED_25 = squeeze(h5read(Hfile, Hdset));

  % NHC Best Track (every six hours) data
  % time step 1 from the simulation is where the NHC data starts
  % NHC wind is in mph, multiply by 0.447 to convert to m/s
  % Aug 22, 2006, 6Z through Aug 24, 2006, 18Z (11 points)
  NHC_WIND  = [ 35 35 35 40 45 50 50 50 50 50 50 ] * 0.447;
  NHC_TIMES = (0:6:60);
  
  % plot
  OutFile = sprintf('%s/TsDebbyWind_%s.jpg', Pdir, Case);

  Fsize = 22;
  LineW = 2;

  LegendFsize = 15;

  Xrange = [ -6 66 ];
  Yrange = [  0 25 ];
  
  FigWind = figure;
  SimST = plot(TIMES, [ SPEED SPEED_25 ], 'LineWidth', LineW);
  set(gca, 'FontSize', Fsize);
  hold on;
  NhcST = plot(NHC_TIMES, NHC_WIND, '+k', 'LineWidth', LineW);
  
%  T = title('(c)');
%  set(T, 'Units', 'Normalized');
%  set(T, 'HorizontalAlignment', 'Left');
%  Tpos = get(T, 'Position');
%  Tpos(1) = 0; % line up with left edge of plot area
%  set(T, 'Position', Tpos);
  title('Maximum Wind Speed');
  
  xlabel('Time');
  xlim(Xrange);
  set(gca,'xtick', (6:24:54));
  set(gca,'xticklabel', { 'Aug22:12Z', 'Aug23:12Z', 'Aug24:12Z' });
  ylabel('Wind Speed (m/s)');
  ylim(Yrange);
  legend([ NhcST SimST' ], 'NHC Best Track', SimLabel, Sim25Label, 'Location', 'SouthEast');
  
  fprintf('Writing: %s\n', OutFile);
  saveas(FigWind, OutFile);
  close(FigWind);

end
