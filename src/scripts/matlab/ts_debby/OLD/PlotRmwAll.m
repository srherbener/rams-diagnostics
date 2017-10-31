function [ ] = PlotRmwAll(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;

  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  Cases = {
   'TSD_DRY_DUST'
   'TSD_DRY_NODUST'
   'TSD_MOIST_DUST'
   'TSD_MOIST_NODUST'
   };

  LegText = {
%   'NHC Best Track'
   'DD'
   'DN'
   'MD'
   'MN'
   };

  % read in the RMW data
  Nc = length(Cases);
  for icase = 1:Nc
    Hfile = sprintf('%s/ts_size_%s.h5', Ddir, Cases{icase});
    Hdset = '/rmw';
    fprintf('Reading: %s (%s)\n', Hfile, Hdset);

    % Read in time series, TS will be a column vector
    TS = squeeze(h5read(Hfile, Hdset));

    if (icase == 1)
      Nt = length(TS);
      TIMES = squeeze(h5read(Hfile, '/time')) - 42.0;
      RMW = zeros([ Nt Nc ]);
    end
    RMW(:,icase) = TS;
  end

  % plot
  OutFile = sprintf('%s/TsDebbyRmwAll.jpg', Pdir);

  Fsize = 22;
  LineW = 2;

  LegendFsize = 15;

  Xrange = [ -6 66 ];
  Yrange = [  0 150 ];
  
  Fig = figure;

  Plot = plot(TIMES, RMW, 'LineWidth', LineW);
  set(gca, 'FontSize', Fsize);
  xlabel('Time');
  xlim (Xrange);
  set(gca,'xtick', (6:24:54));
  set(gca,'xticklabel', { 'Aug22:12Z' 'Aug23:12Z' 'Aug24:12Z' });
  ylabel('RMW (km)');
  ylim(Yrange);

  legend(Plot, LegText, 'Location', 'SouthEast', 'FontSize', LegendFsize);
  
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
