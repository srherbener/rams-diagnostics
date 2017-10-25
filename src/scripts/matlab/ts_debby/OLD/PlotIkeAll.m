function [ ] = PlotAvgIkeAll(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Tdir = Config.TsavgDir;

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

  % read in the Vt data
  Nc = length(Cases);
  for icase = 1:Nc
    Hfile = sprintf('%s/horiz_ke_%s.h5', Tdir, Cases{icase});
    Hdset = '/horiz_ke';
    fprintf('Reading: %s (%s)\n', Hfile, Hdset);

    % Read in time series, TS will be a column vector
    TS = squeeze(h5read(Hfile, Hdset)) .* 1e-7; % make plot look nicer

    if (icase == 1)
      Nt = length(TS);
      TIMES = (squeeze(h5read(Hfile, '/t_coords')) ./3600) - 42.0; % hr
      IKE = zeros([ Nt Nc ]);
    end
    IKE(:,icase) = TS;
  end

  % plot
  OutFile = sprintf('%s/TsDebbyIkeAll.jpg', Pdir);

  Fsize = 22;
  LineW = 2;

  LegendFsize = 15;

  Xrange = [ -6 66 ];
  Yrange = [  0 4 ];
  
  Fig = figure;
  Plot = plot(TIMES, IKE, 'LineWidth', LineW);
  set(gca, 'FontSize', Fsize);
  xlabel('Time');
  xlim (Xrange);
  set(gca,'xtick', (6:24:54));
  set(gca,'xticklabel', { 'Aug22:12Z' 'Aug23:12Z' 'Aug24:12Z' });
  ylabel('IKE (10^7 J)');
  ylim(Yrange);

  legend(Plot, LegText, 'Location', 'SouthEast', 'FontSize', LegendFsize);
  
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
