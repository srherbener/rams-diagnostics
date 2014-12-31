function [ ] = PlotMaxWindAll(ConfigFile)

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
    %Hfile = sprintf('%s/max_speed25m_%s.h5', Tdir, Cases{icase});
    %Hdset = '/max_speed25m';
    Hfile = sprintf('%s/max_speed10m_%s.h5', Tdir, Cases{icase});
    Hdset = '/max_speed10m';
    fprintf('Reading: %s (%s)\n', Hfile, Hdset);

    % Read in time series, S will be a column vector
    S = squeeze(h5read(Hfile, Hdset));

    if (icase == 1)
      Nt = length(S);
      TIMES = (squeeze(h5read(Hfile, '/t_coords')) ./3600) - 42.0; % hr
      SPEED_T = zeros([ Nt Nc ]);
    end
    SPEED_T(:,icase) = S;
  end


  % NHC Best Track (every six hours) data
  % time step 1 from the simulation is where the NHC data starts
  % NHC wind is in mph, multiply by 0.447 to convert to m/s
  %
  NHC_WIND  = [ 35 35 35 40 45 50 50 50 50 50 50 ] * 0.447;
  NHC_TIMES = (0:6:60);

  % plot
  OutFile = sprintf('%s/TsDebbyWindAll.jpg', Pdir);

  Fsize = 22;
  LineW = 2;

  LegendFsize = 15;

  Xrange = [ -6 66 ];
  Yrange = [  0 25 ];
  
  FigWind = figure;
  SimST = plot(TIMES, SPEED_T, 'LineWidth', LineW);
  set(gca, 'FontSize', Fsize);
  hold on;
%  NhcST = plot(NHC_TIMES, NHC_WIND, '+k', 'LineWidth',  LineW);
  title('Maximum 10m Wind Speed');
  xlabel('Time');
  xlim (Xrange);
  set(gca,'xtick', (6:24:54));
  set(gca,'xticklabel', { 'Aug22:12Z' 'Aug23:12Z' 'Aug24:12Z' });
  ylabel('Wind Speed (m s^-^1)');
  ylim(Yrange);

%  legend([ NhcST SimST' ], LegText, 'Location', 'SouthEast', 'FontSize', LegendFsize);
  legend(SimST, LegText, 'Location', 'SouthEast', 'FontSize', LegendFsize);
  
  fprintf('Writing: %s\n', OutFile);
  saveas(FigWind, OutFile);
  close(FigWind);

end
