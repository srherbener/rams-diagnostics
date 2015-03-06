function [ ] = PlotMaxWindAll(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Tdir = Config.TsavgDir;

  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  Cases = {
%   'TSD_DRY_DUST'
%   'TSD_DRY_NODUST'
%   'TSD_MOIST_DUST'
%   'TSD_MOIST_NODUST'

   'TSD_SAL_DUST'
   'TSD_SAL_NODUST'
   };
  Nc = length(Cases);

  LegText = {
   'NHC Best Track'
%   'DD'
%   'DN'
%   'MD'
%   'MN'

   'SAL\_DUST'
   'SAL\_NODUST'
   };

%  InFprefix = 'max_speed10m';
%  InVname = '/max_speed10m';
   InFprefix = 'max_speed25m';
   InVname = '/max_speed25m';

  for icase = 1:Nc
    Hfile = sprintf('%s/%s_%s.h5', Tdir, InFprefix, Cases{icase});
    InFiles{icase} = Hfile;
  end

  % NHC Best Track (every six hours) data
  % time step 1 from the simulation is where the NHC data starts
  % NHC wind is in mph, multiply by 0.447 to convert to m/s
  %
  NHC_WIND  = [ 35 35 35 40 45 50 50 50 50 50 50 ] * 0.447;
  NHC_TIMES = (0:6:60);

  % first find the max number of time steps
  MaxNt = 0;
  for icase = 1:Nc
     Hfile = InFiles{icase};
     Hinfo = h5info(Hfile, '/t_coords');
     Tsize = Hinfo.Dataspace.Size;
     if (Tsize > MaxNt)
       TIMES = (squeeze(h5read(Hfile, '/t_coords')) ./3600) - 42.0; % hr
       MaxNt = Tsize;
     end 
  end

  SPEED_T = nan([ MaxNt Nc ]);

  % read in the Vt data
  for icase = 1:Nc
    Hfile = InFiles{icase};
    fprintf('Reading: %s (%s)\n', Hfile, InVname);

    % Read in time series, S will be a column vector
    S = squeeze(h5read(Hfile, InVname));
    Nt = length(S);
    SPEED_T(1:Nt,icase) = S;
  end

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
  NhcST = plot(NHC_TIMES, NHC_WIND, '+k', 'LineWidth',  LineW);
  %title('Maximum 10m Wind Speed');
  xlabel('Time');
  xlim (Xrange);
  set(gca,'xtick', (6:24:54));
  set(gca,'xticklabel', { 'Aug22:12Z' 'Aug23:12Z' 'Aug24:12Z' });
  ylabel('Wind Speed (m s^-^1)');
  ylim(Yrange);

  legend([ NhcST SimST' ], LegText, 'Location', 'SouthEast', 'FontSize', LegendFsize);
%  legend(SimST, LegText, 'Location', 'SouthEast', 'FontSize', LegendFsize);
  
  fprintf('Writing: %s\n', OutFile);
  saveas(FigWind, OutFile);
  close(FigWind);

end
