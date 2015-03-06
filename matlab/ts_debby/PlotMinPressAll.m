function [ ] = PlotMinPressAll(ConfigFile)

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

%  InFprefix = 'min_sea_press';
%  InVname = '/min_sea_press';
   InFprefix = 'min_press25m';
   InVname = '/min_press25m';

  for icase = 1:Nc
    Hfile = sprintf('%s/%s_%s.h5', Tdir, InFprefix, Cases{icase});
    InFiles{icase} = Hfile;
  end

  % NHC Best Track (every six hours) data
  % time step 1 from the simulation is where the NHC data starts
  %
  NHC_PRESS = [ 1007 1007 1005 1003 1002 1001 1001 1000  999 1000 1000 ];
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

  PRESS = nan([ MaxNt Nc ]);

  % read in the Vt data
  for icase = 1:Nc
    Hfile = InFiles{icase};
    fprintf('Reading: %s (%s)\n', Hfile, InVname);

    % Read in time series, P will be a column vector
    P = squeeze(h5read(Hfile, InVname));
    Nt = length(P);

    PRESS(1:Nt,icase) = P;
  end

  % plot
  OutFile = sprintf('%s/TsDebbyPressAll.jpg', Pdir);

  Fsize = 22;
  LineW = 2;

  LegendFsize = 15;

  Xrange = [  -6   66 ];
  Yrange = [ 995 1010 ];
  
  FigWind = figure;
  SimST = plot(TIMES, PRESS, 'LineWidth', LineW);
  set(gca, 'FontSize', Fsize);
  hold on;
  NhcST = plot(NHC_TIMES, NHC_PRESS, '+k', 'LineWidth',  LineW);
  %title('Minimum SLP');
  xlabel('Time');
  xlim (Xrange);
  set(gca,'xtick', (6:24:54));
  set(gca,'xticklabel', { 'Aug22:12Z' 'Aug23:12Z' 'Aug24:12Z' });
  ylabel('Pressure (mb)');
  ylim(Yrange);

  legend([ NhcST SimST' ], LegText, 'Location', 'SouthEast', 'FontSize', LegendFsize);
%  legend(SimST, LegText, 'Location', 'SouthEast', 'FontSize', LegendFsize);
  
  fprintf('Writing: %s\n', OutFile);
  saveas(FigWind, OutFile);
  close(FigWind);

end
