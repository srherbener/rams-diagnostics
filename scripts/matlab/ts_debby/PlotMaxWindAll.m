function [ ] = PlotMaxWindAll()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  Cases = {
%   'TSD_DRY_DUST'
%   'TSD_DRY_NODUST'
%   'TSD_MOIST_DUST'
%   'TSD_MOIST_NODUST'

   'TSD_SAL_DUST'
%   'TSD_SAL_NODUST'
%   'TSD_NONSAL_DUST'
%   'TSD_NONSAL_NODUST'
   };
  Nc = length(Cases);

  LegText = {
   'NHC Best Track'
   'TAFB Obs'
%   'DD'
%   'DN'
%   'MD'
%   'MN'

   'SAL\_DUST'
%   'SAL\_NODUST'
%   'NONSAL\_DUST'
%   'NONSAL\_NODUST'
   };

   InFprefix = 'DIAGS/storm_meas';
   InVname = '/max_wind';

  for icase = 1:Nc
    InFiles{icase} = sprintf('%s_%s.h5', InFprefix, Cases{icase});
  end

  % NHC Best Track (every six hours) data
  % time step 1 from the simulation is where the NHC data starts
  % NHC wind is in mph, multiply by 0.447 to convert to m/s
  %
  % TAFB - Tropical Analysis and Forecast Branch
  % TAFB wind is in kts, multiply by 0.514 to convert to m/s
  %
  % time == 0 is Aug 22, 2006 at 06Z
  %
  NHC_WIND  = [ 35 35 35 40 45 50 50 50 50 50 50 ] * 0.447;
  TAFB_WIND  = [ 30 30 30 35 35 35 35 35 35 35 35 ] * 0.514;
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

    % Read in time series
    %  If 1D var, (t), then done at this point
    %  If 2D var, (r,t), find max along R dimension
    %  If 3D var, (r,z,t), select k=2 level, find max along R dimension
    S = squeeze(h5read(Hfile, InVname));
%    S = squeeze(max(S, [], 1));
%    S = squeeze(max(squeeze(S(:,2,:)), [], 1));
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
  TafbST = plot(NHC_TIMES, TAFB_WIND, 'ok', 'LineWidth',  LineW);
  %title('Maximum 10m Wind Speed');
  xlabel('Time');
  xlim (Xrange);
  set(gca,'xtick', (6:24:54));
  set(gca,'xticklabel', { 'Aug22:12Z' 'Aug23:12Z' 'Aug24:12Z' });
  ylabel('Wind Speed (m s^-^1)');
  ylim(Yrange);

  legend([ NhcST TafbST SimST' ], LegText, 'Location', 'SouthEast', 'FontSize', LegendFsize);
%  legend(SimST, LegText, 'Location', 'SouthEast', 'FontSize', LegendFsize);
  
  fprintf('Writing: %s\n', OutFile);
  saveas(FigWind, OutFile);
  close(FigWind);

end
