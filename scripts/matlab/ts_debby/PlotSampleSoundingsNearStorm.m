function [ ] = PlotSmapleSoundingsNearStorm()

  OutFile = 'Plots/SampleSoundingsNearStorm.jpg';

  DeltaIndex = 83; % Delta for moving from storm center to region where sample is taken.
                   % 83 grid cells and 3km per cell -> 250km (beginning of outer rainbands)

  T1 = 52;  % time step numbers for selecting out the time period of peak intensity
  T2 = 62;  % for the SD and NSD cases

  fprintf('****************** Creating sample soundings near storm plot ******************\n');

  % ************************** Read in the SAL case ***********************************
  SalTempFile        = 'HDF5/TSD_SAL_DUST/HDF5/tempc-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  SalDewPtFile       = 'HDF5/TSD_SAL_DUST/HDF5/dewptc-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  SalStormCtrFile    = 'HDF5/TSD_SAL_DUST/HDF5/storm_center-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  NonSalTempFile     = 'HDF5/TSD_NONSAL_DUST/HDF5/tempc-TSD_NONSAL_DUST-AS-2006-08-20-120000-g3.h5';
  NonSalDewPtFile    = 'HDF5/TSD_NONSAL_DUST/HDF5/dewptc-TSD_NONSAL_DUST-AS-2006-08-20-120000-g3.h5';
  NonSalStormCtrFile = 'HDF5/TSD_NONSAL_DUST/HDF5/storm_center-TSD_NONSAL_DUST-AS-2006-08-20-120000-g3.h5';

  TempVar = '/tempc';
  DewPtVar = '/dewptc';
  XindVar = '/press_cent_x_index';
  YindVar = '/press_cent_y_index';

  fprintf('  SAL Case:\n');
  fprintf('    Reading: %s (%s)\n', SalTempFile, TempVar);
  fprintf('    Reading: %s (%s)\n', SalDewPtFile, DewPtVar);
  fprintf('    Reading: %s (%s, %s)\n', SalStormCtrFile, XindVar, YindVar);
  fprintf('\n');
  fprintf('  non-SAL Case:\n');
  fprintf('    Reading: %s (%s)\n', NonSalTempFile, TempVar);
  fprintf('    Reading: %s (%s)\n', NonSalDewPtFile, DewPtVar);
  fprintf('    Reading: %s (%s, %s)\n', NonSalStormCtrFile, XindVar, YindVar);
  fprintf('\n');

  % Take the first time step
  X = squeeze(h5read(SalTempFile, '/x_coords'));
  Y = squeeze(h5read(SalTempFile, '/y_coords'));
  Z = squeeze(h5read(SalTempFile, '/z_coords'));
  T = squeeze(h5read(SalTempFile, '/t_coords'));

  Nx = length(X);
  Ny = length(Y);
  Nz = length(Z);
  Nt = length(T);

  % Get the storm centers
  SAL_SCX = squeeze(h5read(SalStormCtrFile, XindVar));
  SAL_SCY = squeeze(h5read(SalStormCtrFile, YindVar));
  NSAL_SCX = squeeze(h5read(NonSalStormCtrFile, XindVar));
  NSAL_SCY = squeeze(h5read(NonSalStormCtrFile, YindVar));

  % Set the z limits
  Z_KM = Z .* 1e-3;
  Z1 = find(Z_KM >=  0, 1, 'first');
  Z2 = find(Z_KM <= 30, 1, 'last');
 
  Z_KM = Z_KM(Z1:Z2);

  % Walk through the selected time steps and grab the soundings from the same
  % relative distance (DeltaIndex) from the storm center.
  SampNx = 10;
  SampNy = 10;
  SampNz = (Z2-Z1) + 1;
  SampNt = (T2-T1) + 1;

  SAL_T   = zeros([ SampNx SampNy SampNz SampNt ]);
  SAL_TD  = zeros([ SampNx SampNy SampNz SampNt ]);
  NSAL_T  = zeros([ SampNx SampNy SampNz SampNt ]);
  NSAL_TD = zeros([ SampNx SampNy SampNz SampNt ]);

  for it = T1:T2
     ix = SAL_SCX(it);               % sample is 250 km north of storm center
     iy = SAL_SCY(it) + DeltaIndex;
     Start = double([ ix iy Z1 it ]);
     Count = double([ SampNx SampNy SampNz 1 ]);
     SAL_T(:,:,:,(it-T1)+1)   = squeeze(h5read(SalTempFile, TempVar, Start, Count));
     SAL_TD(:,:,:,(it-T1)+1)  = squeeze(h5read(SalDewPtFile, DewPtVar, Start, Count));
    
     ix = NSAL_SCX(it);               % sample is 250 km north of storm center
     iy = NSAL_SCY(it) + DeltaIndex;
     Start = double([ ix iy Z1 it ]);
     Count = double([ SampNx SampNy SampNz 1 ]);
     NSAL_T(:,:,:,(it-T1)+1)   = squeeze(h5read(NonSalTempFile, TempVar, Start, Count));
     NSAL_TD(:,:,:,(it-T1)+1)  = squeeze(h5read(NonSalDewPtFile, DewPtVar, Start, Count));
  end

  % Vars are organized as: (x,y,z,t)
  % Construct a spatial mean to produce a single sounding
  SAL_T_AVG   = squeeze(mean(mean(mean(SAL_T, 4), 2), 1));
  SAL_TD_AVG  = squeeze(mean(mean(mean(SAL_TD, 4), 2), 1));
  NSAL_T_AVG  = squeeze(mean(mean(mean(NSAL_T, 4), 2), 1));
  NSAL_TD_AVG = squeeze(mean(mean(mean(NSAL_TD, 4), 2), 1));

  %******************************************************************************
  % single panel plot
  % put SAL lines in solid, NONSAL lines in dashed
  % put T in red, Td in green
  Fig = figure;

  Xlab = 'T (^{\circ}C)';
  Ylab = 'Height (km)';

  Xlimits = [ -100 50 ];
  Xticks = [ -50 0 50 ];

  LegText = {
    'SAL T'
    'SAL T_d'
    'non-SAL T'
    'non-SAL T_d'
    };

  % patch inputs for drawing a rectangle
  SalX = [ Xlimits(1) Xlimits(2) Xlimits(2) Xlimits(1) ]; % full x-axis range
  SalY = [ 1 1 5 5 ]; % height in km
  SalColor = str2rgb('sandybrown');

  SalLabX = Xlimits(1) + 5;
  SalLabY = 2;
  SalFsize = 30;

  Fsize = 25;
  LegFsize = 15;
  LineWidth = 2;

  hold on;
  line(SAL_T_AVG, Z_KM, 'LineWidth', LineWidth, 'Color', str2rgb('red'), 'LineStyle', '-');
  line(SAL_TD_AVG, Z_KM, 'LineWidth', LineWidth, 'Color', str2rgb('green'), 'LineStyle', '-');
  line(NSAL_T_AVG, Z_KM, 'LineWidth', LineWidth, 'Color', str2rgb('red'), 'LineStyle', '--');
  line(NSAL_TD_AVG, Z_KM, 'LineWidth', LineWidth, 'Color', str2rgb('green'), 'LineStyle', '--');

  set (gca, 'FontSize', Fsize);
  xlim(Xlimits);
  set (gca, 'XTick', Xticks);

  title('Near TS Debby');
  xlabel(Xlab);
  ylabel(Ylab);

  legend(LegText, 'FontSize', LegFsize);

  % mark sal location
  patch(SalX, SalY, SalColor, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
  text(SalLabX, SalLabY, 'SAL', 'FontSize', SalFsize);

  fprintf('  Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
