function [ ] = PlotSmapleSoundings(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  OutFile = sprintf('%s/SampleSoundings.jpg', Pdir);

  % Sample region, in the SAL location
  X1 = 290;
  X2 = 300;

  Y1 = 385;
  Y2 = 395;

  fprintf('****************** Creating sample soundings plot ******************\n');

  % ************************** Read in the SAL case ***********************************
  SalTempFile     = 'HDF5/tempc_init-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  SalDewPtFile    = 'HDF5/dewptc_init-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  NonSalTempFile  = 'HDF5/tempc_init-TSD_NONSAL_DUST-AS-2006-08-20-120000-g3.h5';
  NonSalDewPtFile = 'HDF5/dewptc_init-TSD_NONSAL_DUST-AS-2006-08-20-120000-g3.h5';

  TempVar = '/tempc';
  DewPtVar = '/dewptc';

  fprintf('  SAL Case:\n');
  fprintf('    Reading: %s (%s)\n', SalTempFile, TempVar);
  fprintf('    Reading: %s (%s)\n', SalDewPtFile, DewPtVar);
  fprintf('\n');
  fprintf('  non-SAL Case:\n');
  fprintf('    Reading: %s (%s)\n', NonSalTempFile, TempVar);
  fprintf('    Reading: %s (%s)\n', NonSalDewPtFile, DewPtVar);
  fprintf('\n');

  SAL_T   = squeeze(h5read(SalTempFile, TempVar));
  SAL_TD  = squeeze(h5read(SalDewPtFile, DewPtVar));
  NSAL_T  = squeeze(h5read(NonSalTempFile, TempVar));
  NSAL_TD = squeeze(h5read(NonSalDewPtFile, DewPtVar));

  % Grab the z values, and set the z limits
  Z = squeeze(h5read(SalTempFile, '/z_coords')) ./1000;  % km

  Z1 = find(Z >= 0, 1, 'first');
  Z2 = find(Z <= 30, 1, 'last');
 
  Z = Z(Z1:Z2);

  % Vars are organized as: (x,y,z)
  % Construct a spatial mean to produce a single sounding
  SAL_T_AVG   = squeeze(mean(mean(SAL_T(X1:X2,Y1:Y2,Z1:Z2), 2), 1));
  SAL_TD_AVG  = squeeze(mean(mean(SAL_TD(X1:X2,Y1:Y2,Z1:Z2), 2), 1));
  NSAL_T_AVG  = squeeze(mean(mean(NSAL_T(X1:X2,Y1:Y2,Z1:Z2), 2), 1));
  NSAL_TD_AVG = squeeze(mean(mean(NSAL_TD(X1:X2,Y1:Y2,Z1:Z2), 2), 1));

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

  plot(SAL_T_AVG, Z, 'LineWidth', LineWidth, 'Color', str2rgb('red'), 'LineStyle', '-');
  hold on
  plot(SAL_TD_AVG, Z, 'LineWidth', LineWidth, 'Color', str2rgb('green'), 'LineStyle', '-');
  plot(NSAL_T_AVG, Z, 'LineWidth', LineWidth, 'Color', str2rgb('red'), 'LineStyle', '--');
  plot(NSAL_TD_AVG, Z, 'LineWidth', LineWidth, 'Color', str2rgb('green'), 'LineStyle', '--');

  set (gca, 'FontSize', Fsize);
  xlim(Xlimits);
  set (gca, 'XTick', Xticks);

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
