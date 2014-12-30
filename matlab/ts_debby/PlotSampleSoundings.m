function [ ] = PlotSmapleSoundings(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  OutFile = sprintf('%s/SampleSoundings.jpg', Pdir);

  % Sample region, in the SAL location
  X1 = 625;
  X2 = 635;
  %X1 = 630;
  %X2 = 630;

  Y1 = 535;
  Y2 = 545;
  %Y1 = 540;
  %Y2 = 540;

  %T1 = 30;
  %T2 = 35;
  T1 = 1;
  T2 = 39;

  fprintf('****************** Creating sample soundings plot ******************\n');

  % ************************** Read in the DRY case ***********************************
  TempFile = 'HDF5/tempc-TSD_DRY_DUST-AS-2006-08-20-120000-g3.h5';
  DewPtFile = 'HDF5/dewptc-TSD_DRY_DUST-AS-2006-08-20-120000-g3.h5';

  TempVar = 'tempc';
  DewPtVar = 'dewptc';

  fprintf('  DRY Case:\n');
  fprintf('    Reading: %s (%s)\n', TempFile, TempVar);
  fprintf('    Reading: %s (%s)\n', DewPtFile, DewPtVar);
  fprintf('\n');

  T_DS  = ncgeodataset(TempFile);
  TD_DS = ncgeodataset(DewPtFile);

  T_VAR  = T_DS.geovariable(TempVar);
  TD_VAR = TD_DS.geovariable(DewPtVar);

  Z_VAR = T_DS.geovariable('z_coords');

  % Grab the z values, and set the z limits
  Z = Z_VAR.data(:) ./1000;  % km

  Z1 = find(Z >= 0, 1, 'first');
  Z2 = find(Z <= 30, 1, 'last');
 
  Z = Z(Z1:Z2);

  % Vars are organized as: (t,z,y,x)
  TEMP  = T_VAR.data(1:end,Z1:Z2,Y1:Y2,X1:X2);
  DEWPT = TD_VAR.data(1:end,Z1:Z2,Y1:Y2,X1:X2);

  % Construct a spatial and temporal mean to produce a single sounding
  TEMP_AVG = squeeze(mean(mean(mean(TEMP, 4), 3), 1))';
  DEWPT_AVG = squeeze(mean(mean(mean(DEWPT, 4), 3), 1))';

  DRY_PDATA = [ TEMP_AVG DEWPT_AVG ];

  % ************************** Read in the MOIST case ***********************************
  TempFile = 'HDF5/tempc-TSD_MOIST_NODUST-AS-2006-08-20-120000-g3.h5';
  DewPtFile = 'HDF5/dewptc-TSD_MOIST_NODUST-AS-2006-08-20-120000-g3.h5';

  TempVar = 'tempc';
  DewPtVar = 'dewptc';

  fprintf('  MOIST Case:\n');
  fprintf('    Reading: %s (%s)\n', TempFile, TempVar);
  fprintf('    Reading: %s (%s)\n', DewPtFile, DewPtVar);
  fprintf('\n');

  T_DS  = ncgeodataset(TempFile);
  TD_DS = ncgeodataset(DewPtFile);

  T_VAR  = T_DS.geovariable(TempVar);
  TD_VAR = TD_DS.geovariable(DewPtVar);

  % Vars are organized as: (t,z,y,x)
  TEMP  = T_VAR.data(1:end,Z1:Z2,Y1:Y2,X1:X2);
  DEWPT = TD_VAR.data(1:end,Z1:Z2,Y1:Y2,X1:X2);

  % Construct a spatial and temporal mean to produce a single sounding
  TEMP_AVG = squeeze(mean(mean(mean(TEMP, 4), 3), 1))';
  DEWPT_AVG = squeeze(mean(mean(mean(DEWPT, 4), 3), 1))';

  MOIST_PDATA = [ TEMP_AVG DEWPT_AVG ];


  %******************************************************************************
  % Two panel plot
  Fig = figure;

  Xlab = 'T (^{\circ}C)';
  Ylab = 'Height (km)';

  Xlimits = [ -100 50 ];
  Xticks = [ -50 0 50 ];

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

  % *** DRY ***
  subplot(1,2,1);
  plot(DRY_PDATA, Z, 'LineWidth', LineWidth);
  set (gca, 'FontSize', Fsize);
  xlim(Xlimits);
  set (gca, 'XTick', Xticks);

  title('DRY');
  xlabel(Xlab);
  ylabel(Ylab);

  % mark sal location
  patch(SalX, SalY, SalColor, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
  text(SalLabX, SalLabY, 'SAL', 'FontSize', SalFsize);

  % *** MOIST ***
  subplot(1,2,2);
  plot(MOIST_PDATA, Z, 'LineWidth', LineWidth);
  set (gca, 'FontSize', Fsize);
  xlim(Xlimits);
  set (gca, 'XTick', Xticks);

  % shut off labels on y-axis, will share with left panel
  title('MOIST');
  xlabel(Xlab);
  set (gca, 'YTickLabel', []);

  legend({ 'Temp.' 'Dew Pt.' }, 'FontSize', LegFsize);

  % mark sal location
  patch(SalX, SalY, SalColor, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
  text(SalLabX, SalLabY, 'SAL', 'FontSize', SalFsize);
  

  fprintf('  Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
