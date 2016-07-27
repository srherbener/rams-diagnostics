function [ ] = PlotSmapleSoundings()

  OutFile = 'Plots/SampleSoundings.jpg';

  % Sample region, in the SAL location
  X1 = 290;
  X2 = 300;

  Y1 = 385;
  Y2 = 395;

  fprintf('****************** Creating sample soundings plot ******************\n');

  % ************************** Read in the SAL case ***********************************
  SalTempFile     = 'HDF5/TSD_SAL_DUST/HDF5/tempc-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  SalDewPtFile    = 'HDF5/TSD_SAL_DUST/HDF5/dewptc-TSD_SAL_DUST-AS-2006-08-20-120000-g3.h5';
  NonSalTempFile  = 'HDF5/TSD_NONSAL_DUST/HDF5/tempc-TSD_NONSAL_DUST-AS-2006-08-20-120000-g3.h5';
  NonSalDewPtFile = 'HDF5/TSD_NONSAL_DUST/HDF5/dewptc-TSD_NONSAL_DUST-AS-2006-08-20-120000-g3.h5';

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

  % Take the first time step
  X = squeeze(h5read(SalTempFile, '/x_coords'));
  Y = squeeze(h5read(SalTempFile, '/y_coords'));
  Z = squeeze(h5read(SalTempFile, '/z_coords'));
  T = squeeze(h5read(SalTempFile, '/t_coords'));

  Nx = length(X);
  Ny = length(Y);
  Nz = length(Z);
  Nt = length(T);
  
  Start = [ 1 1 1 2 ]; % use time step 2 since non-SAL env is not applied
                       % until after time step 1 file is written out.
  Count = [ Nx Ny Nz 1 ];

  SAL_T   = squeeze(h5read(SalTempFile, TempVar, Start, Count));
  SAL_TD  = squeeze(h5read(SalDewPtFile, DewPtVar, Start, Count));
  NSAL_T  = squeeze(h5read(NonSalTempFile, TempVar, Start, Count));
  NSAL_TD = squeeze(h5read(NonSalDewPtFile, DewPtVar, Start, Count));

  % Set the z limits
  Z_KM = Z .* 1e-3;
  Z1 = find(Z_KM >=  0, 1, 'first');
  Z2 = find(Z_KM <= 30, 1, 'last');
 
  Z_KM = Z_KM(Z1:Z2);

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

  hold on;
  line(SAL_T_AVG, Z_KM, 'LineWidth', LineWidth, 'Color', str2rgb('red'), 'LineStyle', '-');
  line(SAL_TD_AVG, Z_KM, 'LineWidth', LineWidth, 'Color', str2rgb('green'), 'LineStyle', '-');
  line(NSAL_T_AVG, Z_KM, 'LineWidth', LineWidth, 'Color', str2rgb('red'), 'LineStyle', '--');
  line(NSAL_TD_AVG, Z_KM, 'LineWidth', LineWidth, 'Color', str2rgb('green'), 'LineStyle', '--');

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
