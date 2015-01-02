function [ ] = PlotSmapleRH(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  OutFileDry   = sprintf('%s/SampleRhDry.jpg', Pdir);
  OutFileMoist = sprintf('%s/SampleRhMoist.jpg', Pdir);

  % Sample: 2nd timestep, 2400m AGL -> Z = 22
  T1 = 2;
  Z1 = 22; % 2400 m AGL

  fprintf('****************** Creating sample soundings plot ******************\n');

  % ************************** Read in the DRY case ***********************************
  RhFile = 'HDF5/relhum-TSD_DRY_NODUST-AS-2006-08-20-120000-g3.h5';
  RhVar = 'relhum';

  fprintf('  DRY Case:\n');
  fprintf('    Reading: %s (%s)\n', RhFile, RhVar);
  fprintf('\n');

  RH_DS  = ncgeodataset(RhFile);
  RH_VAR  = RH_DS.geovariable(RhVar);
  LON_VAR = RH_DS.geovariable('x_coords');
  LAT_VAR = RH_DS.geovariable('y_coords');

  % Vars are organized as: (t,z,y,x)
  RH_DRY  = squeeze(RH_VAR.data(T1,Z1,:,:));
  LON     = squeeze(LON_VAR.data(:));
  LAT     = squeeze(LAT_VAR.data(:));

  % ************************** Read in the MOIST case ***********************************
  RhFile = 'HDF5/relhum-TSD_MOIST_NODUST-AS-2006-08-20-120000-g3.h5';
  RhVar = 'relhum';

  fprintf('  MOIST Case:\n');
  fprintf('    Reading: %s (%s)\n', RhFile, RhVar);
  fprintf('\n');

  RH_DS  = ncgeodataset(RhFile);
  RH_VAR  = RH_DS.geovariable(RhVar);

  % Vars are organized as: (t,z,y,x)
  RH_MOIST  = squeeze(RH_VAR.data(T1,Z1,:,:));

  %******************************************************************************
  Fsize = 25;
  LegFsize = 15;

  % for Grid3
  LatBounds = [ 7 24 ];
  LonBounds = [ -40 -14 ];

  Crange = [ 0 100 ];
  Clevs = [ 10 20 30 40 50 60 70 80 90 ];

  CoastColor = str2rgb('WhiteSmoke');

  % *** DRY ***
  Fig = figure;
  set (gca, 'FontSize', Fsize);
  m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
  m_coast('color', CoastColor, 'linestyle', '-', 'linewidth', 3);
  m_grid('linestyle','none','box','fancy','tickdir','out');
  hold on;
  m_contourf(LON,LAT,RH_DRY,Clevs);
  caxis(Crange);
  colorbar;
  title('DRY');

  fprintf('  Writing: %s\n', OutFileDry);
  saveas(Fig, OutFileDry);
  close(Fig);

  % *** MOIST ***
  Fig = figure;
  set (gca, 'FontSize', Fsize);
  m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
  m_coast('color', CoastColor, 'linestyle', '-', 'linewidth', 3);
  m_grid('linestyle','none','box','fancy','tickdir','out');
  hold on;
  m_contourf(LON,LAT,RH_MOIST,Clevs);
  caxis(Crange);
  colorbar;
  title('MOIST');

  fprintf('  Writing: %s\n', OutFileMoist);
  saveas(Fig, OutFileMoist);
  close(Fig);

end
