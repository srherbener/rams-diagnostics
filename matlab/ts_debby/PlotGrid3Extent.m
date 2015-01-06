function [ ] = PlotGrid3Extent(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  OutFile   = sprintf('%s/Grid3Extent.jpg', Pdir);

  %******************************************************************************
  Fsize = 22;
  LegFsize = 15;

  % for Grid3
  LatBounds = [ 7 24 ];
  LonBounds = [ -40 -14 ];

  CoastColor = str2rgb('Black');

  % *** DRY ***
  Fig = figure;
  set (gca, 'FontSize', Fsize);
  m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
  m_coast('color', CoastColor, 'linestyle', '-', 'linewidth', 3);
  m_grid('linestyle','none','box','fancy','tickdir','out');
  title('Grid3');

  fprintf('  Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
