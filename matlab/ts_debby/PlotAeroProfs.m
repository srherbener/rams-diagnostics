function [ ] = PlotAeroProfs(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  InFile = 'z.NAMMA_RAMS_LEVELS_aerosols.txt';
  OutFile = sprintf('%s/NAMMA_aerosol_profiles.jpg', Pdir);

  % Read in the aerosol profiles
  % After the text scan, InData will contain:
  %    InData{1}  ->  Z
  %    InData{2}  ->  CCN # conc
  %    InData{3}  ->  Dust 1 (small) # conc
  %    InData{4}  ->  Dust 2 (large) # conc
  %    InData{5}  ->  IN # conc
  %
  %    all # conc in #/cc
  %
  fprintf('****************** Creating aerosol profile plot ******************\n', InFile);
  fprintf('  Reading: %s\n', InFile);
  InFileId = fopen(InFile);
  InData = textscan(InFileId, '%f %f %f %f %f');
  fclose(InFileId);

  Z        = InData{1} ./ 1000; % km
  CCN_CONC = InData{2};
  D1_CONC  = InData{3};
  D2_CONC  = InData{4};
  IN_CONC  = InData{5};

  % Plot from 0 to 10km AGL
  Z1 = find(Z >= 0, 1, 'first');
  Z2 = find(Z <= 10, 1, 'last');

  Z        = Z(Z1:Z2);
  CCN_CONC = CCN_CONC(Z1:Z2);
  D1_CONC  = D1_CONC(Z1:Z2);
  D2_CONC  = D2_CONC(Z1:Z2);
  IN_CONC  = IN_CONC(Z1:Z2);

  % At this point, all quantities are column vectors with 55 entries
  PDATA = [ CCN_CONC D1_CONC D2_CONC ];

  Fig = figure;

  plot(PDATA, Z, 'LineWidth', 2);
  set (gca, 'FontSize', 25);

  T = title('(b)');
  set(T, 'Units', 'Normalized');
  set(T, 'HorizontalAlignment', 'Left');
  Tpos = get(T, 'Position');
  Tpos(1) = 0; % line up with left edge of plot area
  set(T, 'Position', Tpos);

  xlabel('Number Concentration (cm^-^3)');
  ylabel('Height (km)');

  legend({ 'Sulfate' 'Small Dust' 'Large Dust' });
  
  fprintf('  Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
