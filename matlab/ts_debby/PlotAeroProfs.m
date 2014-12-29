function [ ] = PlotAeroProfs(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  DustInFile = 'z.NAMMA_RAMS_LEVELS_aerosols.txt';
  NoDustInFile = 'z.NO_DUST_aerosols.txt';
  OutFile = sprintf('%s/NAMMA_aerosol_profiles.jpg', Pdir);

  fprintf('****************** Creating aerosol profile plot ******************\n');

  %*************************************************************************
  % Read in and assemble the aerosol profiles for the CLEAN case
  % After the text scan, InData will contain:
  %    InData{1}  ->  'CCN-init'
  %    InData{2}  ->  level number
  %    InData{3}  ->  Z
  %    InData{4}  ->  CCN # conc in #/mg
  %    InData{5}  ->  CCN # conc in #/cc
  %
  % This case has both dust modes set to zero
  %
  fprintf('  Reading: %s\n', NoDustInFile);
  InFileId = fopen(NoDustInFile);
  InData = textscan(InFileId, '%s %d %f %f %f');
  fclose(InFileId);

  Z        = InData{3} ./ 1000; % km
  CCN_CONC = InData{5}; % #/cc

  % Plot from 0 to 10km AGL
  Z1 = find(Z >= 0, 1, 'first');
  Z2 = find(Z <= 10, 1, 'last');

  Z        = Z(Z1:Z2);
  Nz       = length(Z);
  CCN_CONC = CCN_CONC(Z1:Z2);
  D1_CONC  = zeros([ Nz 1 ]);
  D2_CONC  = zeros([ Nz 1 ]);

  % At this point, all quantities are column vectors
  CLEAN_PDATA = [ CCN_CONC D1_CONC D2_CONC ];


  %*************************************************************************
  % Read in and assemble the aerosol profiles for the DUST case
  % After the text scan, InData will contain:
  %    InData{1}  ->  Z
  %    InData{2}  ->  CCN # conc
  %    InData{3}  ->  Dust 1 (small) # conc
  %    InData{4}  ->  Dust 2 (large) # conc
  %    InData{5}  ->  IN # conc
  %
  %    all # conc in #/cc
  %
  fprintf('****************** Creating aerosol profile plot ******************\n');
  fprintf('  Reading: %s\n', DustInFile);
  InFileId = fopen(DustInFile);
  InData = textscan(InFileId, '%f %f %f %f %f');
  fclose(InFileId);

  Z        = InData{1} ./ 1000; % km
  CCN_CONC = InData{2};
  D1_CONC  = InData{3};
  D2_CONC  = InData{4};

  % Plot from 0 to 10km AGL
  Z1 = find(Z >= 0, 1, 'first');
  Z2 = find(Z <= 10, 1, 'last');

  Z        = Z(Z1:Z2);
  CCN_CONC = CCN_CONC(Z1:Z2);
  D1_CONC  = D1_CONC(Z1:Z2);
  D2_CONC  = D2_CONC(Z1:Z2);

  % At this point, all quantities are column vectors
  DUST_PDATA = [ CCN_CONC D1_CONC D2_CONC ];


  %******************************************************************************
  % Two panel plot
  Fig = figure;

  Xlab = 'N_a (cm^-^3)';
  Ylab = 'Height (km)';

  Xlimits = [ -10 300 ];
  Xticks = [ 0 150 300 ];

  Fsize = 25;
  LegFsize = 15;
  LineWidth = 2;

  % *** DUST ***
  subplot(1,2,1);
  plot(DUST_PDATA, Z, 'LineWidth', LineWidth);
  set (gca, 'FontSize', Fsize);
  xlim(Xlimits);
  set (gca, 'XTick', Xticks);

  title('DUSTY');
  xlabel(Xlab);
  ylabel(Ylab);


  % *** CLEAN ***
  subplot(1,2,2);
  plot(CLEAN_PDATA, Z, 'LineWidth', LineWidth);
  set (gca, 'FontSize', Fsize);
  xlim(Xlimits);
  set (gca, 'XTick', Xticks);

  % shut off labels on y-axis, will share with left panel
  title('CLEAN');
  xlabel(Xlab);
  set (gca, 'YTickLabel', []);

  legend({ 'Sulfate' 'Small Dust' 'Large Dust' }, 'FontSize', LegFsize);
  
  fprintf('  Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
