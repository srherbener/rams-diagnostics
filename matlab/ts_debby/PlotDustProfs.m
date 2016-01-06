function [ ] = PlotDustProfs()


  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  DustInFile = 'z.NAMMA_RAMS_LEVELS_aerosols.txt';
  OutFile = sprintf('%s/NAMMA_dust_profiles.jpg', Pdir);

  Zmin =  0; %km
  Zmax = 15; %km

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
  D1_CONC  = InData{3};
  D2_CONC  = InData{4};

  % Plot from 0 to 20km AGL
  Z1 = find(Z >= Zmin, 1, 'first');
  Z2 = find(Z <= Zmax, 1, 'last');

  Z        = Z(Z1:Z2);
  D1_CONC  = D1_CONC(Z1:Z2);
  D2_CONC  = D2_CONC(Z1:Z2);

  % At this point, all quantities are column vectors
  DUST_PDATA = [ D1_CONC D2_CONC ];


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

  title('DUST');
  xlabel(Xlab);
  ylabel(Ylab);

  legend({ 'Small Dust' 'Large Dust' }, 'FontSize', LegFsize);

  fprintf('  Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
