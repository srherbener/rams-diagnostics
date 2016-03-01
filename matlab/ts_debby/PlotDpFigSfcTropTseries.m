function [ ] = PlotDpFigSfcTropTseries()
% PlotDpFigDustTransport 

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  Fsize = 13;
  
  % Time series of upper level Mdhy
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_dust_hydro_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDHY = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MDHY = TS_SAL_MDHY .* 1e-9;  % convert 1e-3 Tg

  T    = squeeze(h5read(InFile, '/t_coords'))./3600 -42;

  % Time series of upper level Mdregen
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ra_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDREGEN = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MDREGEN = TS_SAL_MDREGEN .* 1e-12;  % convert Tg

  % Time series of upper level Md1
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_d1_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MD1 = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MD1 = TS_SAL_MD1 .* 1e-9;  % convert 1e-3 Tg

  % Time series of upper level Md2
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_d2_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MD2 = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MD2 = TS_SAL_MD2 .* 1e-12;  % convert Tg

  % all on one plot
  OutFile = sprintf('%s/DpFigDustSfcTrop.jpg', Pdir);
  Fig = figure;

  Paxes = subplot(4,1,1);
  PlotDpFigTseries(Paxes, T, TS_SAL_MD1, 'a', 'SAL\_AR', 'M_d_1 (10^-^3 Tg)', Fsize, 0, 1, { }, 'none');

  Paxes = subplot(4,1,2);
  PlotDpFigTseries(Paxes, T, TS_SAL_MD2, 'b', 'SAL\_AR', 'M_d_2 (Tg)', Fsize, 0, 1, { }, 'none');

  Paxes = subplot(4,1,3);
  PlotDpFigTseries(Paxes, T, TS_SAL_MDHY, 'c', 'SAL\_AR', 'M_d_h_y (10^-^3 Tg)', Fsize, 0, 1, { }, 'none');

  Paxes = subplot(4,1,4);
  PlotDpFigTseries(Paxes, T, TS_SAL_MDREGEN, 'd', 'SAL\_AR', 'M_d_r_g (Tg)', Fsize, 1, 1, { }, 'none');

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
