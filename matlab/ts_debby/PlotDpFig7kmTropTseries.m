function [ ] = PlotDpFig7kmTropTseries()
% PlotDpFigDustTransport 

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  Fsize = 13;
  
  % Time series of upper level Mdhy
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_dust_hydro_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDHY = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MDHY = TS_SAL_MDHY .* 1e-9;  % convert 1e-3 Tg

  T    = squeeze(h5read(InFile, '/t_coords'))./3600 -42;

  % Time series of upper level Mdregen
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ra_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDREGEN = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MDREGEN = TS_SAL_MDREGEN .* 1e-9;  % convert 1e-3 Tg

  % Time series of upper level Md1
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_d1_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MD1 = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MD1 = TS_SAL_MD1 .* 1e-9;  % convert 1e-3 Tg

  % Time series of upper level Md2
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_d2_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MD2 = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MD2 = TS_SAL_MD2 .* 1e-9;  % convert 1e-3 Tg

  % Time series of upper level Mdadv (residual)
  InFile = 'DIAGS/residual_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_residual_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDRES = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MDRES = TS_SAL_MDRES .* 1e-9;  % convert 1e-3 Tg

  % Combine all of the time series into one panel
  TS_SAL_DUST = [ TS_SAL_MD1 TS_SAL_MD2 TS_SAL_MDHY TS_SAL_MDREGEN ];
  LegText = { 'M_d_1' 'M_d_2' 'M_d_h_y' 'M_d_r_g' };

  % dust in hydrometeors
  OutFile = sprintf('%s/DpFig7kmTropMdhy.jpg', Pdir);
  Fig = figure;

  Paxes = subplot(4,1,1);
  PlotDpFigTseries(Paxes, T, TS_SAL_MDHY, 'a', 'SAL\_AR (Upper Levels)', 'M_d_h_y (10^-^3 Tg)', Fsize, 1, 1, { }, 'none');

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

  % regenerated dust
  OutFile = sprintf('%s/DpFig7kmTropMdrg.jpg', Pdir);
  Fig = figure;

  Paxes = subplot(4,1,1);
  PlotDpFigTseries(Paxes, T, TS_SAL_MDREGEN, 'b', 'SAL\_AR (Upper Levels)', 'M_d_r_g (10^-^3 Tg)', Fsize, 1, 1, { }, 'none');

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

  % regenerated dust
  OutFile = sprintf('%s/DpFig7kmTropMdadv.jpg', Pdir);
  Fig = figure;

  Paxes = subplot(4,1,1);
  PlotDpFigTseries(Paxes, T, TS_SAL_MDRES, 'c', 'SAL\_AR (Upper Levels)', 'M_d_a_d_v (10^-^3 Tg)', Fsize, 1, 1, { }, 'none');

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

  % all on one plot
  OutFile = sprintf('%s/DpFigDustHlev.jpg', Pdir);
  Fig = figure;

  % Total aero mass time series
  PlotDpFigTseries(gca, T, TS_SAL_DUST, 'e', 'SAL\_AR (Upper Levels)', 'Total Mass (10^-^3 Tg)', Fsize, 1, 1, LegText, 'NorthEast');
  
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
