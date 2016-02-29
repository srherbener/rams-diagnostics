function [ ] = PlotDpFigHlevDust()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

  % Time series of Md, at high levels
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_aero_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MD = squeeze(h5read(InFile, InVname)); % g

  T = squeeze(h5read(InFile, '/t_coords'))./3600 - 42; % sim time in hours

  % Time series of Md1, at high levels
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_d1_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MD1 = squeeze(h5read(InFile, InVname)); % g

  % Time series of Md2, at high levels
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_d2_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MD2 = squeeze(h5read(InFile, InVname)); % g

  % Time series of Mdhy, at high levels
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_dust_hydro_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDHY = squeeze(h5read(InFile, InVname)); % g

  % Time series of Mdra, at high levels
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ra_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDRA = squeeze(h5read(InFile, InVname)); % g

  % Combine all of the time series into one panel
  TS_SAL_DUST = [ TS_SAL_MD1 TS_SAL_MD2 TS_SAL_MDHY TS_SAL_MDRA ];
  TS_SAL_DUST = TS_SAL_DUST .* 1e-9; % convert to 1e-3 Tg

  LegText = { 'M_d_1' 'M_d_2' 'M_d_h_y' 'M_d_r_e_g_e_n' };

  % plot
  Fig = figure;

  % Total aero mass time series
  PlotDpFigTseries(gca, T, TS_SAL_DUST, 'e', 'SAL\_AR (Upper Levels)', 'Total Mass (10^-^3 Tg)', Fsize, 1, 1, LegText, 'NorthEast');
  
  OutFile = sprintf('%s/DpFigDustHlev.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end

