function [ ] = PlotDpFigDustMassTseries()
% PlotDpFigDustRemoval 

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  Fsize = 13;
  

  %%%%%%%%%%%%%%%% Volumn integrated dust - surface to tropopause %%%%%%%%%%%%%
  % Time series of Md
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MD = squeeze(h5read(InFile, InVname)); % g
  T    = squeeze(h5read(InFile, '/t_coords'))./3600 - 42;

  % Time series of Mdhy
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_hydro_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDHY = squeeze(h5read(InFile, InVname)); % g

  % Time series of Mdrgn
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_ra_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDRGN = squeeze(h5read(InFile, InVname)); % g

  % Time series of Mdsfc
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_sfc_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDSFC = squeeze(h5read(InFile, InVname)); % g

  % Time series of Mdadv
  InFile = 'DIAGS/residual_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_residual_rmvd_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDADV = squeeze(h5read(InFile, InVname)); % g

  TS_DUST_AL = [ TS_SAL_MD TS_SAL_MDHY TS_SAL_MDRGN TS_SAL_MDADV TS_SAL_MDSFC ] .* 1e-12; % convert to Tg
  LegTextAl = { 'M_d' 'M_d_h_y' 'M_d_r_g_n' 'M_d_a_d_v' 'M_d_s_f_c' };
  LegLocAl = 'NorthEastOutside';

  fprintf('\n');

  %%%%%%%%%%%%%%%% Volumn integrated dust - 7 km AGL to tropopause %%%%%%%%%%%%%
  % Time series of upper level Md
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MD = squeeze(h5read(InFile, InVname)); % g

  % Time series of upper level Mdhy
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_hydro_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDHY = squeeze(h5read(InFile, InVname)); % g

  % Time series of upper level Mdrgn
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_ra_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDRGN = squeeze(h5read(InFile, InVname)); % g

  % Time series of upper level Mdadv
  InFile = 'DIAGS/residual_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_residual_rmvd_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDADV = squeeze(h5read(InFile, InVname)); % g

  TS_DUST_UL = [ TS_SAL_MD TS_SAL_MDHY TS_SAL_MDRGN TS_SAL_MDADV ] .* 1e-9; % convert to 1e-3 Tg
  LegTextUl = { 'M_d' 'M_d_h_y' 'M_d_r_g_n' 'M_d_a_d_v' };
  LegLocUl = 'NorthEastOutside';

  fprintf('\n');

  % plot, 4 panels (2x2)
  OutFile = sprintf('%s/DpFig3_DustMassTseries.jpg', Pdir);
  Lcolors = { 'blue' 'orangered' 'goldenrod' 'purple' 'forestgreen' };
  
  Fig = figure;

  Paxes = subplot(2,1,1);
  PlotDpFigTseries(Paxes, T, TS_DUST_AL, Lcolors, 'a', 'All Levels', 'Mass (Tg)', 'log', [ 1e-4 1e1 ], Fsize, 0, 1, LegTextAl, LegLocAl);

  Paxes = subplot(2,1,2);
  PlotDpFigTseries(Paxes, T, TS_DUST_UL, Lcolors, 'b', 'Upper Levels', 'Mass (10^-^3 Tg)', 'linear', [ -0.5 3.5 ], Fsize, 1, 1, LegTextUl, LegLocUl);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
