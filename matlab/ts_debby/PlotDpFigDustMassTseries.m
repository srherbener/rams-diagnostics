function [ ] = PlotDpFigDustMassTseries()
% PlotDpFigDustRemoval 

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  Fsize = 13;

  % color scheme
  ColorMd     = 'blue';
  ColorMdhy   = 'orangered';
  ColorMdrgn  = 'goldenrod';
  ColorMdadv  = 'purple';
  ColorMdsfc  = 'forestgreen';
  
  Lcolors = { 'blue' 'orangered' 'goldenrod' 'purple' 'forestgreen' };

  %%%%%%%%%%%%%%%% Volumn integrated dust - surface to tropopause %%%%%%%%%%%%%
  % Time series of Md
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MD = squeeze(h5read(InFile, InVname)); % g
  T    = squeeze(h5read(InFile, '/t_coords'))./3600 - 42;

  % normalization of all quantities are relative to the total initial mass at all levels
  NORM_FAC = 1 / TS_SAL_MD(1);

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

  % Time series of Md rate
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_total_mass_rate';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MD_RATE = squeeze(h5read(InFile, InVname)); % g/s

  % Time series of Mdadv rate
  InFile = 'DIAGS/residual_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_residual_rmvd_total_mass_rate';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDADV_RATE = squeeze(h5read(InFile, InVname)); % g/s

  % Rates have Nt-1 points
  T_RATE = T(1:end-1);

  % normalize to initial dust amount
  TS_DUST_AL = [ TS_SAL_MD TS_SAL_MDHY TS_SAL_MDRGN TS_SAL_MDSFC ] .* NORM_FAC;
  LegTextAl = { 'M_D' 'M_D_H_Y' 'M_D_R_G_N' 'M_D_S_F_C' };
  LegLocAl = 'NorthEastOutside';
  LcolorsAl = { ColorMd ColorMdhy ColorMdrgn ColorMdsfc };

  % rates
  % Note Mdadv rate is a removal term meaning that a positive rate is removal from the SAL_AR
  TS_DUST_AL_RATE = [ TS_SAL_MD_RATE TS_SAL_MDADV_RATE ] .* 1e-12 .* 3600 .* 24; % Tg/d
  LegTextAlRate = { 'M_D' 'M_D_A_D_V' };
  LegLocAlRate = 'NorthEastOutside';
  LcolorsAlRate = { ColorMd ColorMdadv };

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

  % Time series of upper level Md rate
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_total_mass_hlev_rate';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MD_RATE = squeeze(h5read(InFile, InVname)); % g/s

  % Time series of upper level Mdadv
  InFile = 'DIAGS/residual_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_residual_rmvd_total_mass_hlev_rate';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDADV_RATE = squeeze(h5read(InFile, InVname)); % g/s

  % normalize to initial dust amount, ie [ TS_SAL_MD(1) ]
  TS_DUST_UL = [ TS_SAL_MD TS_SAL_MDHY TS_SAL_MDRGN ] .* NORM_FAC;
  LegTextUl = { 'M_D' 'M_D_H_Y' 'M_D_R_G_N' };
  LegLocUl = 'NorthEastOutside';
  LcolorsUl = { ColorMd ColorMdhy ColorMdrgn };

  % rates
  TS_DUST_UL_RATE = [ TS_SAL_MD_RATE TS_SAL_MDADV_RATE ] .* 1e-9 .* 3600 .* 24; % convert to 1e-3 Tg/d
  LegTextUlRate = { 'M_D' 'M_D_A_D_V' };
  LegLocUlRate = 'NorthEastOutside';
  LcolorsUlRate = { ColorMd ColorMdadv };

  fprintf('\n');

  %%%%%%%%%%%%%%%% precip rate %%%%%%%%%%%%%
  % Time series of precip rate
  InFile = 'DIAGS/hist_meas_ts_pcprate_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_pcprate_ts';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_PR = squeeze(h5read(InFile, InVname)); % mm/h


  % plot, 4 panels (2x2)
  OutFile = sprintf('%s/DpFig3_DustMassTseries.jpg', Pdir);
  Lcolors = { 'blue' 'orangered' 'goldenrod' 'purple' 'forestgreen' };
  
  Fig = figure;

  Paxes = subplot(4,1,1);
  PlotDpFigTseries(Paxes, T_RATE, TS_DUST_AL_RATE, LcolorsAlRate, 'a', 'All Levels', '\DeltaM (Tg d^-^1)', 'linear', [ -2 2 ], Fsize, 0, 1, LegTextAlRate, LegLocAlRate);
  % mark y = zero
  hold on
  line([ 0 60 ], [ 0 0 ], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
  hold off

  Paxes = subplot(4,1,2);
  % adjust size so that x-axis lines up with the other panels
  Ploc = get(Paxes, 'Position');
  Ploc(3) = Ploc(3) * 0.805;
  Ploc = set(Paxes, 'Position', Ploc);
  PlotDpFigTseries(Paxes, T, TS_SAL_PR, { 'blue' }, 'b', '', 'PR (mm h^-^1)', 'linear', [ 0 0.45 ], Fsize, 0, 1, { '' }, 'none');

  Paxes = subplot(4,1,3);
  PlotDpFigTseries(Paxes, T, TS_DUST_AL, LcolorsAl, 'c', 'All Levels', 'Norm. Mass', 'log', [ 1e-4 2 ], Fsize, 0, 1, LegTextAl, LegLocAl);

  Paxes = subplot(4,1,4);
  PlotDpFigTseries(Paxes, T, TS_DUST_UL, LcolorsUl, 'd', 'Upper Levels', 'Norm. Mass', 'log', [ 1e-6  1e-2 ], Fsize, 1, 1, LegTextUl, LegLocUl);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
