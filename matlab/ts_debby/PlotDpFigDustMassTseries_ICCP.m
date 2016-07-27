function [ ] = PlotDpFigDustMassTseries_ICCP()
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
  InVname = '/sal_ar_residual_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDADV = squeeze(h5read(InFile, InVname)); % g

  % Absolute values of Md and Mdadv
  % NanSeries is there to suppress plotting of the lines, while allowing the legend
  % to take on the colors of the lines in the other panels.
  NanSeries = TS_SAL_MD .* nan;
  TS_DUST_AL_MD = [ TS_SAL_MD NanSeries NanSeries NanSeries NanSeries ] .* 1e-12;  % Tg
  LegTextAlMdMdadv = { 'M_D' 'M_D_A_D_V' 'M_D_H_Y' 'M_D_R_G_N' 'M_D_S_F_C' };
  LegLocAlMdMdadv  = 'none';
  LcolorsAlMdMdadv = { ColorMd ColorMdadv ColorMdhy ColorMdrgn ColorMdsfc };

  % normalize to initial dust amount
  TS_DUST_AL = [ TS_SAL_MDADV TS_SAL_MDHY TS_SAL_MDRGN TS_SAL_MDSFC ] .* NORM_FAC;
  LegTextAl = { 'M_D_A_D_V' 'M_D_H_Y' 'M_D_R_G_N' 'M_D_S_F_C' };
  LegLocAl = 'none';
  LcolorsAl = { ColorMdadv ColorMdhy ColorMdrgn ColorMdsfc };

  fprintf('\n');

  %%%%%%%%%%%%%%%% Volumn integrated dust - 4 km AGL to 7 km AGL %%%%%%%%%%%%%
  % Time series of mid level Mdhy
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_hydro_total_mass_mlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDHY = squeeze(h5read(InFile, InVname)); % g

  % Time series of mid level Mdrgn
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_ra_total_mass_mlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDRGN = squeeze(h5read(InFile, InVname)); % g

  % normalize to initial dust amount, ie [ TS_SAL_MD(1) ]
  TS_DUST_ML = [ TS_SAL_MDHY TS_SAL_MDRGN ] .* NORM_FAC;
  LegTextMl = { 'M_D_H_Y' 'M_D_R_G_N' };
  LegLocMl = 'none';
  LcolorsMl = { ColorMdhy ColorMdrgn };

  fprintf('\n');

  %%%%%%%%%%%%%%%% Volumn integrated dust - 7 km AGL to tropopause %%%%%%%%%%%%%
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

  % normalize to initial dust amount, ie [ TS_SAL_MD(1) ]
  TS_DUST_UL = [ TS_SAL_MDHY TS_SAL_MDRGN ] .* NORM_FAC;
  LegTextUl = { 'M_D_H_Y' 'M_D_R_G_N' };
  LegLocUl = 'none';
  LcolorsUl = { ColorMdhy ColorMdrgn };

  fprintf('\n');

  % plot, 4 panels (2x2)
  OutFile = sprintf('%s/DpFig3_DustMassTseries_ICCP.jpg', Pdir);
  
  Fig = figure;

  Paxes = subplot(4,1,1);
  % adjust size so that x-axis lines up with the other panels
  PlotDpFigTseries(Paxes, T, TS_DUST_AL_MD, LcolorsAlMdMdadv, 'a', 'All Levels', 'Mass (Tg)', 'linear', [ 0 3 ], Fsize, 0, 1, LegTextAlMdMdadv, LegLocAlMdMdadv);
  columnlegend(3, LegTextAlMdMdadv, 'Location', 'NorthEast');

  Paxes = subplot(4,1,2);
  PlotDpFigTseries(Paxes, T, TS_DUST_AL, LcolorsAl, 'b', 'All Levels', 'Norm. Mass', 'log', [ 1e-4 2 ], Fsize, 0, 1, LegTextAl, LegLocAl);

  Paxes = subplot(4,1,3);
  PlotDpFigTseries(Paxes, T, TS_DUST_ML, LcolorsMl, 'c', 'Mid Levels', 'Norm. Mass', 'log', [ 1e-5  1e-1 ], Fsize, 0, 1, LegTextMl, LegLocMl);

  Paxes = subplot(4,1,4);
  PlotDpFigTseries(Paxes, T, TS_DUST_UL, LcolorsUl, 'd', 'Upper Levels', 'Norm. Mass', 'log', [ 1e-6  1e-2 ], Fsize, 1, 1, LegTextUl, LegLocUl);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
