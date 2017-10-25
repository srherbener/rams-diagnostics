function [ ] = PlotDpFigDustMassTseries_ICCP_2()

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

  ColorMdrgn1  = 'orangered';
  ColorMdrgn2  = 'forestgreen';
  
  %%%%%%%%%%%%%%%% Volumn integrated dust - surface to tropopause %%%%%%%%%%%%%
  % Time series of Md
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MD = squeeze(h5read(InFile, InVname)); % g
  T    = squeeze(h5read(InFile, '/t_coords'))./3600 - 42;

  % Absolute values of Md and Mdadv
  TS_DUST_AL_MD = TS_SAL_MD .* 1e-12;  % Tg
  LegTextAlMdMdadv = { 'M_D' };
  LegLocAlMdMdadv  = 'NorthEastOutside';
  LcolorsAlMdMdadv = { ColorMd };

  % normalization of all quantities are relative to the total initial mass at all levels
  NORM_FAC = 1 / TS_SAL_MD(1);

  %%%%%%%%%%%%%%%% Volumn integrated dust - 4 km AGL to 7 km AGL %%%%%%%%%%%%%
  % Time series of mid level Mdrgn
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_ra_total_mass_mlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDRGN = squeeze(h5read(InFile, InVname)); % g

  % Time series of mid level Mdrgn1
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_ra1_total_mass_mlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDRGN1 = squeeze(h5read(InFile, InVname)); % g

  % Time series of mid level Mdrgn2
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_ra2_total_mass_mlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDRGN2 = squeeze(h5read(InFile, InVname)); % g

  % normalize to initial dust amount, ie [ TS_SAL_MD(1) ]
  TS_DUST_ML = TS_SAL_MDRGN .* NORM_FAC;
  LegTextMl = { 'M_D_R_G_N' };
  LegLocMl = 'NorthEastOutside';
  LcolorsMl = { ColorMdrgn };

  % normalize to initial dust amount, ie [ TS_SAL_MD(1) ]
  TS_DUST_ML_SEP = [ TS_SAL_MDRGN1 TS_SAL_MDRGN2 ] .* NORM_FAC;
  LegTextMlSep = { 'M_D_R_G_N_1' 'M_D_R_G_N_2' };
  LegLocMlSep = 'NorthEastOutside';
  LcolorsMlSep = { ColorMdrgn1 ColorMdrgn2 };

  fprintf('\n');

  % plot, 4 panels (2x2)
  OutFile = sprintf('%s/DpFig3_DustMassTseries_ICCP_2.jpg', Pdir);
  
  Fig = figure;

  Paxes = subplot(4,1,1);
  % adjust size so that x-axis lines up with the other panels
  Ploc = get(Paxes, 'Position');
  Ploc(3) = Ploc(3) * 0.958;
  Ploc = set(Paxes, 'Position', Ploc);
  PlotDpFigTseries(Paxes, T, TS_DUST_AL_MD, LcolorsAlMdMdadv, 'a', 'All Levels', 'Mass (Tg)', 'linear', [ 0 2 ], Fsize, 0, 1, LegTextAlMdMdadv, LegLocAlMdMdadv);

  Paxes = subplot(4,1,2);
  PlotDpFigTseries(Paxes, T, TS_DUST_ML, LcolorsMl, 'b', 'Mid Levels', 'Norm. Mass', 'log', [ 1e-5  1e-1 ], Fsize, 0, 1, LegTextMl, LegLocMl);

  Paxes = subplot(4,1,3);
  PlotDpFigTseries(Paxes, T, TS_DUST_ML_SEP, LcolorsMlSep, 'c', 'Mid Levels', 'Norm. Mass', 'log', [ 1e-5  1e-1 ], Fsize, 1, 1, LegTextMlSep, LegLocMlSep);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
