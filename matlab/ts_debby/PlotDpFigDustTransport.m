function [ ] = PlotDpFigDustTransport()
% PlotDpFigDustTransport 

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  Fsize = 13;
  

  % Time series of upper level Md
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MD = squeeze(h5read(InFile, InVname)); % g
  T    = squeeze(h5read(InFile, '/t_coords'))./3600 -42;

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
  InVname = '/sal_ar_residual_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  TS_SAL_MDADV = squeeze(h5read(InFile, InVname)); % g

  TS_DUST = [ TS_SAL_MD TS_SAL_MDHY TS_SAL_MDRGN TS_SAL_MDADV ] .* 1e-9; % convert to 1e-3 Tg
  LegText = { 'M_d' 'M_d_h_y' 'M_d_r_g_n' 'M_d_a_d_v' };
  LegLoc = 'NorthEast';

  % plot, 4 panels (2x2)
  OutFile = sprintf('%s/DpFig3_DustTransport.jpg', Pdir);
  
  Fig = figure;

  PlotDpFigTseries(gca, T, TS_DUST, 'a', 'Upper Levels', 'Mass (10^-^3 Tg)', 'linear', [ -0.5 3.5 ], Fsize, 1, 1, LegText, LegLoc);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
