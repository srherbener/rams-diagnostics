function [ ] = PlotDpFigHovmollers()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

  % Hovmoller data for Md
  InFile = 'DIAGS/storm_hovmollers_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  HOV_DUST = squeeze(h5read(InFile, InVname));
  Z    = squeeze(h5read(InFile, '/z_coords'))./1000;         % convert to km
  T    = squeeze(h5read(InFile, '/t_coords'))./3600 - 42;    % convert to sim time in hours

  % Hovmoller data for Mdhy
  InFile = 'DIAGS/storm_hovmollers_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_dust_hydro';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  HOV_DHY = squeeze(h5read(InFile, InVname));

  % Hovmoller data for Mdrgn
  InFile = 'DIAGS/storm_hovmollers_TSD_SAL_DUST.h5';
  InVname = '/sal_ar_ra_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);
  HOV_DRGN = squeeze(h5read(InFile, InVname));

  % plot
  Fig = figure;

  % Md
  Paxes = subplot(3,1,1);
  PlotDpFigDustHov(Paxes, T, Z, HOV_DUST, 'a', 'M_D', Fsize, 0, 1, 0, 2);

  % Mdhy
  Paxes = subplot(3,1,2);
  PlotDpFigDustHov(Paxes, T, Z, HOV_DHY, 'b', 'M_D_H_Y', Fsize, 0, 1, 0, 0);

  % Mdrgn
  Paxes = subplot(3,1,3);
  PlotDpFigDustHov(Paxes, T, Z, HOV_DRGN, 'c', 'M_D_R_G_N', Fsize, 1, 1, 0, 1);
  
  OutFile = sprintf('%s/DpFig2_Hovmollers.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
