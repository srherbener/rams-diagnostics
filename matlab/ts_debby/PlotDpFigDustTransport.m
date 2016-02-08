function [ ] = PlotDpFigDustTransport()
% PlotDpFigDustTransport 

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  Fsize = 13;
  
  % Hovmoller data for dust-in-hydrometeors
  InFile = 'DIAGS/storm_hovmollers_TSD_SAL_DUST.h5';
  InVname = '/storm_dust_hydro';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  HOV_MDHY = squeeze(h5read(InFile, InVname));
  Z    = squeeze(h5read(InFile, '/z_coords'))./1000;
  T    = squeeze(h5read(InFile, '/t_coords'))./3600 -42;

  % Hovmoller data for dust
  InFile = 'DIAGS/storm_hovmollers_TSD_SAL_DUST.h5';
  InVname = '/storm_aero_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  HOV_MD = squeeze(h5read(InFile, InVname));

  % Time series of upper level Mdhy
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/storm_aero_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_STORM_MD = squeeze(h5read(InFile, InVname)); % g
  TS_STORM_MD = TS_STORM_MD .* 1e-9;  % convert to 1e-3 Tg

  % Time series of upper level Mdhy
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/storm_dust_hydro_total_mass_hlev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_STORM_MDHY = squeeze(h5read(InFile, InVname)); % g
  TS_STORM_MDHY = TS_STORM_MDHY .* 1e-9;  % convert 1e-3 Tg

  % plot, 4 panels (2x2)
  OutFile = sprintf('%s/DpFigDustTransport.jpg', Pdir);
  
  Fig = figure;

  % create the Hovmoller of storm_dust_hydro
  Paxes = subplot(4,1,1);
  PlotDpFigDustHov(Paxes, T, Z, HOV_MD, 'a', 'STORM', Fsize, 0, 1, 0, 0);
  
  Paxes = subplot(4,1,2);
  PlotDpFigTseries(Paxes, T, TS_STORM_MD, 'b', 'STORM', 'M_d (10^-^3 Tg)', Fsize, 0, 1);

  Paxes = subplot(4,1,3);
  PlotDpFigDustHov(Paxes, T, Z, HOV_MDHY, 'c', 'STORM', Fsize, 0, 1, 0, 0);
  
  Paxes = subplot(4,1,4);
  PlotDpFigTseries(Paxes, T, TS_STORM_MDHY, 'd', 'STORM', 'M_d_h_y (10^-^3 Tg)', Fsize, 1, 1);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
