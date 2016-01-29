function [ ] = PlotDpFigDustRainout()
% PlotDpFigDustRainout 

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  Fsize = 13;
  
  % Hovmoller data for dust-in-hydrometeors
  InFile = 'DIAGS/storm_hovmollers_TSD_SAL_DUST.h5';
  InVname = '/sal_dust_hydro';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  HOV_DUST = squeeze(h5read(InFile, InVname));
  Z    = squeeze(h5read(InFile, '/z_coords'))./1000;
  T    = squeeze(h5read(InFile, '/t_coords'))./3600 - 42;

  % Time series of lower level Mdhy
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_dust_hydro_total_mass_llev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDHY = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MDHY = TS_SAL_MDHY .* 1e-9;  % convert to millions of kg

  % Time series of precip rate
  InFile = 'DIAGS/hist_meas_ts_pcprate_TSD_SAL_DUST.h5';
  InVname = '/sal_pcprate_ts';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_PR = squeeze(h5read(InFile, InVname)); % mm/h

  % Time series of surface dust accumulation
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_dust_sfc_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_DSFC = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_DSFC = TS_SAL_DSFC .* 1e-12;  % convert to billions of kg

  % plot, 4 panels (2x2)
  OutFile = sprintf('%s/DpFigDustRainout.jpg', Pdir);
  
  Fig = figure;

  % create the Hovmoller of sal_dust_hydro
  Paxes = subplot(4,1,1);
  PlotDpFigDustHov(Paxes, T, Z, HOV_DUST, 'a', '', Fsize, 0, 1, -1, 0);
  
  Paxes = subplot(4,1,2);
  PlotDpFigTseries(Paxes, T, TS_SAL_MDHY, 'b', '', 'M_d_h_y (10^6 kg)', Fsize, 0, 1);

  Paxes = subplot(4,1,3);
  PlotDpFigTseries(Paxes, T, TS_SAL_PR, 'c', '', 'PR (mm h^-^1)', Fsize, 0, 1);

  Paxes = subplot(4,1,4);
  PlotDpFigTseries(Paxes, T, TS_SAL_DSFC, 'd', '', 'M_d (10^9 kg)', Fsize, 1, 1);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
