function [ ] = PlotFsFigEvapThetaeXsection()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  % SAP, SD, cloud evaporation (mix ratio)
  InFname = 'DIAGS/storm_xsections_TSD_SAL_DUST.h5';
  InVname = '/all_s_cloud_evap';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  S_EVAP_SD = squeeze(h5read(InFname, InVname));

  X = squeeze(h5read(InFname, '/x_coords')) ./ 1000; % km
  Z = squeeze(h5read(InFname, '/z_coords')) ./ 1000; % km

  % SAP, SND, cloud evaporation (mix ratio)
  InFname = 'DIAGS/storm_xsections_TSD_SAL_NODUST.h5';
  InVname = '/all_s_cloud_evap';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  S_EVAP_SND = squeeze(h5read(InFname, InVname));

  % SAP, SD, theta-e
  InFname = 'DIAGS/storm_xsections_TSD_SAL_DUST.h5';
  InVname = '/all_s_theta_e';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  S_THETAE_SD = squeeze(h5read(InFname, InVname));

  % SAP, SND, theta-e
  InFname = 'DIAGS/storm_xsections_TSD_SAL_NODUST.h5';
  InVname = '/all_s_theta_e';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  S_THETAE_SND = squeeze(h5read(InFname, InVname));


  % Differences
  S_EVAP_DIFF = S_EVAP_SND - S_EVAP_SD;
  S_THETAE_DIFF = S_THETAE_SND - S_THETAE_SD;

  % Plot: 6 panels (3x2)
  OutFile = sprintf('%s/FsFigEvapThetaeXsection.jpg', Pdir);
  Fig = figure;

  Fsize = 13;
  Xlim = [ 0 500 ];
  Ylim = [ 0   6 ];

  % Evap settings
  EvapClim = [ -2.5 0 ];
  EvapClimDiff = [ -0.5 0.5 ];
  EvapClevs = -2.5:0.2:0;
  EvapClevsDiff = -0.5:0.01:0.5;
  EvapCmap = 'parula';
  EvapCmapDiff = 'redblue';

  % Theta-e settings
  TempClim = [ 335 355 ];
  TempClimDiff = [ -5 5 ];
  TempClevs = 335:2:355;
  TempClevsDiff = -5:0.1:5;
  TempCmap = 'parula';
  TempCmapDiff = 'redblue';

  PlocInc = 0.055;

  % Liq. Evap., left column
  Paxes = subplot(3, 2, 1);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_EVAP_SD', 'a', 'SAP (SD)', 'Radius (km)', Xlim, 'Height (km)', Ylim, EvapCmap, EvapClim, EvapClevs, Fsize, 0, 1);

  Paxes = subplot(3, 2, 3);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_EVAP_SND', 'c', 'SAP (SND)', 'Radius (km)', Xlim, 'Height (km)', Ylim, EvapCmap, EvapClim, EvapClevs, Fsize, 0, 1);

  Paxes = subplot(3, 2, 5);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_EVAP_DIFF', 'e', 'SAP (SND-SD)', 'Radius (km)', Xlim, 'Height (km)', Ylim, EvapCmapDiff, EvapClimDiff, EvapClevsDiff, Fsize, 1, 1);
%  text(0, -1, 'C', 'Color', 'b', 'FontSize', Fsize);
%  text(1700, -1, 'D', 'Color', 'b', 'FontSize', Fsize);

  % Tangential velocity in the right column
  Paxes = subplot(3, 2, 2);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_THETAE_SD', 'b', 'SAP (SD)', 'Radius (km)', Xlim, 'Height (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 0);

  Paxes = subplot(3, 2, 4);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_THETAE_SND', 'd', 'SAP (SND)', 'Radius (km)', Xlim, 'Height (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 0);

  Paxes = subplot(3, 2, 6);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_THETAE_DIFF', 'f', 'SAP (SND-SD)', 'Radius (km)', Xlim, 'Height (km)', Ylim, TempCmapDiff, TempClimDiff, TempClevsDiff, Fsize, 1, 0);
%  text(0, -1, 'C', 'Color', 'b', 'FontSize', Fsize);
%  text(1700, -1, 'D', 'Color', 'b', 'FontSize', Fsize);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
