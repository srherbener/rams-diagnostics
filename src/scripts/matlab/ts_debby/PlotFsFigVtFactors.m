function [ ] = PlotFsFigVtFactors()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  % Factors from factor separation method:
  %    1: SAL Env
  %    2: Dust
  %
  %   F0:  independent of SAL, Dust:        NSND
  %   F1:  SAL Env:                         SND - NSND
  %   F2:  Dust:                            NSD - NSND
  %   F12: Interaction of SAL Env and Dust: SD - (SND+NSD) + NSND

  % Tangential velocity relative to storm center, Vt
  % Pre-SAL time period
  InFname = 'DIAGS/storm_xsections_TSD_SAL_DUST.h5';
  InVname = '/all_ps_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_VT_SD = squeeze(h5read(InFname, InVname))';

  X = squeeze(h5read(InFname, '/x_coords')) ./ 1000; % km
  Z = squeeze(h5read(InFname, '/z_coords')) ./ 1000; % km

  InFname = 'DIAGS/storm_xsections_TSD_SAL_NODUST.h5';
  InVname = '/all_ps_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_VT_SND = squeeze(h5read(InFname, InVname))';

  InFname = 'DIAGS/storm_xsections_TSD_NONSAL_DUST.h5';
  InVname = '/all_ps_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_VT_NSD = squeeze(h5read(InFname, InVname))';

  InFname = 'DIAGS/storm_xsections_TSD_NONSAL_NODUST.h5';
  InVname = '/all_ps_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_VT_NSND = squeeze(h5read(InFname, InVname))';

  % SAL time period
  InFname = 'DIAGS/storm_xsections_TSD_SAL_DUST.h5';
  InVname = '/all_s_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  S_VT_SD = squeeze(h5read(InFname, InVname))';

  X = squeeze(h5read(InFname, '/x_coords')) ./ 1000; % km
  Z = squeeze(h5read(InFname, '/z_coords')) ./ 1000; % km

  InFname = 'DIAGS/storm_xsections_TSD_SAL_NODUST.h5';
  InVname = '/all_s_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  S_VT_SND = squeeze(h5read(InFname, InVname))';

  InFname = 'DIAGS/storm_xsections_TSD_NONSAL_DUST.h5';
  InVname = '/all_s_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  S_VT_NSD = squeeze(h5read(InFname, InVname))';

  InFname = 'DIAGS/storm_xsections_TSD_NONSAL_NODUST.h5';
  InVname = '/all_s_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  S_VT_NSND = squeeze(h5read(InFname, InVname))';


  % Difference showing impact due to SAL Env (F1)
  PS_VT_F1 = PS_VT_SND - PS_VT_NSND;
  S_VT_F1 = S_VT_SND - S_VT_NSND;

  % Difference showing of impact due to Dust (F2)
  PS_VT_F2 = PS_VT_NSD - PS_VT_NSND;
  S_VT_F2 = S_VT_NSD - S_VT_NSND;

  % Difference showing of impact due to interaction
  % of SAL Env Dust (F12)
  PS_VT_F12 = PS_VT_SD - (PS_VT_SND + PS_VT_NSD) + PS_VT_NSND;
  S_VT_F12 = S_VT_SD - (S_VT_SND + S_VT_NSD) + S_VT_NSND;

  % Plot: 8 panels (4x2) - Pre-SAL time period`
  OutFile = sprintf('%s/FsFigVtPsapFactors.jpg', Pdir);
  Fig = figure;

  Fsize = 13;
  Xlim = [ 0 500 ];
  Ylim = [ 0 8 ];

  VelClim = [ 0 20 ];
  VelClimDiff = [ -2 2 ];
  VelClevs = 0:2:20;
  VelClevsDiff = -2:0.2:2;
  VelCmap = 'jet';
  VelCmapDiff = 'redblue';

  PlocInc = 0.055;

  % Sim data in left column
  Paxes = subplot(4, 2, 1);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_VT_NSND, 'a', 'PSAP: V_t (NSND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 1);

  Paxes = subplot(4, 2, 3);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_VT_SND, 'c', 'V_t (SND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 1);

  Paxes = subplot(4, 2, 5);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_VT_NSD, 'e', 'V_t (NSD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 1);

  Paxes = subplot(4, 2, 7);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_VT_SD, 'g', 'V_t (SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 1, 1);


  % Factors in right column, skip the top column
  Paxes = subplot(4, 2, 4);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_VT_F1, 'd', 'V_t (F1)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, VelClimDiff, VelClevsDiff, Fsize, 0, 0);

  Paxes = subplot(4, 2, 6);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_VT_F2, 'f', 'V_t (F2)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, VelClimDiff, VelClevsDiff, Fsize, 0, 0);

  Paxes = subplot(4, 2, 8);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_VT_F12, 'f', 'V_t (F12)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, VelClimDiff, VelClevsDiff, Fsize, 1, 0);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);


  % Plot: 8 panels (4x2) - SAL time period`
  OutFile = sprintf('%s/FsFigVtSapFactors.jpg', Pdir);
  Fig = figure;

  Fsize = 13;
  Xlim = [ 0 500 ];
  Ylim = [ 0 8 ];

  VelClim = [ 0 20 ];
  VelClimDiff = [ -5 5 ];
  VelClevs = 0:2:20;
  VelClevsDiff = -5:0.5:5;
  VelCmap = 'jet';
  VelCmapDiff = 'redblue';

  PlocInc = 0.055;

  % Sim data in left column
  Paxes = subplot(4, 2, 1);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_VT_NSND, 'a', 'SAP: V_t (NSND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 1);

  Paxes = subplot(4, 2, 3);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_VT_SND, 'c', 'V_t (SND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 1);

  Paxes = subplot(4, 2, 5);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_VT_NSD, 'e', 'V_t (NSD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 1);

  Paxes = subplot(4, 2, 7);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_VT_SD, 'g', 'V_t (SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 1, 1);


  % Factors in right column, skip the top column
  Paxes = subplot(4, 2, 4);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_VT_F1, 'd', 'V_t (F1)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, VelClimDiff, VelClevsDiff, Fsize, 0, 0);

  Paxes = subplot(4, 2, 6);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_VT_F2, 'f', 'V_t (F2)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, VelClimDiff, VelClevsDiff, Fsize, 0, 0);

  Paxes = subplot(4, 2, 8);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, S_VT_F12, 'f', 'V_t (F12)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, VelClimDiff, VelClevsDiff, Fsize, 1, 0);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
