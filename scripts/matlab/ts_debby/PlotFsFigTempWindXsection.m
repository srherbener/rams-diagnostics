function [ ] = PlotFsFigTempWindXsection()

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

  % Velocity normal to ptrack, and along ptrack
  InFname = 'DIAGS/ptrack_hvelocity_TSD_SAL_NODUST.h5';
  InVname = '/ps_v';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_V_SND = squeeze(h5read(InFname, InVname))';

  X = squeeze(h5read(InFname, '/x_coords'));         % km
  Z = squeeze(h5read(InFname, '/z_coords')) ./ 1000; % km

  InFname = 'DIAGS/ptrack_hvelocity_TSD_NONSAL_DUST.h5';
  InVname = '/ps_v';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_V_NSD = squeeze(h5read(InFname, InVname))';

  InFname = 'DIAGS/ptrack_hvelocity_TSD_NONSAL_NODUST.h5';
  InVname = '/ps_v';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_V_NSND = squeeze(h5read(InFname, InVname))';


  % Theta along ptrack
  InFname = 'DIAGS/ptrack_avgs_TSD_SAL_NODUST.h5';
  InVname = '/ps_theta';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_T_SND = squeeze(h5read(InFname, InVname))';

  InFname = 'DIAGS/ptrack_avgs_TSD_NONSAL_DUST.h5';
  InVname = '/ps_theta';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_T_NSD = squeeze(h5read(InFname, InVname))';

  % pre-SAL, NSND, theta
  InFname = 'DIAGS/ptrack_avgs_TSD_NONSAL_NODUST.h5';
  InVname = '/ps_theta';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_T_NSND = squeeze(h5read(InFname, InVname))';


  % Difference showing impact due to SAL Env (F1)
  PS_V_F1 = PS_V_SND - PS_V_NSND;
  PS_T_F1 = PS_T_SND - PS_T_NSND;

  % Difference showing of impact due to Dust (F2)
  PS_V_F2 = PS_V_NSD - PS_V_NSND;
  PS_T_F2 = PS_T_NSD - PS_T_NSND;

  % Plot: 6 panels (3x2)
  OutFile = sprintf('%s/FsFigTempWindXsection.jpg', Pdir);
  Fig = figure;

  Fsize = 13;
  Xlim = [ 0 1850 ];
  Ylim = [ 0 5.5 ];

%  Theta settings
  TempClim = [ 290 320 ];
  TempClimDiff = [ -4 4 ];
  TempClevs = 290:4:320;
  TempClevsDiff = -4:0.1:4;

%%  Temp C settings
%  TempClim = [ -10 25 ];
%  TempClimDiff = [ -4 4 ];
%  TempClevs = -10:5:25;
%  TempClevsDiff = -4:0.1:4;

  TempCmap = 'parula';
  TempCmapDiff = 'redblue';

  VelClim = [ -10 20 ];
  VelClimDiff = [ -4 4 ];
  VelClevs = -10:2:20;
  VelClevsDiff = -4:0.1:4;
  VelCmap = 'jet';
  VelCmapDiff = 'redblue';

  PlocInc = 0.055;

  % Theta, left column
  Paxes = subplot(5, 2, 1);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_NSND, 'a', 'PTRACK: PSAP (NSND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 1);

  Paxes = subplot(5, 2, 3);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_SND, 'c', 'PSAP (SND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 1);

  Paxes = subplot(5, 2, 5);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_NSD, 'g', 'PSAP (NSD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 1);

  Paxes = subplot(5, 2, 7);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_F1, 'e', 'PSAP (F1)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, TempClimDiff, TempClevsDiff, Fsize, 0, 1);

  Paxes = subplot(5, 2, 9);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_F2, 'i', 'PSAP (F2)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, TempClimDiff, TempClevsDiff, Fsize, 1, 1);
  text(0, -2, 'C', 'Color', 'b', 'FontSize', Fsize);
  text(1700, -2, 'D', 'Color', 'b', 'FontSize', Fsize);

  % Tangential velocity in the right column
  Paxes = subplot(5, 2, 2);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_V_NSND, 'b', 'PTRACK: PSAP (NSND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 0);

  Paxes = subplot(5, 2, 4);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_V_SND, 'd', 'PSAP (SND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 0);

  Paxes = subplot(5, 2, 6);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_V_NSD, 'h', 'PSAP (NSD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 0);

  Paxes = subplot(5, 2, 8);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_V_F1, 'f', 'PSAP (F1)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, VelClimDiff, VelClevsDiff, Fsize, 0, 0);

  Paxes = subplot(5, 2, 10);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_V_F2, 'j', 'PSAP (F2)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, VelClimDiff, VelClevsDiff, Fsize, 1, 0);
  text(0, -2, 'C', 'Color', 'b', 'FontSize', Fsize);
  text(1700, -2, 'D', 'Color', 'b', 'FontSize', Fsize);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
