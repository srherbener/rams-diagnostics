function [ ] = PlotFsFigTempWindByDustXsection()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  % pre-SAL, SD, tangential to ptrack velocity
  InFname = 'DIAGS/ptrack_hvelocity_TSD_SAL_DUST.h5';
  InVname = '/ps_v';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_V_SD = squeeze(h5read(InFname, InVname));

  X = squeeze(h5read(InFname, '/x_coords'));         % km
  Z = squeeze(h5read(InFname, '/z_coords')) ./ 1000; % km

  % pre-SAL, SND, tangential to ptrack velocity
  InFname = 'DIAGS/ptrack_hvelocity_TSD_SAL_NODUST.h5';
  InVname = '/ps_v';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_V_SND = squeeze(h5read(InFname, InVname));

  % pre-SAL, SD, theta
  InFname = 'DIAGS/ptrack_avgs_TSD_SAL_DUST.h5';
  InVname = '/ps_theta';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_T_SD = squeeze(h5read(InFname, InVname));

  % pre-SAL, SND, theta
  InFname = 'DIAGS/ptrack_avgs_TSD_SAL_NODUST.h5';
  InVname = '/ps_theta';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_T_SND = squeeze(h5read(InFname, InVname));


  % Differences
  PS_V_DIFF = PS_V_SND - PS_V_SD;
  PS_T_DIFF = PS_T_SND - PS_T_SD;

  % Plot: 6 panels (3x2)
  OutFile = sprintf('%s/FsFigTempWindByDustXsection.jpg', Pdir);
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
%  TempClevs = -10:1:25;
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
  Paxes = subplot(3, 2, 1);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_SD', 'a', 'PTRACK: PSAP (SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 1);

  Paxes = subplot(3, 2, 3);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_SND', 'c', 'PSAP (SND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 1);

  Paxes = subplot(3, 2, 5);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_DIFF', 'e', 'PSAP (SND-SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, TempClimDiff, TempClevsDiff, Fsize, 1, 1);
  text(0, -1, 'C', 'Color', 'b', 'FontSize', Fsize);
  text(1700, -1, 'D', 'Color', 'b', 'FontSize', Fsize);

  % Tangential velocity in the right column
  Paxes = subplot(3, 2, 2);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_V_SD', 'b', 'PTRACK: PSAP (SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 0);

  Paxes = subplot(3, 2, 4);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_V_SND', 'd', 'PSAP (SND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 0);

  Paxes = subplot(3, 2, 6);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_V_DIFF', 'f', 'PSAP (SND-SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, VelClimDiff, VelClevsDiff, Fsize, 1, 0);
  text(0, -1, 'C', 'Color', 'b', 'FontSize', Fsize);
  text(1700, -1, 'D', 'Color', 'b', 'FontSize', Fsize);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
