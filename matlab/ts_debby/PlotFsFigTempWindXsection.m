function [ ] = PlotFsFigTempWindXsection()

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

  % pre-SAL, NSD, tangential to ptrack velocity
  InFname = 'DIAGS/ptrack_hvelocity_TSD_NONSAL_DUST.h5';
  InVname = '/ps_v';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_V_NSD = squeeze(h5read(InFname, InVname));

  % pre-SAL, SD, theta
  InFname = 'DIAGS/ptrack_avgs_TSD_SAL_DUST.h5';
  InVname = '/ps_theta';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_T_SD = squeeze(h5read(InFname, InVname));

  X = squeeze(h5read(InFname, '/x_coords'));         % km
  Z = squeeze(h5read(InFname, '/z_coords')) ./ 1000; % km

  % pre-SAL, NSD, theta
  InFname = 'DIAGS/ptrack_avgs_TSD_NONSAL_DUST.h5';
  InVname = '/ps_theta';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_T_NSD = squeeze(h5read(InFname, InVname));


  % Differences
  PS_V_DIFF = PS_V_NSD - PS_V_SD;
  PS_T_DIFF = PS_T_NSD - PS_T_SD;

  % Plot: 6 panels (3x2)
  OutFile = sprintf('%s/FsFigTempWindXsection.jpg', Pdir);
  Fig = figure;

  Fsize = 13;
  Xlim = [ 0 1850 ];
  Ylim = [ 0 5.5 ];

%  Theta settings
  TempClim = [ 290 320 ];
  TempClimDiff = [ -4 4 ];
  TempClevs = 290:2:320;
  TempClevsDiff = -4:0.1:4;

%%  Temp C settings
%  TempClim = [ -10 25 ];
%  TempClimDiff = [ -4 4 ];
%  TempClevs = -10:1:25;
%  TempClevsDiff = -4:0.1:4;

  TempCmap = 'default';
  TempCmapDiff = 'redblue';

  VelClim = [ -10 20 ];
  VelClimDiff = [ -4 4 ];
  VelClevs = -10:1:20;
  VelClevsDiff = -4:0.1:4;
  VelCmap = 'default';
  VelCmapDiff = 'redblue';

  PlocInc = 0.055;

  % Theta, left column
  Paxes = subplot(3, 2, 1);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_SD', 'a', 'Pre-SAL (SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 1);

  Paxes = subplot(3, 2, 3);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_NSD', 'c', 'Pre-SAL (NSD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 1);

  Paxes = subplot(3, 2, 5);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_DIFF', 'e', 'Pre-SAL (NSD-SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, TempClimDiff, TempClevsDiff, Fsize, 1, 1);

  % Tangential velocity in the right column
  Paxes = subplot(3, 2, 2);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_V_SD', 'b', 'Pre-SAL (SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 0);

  Paxes = subplot(3, 2, 4);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_V_NSD', 'd', 'Pre-SAL (NSD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmap, VelClim, VelClevs, Fsize, 0, 0);

  Paxes = subplot(3, 2, 6);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_V_DIFF', 'f', 'Pre-SAL (NSD-SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, VelCmapDiff, VelClimDiff, VelClevsDiff, Fsize, 1, 0);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
