function [ ] = PlotFsFigTempPtrack()

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

  % Theta along ptrack
  InFname = 'DIAGS/ptrack_avgs_TSD_SAL_DUST.h5';
  InVname = '/ps_theta';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_T_SD = squeeze(h5read(InFname, InVname))';

  X = squeeze(h5read(InFname, '/x_coords'));         % km
  Z = squeeze(h5read(InFname, '/z_coords')) ./ 1000; % km

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
  PS_T_F1 = PS_T_SND - PS_T_NSND;

  % Difference showing impact due to Dust (F2)
  PS_T_F2 = PS_T_NSD - PS_T_NSND;

  % Difference showing impact due to interaction 
  % of SAL Env and Dust (F12)
  PS_T_F12 = PS_T_SD - (PS_T_SND + PS_T_NSD) + PS_T_NSND;

  % Plot: 8 panels (4x2)
  OutFile = sprintf('%s/FsFigTempPtrack.jpg', Pdir);
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

  PlocInc = 0.055;

  % Sim data in left column
  Paxes = subplot(4, 2, 1);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_NSND, 'a', 'PTRACK: PSAP (NSND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 1);

  Paxes = subplot(4, 2, 3);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_SND, 'c', 'PSAP (SND)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 1);

  Paxes = subplot(4, 2, 5);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_NSD, 'e', 'PSAP (NSD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 0, 1);

  Paxes = subplot(4, 2, 7);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_SD, 'g', 'PSAP (SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmap, TempClim, TempClevs, Fsize, 1, 1);
  % mark PTRACK reference points near x-axis of bottom panel
  text(0, -2, 'C', 'Color', 'b', 'FontSize', Fsize);
  text(1700, -2, 'D', 'Color', 'b', 'FontSize', Fsize);


  % Factors in right column (skip top slot)
  Paxes = subplot(4, 2, 4);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_F1, 'd', 'PSAP (F1)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmapDiff, TempClimDiff, TempClevsDiff, Fsize, 0, 0);

  Paxes = subplot(4, 2, 6);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_F2, 'f', 'PSAP (F2)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmapDiff, TempClimDiff, TempClevsDiff, Fsize, 0, 0);

  Paxes = subplot(4, 2, 8);
  Ploc = get(Paxes, 'Position'); % cover a wider portion of the subplot region
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PlotFsFigXsection(Paxes, X, Z, PS_T_F12, 'h', 'PSAP (F12)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, TempCmapDiff, TempClimDiff, TempClevsDiff, Fsize, 1, 0);

  % mark PTRACK reference points near x-axis of bottom panel
  text(0, -2, 'C', 'Color', 'b', 'FontSize', Fsize);
  text(1700, -2, 'D', 'Color', 'b', 'FontSize', Fsize);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
