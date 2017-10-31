function [ ] = PlotFsFigStormVt()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  % pre-SAL, SD, u
  % Use sim time 20 h, which is half-way through the Pre-SAL period
  InFname = 'DIAGS/sample_hwind_TSD_SAL_DUST.h5';
  InVname = '/u_t20_z35';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  SD_U = squeeze(h5read(InFname, InVname));

  X = squeeze(h5read(InFname, '/x_coords'));
  Y = squeeze(h5read(InFname, '/y_coords'));

  % pre-SAL, SD, v
  InFname = 'DIAGS/sample_hwind_TSD_SAL_DUST.h5';
  InVname = '/v_t20_z35';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  SD_V = squeeze(h5read(InFname, InVname));

  % pre-SAL, NSD, u
  InFname = 'DIAGS/sample_hwind_TSD_NONSAL_DUST.h5';
  InVname = '/u_t20_z35';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  NSD_U = squeeze(h5read(InFname, InVname));

  % pre-SAL, NSD, v
  InFname = 'DIAGS/sample_hwind_TSD_NONSAL_DUST.h5';
  InVname = '/v_t20_z35';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  NSD_V = squeeze(h5read(InFname, InVname));

  % pre-SAL, SND, u
  InFname = 'DIAGS/sample_hwind_TSD_SAL_NODUST.h5';
  InVname = '/u_t20_z35';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  SND_U = squeeze(h5read(InFname, InVname));

  % pre-SAL, SND, v
  InFname = 'DIAGS/sample_hwind_TSD_SAL_NODUST.h5';
  InVname = '/v_t20_z35';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  SND_V = squeeze(h5read(InFname, InVname));

  % X,Y,U,V have too many points which muddies up the quiver plot
  % Take every nth point in order to make the plots more readable
  WvectInc = 60;
  WV_X = X(1:WvectInc:end); 
  WV_Y = Y(1:WvectInc:end); 
  WV_SD_U = SD_U(1:WvectInc:end, 1:WvectInc:end);
  WV_SD_V = SD_V(1:WvectInc:end, 1:WvectInc:end);
  WV_NSD_U = NSD_U(1:WvectInc:end, 1:WvectInc:end);
  WV_NSD_V = NSD_V(1:WvectInc:end, 1:WvectInc:end);
  WV_SND_U = SND_U(1:WvectInc:end, 1:WvectInc:end);
  WV_SND_V = SND_V(1:WvectInc:end, 1:WvectInc:end);

  % pre-SAL, SD, Vt
  InFname = 'DIAGS/storm_xsections_TSD_SAL_DUST.h5';
  InVname = '/all_ps_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_SD_VT = squeeze(h5read(InFname, InVname))';

  R = squeeze(h5read(InFname, '/x_coords')) ./ 1000; % km
  Z = squeeze(h5read(InFname, '/z_coords')) ./ 1000; % km

  % pre-SAL, SD, Updrafts
  InFname = 'DIAGS/storm_xsections_TSD_SAL_DUST.h5';
  InVname = '/all_ps_updraft';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_SD_UP = squeeze(h5read(InFname, InVname))' .* 100; % cm/s

  % pre-SAL, NSD, Vt
  InFname = 'DIAGS/storm_xsections_TSD_NONSAL_DUST.h5';
  InVname = '/all_ps_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_NSD_VT = squeeze(h5read(InFname, InVname))';

  % pre-SAL, NSD, Updrafts
  InFname = 'DIAGS/storm_xsections_TSD_NONSAL_DUST.h5';
  InVname = '/all_ps_updraft';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_NSD_UP = squeeze(h5read(InFname, InVname))' .* 100; % cm/s

  % pre-SAL, SND, Vt
  InFname = 'DIAGS/storm_xsections_TSD_SAL_NODUST.h5';
  InVname = '/all_ps_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_SND_VT = squeeze(h5read(InFname, InVname))';

  % pre-SAL, SND, Updrafts
  InFname = 'DIAGS/storm_xsections_TSD_SAL_NODUST.h5';
  InVname = '/all_ps_updraft';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_SND_UP = squeeze(h5read(InFname, InVname))' .* 100; % cm/s

  % Vt diffs
  PS_VT_DIFF_S = PS_SD_VT - PS_NSD_VT;
  PS_VT_DIFF_D = PS_SD_VT - PS_SND_VT;

  Rlim = [ 0 450 ];
  Zlim = [ 0   7 ];
  VtCmap = 'default';
  VtClim = [ 3 15 ];
  VtClevs = 3:1:15;

  UpClevs = 0:4:20;

  VtDiffCmap  = 'redblue';
  VtDiffClim  = [ -3 3 ];
  VtDiffClevs = -3:0.1:3;

  % Plot: 4 panels (2x2)
  OutFile = sprintf('%s/FsFigStormVt.jpg', Pdir);
  Fig = figure;

  Fsize = 12;
  Vscale = 1.5; % scale for vector size on quiver plot

  % Vt, left column
  Paxes = subplot(3, 3, 1);
  PlotFsFigStormStruct(Paxes, R, Z, PS_SD_VT, PS_SD_UP, 'a', 'PSAP (SD)', 'Radius (km)', Rlim, 'Z (km)', Zlim, VtCmap, VtClim, VtClevs, UpClevs, Fsize, 0, 1);

  Paxes = subplot(3, 3, 4);
  PlotFsFigStormStruct(Paxes, R, Z, PS_NSD_VT, PS_NSD_UP, 'd', '(NSD)', 'Radius (km)', Rlim, 'Z (km)', Zlim, VtCmap, VtClim, VtClevs, UpClevs, Fsize, 0, 1);

  Paxes = subplot(3, 3, 7);
  PlotFsFigStormStruct(Paxes, R, Z, PS_SND_VT, PS_SND_UP, 'g', '(SND)', 'Radius (km)', Rlim, 'Z (km)', Zlim, VtCmap, VtClim, VtClevs, UpClevs, Fsize, 1, 1);

  % Vt diffs, middle column
  Paxes = subplot(3, 3, 5);
  PlotFsFigXsection(Paxes, R, Z, PS_VT_DIFF_S, 'e', '(SD-NSD)', 'Radius (km)', Rlim, 'Z (km)', Zlim, VtDiffCmap, VtDiffClim, VtDiffClevs, Fsize, 0, 0);

  Paxes = subplot(3, 3, 8);
  PlotFsFigXsection(Paxes, R, Z, PS_VT_DIFF_D, 'h', '(SD-SND)', 'Radius (km)', Rlim, 'Z (km)', Zlim, VtDiffCmap, VtDiffClim, VtDiffClevs, Fsize, 1, 0);


  % Wind vectors, right column
  % Using sim time 20 h which correposonds to 2Z on Aug 23.
  Pinc1 = 0.02;
  Pinc3 = 0.02;

  Pinc2 = 0.022;
  Pinc4 = 0.02;

  Paxes = subplot(3, 3, 3);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - Pinc1;
  Ploc(3) = Ploc(3) + Pinc3;
  Ploc(2) = Ploc(2) - Pinc2;
  Ploc(4) = Ploc(4) + Pinc4;
  set(Paxes, 'Position', Ploc);
  PlotFsFigHwindVectors(Paxes, WV_X, WV_Y, WV_SD_U, WV_SD_V, 'c', '02Z, 23Aug (SD)', Fsize, Vscale);

  Paxes = subplot(3, 3, 6);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - Pinc1;
  Ploc(3) = Ploc(3) + Pinc3;
  Ploc(2) = Ploc(2) - Pinc2;
  Ploc(4) = Ploc(4) + Pinc4;
  set(Paxes, 'Position', Ploc);
  PlotFsFigHwindVectors(Paxes, WV_X, WV_Y, WV_NSD_U, WV_NSD_V, 'f', '02Z, 23Aug (NSD)', Fsize, Vscale);

  Paxes = subplot(3, 3, 9);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - Pinc1;
  Ploc(3) = Ploc(3) + Pinc3;
  Ploc(2) = Ploc(2) - Pinc2;
  Ploc(4) = Ploc(4) + Pinc4;
  set(Paxes, 'Position', Ploc);
  PlotFsFigHwindVectors(Paxes, WV_X, WV_Y, WV_SND_U, WV_SND_V, 'i', '02Z, 23Aug (SND)', Fsize, Vscale);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
