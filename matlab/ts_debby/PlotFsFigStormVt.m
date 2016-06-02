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

  % X,Y,U,V have too many points which muddies up the quiver plot
  % Take every nth point in order to make the plots more readable
  WvectInc = 60;
  WV_X = X(1:WvectInc:end); 
  WV_Y = Y(1:WvectInc:end); 
  WV_SD_U = SD_U(1:WvectInc:end, 1:WvectInc:end);
  WV_SD_V = SD_V(1:WvectInc:end, 1:WvectInc:end);
  WV_NSD_U = NSD_U(1:WvectInc:end, 1:WvectInc:end);
  WV_NSD_V = NSD_V(1:WvectInc:end, 1:WvectInc:end);

  % pre-SAL, SD, Vt
  InFname = 'DIAGS/storm_xsections_TSD_SAL_DUST.h5';
  InVname = '/all_ps_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_SD_VT = squeeze(h5read(InFname, InVname));

  R = squeeze(h5read(InFname, '/x_coords')) ./ 1000; % km
  Z = squeeze(h5read(InFname, '/z_coords')) ./ 1000; % km

  % pre-SAL, SD, Updrafts
  InFname = 'DIAGS/storm_xsections_TSD_SAL_DUST.h5';
  InVname = '/all_ps_updraft';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_SD_UP = squeeze(h5read(InFname, InVname)) .* 100; % cm/s

  % pre-SAL, NSD, Vt
  InFname = 'DIAGS/storm_xsections_TSD_NONSAL_DUST.h5';
  InVname = '/all_ps_speed_t';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_NSD_VT = squeeze(h5read(InFname, InVname));

  % pre-SAL, NSD, Updrafts
  InFname = 'DIAGS/storm_xsections_TSD_NONSAL_DUST.h5';
  InVname = '/all_ps_updraft';
  fprintf('Reading: %s (%s)\n', InFname, InVname);
  PS_NSD_UP = squeeze(h5read(InFname, InVname)) .* 100; % cm/s

  Rlim = [ 0 450 ];
  Zlim = [ 0   7 ];
  VtCmap = 'default';
  VtClim = [ 3 15 ];
  VtClevs = 3:1:15;

  UpClevs = 0:4:20;

  % Plot: 4 panels (2x2)
  OutFile = sprintf('%s/FsFigStormVt.jpg', Pdir);
  Fig = figure;

  Fsize = 13;
  Vscale = 1.5; % scale for vector size on quiver plot

  % SD, left column
  % Using sim time 20 h which correposonds to 2Z on Aug 23.
  Paxes = subplot(2, 2, 1);
  PlotFsFigStormStruct(Paxes, R, Z, PS_SD_VT', PS_SD_UP', 'a', 'PSAP (SD)', 'Radius (km)', Rlim, 'Z (km)', Zlim, VtCmap, VtClim, VtClevs, UpClevs, Fsize, 1, 1);

  Paxes = subplot(2, 2, 3);
  PlotFsFigHwindVectors(Paxes, WV_X, WV_Y, WV_SD_U, WV_SD_V, 'c', '02Z, 23Aug (SD)', Fsize, Vscale);

  % NSD, right column
  Paxes = subplot(2, 2, 2);
  PlotFsFigStormStruct(Paxes, R, Z, PS_NSD_VT', PS_NSD_UP', 'b', 'PSAP (NSD)', 'Radius (km)', Rlim, 'Z (km)', Zlim, VtCmap, VtClim, VtClevs, UpClevs, Fsize, 1, 0);

  Paxes = subplot(2, 2, 4);
  PlotFsFigHwindVectors(Paxes, WV_X, WV_Y, WV_NSD_U, WV_NSD_V, 'd', '02Z, 23Aug (NSD)', Fsize, Vscale);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
