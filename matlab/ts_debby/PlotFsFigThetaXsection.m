function [ ] = PlotFsFigThetaXsection()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  CaseList = {
    { 'TSD_SAL_DUST'      'SD'   }
%    { 'TSD_SAL_NODUST'    'SND'  }
    { 'TSD_NONSAL_DUST'   'NSD'  }
%    { 'TSD_NONSAL_NODUST' 'NSND' }
    };
  Ncases = length(CaseList);

  % Want to use 0 for the first entry, but this file got written during initialization
  % before the non-SAL environment was introduced.
  SimTimes = [ 0.5 15 30 45 60 ];
  Nst = length(SimTimes);

  Zbot = 0; % km
  Ztop = 6; 
  
  Fsize = 13;

  % Read in the theta slices
  %   x-z planes where x is along the x-section line
  for icase = 1:Ncases
    Case  = CaseList{icase}{1};
    Label = CaseList{icase}{2};

    InFname = sprintf('XsectionData/strack_theta_%s.h5', Case);
    InVname = '/theta';
    fprintf('Reading %s (%s)\n', InFname, InVname);

    % 
    if (icase == 1)
      X = squeeze(h5read(InFname, '/x_coords')); % km
      Z = squeeze(h5read(InFname, '/z_coords'))./1000; % convert to km
      T = squeeze(h5read(InFname, '/t_coords'))./3600 - 42; % convert to sim time in hours starting with zero

      Nx = length(X);

      Z1 = find(Z >= Zbot, 1, 'first');
      Z2 = find(Z <= Ztop, 1, 'last');
      Z = Z(Z1:Z2);
      Nz = length(Z);
    end
     
    % for each time, extract the x-z slice
    for it = 1:Nst
      T1 = find(T >= SimTimes(it), 1, 'first');
      THETA(icase,it,:,:) = squeeze(h5read(InFname, InVname, [ 1 1 Z1 T1 ], [ Nx 1 Nz 1 ]));
    end
  end
  fprintf('\n');

  % Make difference plots: SD - NSD
  THETA_DIFF = squeeze(THETA(1,:,:,:) - THETA(2,:,:,:));

  % Plot: 10 panels (5x2)
  OutFile = sprintf('%s/FsFigThetaXsection.jpg', Pdir);
  Fig = figure;

  Xlim = [ 0 2150 ];
  Ylim = [ 0 5.5 ];
  Clim = [ 295 325 ];
  Clevs = 295:6:325;
  Cmap = 'parula';

  DiffClim = [ -5 5 ];
  DiffClevs = -5:1:5;
  DiffCmap = 'redblue';

  PlocInc = 0.05;

  % TSD_SAL_DUST
  Paxes = subplot(5, 3, 1);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA(1,1,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'a', '06Z, 22Aug (SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Cmap, Clim, Clevs, Fsize, 0, 1);

  Paxes = subplot(5, 3, 4);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA(1,2,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'd', '21Z, 22Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Cmap, Clim, Clevs, Fsize, 0, 1);

  Paxes = subplot(5, 3, 7);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA(1,3,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'g', '12Z, 23Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Cmap, Clim, Clevs, Fsize, 0, 1);

  Paxes = subplot(5, 3, 10);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA(1,4,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'j', '03Z, 24Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Cmap, Clim, Clevs, Fsize, 0, 1);

  Paxes = subplot(5, 3, 13);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA(1,5,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'm', '18Z, 24Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Cmap, Clim, Clevs, Fsize, 1, 1);
  text(0, -3.6, 'A', 'Color', 'r', 'FontSize', Fsize);
  text(1900, -3.6, 'B', 'Color', 'r', 'FontSize', Fsize);

  % TSD_NONSAL_DUST
  Paxes = subplot(5, 3, 2);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA(2,1,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'b', '06Z, 22Aug (NSD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Cmap, Clim, Clevs, Fsize, 0, 0);

  Paxes = subplot(5, 3, 5);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA(2,2,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'e', '21Z, 22Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Cmap, Clim, Clevs, Fsize, 0, 0);

  Paxes = subplot(5, 3, 8);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA(2,3,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'h', '12Z, 23Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Cmap, Clim, Clevs, Fsize, 0, 0);

  Paxes = subplot(5, 3, 11);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA(2,4,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'k', '03Z, 24Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Cmap, Clim, Clevs, Fsize, 0, 0);

  Paxes = subplot(5, 3, 14);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA(2,5,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'n', '18Z, 24Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Cmap, Clim, Clevs, Fsize, 1, 0);
  text(0, -3.6, 'A', 'Color', 'r', 'FontSize', Fsize);
  text(1900, -3.6, 'B', 'Color', 'r', 'FontSize', Fsize);

  % DIFFERENCE
  Paxes = subplot(5, 3, 3);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_DIFF(1,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'c', '06Z, 22Aug (SD-NSD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 3, 6);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_DIFF(2,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'f', '21Z, 22Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 3, 9);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_DIFF(3,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'i', '12Z, 23Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 3, 12);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_DIFF(4,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'l', '03Z, 24Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 3, 15);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_DIFF(5,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'o', '18Z, 24Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 1, 0);
  text(0, -3.6, 'A', 'Color', 'r', 'FontSize', Fsize);
  text(1900, -3.6, 'B', 'Color', 'r', 'FontSize', Fsize);


  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
