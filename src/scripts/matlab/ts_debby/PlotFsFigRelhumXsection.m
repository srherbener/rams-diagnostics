function [ ] = PlotFsFigRelhumXsection()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  CaseList = {
    { 'TSD_NONSAL_NODUST' 'NSND' }
    { 'TSD_SAL_NODUST'    'SND'  }
    { 'TSD_NONSAL_DUST'   'NSD'  }
    { 'TSD_SAL_DUST'      'SD'   }
    };
  Ncases = length(CaseList);

  % Want to use 0 for the first entry, but this file got written during initialization
  % before the non-SAL environment was introduced.
  SimTimes = [ 0.5 15 30 45 60 ];
  Nst = length(SimTimes);

  Zbot = 0; % km
  Ztop = 6; 
  
  Fsize = 8;

  % Read in the relhum slices
  %   x-z planes where x is along the x-section line
  for icase = 1:Ncases
    Case  = CaseList{icase}{1};
    Labels{icase} = CaseList{icase}{2};

    InFname = sprintf('XsectionData/strack_relhum_%s.h5', Case);
    InVname = '/relhum';
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
      RELHUM(icase,it,:,:) = squeeze(h5read(InFname, InVname, [ 1 1 Z1 T1 ], [ Nx 1 Nz 1 ]));
    end
  end
  fprintf('\n');

  % Calculate factors:
  %   F1 = SND - NSND
  %   F2 = NSD - NSND
  %   F12 = SD - (SND+NSD) + NSND
  F1  = squeeze(RELHUM(2,:,:,:) - RELHUM(1,:,:,:));
  F2  = squeeze(RELHUM(3,:,:,:) - RELHUM(1,:,:,:));
  F12 = squeeze(RELHUM(4,:,:,:) -(RELHUM(2,:,:,:)+RELHUM(3,:,:,:)) + RELHUM(1,:,:,:));

  % Plot: 10 panels (5x2)
  OutFile = sprintf('%s/FsFigRelhumXsection.jpg', Pdir);
  Fig = figure;

  Xlim = [ 0 2150 ];
  Ylim = [ 0 5.5 ];
  Clim = [ 0 100 ];
  Clevs = 0:5:100;
  Cmap = 'parula';

  DiffClim  = [ -50 50 ];
  DiffClevs = -50:10:50;
  DiffCmap = 'redblue';

  Xlab = 'Distance (km)';
  Ylab = 'Z (km)';

  PlocInc = 0.05;

  % CONTROL (NSND)
  Paxes = subplot(5, 4, 1);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(RELHUM(1,1,:,:))';
  Ptitle = sprintf('06Z, 22Aug (%s)', Labels{1});
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'a', Ptitle, Xlab, Xlim, Ylab, Ylim, Cmap, Clim, Clevs, Fsize, 0, 1);

  Paxes = subplot(5, 4, 5);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(RELHUM(1,2,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'e', '21Z, 22Aug', Xlab, Xlim, Ylab, Ylim, Cmap, Clim, Clevs, Fsize, 0, 1);

  Paxes = subplot(5, 4, 9);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(RELHUM(1,3,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'i', '12Z, 23Aug', Xlab, Xlim, Ylab, Ylim, Cmap, Clim, Clevs, Fsize, 0, 1);

  Paxes = subplot(5, 4, 13);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(RELHUM(1,4,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'm', '03Z, 24Aug', Xlab, Xlim, Ylab, Ylim, Cmap, Clim, Clevs, Fsize, 0, 1);

  Paxes = subplot(5, 4, 17);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(RELHUM(1,5,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'q', '18Z, 24Aug', Xlab, Xlim, Ylab, Ylim, Cmap, Clim, Clevs, Fsize, 1, 1);
  text(0, -3.6, 'A', 'Color', 'r', 'FontSize', Fsize);
  text(1900, -3.6, 'B', 'Color', 'r', 'FontSize', Fsize);

  % F1
  Paxes = subplot(5, 4, 2);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F1(1,:,:))';
  Ptitle = '06Z, 22Aug (F1)';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'b', Ptitle, Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 6);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F1(2,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'f', '21Z, 22Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 10);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F1(3,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'j', '12Z, 23Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 14);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F1(4,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'n', '03Z, 24Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 18);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F1(5,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'r', '18Z, 24Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 1, 0);
  text(0, -3.6, 'A', 'Color', 'r', 'FontSize', Fsize);
  text(1900, -3.6, 'B', 'Color', 'r', 'FontSize', Fsize);

  % F2
  Paxes = subplot(5, 4, 3);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F2(1,:,:))';
  Ptitle = '06Z, 22Aug (F2)';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'c', Ptitle, Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 7);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F2(2,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'g', '21Z, 22Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 11);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F2(3,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'k', '12Z, 23Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 15);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F2(4,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'o', '03Z, 24Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 19);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F2(5,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 's', '18Z, 24Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 1, 0);
  text(0, -3.6, 'A', 'Color', 'r', 'FontSize', Fsize);
  text(1900, -3.6, 'B', 'Color', 'r', 'FontSize', Fsize);

  % F12
  Paxes = subplot(5, 4, 4);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F12(1,:,:))';
  Ptitle = '06Z, 22Aug (F12)';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'd', Ptitle, Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 8);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F12(2,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'h', '21Z, 22Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 12);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F12(3,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'l', '12Z, 23Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 16);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F12(4,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'p', '03Z, 24Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 0, 0);

  Paxes = subplot(5, 4, 20);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(F12(5,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 't', '18Z, 24Aug', Xlab, Xlim, Ylab, Ylim, DiffCmap, DiffClim, DiffClevs, Fsize, 1, 0);
  text(0, -3.6, 'A', 'Color', 'r', 'FontSize', Fsize);
  text(1900, -3.6, 'B', 'Color', 'r', 'FontSize', Fsize);


  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
