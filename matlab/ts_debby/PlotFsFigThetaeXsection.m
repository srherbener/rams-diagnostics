function [ ] = PlotFsFigThetaeXsection()

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

  % Read in the theta-e slices
  %   x-z planes where x is along the x-section line
  for icase = 1:Ncases
    Case  = CaseList{icase}{1};
    Label = CaseList{icase}{2};

    InFname = sprintf('XsectionData/strack_theta_e_%s.h5', Case);
    InVname = '/theta_e';
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
      THETA_E(icase,it,:,:) = squeeze(h5read(InFname, InVname, [ 1 1 Z1 T1 ], [ Nx 1 Nz 1 ]));
    end
  end
  fprintf('\n');

  % Plot: 10 panels (5x2)
  OutFile = sprintf('%s/FsFigThetaeXsection.jpg', Pdir);
  Fig = figure;

  Xlim = [ 0 2150 ];
  Ylim = [ 0 5.5 ];
  Clim = [ 330 357 ];

  PlocInc = 0.05;

  % TSD_SAL_DUST
  Paxes = subplot(5, 2, 1);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_E(1,1,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'a', '06Z, 22Aug (SD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Clim, Fsize, 0, 1);

  Paxes = subplot(5, 2, 3);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_E(1,2,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'c', '21Z, 22Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Clim, Fsize, 0, 1);

  Paxes = subplot(5, 2, 5);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_E(1,3,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'e', '12Z, 23Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Clim, Fsize, 0, 1);

  Paxes = subplot(5, 2, 7);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_E(1,4,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'g', '03Z, 24Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Clim, Fsize, 0, 1);

  Paxes = subplot(5, 2, 9);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_E(1,5,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'i', '18Z, 24Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Clim, Fsize, 1, 1);

  % TSD_NONSAL_DUST
  Paxes = subplot(5, 2, 2);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_E(2,1,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'b', '06Z, 22Aug (NSD)', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Clim, Fsize, 0, 0);

  Paxes = subplot(5, 2, 4);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_E(2,2,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'd', '21Z, 22Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Clim, Fsize, 0, 0);

  Paxes = subplot(5, 2, 6);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_E(2,3,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'f', '12Z, 23Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Clim, Fsize, 0, 0);

  Paxes = subplot(5, 2, 8);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_E(2,4,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'h', '03Z, 24Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Clim, Fsize, 0, 0);

  Paxes = subplot(5, 2, 10);
  Ploc = get(Paxes, 'Position');
  Ploc(1) = Ploc(1) - PlocInc;
  Ploc(3) = Ploc(3) + PlocInc;
  set(Paxes, 'Position', Ploc);
  PDATA = squeeze(THETA_E(2,5,:,:))';
  PlotFsFigXsection(Paxes, X, Z, PDATA, 'j', '18Z, 24Aug', 'Linear Distance (km)', Xlim, 'Z (km)', Ylim, Clim, Fsize, 1, 0);


  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
