function [ ] = PlotFsFigPcprateXsection()

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

  % Read in the pcprate slices
  %   x-z planes where x is along the x-section line
  for icase = 1:Ncases
    Case  = CaseList{icase}{1};
    Labels{icase} = CaseList{icase}{2};

    InFname = sprintf('XsectionData/strack_pcprate_%s.h5', Case);
    InVname = '/pcprate';
    fprintf('Reading %s (%s)\n', InFname, InVname);

    % 
    if (icase == 1)
      X = squeeze(h5read(InFname, '/x_coords')); % km
      Z = squeeze(h5read(InFname, '/z_coords'))./1000; % convert to km
      T = squeeze(h5read(InFname, '/t_coords'))./3600 - 42; % convert to sim time in hours starting with zero

      Nx = length(X);
    end
     
    % for each time, extract the x-z slice
    for it = 1:Nst
      T1 = find(T >= SimTimes(it), 1, 'first');
      PCPRATE(icase,it,:) = squeeze(h5read(InFname, InVname, [ 1 1 1 T1 ], [ Nx 1 1 1 ]));
    end
  end
  fprintf('\n');

  % Calculate factors:
  %   F1 = SND - NSND
  %   F2 = NSD - NSND
  %   F12 = SD - (SND+NSD) + NSND
  F1  = squeeze(PCPRATE(2,:,:) - PCPRATE(1,:,:));
  F2  = squeeze(PCPRATE(3,:,:) - PCPRATE(1,:,:));
  F12 = squeeze(PCPRATE(4,:,:) -(PCPRATE(2,:,:)+PCPRATE(3,:,:)) + PCPRATE(1,:,:));

  % Plot: 10 panels (5x2)
  OutFile = sprintf('%s/FsFigPcprateXsection.jpg', Pdir);
  Fig = figure;

  Xlim = [ 0 2150 ];
  Ylim = [ 0 11 ];

  DiffYlim = [ -4 4 ];

  Xlab = 'Distance (km)';
  Ylab = 'PR (mm h^-^1)';

  PlocInc = 0.05;

  Xscale = 'linear';
  Yscale = 'linear';

  LegText = { '' };
  LegLoc = 'none';
  Pcolors = { 'black' };

  % CONTROL (NSND)
  Paxes = subplot(5, 4, 1);
  PDATA = squeeze(PCPRATE(1,1,:));
  Ptitle = sprintf('06Z, 22Aug (%s)', Labels{1});
  PlotFsFigLine(Paxes, X, PDATA, 'a', Ptitle, Xlab, Xlim, Xscale, 0, Ylab, Ylim, Yscale, 1, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 5);
  PDATA = squeeze(PCPRATE(1,2,:));
  PlotFsFigLine(Paxes, X, PDATA, 'e', '21Z, 22Aug', Xlab, Xlim, Xscale, 0, Ylab, Ylim, Yscale, 1, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 9);
  PDATA = squeeze(PCPRATE(1,3,:));
  PlotFsFigLine(Paxes, X, PDATA, 'i', '12Z, 23Aug', Xlab, Xlim, Xscale, 0, Ylab, Ylim, Yscale, 1, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 13);
  PDATA = squeeze(PCPRATE(1,4,:));
  PlotFsFigLine(Paxes, X, PDATA, 'm', '03Z, 24Aug', Xlab, Xlim, Xscale, 0, Ylab, Ylim, Yscale, 1, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 17);
  PDATA = squeeze(PCPRATE(1,5,:));
  PlotFsFigLine(Paxes, X, PDATA, 'q', '18Z, 24Aug', Xlab, Xlim, Xscale, 1, Ylab, Ylim, Yscale, 1, Fsize, LegText, LegLoc, Pcolors);
  text(0, -3.6, 'A', 'Color', 'r', 'FontSize', Fsize);
  text(1900, -3.6, 'B', 'Color', 'r', 'FontSize', Fsize);

  % F1
  Paxes = subplot(5, 4, 2);
  PDATA = squeeze(F1(1,:))';
  Ptitle = '06Z, 22Aug (F1)';
  PlotFsFigLine(Paxes, X, PDATA, 'b', Ptitle, Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 1, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 6);
  PDATA = squeeze(F1(2,:))';
  PlotFsFigLine(Paxes, X, PDATA, 'f', '21Z, 22Aug', Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 1, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 10);
  PDATA = squeeze(F1(3,:))';
  PlotFsFigLine(Paxes, X, PDATA, 'j', '12Z, 23Aug', Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 1, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 14);
  PDATA = squeeze(F1(4,:))';
  PlotFsFigLine(Paxes, X, PDATA, 'n', '03Z, 24Aug', Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 1, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 18);
  PDATA = squeeze(F1(5,:)');
  PlotFsFigLine(Paxes, X, PDATA, 'r', '18Z, 24Aug', Xlab, Xlim, Xscale, 1, Ylab, DiffYlim, Yscale, 1, Fsize, LegText, LegLoc, Pcolors);
  text(0, -3.6, 'A', 'Color', 'r', 'FontSize', Fsize);
  text(1900, -3.6, 'B', 'Color', 'r', 'FontSize', Fsize);

  % F2
  Paxes = subplot(5, 4, 3);
  PDATA = squeeze(F2(1,:))';
  Ptitle = '06Z, 22Aug (F2)';
  PlotFsFigLine(Paxes, X, PDATA, 'c', Ptitle, Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 0, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 7);
  PDATA = squeeze(F2(2,:))';
  PlotFsFigLine(Paxes, X, PDATA, 'g', '21Z, 22Aug', Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 0, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 11);
  PDATA = squeeze(F2(3,:))';
  PlotFsFigLine(Paxes, X, PDATA, 'k', '12Z, 23Aug', Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 0, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 15);
  PDATA = squeeze(F2(4,:))';
  PlotFsFigLine(Paxes, X, PDATA, 'o', '03Z, 24Aug', Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 0, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 19);
  PDATA = squeeze(F2(5,:)');
  PlotFsFigLine(Paxes, X, PDATA, 's', '18Z, 24Aug', Xlab, Xlim, Xscale, 1, Ylab, DiffYlim, Yscale, 0, Fsize, LegText, LegLoc, Pcolors);
  text(0, -3.6, 'A', 'Color', 'r', 'FontSize', Fsize);
  text(1900, -3.6, 'B', 'Color', 'r', 'FontSize', Fsize);

  % F12
  Paxes = subplot(5, 4, 4);
  PDATA = squeeze(F12(1,:))';
  Ptitle = '06Z, 22Aug (F12)';
  PlotFsFigLine(Paxes, X, PDATA, 'd', Ptitle, Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 0, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 8);
  PDATA = squeeze(F12(2,:))';
  PlotFsFigLine(Paxes, X, PDATA, 'h', '21Z, 22Aug', Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 0, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 12);
  PDATA = squeeze(F12(3,:))';
  PlotFsFigLine(Paxes, X, PDATA, 'l', '12Z, 23Aug', Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 0, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 16);
  PDATA = squeeze(F12(4,:))';
  PlotFsFigLine(Paxes, X, PDATA, 'p', '03Z, 24Aug', Xlab, Xlim, Xscale, 0, Ylab, DiffYlim, Yscale, 0, Fsize, LegText, LegLoc, Pcolors);

  Paxes = subplot(5, 4, 20);
  PDATA = squeeze(F12(5,:)');
  PlotFsFigLine(Paxes, X, PDATA, 't', '18Z, 24Aug', Xlab, Xlim, Xscale, 1, Ylab, DiffYlim, Yscale, 0, Fsize, LegText, LegLoc, Pcolors);
  text(0, -3.6, 'A', 'Color', 'r', 'FontSize', Fsize);
  text(1900, -3.6, 'B', 'Color', 'r', 'FontSize', Fsize);


  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);

end
