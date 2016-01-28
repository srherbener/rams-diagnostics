function [ ] = PlotDpFigDustRainout()
% PlotDpFigDustRainout 

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end

  Fsize = 13;
  
  % Hovmoller data for dust-in-hydrometeors
  InFile = 'DIAGS/storm_hovmollers_TSD_SAL_DUST.h5';
  InVname = '/sal_dust_hydro';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  HOV_DUST = squeeze(h5read(InFile, InVname));
  Z    = squeeze(h5read(InFile, '/z_coords'));
  T    = squeeze(h5read(InFile, '/t_coords'));

  % Do Hovmoller on log scale. Declaring color axis as logrithmic doesn't
  % work right with colorbar markings so instead apply log scaling to
  % data and force the appropriate tick labels on the colorbar.
  HOV_DUST = log10(HOV_DUST);   % logarithmic scaling
  Z = Z ./ 1000;        % km
  T = (T ./ 3600) - 42; % convert to simulation hours

  % Time series of lower level Mdhy
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_dust_hydro_total_mass_llev';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_MDHY = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_MDHY = TS_SAL_MDHY .* 1e-9;  % convert to millions of kg

  % Time series of precip rate
  InFile = 'DIAGS/hist_meas_ts_pcprate_TSD_SAL_DUST.h5';
  InVname = '/sal_pcprate_ts';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_PR = squeeze(h5read(InFile, InVname)); % mm/h

  % Time series of surface dust accumulation
  InFile = 'DIAGS/total_mass_TSD_SAL_DUST.h5';
  InVname = '/sal_dust_sfc_total_mass';
  fprintf('Reading: %s (%s)\n', InFile, InVname);

  TS_SAL_DSFC = squeeze(h5read(InFile, InVname)); % g
  TS_SAL_DSFC = TS_SAL_DSFC .* 1e-12;  % convert to billions of kg

  % plot, 4 panels (2x2)
  OutFile = sprintf('%s/DpFigDustRainout.jpg', Pdir);
  
  Fig = figure;

  % create the Hovmoller of sal_dust_hydro
  Paxes = subplot(4,1,1);
  PlotDustHov(Paxes, T, Z, HOV_DUST, 'a', '', Fsize, 0, 1);
  
  Paxes = subplot(4,1,2);
  PlotDpTseries(Paxes, T, TS_SAL_MDHY, 'b', '', 'M_d_h_y (10^6 kg)', Fsize, 0, 1);

  Paxes = subplot(4,1,3);
  PlotDpTseries(Paxes, T, TS_SAL_PR, 'c', '', 'PR (mm h^-^1)', Fsize, 0, 1);

  Paxes = subplot(4,1,4);
  PlotDpTseries(Paxes, T, TS_SAL_DSFC, 'd', '', 'M_d (10^9 kg)', Fsize, 1, 1);

  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlotDustHov(Paxes, X, Y, Z, Pmarker, Ptitle, Fsize, ShowX, ShowY)

  axes(Paxes);

  contourf(X, Y, Z, 20, 'LineStyle', 'none');

  % place colorbar beneath contour plot, but don't allow matlab to put in large gaps
  GapSize = 0.05;
  PaxesLoc = get(Paxes, 'Position'); % open up a small gap below plot
  CbarLoc = PaxesLoc;
  PaxesLoc(2) = PaxesLoc(2) + GapSize;
  PaxesLoc(4) = PaxesLoc(4) - GapSize;
  set(Paxes, 'Position', PaxesLoc);
  
  CbarLoc(4) = GapSize * 0.8;  % Force colorbar into the above gap
  Cbar = colorbar('Location', 'SouthOutside', 'Position', CbarLoc);
  caxis([ -3 0 ]);
  ylim([ 0 6 ]);

  set(Cbar, 'Ticks', [ -2 -1 0 ]);
  set(Cbar, 'TickLabels', { '10^-^2' '10^-^1' '1' });

  set(Paxes, 'FontSize', Fsize);
  set(Paxes, 'LineWidth', 2);
  set(Paxes, 'TickLength', [ 0.025 0.025 ]);

  set(Paxes, 'XTick', [ 6 30 54 ]);
  if (ShowX > 0)
    set(Paxes, 'XTickLabel', { ' 12Z\newline22Aug' ' 12Z\newline23Aug' ' 12Z\newline24Aug' });
  else
    set(Paxes, 'XTickLabel', {});
  end

  set(Paxes, 'YTick', [ 0 5 10 15 ]);
  if (ShowY > 0)
    ylabel('Height (km)');
  else
    set(Paxes, 'YTickLabel', {});
  end

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlotDpTseries(Paxes, X, Y, Pmarker, Ptitle, Ylabel, Fsize, ShowX, ShowY)

  axes(Paxes);

  plot(X, Y, 'LineWidth', 2);

  set(Paxes, 'FontSize', Fsize);
  set(Paxes, 'LineWidth', 2);
  set(Paxes, 'TickLength', [ 0.025 0.025 ]);

  set(Paxes, 'XTick', [ 6 30 54 ]);
  if (ShowX > 0)
    set(Paxes, 'XTickLabel', { ' 12Z\newline22Aug' ' 12Z\newline23Aug' ' 12Z\newline24Aug' });
  else
    set(Paxes, 'XTickLabel', {});
  end

  ylabel(Ylabel);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

end
