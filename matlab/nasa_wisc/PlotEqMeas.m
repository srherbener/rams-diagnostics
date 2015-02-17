function [ ] = PlotEqMeas(ConfigFile)
% PlotEqMeas plot RCE measurements

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;
  Pdir = Config.PlotDir;

  % make sure output directory exists
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  % case_name temp_var_name cloud_var_name vapor_var_name select_time (hrs) plot_title
  PlotSets = {
      { 'RCE50_RECT'      'eq_meas' 'therm_heat_flux' 'rad_flux_div' 'SZA 50'       'EqMeas' }
      { 'RCE50_RECT_S300' 'eq_meas' 'therm_heat_flux' 'rad_flux_div' 'SZA 50, S300' 'EqMeas' }
      { 'RCE50_RECT_S303' 'eq_meas' 'therm_heat_flux' 'rad_flux_div' 'SZA 50, S303' 'EqMeas' }
      { 'RCE50_SQ'        'eq_meas' 'therm_heat_flux' 'rad_flux_div' 'SZA 50, SQ'   'EqMeas' }
      { 'RCE_MATT'        'eq_meas' 'therm_heat_flux' 'rad_flux_div' 'MATT'         'EqMeas' }
    };
  Nset = length(PlotSets);


  fprintf('***************************************************************\n');
  fprintf('Generating equillibrium measurement plots:\n');

  for iset = 1:Nset
    Case       = PlotSets{iset}{1};
    Fprefix    = PlotSets{iset}{2};
    ThfVname   = PlotSets{iset}{3};
    QradVname  = PlotSets{iset}{4};
    Ptitle     = PlotSets{iset}{5};
    OutFprefix = PlotSets{iset}{6};

    InFile = sprintf('%s/%s_%s.h5', Ddir, Fprefix, Case);
    OutFile = sprintf('%s/%s_%s.jpg', Pdir, OutFprefix, Case);
    ThfVname = sprintf('/%s', ThfVname);
    QradVname = sprintf('/%s', QradVname);

    fprintf('  Case: %s\n', Case);

    % Read in TEMP, CLOUD, and VAPOR data
    % After squeezing, vars will be: (z,t)
    fprintf('    Reading: %s (%s, %s)\n', InFile, ThfVname, QradVname);

    THF  = squeeze(h5read(InFile, ThfVname));
    QRAD = squeeze(h5read(InFile, QradVname));
    T = squeeze(hdf5read(InFile, '/t_coords'));  % days

    % Create the plot
    PDATA = [ THF QRAD ];
    LegendText = { 'THF' 'Q_r_a_d' };
    Fsize = 20;

    Fig = figure;

    plot(T, PDATA, 'LineWidth', 1.5);
    set(gca, 'FontSize', Fsize);
    legend(LegendText);
    title(Ptitle);
    xlabel('Simulation Time (d)');
    ylabel('Flux (W m^2)');

    fprintf('    Writing: %s\n', OutFile);
    fprintf('\n');

    saveas(Fig, OutFile);
    close(Fig);

    % Generate a plot of the components of the Qrad measurement
    fprintf('    Reading Qrad measurement components\n');
    SFC_ALB  = squeeze(h5read(InFile, '/avg_sfc_alb'));
    SFC_LWDN = squeeze(h5read(InFile, '/avg_sfc_lwdn'));
    SFC_LWUP = squeeze(h5read(InFile, '/avg_sfc_lwup'));
    SFC_SWDN = squeeze(h5read(InFile, '/avg_sfc_swdn'));
    TOP_LWUP = squeeze(h5read(InFile, '/avg_top_lwup'));
    TOP_SWDN = squeeze(h5read(InFile, '/avg_top_swdn'));
    TOP_SWUP = squeeze(h5read(InFile, '/avg_top_swup'));

    SFC_SWUP = SFC_ALB .* SFC_SWDN;

    PDATA = [ SFC_SWDN SFC_LWDN SFC_SWUP SFC_LWUP TOP_SWDN TOP_SWUP TOP_LWUP ];
    LegendText = { 'SFC\_SWDN' 'SFC\_LWDN' 'SFC\_SWUP' 'SFC\_LWUP' 'TOP\_SWDN' 'TOP\_SWUP' 'TOP\_LWUP' };

    Fig = figure;

    PtitleComp = sprintf('%s Qrad Components', Ptitle);
    OutFileComp = sprintf('%s/%s_%s_QradComp.jpg', Pdir, OutFprefix,Case);

    plot(T, PDATA, 'LineWidth', 1.5);
    set(gca, 'FontSize', Fsize);
    legend(LegendText);
    title(PtitleComp);
    xlabel('Simulation Time (d)');
    ylabel('Flux (W m^2)');

    fprintf('    Writing: %s\n', OutFileComp);
    fprintf('\n');

    saveas(Fig, OutFileComp);
    close(Fig);

    % Generate a plot of the components of the THF measurement
    fprintf('    Reading THF measurement components\n');
    SFC_LAT  = squeeze(h5read(InFile, '/avg_sfc_lat'));
    SFC_SENS = squeeze(h5read(InFile, '/avg_sfc_sens'));

    PDATA = [ SFC_LAT SFC_SENS ];
    LegendText = { 'SFC\_LAT' 'SFC\_SENS' };

    Fig = figure;

    PtitleComp = sprintf('%s THF Components', Ptitle);
    OutFileComp = sprintf('%s/%s_%s_ThfComp.jpg', Pdir, OutFprefix,Case);

    plot(T, PDATA, 'LineWidth', 1.5);
    set(gca, 'FontSize', Fsize);
    legend(LegendText);
    title(PtitleComp);
    xlabel('Simulation Time (d)');
    ylabel('Flux (W m^2)');

    fprintf('    Writing: %s\n', OutFileComp);
    fprintf('\n');

    saveas(Fig, OutFileComp);
    close(Fig);
  end
end
