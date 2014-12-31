function [ ] = PlotFactSep(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);
  
  Ddir = Config.DiagDir;

  Pdir = Config.PlotDir;
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 22;
  LegendFsize = 15;

  InFile = sprintf('%s/factor_sep.h5', Ddir);

  fprintf('Generating plots for factor separation analysis\n');
  fprintf('\n');

  SimList = {
    { 'DD' 'TSD_DRY_DUST' }
    { 'DN' 'TSD_DRY_NODUST' }
    { 'MD' 'TSD_MOIST_DUST' }
    { 'MN' 'TSD_MOIST_NODUST' }
    };
  Ns = length(SimList);
  
  VarList = { 'avg_wind' 'avg_press' 'avg_ike' 'avg_rmw' };
  Nv = length(VarList);

  FACT_LABELS = { 'NO' 'DA' 'DU' 'NC' };

  for ivar = 1:Nv
    VarName = VarList{ivar};

    % read in the sim averages
    for isim = 1:Ns
      SimLabel = SimList{isim}{1}; 
      SimCase  = SimList{isim}{2}; 
      
      BPLOT_LABELS{isim} = SimLabel;
      
      InVar  = sprintf('/%s/%s', SimCase, VarName);
      fprintf('  Reading: %s (%s)\n', InFile, InVar);
      SIM_AVGS(isim) = h5read(InFile, InVar);
    end
    fprintf('\n');

    % The MN (4th) entry is the simulation without the SAL features.
    % Covert the entries to percent changes from the
    % fourth (MN) entry for making the plot easier to read.
    SIM_PCT_CHNG = ((SIM_AVGS - SIM_AVGS(4)) ./ SIM_AVGS(4)) .* 100;

    % read in the factor separations
    InVar  = sprintf('/FACTOR_SEP/%s_f0', VarName);
    fprintf('  Reading: %s (%s)\n', InFile, InVar);
    F0 = h5read(InFile, InVar);

    InVar  = sprintf('/FACTOR_SEP/%s_f1', VarName);
    fprintf('  Reading: %s (%s)\n', InFile, InVar);
    F1 = h5read(InFile, InVar);

    InVar  = sprintf('/FACTOR_SEP/%s_f2', VarName);
    fprintf('  Reading: %s (%s)\n', InFile, InVar);
    F2 = h5read(InFile, InVar);

    InVar  = sprintf('/FACTOR_SEP/%s_f12', VarName);
    fprintf('  Reading: %s (%s)\n', InFile, InVar);
    F12 = h5read(InFile, InVar);

    FS_DATA = [ F0 F1 F2 F12 ];
    % Express as percent changes from MN, note that
    % factors are already differences.
    FACT_PCT_CHNG = (FS_DATA ./ SIM_AVGS(4)) .* 100;

    fprintf('\n');

    % Create the plots, 2 panels, left - sim data, right - factor separation
    Fig = figure;

    OutFile = sprintf('%s/fs_%s.jpg', Pdir, VarName);

    % Put sims and factors in the same order: dry air, dust, dry + dust
    SimOrder = [ 2 3 1 ];
    FactOrder = [ 2 3 4 ];

    subplot(1,2,1);
    bar(SIM_PCT_CHNG(SimOrder));
    set(gca, 'FontSize', Fsize);
    set(gca, 'XTickLabel', BPLOT_LABELS(SimOrder));
    title ('Change from MN');
    ylabel('Percent');

    subplot(1,2,2);
    bar(FACT_PCT_CHNG(FactOrder));
    set(gca, 'FontSize', Fsize);
    set(gca, 'XTickLabel', FACT_LABELS(FactOrder));
    title('SAL Factors');

    saveas(Fig, OutFile);
    close(Fig);
    fprintf('\n');
  end

end
