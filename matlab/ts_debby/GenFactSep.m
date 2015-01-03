function [ ] = GenFactSep(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;

  % Ignore the first 10 hours of the time series since this contains the
  % end of the rapid intensification phases of the simulated storms.
  T1 = 10;

  VarList = {
    %  Avg Var           File Prefix                 File Var Name
    { 'avg_wind'    'TsAveragedData/max_speed10m'  '/max_speed10m'  }
    { 'avg_press'   'TsAveragedData/min_sea_press' '/min_sea_press' }
    { 'avg_ike'     'TsAveragedData/horiz_ke'      '/horiz_ke'      }
    { 'avg_rmw'     'DIAGS/ts_size'                '/rmw'           }
    { 'avg_pcprate' 'DIAGS/ts_avg_pcprate'         '/avg_pcprate_1' }  % use filtered rates >= 1 mm/h
    { 'hda_pcprate' 'DIAGS/ts_avg_pcprate'         '/hda_pcprate_1' }  % use filtered rates >= 1 mm/h
    };

  Nv = length(VarList);

  % Put all of the results into one output file.
  % If the file exists, remove it so that the HDF5 commands
  % can create a new datasets.
  OutFile = sprintf('%s/factor_sep.h5', Ddir);
  if (exist(OutFile, 'file') == 2)
    delete(OutFile);
  end

  % First create time averages of metrics that will be used for the factor
  % separation. Store these averages into the output file so they can be
  % plotted, etc later on.
  
  for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;

    fprintf('*****************************************************************\n');
    fprintf('Generating averages of metrics for factor separation analysis:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for ivar = 1:Nv
      AvgVar   = VarList{ivar}{1};
      InPrefix = VarList{ivar}{2};
      InVar    = VarList{ivar}{3};

      InFile = sprintf('%s_%s.h5', InPrefix , Case);

      % read in var, form time average, write out result
      fprintf('  Reading: %s (%s)\n', InFile, InVar);
      VAR = squeeze(h5read(InFile, InVar));

      VAR_AVG = nanmean(VAR(T1:end));

      OutVar = sprintf('/%s/%s', Case, AvgVar);
      fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
      h5create(OutFile, OutVar,  size(VAR_AVG));
      h5write( OutFile, OutVar,  VAR_AVG);
    end
    fprintf('\n');

  end

  % Do the factor separation on the four quantities: wind, press, ike, rmw
  %
  %      Factor     Number    Included      Excluded
  %
  %    dry air        1         DRY           MOIST
  %    dusty air      2         DUST         NO_DUST
  %
  %
  %      Simulations  Dry Air    Dusty Air      Case              Label
  %
  %          S12        Yes        Yes       TSD_DRY_DUST           DD
  %          S1         Yes        No        TSD_DRY_NODUST         DN
  %          S2         No         Yes       TSD_MOIST_DUST         MD
  %          S0         No         No        TSD_MOIST_NODUST       MN
  %
  %
  %      Factors              Descrip                     Formula
  %
  %          F0       Part independent of SAL chars.        F0 = S0 = MN
  %          F1       Part due to dry air                   F1 = S1 - S0 = DN - MN
  %          F2       Part due to dusty air                 F2 = S2 - S0 = MD - MN
  %          F12      Part due to combo of dry, dusty air   F12 = S12 - (S1 + S2) + S0 = DD - (DN + MD) + MN
  %
  
  fprintf('*****************************************************************\n');
  fprintf('Generating factor separations:\n');
  fprintf('\n');

  for ivar = 1:Nv
    AvgVar = VarList{ivar}{1};

    % S12 - TSD_DRY_DUST
    InVar = sprintf('/TSD_DRY_DUST/%s', AvgVar);
    DD = h5read(OutFile, InVar);

    % S1 - TSD_DRY_NODUST
    InVar = sprintf('/TSD_DRY_NODUST/%s', AvgVar);
    DN = h5read(OutFile, InVar);

    % S2 - TSD_MOIST_DUST
    InVar = sprintf('/TSD_MOIST_DUST/%s', AvgVar);
    MD = h5read(OutFile, InVar);

    % S0 - TSD_MOIST_NODUST
    InVar = sprintf('/TSD_MOIST_NODUST/%s', AvgVar);
    MN = h5read(OutFile, InVar);
  
    F0  = MN;
    F1  = DN - MN;
    F2  = MD - MN;
    F12 = DD - (DN + MD) + MN;
  
    OutVar = sprintf('/FACTOR_SEP/%s_f0', AvgVar);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
    h5create(OutFile, OutVar,  size(F0));
    h5write( OutFile, OutVar,  F0);
    
    OutVar = sprintf('/FACTOR_SEP/%s_f1', AvgVar);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
    h5create(OutFile, OutVar,  size(F1));
    h5write( OutFile, OutVar,  F1);
  
    OutVar = sprintf('/FACTOR_SEP/%s_f2', AvgVar);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
    h5create(OutFile, OutVar,  size(F2));
    h5write( OutFile, OutVar,  F2);
  
    OutVar = sprintf('/FACTOR_SEP/%s_f12', AvgVar);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
    h5create(OutFile, OutVar,  size(F12));
    h5write( OutFile, OutVar,  F12);

    fprintf('\n');
  end

end 
