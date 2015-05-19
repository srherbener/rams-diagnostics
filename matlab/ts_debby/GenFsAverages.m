function [ ] = GenFsAverages()

  Ddir = 'DIAGS';

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
    'TSD_NONSAL_NODUST'
    };

  Ncases = length(CaseList);

  % Two time periods:
  %    Pre SAL: T = 10 to 30 h (after RI, before encounter SAL)
  %        SAL: T = 40 to 60 h (during SAL)

  PreSalTstart = 10;
  PreSalTend   = 30;
  SalTstart    = 40;
  SalTend      = 60;

  VarList = {
    %  Avg Var           File Name                          File Var Name
    { 'max_wind'    'DIAGS/storm_meas_tseries_<CASE>.h5'    '/max_wind_fa'  }
    { 'min_press'   'DIAGS/storm_meas_tseries_<CASE>.h5'    '/min_slp_fa'   }
    { 'ike'         'DIAGS/storm_meas_tseries_<CASE>.h5'    '/ike'          }
    { 'rmw'         'DIAGS/storm_meas_tseries_<CASE>.h5'    '/rmw_fa'       }

%    { 'avg_pcprate' 'DIAGS/ts_avg_pcprate'         '/avg_pcprate_1' }  % use filtered rates >= 1 mm/h
%    { 'hda_pcprate' 'DIAGS/ts_avg_pcprate'         '/hda_pcprate_1' }  % use filtered rates >= 1 mm/h
    };

  Nvars = length(VarList);

  % Put all of the results into one output file.
  % If the file exists, remove it so that the HDF5 commands
  % can create a new datasets.
  OutFile = sprintf('%s/fs_averages.h5', Ddir);
  if (exist(OutFile, 'file') == 2)
    delete(OutFile);
  end

  % Create time averages of metrics that will be used for the factor separation.
  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating averages of metrics for factor separation analysis:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for ivar = 1:Nvars
      AvgVar   = VarList{ivar}{1};
      InFile   = VarList{ivar}{2};
      InVar    = VarList{ivar}{3};

      InFile = regexprep(InFile, '<CASE>' , Case);

      % read in var, form time average, write out result
      fprintf('  Reading: %s (%s)\n', InFile, InVar);
      VAR = squeeze(h5read(InFile, InVar));
 
      % if on the first var, set up the time indices
      if (ivar == 1)
        T   = (squeeze(h5read(InFile, '/t_coords')) ./ 3600) - 42;  % hours, 0 h -> sim time 42 h

        PS_T1 = find(T >= PreSalTstart, 1, 'first');
        PS_T2 = find(T <= PreSalTend,   1, 'last');

        S_T1 = find(T >= SalTstart, 1, 'first');
        S_T2 = find(T <= SalTend,   1, 'last');
      end

      % calculate the time averages
      PS_AVG = nanmean(VAR(PS_T1:PS_T2));
      S_AVG  = nanmean(VAR(S_T1:S_T2));

      % save the averages (raw data for the factor separation)
      OutVar = sprintf('/%s/ps_avg_%s', Case, AvgVar);
      fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
      h5create(OutFile, OutVar,  size(PS_AVG));
      h5write( OutFile, OutVar,  PS_AVG);

      OutVar = sprintf('/%s/s_avg_%s', Case, AvgVar);
      fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
      h5create(OutFile, OutVar,  size(S_AVG));
      h5write( OutFile, OutVar,  S_AVG);
    end
    fprintf('\n');

  end
end 
