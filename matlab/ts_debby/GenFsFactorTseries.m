function [ ] = GenFsFactors()

  Ddir = 'DIAGS';

  % Factor separation method, 2 factors:
  %
  %      Factor     Number    Included      Excluded
  %
  %    SAL-Env        1         SAL          NONSAL
  %    DUST           2         DUST         NODUST
  %
  %    SAL-Env means the existance of a mid-level dry, warm layer
  %
  % Factor separation according to Stein and Alpert (1993). Their
  % notation uses F0,... for simulations and F-hat0,... for factors.
  % We'll use S0,... for simulations and F0,... for factors.
  %
  %      Simulations  SAL-Env     Dust          Case
  %
  %          S0        No          No        TSD_NONSAL_NODUST 
  %          S1        Yes         No        TSD_SAL_NODUST 
  %          S2        No          Yes       TSD_NONSAL_DUST 
  %          S12       Yes         Yes       TSD_SAL_DUST 
  %
  %      Factors              Descrip                     Formula
  %
  %          F0       Part independent of factors           F0 = S0
  %          F1       Part due to SAL-Env                   F1 = S1 - S0
  %          F2       Part due to DUST                      F2 = S2 - S0
  %          F12      Part due to interactions of           F12 = S12 - (S1 + S2) + S0
  %                   SAL-Env and DUST
  
  
  % Cases - list these out in the factor separation order:
  %   S0
  %   S1
  %   S2
  %   S12
  CaseList = {
    'TSD_NONSAL_NODUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
    'TSD_SAL_DUST'
    };

  Ncases = length(CaseList);

  VarList = {
    %     File Name                    Var Name        
    { 'DIAGS/storm_meas_<CASE>.h5'    '/avg_wind'          }
    { 'DIAGS/storm_meas_<CASE>.h5'    '/avg_wind_t'        }
    { 'DIAGS/storm_meas_<CASE>.h5'    '/max_wind'      }
    { 'DIAGS/storm_meas_<CASE>.h5'    '/max_wind_t'    }
    };

  Nvars = length(VarList);

  fprintf('*****************************************************************\n');
  fprintf('Generating factor separation time series:\n');
  fprintf('\n');

  % Create output file. Delete an existing file so that the h5create/h5write
  % statments will create new datasets.
  OutFile = sprintf('%s/fs_factor_time_series.h5', Ddir);
  if (exist(OutFile, 'file') == 2)
    delete(OutFile);
  end

  for ivar = 1:Nvars
    InFilePat = VarList{ivar}{1};
    InVname   = VarList{ivar}{2};

    % S0
    InFile = regexprep(InFilePat, '<CASE>', CaseList{1});
    fprintf('  Reading: %s (%s)\n', InFile, InVname);
    S0 = squeeze(h5read(InFile, InVname));

    % S1
    InFile = regexprep(InFilePat, '<CASE>', CaseList{2});
    fprintf('  Reading: %s (%s)\n', InFile, InVname);
    S1 = squeeze(h5read(InFile, InVname));

    % S2
    InFile = regexprep(InFilePat, '<CASE>', CaseList{3});
    fprintf('  Reading: %s (%s)\n', InFile, InVname);
    S2 = squeeze(h5read(InFile, InVname));

    % S12
    InFile = regexprep(InFilePat, '<CASE>', CaseList{4});
    fprintf('  Reading: %s (%s)\n', InFile, InVname);
    S12 = squeeze(h5read(InFile, InVname));

    % If first iteration, grab and set up the dimension in the output file
    if (ivar == 1)
      Xname = '/x_coords';
      Yname = '/y_coords';
      Zname = '/z_coords';
      Tname = '/t_coords';

      X = squeeze(h5read(InFile, Xname));
      Y = squeeze(h5read(InFile, Yname));
      Z = squeeze(h5read(InFile, Zname));
      T = squeeze(h5read(InFile, Tname));

      CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
      NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
    end

    % Factors
    F0  = S0;
    F1  = S1 - S0;
    F2  = S2 - S0;
    F12 = S12 - (S1 + S2) + S0;

    % Write simulation values and the factors. Place each of these
    % into a vector.
    %   Sim averages: [ S0 S1 S2 S12 ]
    %   Factors:      [ F0 F1 F2 F12 ]
    S_DATA = [ S0 S1 S2 S12 ];
    F_DATA = [ F0 F1 F2 F12 ];

    DimOrder = { 't' };

    OutVname = sprintf('%s_sims', InVname);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVname)
    h5create(OutFile, OutVname,  size(S_DATA));
    h5write( OutFile, OutVname,  S_DATA);
    AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);

    OutVname = sprintf('%s_factors', InVname);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVname)
    h5create(OutFile, OutVname,  size(F_DATA));
    h5write( OutFile, OutVname,  F_DATA);
    AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);

    fprintf('\n');
  end
end 
