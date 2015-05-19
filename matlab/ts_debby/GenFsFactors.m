function [ ] = GenFactSep()

  Ddir = 'DIAGS';

  % Factor separation method, 2 factors:
  %
  %      Factor     Number    Included      Excluded
  %
  %    "sal" air      1         SAL          NONSAL
  %    dusty air      2         DUST         NODUST
  %
  %    sal air means the existance of a mid-level dry, warm layer
  %
  %
  %      Simulations  Sal Air    Dusty Air      Case
  %
  %          S0        Yes         Yes       TSD_SAL_DUST 
  %          S1        No          Yes       TSD_NONSAL_DUST 
  %          S2        Yes         No        TSD_SAL_NODUST 
  %          S12       No          No        TSD_NONSAL_NODUST 
  %
  %      Factors              Descrip                     Formula
  %
  %          F0       Part independent of SAL chars.        F0 = S0
  %          F1       Part due to dry air                   F1 = S1 - S0
  %          F2       Part due to dusty air                 F2 = S2 - S0
  %          F12      Part due to combo of dry, dusty air   F12 = S12 - (S1 + S2) + S0
  
  % Cases - list these out in the factor separation order:
  %   S0
  %   S1
  %   S2
  %   S12
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_NONSAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_NODUST'
    };

  Ncases = length(CaseList);

  % ps_* -> Pre-SAL time period, s_* -> SAL time period
  VarList = {
    %  Avg Var           File Name                          File Var Name
    { 'ps_max_wind'    'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_max_wind'  }
    { 'ps_min_press'   'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_min_press' }
    { 'ps_ike'         'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_ike'       }
    { 'ps_rmw'         'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_rmw'       }

    { 's_max_wind'     'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_max_wind'   }
    { 's_min_press'    'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_min_press'  }
    { 's_ike'          'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_ike'        }
    { 's_rmw'          'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_rmw'        }
    };

  Nvars = length(VarList);

  fprintf('*****************************************************************\n');
  fprintf('Generating factor separations:\n');
  fprintf('\n');

  % Create output file. Delete an existing file so that the h5create/h5write
  % statments will create new datasets.
  OutFile = sprintf('%s/fs_factors.h5', Ddir);
  if (exist(OutFile, 'file') == 2)
    delete(OutFile);
  end

  for ivar = 1:Nvars
    FsVar    = VarList{ivar}{1};
    InFile   = VarList{ivar}{2};
    InVarPat = VarList{ivar}{3};

    % S0
    InVar = regexprep(InVarPat, '<CASE>', CaseList{1});
    fprintf('  Reading: %s (%s)\n', InFile, InVar);
    S0 = h5read(InFile, InVar);

    % S1
    InVar = regexprep(InVarPat, '<CASE>', CaseList{2});
    fprintf('  Reading: %s (%s)\n', InFile, InVar);
    S1 = h5read(InFile, InVar);

    % S2
    InVar = regexprep(InVarPat, '<CASE>', CaseList{3});
    fprintf('  Reading: %s (%s)\n', InFile, InVar);
    S2 = h5read(InFile, InVar);

    % S12
    InVar = regexprep(InVarPat, '<CASE>', CaseList{4});
    fprintf('  Reading: %s (%s)\n', InFile, InVar);
    S12 = h5read(InFile, InVar);
  
    % Factors
    F0  = S0;
    F1  = S1 - S0;
    F2  = S2 - S0;
    F12 = S12 - (S1 + S2) + S0;
  
    % Write factors into output file
    OutVar = sprintf('/%s_f0', FsVar);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
    h5create(OutFile, OutVar,  size(F0));
    h5write( OutFile, OutVar,  F0);
    
    OutVar = sprintf('/%s_f1', FsVar);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
    h5create(OutFile, OutVar,  size(F1));
    h5write( OutFile, OutVar,  F1);
  
    OutVar = sprintf('/%s_f2', FsVar);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
    h5create(OutFile, OutVar,  size(F2));
    h5write( OutFile, OutVar,  F2);
  
    OutVar = sprintf('/%s_f12', FsVar);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
    h5create(OutFile, OutVar,  size(F12));
    h5write( OutFile, OutVar,  F12);

    fprintf('\n');
  end
end 
