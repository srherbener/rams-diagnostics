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

  % ps_* -> Pre-SAL time period, s_* -> SAL time period
  VarList = {
    %  Avg Var           File Name                          File Var Name

    % Vertical coords -> height
    { 'avg_wind_t'     'DIAGS/fs_averages.h5'    '/<CASE>/avg_wind_t'        }
    { 'max_wind_t'     'DIAGS/fs_averages.h5'    '/<CASE>/avg_max_wind_t'    }
    { 'rmw_t'          'DIAGS/fs_averages.h5'    '/<CASE>/avg_rmw_t'         }
    { 'min_press'      'DIAGS/fs_averages.h5'    '/<CASE>/avg_min_press'     }
    { 'ike'            'DIAGS/fs_averages.h5'    '/<CASE>/avg_ike'           }
    { 'pcprate'        'DIAGS/fs_averages.h5'    '/<CASE>/avg_pcprate'       }

    { 'ps_avg_wind_t'  'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_wind_t'     }
    { 'ps_max_wind_t'  'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_max_wind_t' }
    { 'ps_rmw_t'       'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_rmw_t'      }
    { 'ps_min_press'   'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_min_press'  }
    { 'ps_ike'         'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_ike'        }
    { 'ps_pcprate'     'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_pcprate'    }

    { 's_avg_wind_t'   'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_wind_t'      }
    { 's_max_wind_t'   'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_max_wind_t'  }
    { 's_rmw_t'        'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_rmw_t'       }
    { 's_min_press'    'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_min_press'   }
    { 's_ike'          'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_ike'         }
    { 's_pcprate'      'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_pcprate'     }

    % Vertical coords -> pressure
    { 'avg_wind_t_p'       'DIAGS/fs_averages.h5'    '/<CASE>/avg_wind_t_p'          }
    { 'max_wind_t_p'       'DIAGS/fs_averages.h5'    '/<CASE>/avg_max_wind_t_p'      }
    { 'rmw_t_p'            'DIAGS/fs_averages.h5'    '/<CASE>/avg_rmw_t_p'           }

    { 'ps_avg_wind_t_p'    'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_wind_t_p'       }
    { 'ps_max_wind_t_p'    'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_max_wind_t_p'   }
    { 'ps_rmw_t_p'         'DIAGS/fs_averages.h5'    '/<CASE>/ps_avg_rmw_t_p'        }

    { 's_avg_wind_t_p'     'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_wind_t_p'        }
    { 's_max_wind_t_p'     'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_max_wind_t_p'    }
    { 's_rmw_t_p'          'DIAGS/fs_averages.h5'    '/<CASE>/s_avg_rmw_t_p'         }
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

    % Write simulation values and the factors. Place each of these
    % into a vector.
    %   Sim averages: [ S0 S1 S2 S12 ]
    %   Factors:      [ F0 F1 F2 F12 ]
    %
    % Create a third matrix that shows accumulation of the factors.
    % The first row shows start values of each factor and the second
    % row shows end values. Thus, the magnitude and sign (+/-) information
    % is saved. The first four start/end pairs show the accumulation
    % of the factors and the fifth pair shows the resulting sum.

    S_DATA = [ S0 S1 S2 S12 ];
    F_DATA = [ F0 F1 F2 F12 ];

    P1 = 0;
    P2 = S0;
    P3 = P2 + F1;
    P4 = P3 + F2;
    P5 = P4 + F12;

    BF_DATA = [ P1 P2 P3 P4   0 ;
                P2 P3 P4 P5 S12 ];
    
    OutVar = sprintf('/%s_averages', FsVar);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
    h5create(OutFile, OutVar,  length(S_DATA));
    h5write( OutFile, OutVar,  S_DATA);

    OutVar = sprintf('/%s_factors', FsVar);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
    h5create(OutFile, OutVar,  length(F_DATA));
    h5write( OutFile, OutVar,  F_DATA);

    OutVar = sprintf('/%s_bar_factors', FsVar);
    fprintf('  Writing: %s (%s)\n', OutFile, OutVar)
    h5create(OutFile, OutVar,  size(BF_DATA));
    h5write( OutFile, OutVar,  BF_DATA);

    fprintf('\n');
  end
end 
