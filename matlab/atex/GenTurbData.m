function [ ] = GenTurbData(ConfigFile)
% GenTurbData generate w moment and w flux data from turb_* function output
% from tsavg

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;
    Tdir = Config.TsavgDir;

    % These lists are organized as 2D cell arrays. Each row is one complete spec for one variable.
    % Syntax for rows:
    %    { 'input file prefix' 'input file var name' 'input var  number' 'output var name' }
    VarList = {
      % means
      { 'turb_mmts_w'               'turb_mmts_w'       1 'w'                }
      { 'turb_cov_w_theta'          'turb_cov_theta'    2 'theta'            }
      { 'turb_mmts_theta_e'         'turb_mmts_theta_e' 1 'theta_e'          }
      { 'turb_cov_w_theta_v'        'turb_cov_theta_v'  2 'theta_v'          }
      { 'turb_cov_w_vapor'          'turb_cov_vapor'    2 'vapor'            }
      { 'turb_cov_w_speed'          'turb_cov_speed'    2 'speed'            }
    
      { 'turb_mmts_w_ud0p10'        'turb_mmts_w'       1 'w_ud0p10'         }
      { 'turb_cov_w_theta_ud0p10'   'turb_cov_theta'    2 'theta_ud0p10'     }
      { 'turb_mmts_theta_e_ud0p10'  'turb_mmts_theta_e' 1 'theta_e_ud0p10'   }
      { 'turb_cov_w_theta_v_ud0p10' 'turb_cov_theta_v'  2 'theta_v_ud0p10'   }
      { 'turb_cov_w_vapor_ud0p10'   'turb_cov_vapor'    2 'vapor_ud0p10'     }
      { 'turb_cov_w_speed_ud0p10'   'turb_cov_speed'    2 'speed_ud0p10'     }

      % fluxes (covariances)
      { 'turb_cov_w_theta'          'turb_cov_theta'    3 'w-theta'          }
      { 'turb_cov_w_theta_v'        'turb_cov_theta_v'  3 'w-theta_v'        }
      { 'turb_cov_w_vapor'          'turb_cov_vapor'    3 'w-vapor'          }
      { 'turb_cov_w_speed'          'turb_cov_speed'    3 'w-speed'          }
      
      { 'turb_cov_w_theta_ud0p10'   'turb_cov_theta'    3 'w-theta_ud0p10'   }
      { 'turb_cov_w_theta_v_ud0p10' 'turb_cov_theta_v'  3 'w-theta_v_ud0p10' }
      { 'turb_cov_w_vapor_ud0p10'   'turb_cov_vapor'    3 'w-vapor_ud0p10'   }
      { 'turb_cov_w_speed_ud0p10'   'turb_cov_speed'    3 'w-speed_ud0p10'   }

      % variances
      { 'turb_mmts_w'               'turb_mmts_w'       2 'w-w'              }

      { 'turb_mmts_w_ud0p10'        'turb_mmts_w'       2 'w-w_ud0p10'       }

      % skews
      { 'turb_mmts_w'               'turb_mmts_w'       3 'w-w-w'            }

      { 'turb_mmts_w_ud0p10'        'turb_mmts_w'       3 'w-w-w_ud0p10'     }
      
      };

    % do one hour time average at beginning, middle and end of sampling period
    % also do time average over entire sampling inerval
    TSstart = 12;
    TSend = 13;
    TMstart = 23.5;
    TMend = 24.5;
    TEstart = 35;
    TEend = 36;

    TAstart = 12;
    TAend = 36;

    for icase = 1:length(Config.Cases)
      Case = Config.Cases(icase).Cname;
      OutFname = sprintf('%s/turb_stats_%s.h5', Ddir, Case);

      fprintf('***************************************************************\n');
      fprintf('Generating moment/flux profiles:\n');
      fprintf('  Case: %s\n', Case);
      fprintf('  Output file: %s\n', OutFname);
      fprintf('\n');

      % The vars and coords are written into the file in append mode so write 
      % the file here without append mode in order to create a new file each time
      % this script is run.
      hdf5write(OutFname, 'header', Case);

      for ivar = 1:length(VarList)
        InFprefix = VarList{ivar}{1};
        InVname   = VarList{ivar}{2};
        InOrder   = VarList{ivar}{3};
        OutVname  = VarList{ivar}{4};

        InFname = sprintf('%s/%s_%s.h5', Tdir, InFprefix, Case);

        fprintf('  Input file: %s\n', InFname);
        fprintf('    Var name: %s\n', InVname);
        fprintf('    Var order: %d\n', InOrder);
        fprintf('\n');
        
        % the number of points is always the 4th index of the first
        % dimension
        TURB_DATA = squeeze(hdf5read(InFname, InVname));
        SUM = squeeze(TURB_DATA(InOrder,:,:,:));
        N   = squeeze(TURB_DATA(4,:,:,:));
        Z   = squeeze(hdf5read(InFname, 'z_coords'));
        T   = squeeze(hdf5read(InFname, 't_coords')) / 3600;   % hr
        
        Nz = length(Z);
        Nt = length(T);

        % find the indices of the start, mid, end and "all" time intervals
        TS1 = find(T >= TSstart, 1, 'first');
        TS2 = find(T <= TSend, 1, 'last');

        TM1 = find(T >= TMstart, 1, 'first');
        TM2 = find(T <= TMend, 1, 'last');

        TE1 = find(T >= TEstart, 1, 'first');
        TE2 = find(T <= TEend, 1, 'last');

        TA1 = find(T >= TAstart, 1, 'first');
        TA2 = find(T <= TAend, 1, 'last');

        % N and SUM are organized as (z,t)
        % Create the entire time series profiles, then average over time to
        % get the start, middle, end and all profiles.
        PROF = SUM ./ N;
        
        PROF_S = squeeze(mean(PROF(:,TS1:TS2),2));
        PROF_M = squeeze(mean(PROF(:,TM1:TM2),2));
        PROF_E = squeeze(mean(PROF(:,TE1:TE2),2));
        PROF_A = squeeze(mean(PROF(:,TA1:TA2),2));

        % Write out data - put in dummy x, y and t coordinates
        Xdummy = 1;
        Ydummy = 1;
        Tdummy = 1;

        OutVar = reshape(PROF_S, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_start', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(PROF_M, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_mid', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(PROF_E, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_end', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(PROF_A, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_all', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        fprintf('\n');
      end

      % all vars are written out to the file, now write out the coordinates
      hdf5write(OutFname, 'x_coords', Xdummy, 'WriteMode', 'append');
      hdf5write(OutFname, 'y_coords', Ydummy, 'WriteMode', 'append');
      hdf5write(OutFname, 'z_coords', Z,      'WriteMode', 'append');
      hdf5write(OutFname, 't_coords', Tdummy, 'WriteMode', 'append');
    end
