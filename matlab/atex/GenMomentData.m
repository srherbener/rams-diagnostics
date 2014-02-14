function [ ] = GenMomentData(ConfigFile)
% GenMomentData generate w moment and w flux data

    % Read the config file to get the structure of how the data is laid out in
    % the file system.
    [ Config ] = ReadConfig(ConfigFile);

    Ddir = Config.DiagDir;

    % These lists are organized as 2D cell arrays. Each row is one complete spec for one variable.
    % Syntax for rows:
    %    { 'input file prefix' 'input file var name' 'input var term number' 'input var order number' 'output var name' }
    VarList = {
      % means
      { 'w_M3'           'w-w-w'     1 1 'w'           }
      { 'w_theta_flux'   'w-theta'   2 1 'theta'       }
      { 'theta_e_M1'     'theta_e'   1 1 'theta_e'     }
      { 'w_theta_v_flux' 'w-theta_v' 2 1 'theta_v'     }
      { 'w_vapor_flux'   'w-vapor'   2 1 'vapor'       }
      { 'w_speed_flux'   'w-speed'   2 1 'speed'       }
      { 'cloud_M1'       'cloud'     1 1 'cloud'       }
      { 'cloud_M1_c0p01' 'cloud'     1 1 'cloud_c0p01' }
      { 'cloud_M1_c0p10' 'cloud'     1 1 'cloud_c0p10' }

      % fluxes (covariances)
      { 'w_theta_flux'   'w-theta'   1 2 'w-theta'     }
      { 'w_theta_v_flux' 'w-theta_v' 1 2 'w-theta_v'   }
      { 'w_vapor_flux'   'w-vapor'   1 2 'w-vapor'     }
      { 'w_speed_flux'   'w-speed'   1 2 'w-speed'     }

      % variances
      { 'w_M3'           'w-w-w'     1 2 'w-w'         }

      % skews
      { 'w_M3'           'w-w-w'     1 3 'w-w-w'       }
      };

    Tstart = 12;
    Tmid   = 24;
    Tend   = 36;


    for icase = 1:length(Config.Cases)
      Case = Config.Cases(icase).Cname;
      OutFname = sprintf('%s/moments_%s.h5', Ddir, Case);

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
        InTerm    = VarList{ivar}{3};
        InOrder   = VarList{ivar}{4};
        OutVname  = VarList{ivar}{5};

        InFname = sprintf('%s/%s_%s.h5', Ddir, InFprefix, Case);

        fprintf('  Input file: %s\n', InFname);
        fprintf('    Var name: %s\n', InVname);
        fprintf('    Var term: %d\n', InTerm);
        fprintf('    Var order: %d\n', InOrder);
        fprintf('\n');

        % read in and 
        TERMS = squeeze(hdf5read(InFname, InVname));
        NPTS = squeeze(hdf5read(InFname, 'num_points'));
        Z = squeeze(hdf5read(InFname, 'z_coords'));
        T = squeeze(hdf5read(InFname, 't_coords')) / 3600;   % hr

        Nz = length(Z);
        Nt = length(T);

        % *_M1 files will lose their first two dimensions, so put them back
        if (strcmp(InVname, 'cloud') || strcmp(InVname, 'theta_e'))
          TERMS = reshape(TERMS, [ 1 1 Nz Nt ]);
        end

        % find the indices of the start, mid and end times
        % use these to generate 4 profiles: one each at Tstart, Tmid, Tend
        % and one across range from Tstart to Tend
        T1 = find(T >= Tstart, 1, 'first');
        T2 = find(T >= Tmid, 1, 'first');
        T3 = find(T >= Tend, 1, 'first');

        % MOMENTS is organized: (Nterms, Norder, Nz)
        [ MOMENTS_T1    OUT_NPTS ] = GenMoments(TERMS, NPTS, T1, T1);
        [ MOMENTS_T2    OUT_NPTS ] = GenMoments(TERMS, NPTS, T2, T2);
        [ MOMENTS_T3    OUT_NPTS ] = GenMoments(TERMS, NPTS, T3, T3);
        [ MOMENTS_T1_T3 OUT_NPTS ] = GenMoments(TERMS, NPTS, T1, T3);

        OUT_MOMENTS_T1    = squeeze(MOMENTS_T1(InTerm, InOrder, :));
        OUT_MOMENTS_T2    = squeeze(MOMENTS_T2(InTerm, InOrder, :));
        OUT_MOMENTS_T3    = squeeze(MOMENTS_T3(InTerm, InOrder, :));
        OUT_MOMENTS_T1_T3 = squeeze(MOMENTS_T1_T3(InTerm, InOrder, :));

        % GenMoments() will fill levels with zero count (NPTS == 0) with nans. For
        % cloud mass we want these moments to be zero.
        if (strcmp(InVname, 'cloud'))
          OUT_MOMENTS_T1(isnan(OUT_MOMENTS_T1)) = 0;
          OUT_MOMENTS_T2(isnan(OUT_MOMENTS_T2)) = 0;
          OUT_MOMENTS_T3(isnan(OUT_MOMENTS_T3)) = 0;
          OUT_MOMENTS_T1_T3(isnan(OUT_MOMENTS_T1_T3)) = 0;
        end

        % Write out data - put in dummy x, y and t coordinates
        Xdummy = 1;
        Ydummy = 1;
        Tdummy = 1;

        OutVar = reshape(OUT_MOMENTS_T1, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_T%d', OutVname, Tstart);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(OUT_MOMENTS_T2, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_T%d', OutVname, Tmid);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(OUT_MOMENTS_T3, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_T%d', OutVname, Tend);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(OUT_MOMENTS_T1_T3, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_T%d_T%d', OutVname, Tstart, Tend);
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
