function [ ] = GenMomentFiles()
% GenMomentFiles generate w moment and w flux data

    Mdir = 'MomentData';
    Ddir = 'DIAGS';

    CaseList = {
%      { 'z.atex.ccn0050.sst293' 12 36 }
%      { 'z.atex.ccn0100.sst293' 12 36 }
%      { 'z.atex.ccn0200.sst293' 12 36 }
%      { 'z.atex.ccn0400.sst293' 12 36 }
%      { 'z.atex.ccn0800.sst293' 12 36 }
%      { 'z.atex.ccn1600.sst293' 12 36 }

      { 'z.atex.ccn0050.sst298' 12 36 }
%      { 'z.atex.ccn0100.sst298' 12 36 }
%      { 'z.atex.ccn0200.sst298' 12 36 }
%      { 'z.atex.ccn0400.sst298' 12 36 }
%      { 'z.atex.ccn0800.sst298' 12 36 }
      { 'z.atex.ccn1600.sst298' 12 36 }
      };
    Ncases = length(CaseList);

    %Tfirst = 0; % calc moment for each time step, the take average over time steps.
    Tfirst = 1; % take average over time steps first, then calculate moments.

    % These lists are organized as 2D cell arrays. Each row is one complete spec for one variable.
    % Syntax for rows:
    %    { 'input file prefix' 'input file var name' 'input var term number' 'input var order number' 'output var name' }
    VarList = {
      % means
      { 'cloud_M1_c0p01'        'cloud'     1 1 'cloud_c0p01'    }

%      { 'w_M3_up0p01_all_cld'   'w-w-w'     1 1 'w_up0p01_all_cld' }
%      { 'w_M3_dn0p01_all_cld'   'w-w-w'     1 1 'w_dn0p01_all_cld' }

%      { 'w_theta_v_flux_up0p01_all_cld' 'w-theta_v' 2 1 'theta_v_up0p01_all_cld'   }
%      { 'w_theta_v_flux_dn0p01_all_cld' 'w-theta_v' 2 1 'theta_v_dn0p01_all_cld'   }

%      { 'w_vapor_flux_up0p01_all_cld' 'w-vapor'   2 1 'vapor_up0p01_all_cld'   }
%      { 'w_vapor_flux_dn0p01_all_cld' 'w-vapor'   2 1 'vapor_dn0p01_all_cld'   }

      % fluxes (covariances)
%      { 'w_theta_v_flux_up0p01_all_cld' 'w-theta_v' 1 2 'w-theta_v_up0p01_all_cld' }
%      { 'w_theta_v_flux_dn0p01_all_cld' 'w-theta_v' 1 2 'w-theta_v_dn0p01_all_cld' }

%      { 'w_vapor_flux_up0p01_all_cld' 'w-vapor'   1 2 'w-vapor_up0p01_all_cld'   }
%      { 'w_vapor_flux_dn0p01_all_cld' 'w-vapor'   1 2 'w-vapor_dn0p01_all_cld'   }

      % variances
%      { 'w_M3_up0p01_all_cld'      'w-w-w'     1 2 'w-w_up0p01_all_cld'   }
%      { 'w_M3_dn0p01_all_cld'      'w-w-w'     1 2 'w-w_dn0p01_all_cld'   }

      % skews
%      { 'w_M3_up0p01_all_cld'      'w-w-w'     1 3 'w-w-w_up0p01_all_cld'   }
%      { 'w_M3_dn0p01_all_cld'      'w-w-w'     1 3 'w-w-w_dn0p01_all_cld'   }

      };

    TotalN = 158404; % excluding borders --> 398 * 398

    for icase = 1:Ncases;
      Case   = CaseList{icase}{1};
      Tstart = CaseList{icase}{2};
      Tend   = CaseList{icase}{3};

      OutFname = sprintf('%s/moments_%s.h5', Ddir, Case);

      fprintf('***************************************************************\n');
      fprintf('Generating moment/flux profiles:\n');
      fprintf('  Case: %s\n', Case);
      fprintf('  Temporal Averaging:\n');
      fprintf('     Tstart = %f\n', Tstart);
      fprintf('     Tend   = %f\n', Tend);
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

        InFname = sprintf('%s/%s_%s.h5', Mdir, InFprefix, Case);

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

        % find the indices of the start, mid, end and "all" time intervals
        T1 = find(T >= Tstart, 1, 'first');
        T2 = find(T <= Tend, 1, 'last');

        % MOMENTS is organized: (Nterms, Norder, Nz)
        [ MOMENTS OUT_NPTS ] = GenMoments(TERMS, NPTS, T1, T2, Tfirst);
        OUT_MOMENTS = squeeze(MOMENTS(InTerm, InOrder, :));

        % GenMoments() will fill levels with zero count (NPTS == 0) with nans. For
        % cloud mass we want these moments to be zero.
        if (strcmp(InVname, 'cloud'))
          OUT_MOMENTS(isnan(OUT_MOMENTS)) = 0;
        end

        % Generate a fraction statistic - the ratio of number of points selected 
        % to total number of points in domain
        PROF_FRAC = NPTS ./ TotalN;
        PROF_FRAC = squeeze(mean(PROF_FRAC(:,T1:T2), 2));


        % Write out data - put in dummy x, y and t coordinates
        Xdummy = 1;
        Ydummy = 1;
        Tdummy = 1;

        OutVar = reshape(OUT_MOMENTS, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s', OutVname);
        fprintf('  Writing var: %s\n', OutVarName);
        hdf5write(OutFname, OutVarName, OutVar, 'WriteMode', 'append');

        OutVar = reshape(PROF_FRAC, [ 1 1 Nz 1 ]);
        OutVarName = sprintf('%s_frac', OutVname);
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
