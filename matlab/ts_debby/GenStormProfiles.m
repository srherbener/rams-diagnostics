function [ ] = GenStormProfiles()

  Ddir = 'DIAGS';
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  Cases = {
   'TSD_SAL_DUST'
   'TSD_SAL_NODUST'
   'TSD_NONSAL_DUST'
   'TSD_NONSAL_NODUST'
   };
  Ncases = length(Cases);

  ControlCase = 'TSD_SAL_DUST';

  % Description of profiles
  ProfileList = {
    % in_file in_var out_var pre_sal_rb_var sal_core_var sal_rb_var
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_core_avg_speed_t' '/ps_core_speed_t' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_rb_avg_speed_t' '/ps_rb_speed_t' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_core_avg_speed_t' '/s_core_speed_t' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_rb_avg_speed_t' '/s_rb_speed_t' }

    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_core_avg_speed_r' '/ps_core_speed_r' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_rb_avg_speed_r' '/ps_rb_speed_r' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_core_avg_speed_r' '/s_core_speed_r' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_rb_avg_speed_r' '/s_rb_speed_r' }

    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_core_avg_updraft' '/ps_core_updraft' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_rb_avg_updraft' '/ps_rb_updraft' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_core_avg_updraft' '/s_core_updraft' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_rb_avg_updraft' '/s_rb_updraft' }

    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_core_avg_dndraft' '/ps_core_dndraft' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_rb_avg_dndraft' '/ps_rb_dndraft' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_core_avg_dndraft' '/s_core_dndraft' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_rb_avg_dndraft' '/s_rb_dndraft' }

    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/ps_core_avg_theta_e' '/ps_core_theta_e' }
    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/ps_rb_avg_theta_e' '/ps_rb_theta_e' }
    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/s_core_avg_theta_e' '/s_core_theta_e' }
    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/s_rb_avg_theta_e' '/s_rb_theta_e' }

    { 'DIAGS/hist_meas_pcprate_<CASE>.h5' '/ps_core_avg_pcprate' '/ps_core_pcprate' }
    { 'DIAGS/hist_meas_pcprate_<CASE>.h5' '/ps_rb_avg_pcprate' '/ps_rb_pcprate' }
    { 'DIAGS/hist_meas_pcprate_<CASE>.h5' '/s_core_avg_pcprate' '/s_core_pcprate' }
    { 'DIAGS/hist_meas_pcprate_<CASE>.h5' '/s_rb_avg_pcprate' '/s_rb_pcprate' }

    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/ps_core_avg_vint_cond' '/ps_core_vint_cond' }
    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/ps_rb_avg_vint_cond' '/ps_rb_vint_cond' }
    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/s_core_avg_vint_cond' '/s_core_vint_cond' }
    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/s_rb_avg_vint_cond' '/s_rb_vint_cond' }

    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_core_avg_d1_mass' '/ps_core_d1_mass' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_rb_avg_d1_mass' '/ps_rb_d1_mass' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_core_avg_d1_mass' '/s_core_d1_mass' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_rb_avg_d1_mass' '/s_rb_d1_mass' }

    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_core_avg_d1_num'  '/ps_core_d1_num'  }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_rb_avg_d1_num'  '/ps_rb_d1_num'  }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_core_avg_d1_num'  '/s_core_d1_num'  }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_rb_avg_d1_num'  '/s_rb_d1_num'  }

    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_core_avg_d2_mass' '/ps_core_d2_mass' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_rb_avg_d2_mass' '/ps_rb_d2_mass' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_core_avg_d2_mass' '/s_core_d2_mass' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_rb_avg_d2_mass' '/s_rb_d2_mass' }

    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_core_avg_d2_num'  '/ps_core_d2_num'  }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_rb_avg_d2_num'  '/ps_rb_d2_num'  }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_core_avg_d2_num'  '/s_core_d2_num'  }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_rb_avg_d2_num'  '/s_rb_d2_num'  }

    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/ps_core_avg_cloud'  '/ps_core_cloud' }
    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/ps_rb_avg_cloud'  '/ps_rb_cloud' }
    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/s_core_avg_cloud'  '/s_core_cloud' }
    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/s_rb_avg_cloud'  '/s_rb_cloud' }

    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_core_avg_rain'  '/ps_core_rain'  }
    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_rb_avg_rain'  '/ps_rb_rain'  }
    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_core_avg_rain'  '/s_core_rain'  }
    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_rb_avg_rain'  '/s_rb_rain'  }

    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_core_avg_pris'  '/ps_core_pris'  }
    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_rb_avg_pris'  '/ps_rb_pris'  }
    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_core_avg_pris'  '/s_core_pris'  }
    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_rb_avg_pris'  '/s_rb_pris'  }

    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_core_avg_snow'  '/ps_core_snow'  }
    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_rb_avg_snow'  '/ps_rb_snow'  }
    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_core_avg_snow'  '/s_core_snow'  }
    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_rb_avg_snow'  '/s_rb_snow'  }

    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_core_avg_aggr'  '/ps_core_aggr'  }
    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_rb_avg_aggr'  '/ps_rb_aggr'  }
    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_core_avg_aggr'  '/s_core_aggr'  }
    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_rb_avg_aggr'  '/s_rb_aggr'  }

    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_core_avg_graup'  '/ps_core_graup'  }
    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_rb_avg_graup'  '/ps_rb_graup'  }
    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_core_avg_graup'  '/s_core_graup'  }
    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_rb_avg_graup'  '/s_rb_graup'  }

    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_core_avg_hail'  '/ps_core_hail'  }
    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_rb_avg_hail'  '/ps_rb_hail'  }
    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_core_avg_hail'  '/s_core_hail'  }
    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_rb_avg_hail'  '/s_rb_hail'  }

    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/ps_core_avg_lhf_cool'  '/ps_core_lhf_cool'  }
    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/ps_rb_avg_lhf_cool'  '/ps_rb_lhf_cool'  }
    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/s_core_avg_lhf_cool'  '/s_core_lhf_cool'  }
    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/s_rb_avg_lhf_cool'  '/s_rb_lhf_cool'  }

    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/ps_core_avg_lhf_heat'  '/ps_core_lhf_heat'  }
    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/ps_rb_avg_lhf_heat'  '/ps_rb_lhf_heat'  }
    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/s_core_avg_lhf_heat'  '/s_core_lhf_heat'  }
    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/s_rb_avg_lhf_heat'  '/s_rb_lhf_heat'  }

    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/ps_core_avg_lhv_cool'  '/ps_core_lhv_cool'  }
    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/ps_rb_avg_lhv_cool'  '/ps_rb_lhv_cool'  }
    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/s_core_avg_lhv_cool'  '/s_core_lhv_cool'  }
    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/s_rb_avg_lhv_cool'  '/s_rb_lhv_cool'  }

    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/ps_core_avg_lhv_heat'  '/ps_core_lhv_heat'  }
    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/ps_rb_avg_lhv_heat'  '/ps_rb_lhv_heat'  }
    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/s_core_avg_lhv_heat'  '/s_core_lhv_heat'  }
    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/s_rb_avg_lhv_heat'  '/s_rb_lhv_heat'  }

    };
  Nsets = length(ProfileList);


  for icase = 1:Ncases
    Case = Cases{icase};

    fprintf('**********************************************************\n');
    fprintf('Generating storm profiles for case: %s\n', Case);
    fprintf('\n');

    % Place all vars into one case specific output file.
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('DIAGS/storm_profiles_%s.h5', Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for iset = 1:Nsets
      InFile   = regexprep(ProfileList{iset}{1}, '<CASE>', Case);
      InVname  = ProfileList{iset}{2};
      OutVname = ProfileList{iset}{3};

      ControlInFile = regexprep(ProfileList{iset}{1}, '<CASE>', ControlCase);

      OutDiffVname = sprintf('%s_diff', OutVname);

      % Read in the variables, split into pre-SAL and SAL sections
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      fprintf('  Reading: %s (%s)\n', ControlInFile, InVname);

      VAR = squeeze(h5read(InFile, InVname));
      X   = squeeze(h5read(InFile, '/x_coords'));
      Y   = squeeze(h5read(InFile, '/y_coords'));
      Z   = squeeze(h5read(InFile, '/z_coords'));
      T   = squeeze(h5read(InFile, '/t_coords'));

      CNTL_VAR = squeeze(h5read(ControlInFile, InVname));

      % Change nans to zeros to make the profiles look better on plots.
      VAR(isnan(VAR)) = 0;
      CNTL_VAR(isnan(CNTL_VAR)) = 0;

      DIFF_VAR = CNTL_VAR - VAR;

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % If this is the first set, write out the coordinates into the output file
      % so that subsequent variables can have these coordinates attached.
      if (iset == 1)
        Xname = '/x_coords';
        Yname = '/y_coords';
        Zname = '/z_coords';
        Tname = '/t_coords';
  
        CreateDimensionsXyzt(OutFile, [ 1 ], Y, Z, [ 1 ], Xname, Yname, Zname, Tname);
        % Add COARDS annotations
        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
      end

      % Determine the size and outdims spec for the output using size():
      %     [ 1 1 ]             --> Vsize = 1,         DimOrder = { }
      %     [ n 1 ] or [ 1 n ]  --> Vsize = n,         DimOrder = { 'z' }
      %
      % size() always returns at least two values. At this point, only
      % accommodating up to scalar or vector.
      Vsize = size(VAR);
      if ((Vsize(1) == 1) && (Vsize(2) == 1))
        % VAR is [ 1 1 ], ie. a scalar value
        Vsize = 1;
        DimOrder = { };
      elseif ((Vsize(1) == 1) || (Vsize(2) == 1))
        % VAR is [ 1 n ] or [ n 1 ], ie. a vector value
        % This implies 1D in height which is the z dimension
        Vsize = Nz;
        DimOrder = { 'z' };
      end

      % Collect output vars into one case specific file.
      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
      h5create(OutFile, OutVname, Vsize);
      h5write(OutFile, OutVname, VAR);

      fprintf('  Writing: %s (%s)\n', OutFile, OutDiffVname);
      h5create(OutFile, OutDiffVname, Vsize);
      h5write(OutFile, OutDiffVname, DIFF_VAR);

      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, OutDiffVname, DimOrder, Xname, Yname, Zname, Tname);

      fprintf('\n');
    end
  end
end
