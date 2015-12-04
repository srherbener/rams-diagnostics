function [ ] = GenStormXsections()

  Ddir = 'DIAGS';
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % From TSD_SAL_DUST RAMS output (Reference density)
  RhoAir = [
    1.196
    1.191
    1.185
    1.179
    1.172
    1.165
    1.157
    1.147
    1.136
    1.124
    1.111
    1.097
    1.082
    1.066
    1.051
    1.034
    1.017
    1.000
    0.982
    0.963
    0.945
    0.925
    0.905
    0.883
    0.860
    0.833
    0.804
    0.773
    0.741
    0.711
    0.679
    0.647
    0.612
    0.577
    0.541
    0.505
    0.469
    0.432
    0.394
    0.355
    0.316
    0.279
    0.244
    0.210
    0.179
    0.150
    0.126
    0.105
    0.087
    0.073
    0.062
    0.052
    0.044
    0.038
    0.032
    0.027
    ];

  Cases = {
   'TSD_SAL_DUST'
   'TSD_SAL_NODUST'
   'TSD_NONSAL_DUST'
   'TSD_NONSAL_NODUST'
   };
  Ncases = length(Cases);

  % Description of cross sections
  XsectionList = {
    % in_file in_var pre_sal_out_var sal_out_var
    { 'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_ps_speed_t' '/all_ps_speed_t' }
    { 'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_ps_speed_r' '/all_ps_speed_r' }
    { 'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_s_speed_t' '/all_s_speed_t' }
    { 'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_s_speed_r' '/all_s_speed_r' }

    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_ps_updraft' '/all_ps_updraft' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_ps_dndraft' '/all_ps_dndraft' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_s_updraft' '/all_s_updraft' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_s_dndraft' '/all_s_dndraft' }

    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/all_ps_theta_e' '/all_ps_theta_e' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/all_s_theta_e' '/all_s_theta_e' }

    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/all_ps_theta'  '/all_ps_theta' }
    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/all_s_theta'  '/all_s_theta' }

    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/all_ps_tempc'  '/all_ps_tempc' }
    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/all_s_tempc'  '/all_s_tempc' }

    { 'DIAGS/hist_meas_az_relhum_<CASE>.h5' '/all_ps_relhum'  '/all_ps_relhum' }
    { 'DIAGS/hist_meas_az_relhum_<CASE>.h5' '/all_s_relhum'  '/all_s_relhum' }

    { 'DIAGS/hist_meas_az_pcprate_<CASE>.h5' '/all_ps_pcprate' '/all_ps_pcprate' }
    { 'DIAGS/hist_meas_az_pcprate_<CASE>.h5' '/all_s_pcprate' '/all_s_pcprate' }

    { 'DIAGS/hist_meas_az_vint_cond_<CASE>.h5' '/all_ps_vint_cond' '/all_ps_vint_cond' }
    { 'DIAGS/hist_meas_az_vint_cond_<CASE>.h5' '/all_s_vint_cond' '/all_s_vint_cond' }

    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/all_ps_vapor' '/all_ps_vapor' }
    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/all_s_vapor' '/all_s_vapor' }

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_ps_d1_num'  '/all_ps_d1_num'  }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_ps_d2_num'  '/all_ps_d2_num'  }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_s_d1_num'  '/all_s_d1_num'  }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_s_d2_num'  '/all_s_d2_num'  }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_cloud' '/all_ps_dust_cloud' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_rain'  '/all_ps_dust_rain'  }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_pris'  '/all_ps_dust_pris'  }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_snow'  '/all_ps_dust_snow'  }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_aggr'  '/all_ps_dust_aggr'  }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_graup' '/all_ps_dust_graup' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_hail'  '/all_ps_dust_hail'  }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_hydro' '/all_ps_dust_hydro' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_cloud'  '/all_s_dust_cloud'  }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_rain'   '/all_s_dust_rain'   }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_pris'   '/all_s_dust_pris'   }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_snow'   '/all_s_dust_snow'   }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_aggr'   '/all_s_dust_aggr'   }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_graup'  '/all_s_dust_graup'  }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_hail'   '/all_s_dust_hail'   }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_hydro'  '/all_s_dust_hydro'  }

    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_ps_cloud_mass' '/all_ps_cloud_mass' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_s_cloud_mass'  '/all_s_cloud_mass' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_ps_cloud_num'  '/all_ps_cloud_num' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_s_cloud_num'   '/all_s_cloud_num' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_ps_cloud_diam' '/all_ps_cloud_diam' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_s_cloud_diam'  '/all_s_cloud_diam' }

    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_ps_rain_mass' '/all_ps_rain_mass' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_s_rain_mass'  '/all_s_rain_mass' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_ps_rain_num'  '/all_ps_rain_num' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_s_rain_num'   '/all_s_rain_num' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_ps_rain_diam' '/all_ps_rain_diam' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_s_rain_diam'  '/all_s_rain_diam' }

    { 'DIAGS/hist_meas_az_pris_<CASE>.h5' '/all_ps_pris_mass'      '/all_ps_pris_mass' }
    { 'DIAGS/hist_meas_az_pris_<CASE>.h5' '/all_s_pris_mass'       '/all_s_pris_mass' }

    { 'DIAGS/hist_meas_az_snow_<CASE>.h5' '/all_ps_snow_mass'      '/all_ps_snow_mass' }
    { 'DIAGS/hist_meas_az_snow_<CASE>.h5' '/all_s_snow_mass'       '/all_s_snow_mass' }

    { 'DIAGS/hist_meas_az_aggr_<CASE>.h5' '/all_ps_aggr_mass'      '/all_ps_aggr_mass' }
    { 'DIAGS/hist_meas_az_aggr_<CASE>.h5' '/all_s_aggr_mass'       '/all_s_aggr_mass' }

    { 'DIAGS/hist_meas_az_graup_<CASE>.h5' '/all_ps_graup_mass'      '/all_ps_graup_mass' }
    { 'DIAGS/hist_meas_az_graup_<CASE>.h5' '/all_s_graup_mass'       '/all_s_graup_mass' }

    { 'DIAGS/hist_meas_az_hail_<CASE>.h5' '/all_ps_hail_mass'      '/all_ps_hail_mass' }
    { 'DIAGS/hist_meas_az_hail_<CASE>.h5' '/all_s_hail_mass'       '/all_s_hail_mass' }

    { 'DIAGS/hist_meas_az_tcond_<CASE>.h5' '/all_ps_tcond_mass' '/all_ps_tcond_mass' }
    { 'DIAGS/hist_meas_az_tcond_<CASE>.h5' '/all_s_tcond_mass'  '/all_s_tcond_mass' }

    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/all_ps_lhf_cool' '/all_ps_lhf_cool' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/all_s_lhf_cool'  '/all_s_lhf_cool' }

    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/all_ps_lhf_heat' '/all_ps_lhf_heat' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/all_s_lhf_heat'  '/all_s_lhf_heat' }

    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/all_ps_lhv_cool' '/all_ps_lhv_cool' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/all_s_lhv_cool'  '/all_s_lhv_cool' }

    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/all_ps_lhv_heat' '/all_ps_lhv_heat' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/all_s_lhv_heat'  '/all_s_lhv_heat' }

    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/all_ps_liq_evap' '/all_ps_liq_evap' }
    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/all_s_liq_evap'  '/all_s_liq_evap' }

    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/all_ps_liq_cond' '/all_ps_liq_cond' }
    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/all_s_liq_cond'  '/all_s_liq_cond' }

    { 'DIAGS/hist_meas_az_cloud_evap_<CASE>.h5' '/all_ps_cloud_evap' '/all_ps_cloud_evap' }
    { 'DIAGS/hist_meas_az_cloud_evap_<CASE>.h5' '/all_s_cloud_evap'  '/all_s_cloud_evap' }

    { 'DIAGS/hist_meas_az_cloud_cond_<CASE>.h5' '/all_ps_cloud_cond' '/all_ps_cloud_cond' }
    { 'DIAGS/hist_meas_az_cloud_cond_<CASE>.h5' '/all_s_cloud_cond'  '/all_s_cloud_cond' }

    { 'DIAGS/hist_meas_az_rain_evap_<CASE>.h5' '/all_ps_rain_evap' '/all_ps_rain_evap' }
    { 'DIAGS/hist_meas_az_rain_evap_<CASE>.h5' '/all_s_rain_evap'  '/all_s_rain_evap' }

    { 'DIAGS/hist_meas_az_rain_cond_<CASE>.h5' '/all_ps_rain_cond' '/all_ps_rain_cond' }
    { 'DIAGS/hist_meas_az_rain_cond_<CASE>.h5' '/all_s_rain_cond'  '/all_s_rain_cond' }

    };
  Nsets = length(XsectionList);


  for icase = 1:Ncases
    Case = Cases{icase};

    fprintf('**********************************************************\n');
    fprintf('Generating storm cross-sections for case: %s\n', Case);
    fprintf('\n');

    % Place all vars into one case specific output file.
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('DIAGS/storm_xsections_%s.h5', Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    icount = 0;
    for iset = 1:Nsets
      InFile   = regexprep(XsectionList{iset}{1}, '<CASE>', Case);
      InVname  = XsectionList{iset}{2};
      OutVname = XsectionList{iset}{3};

      % skip this profile set if doing dust and on a NODUST case
      if ((~isempty(regexp(Case, 'NODUST'))) && ...
          ((~isempty(regexp(InVname, 'd[12]_num'))) || (~isempty(regexp(InVname, 'dust_'))) || (~isempty(regexp(InVname, 'dustifn_')))))
        continue
      else
        icount = icount + 1;
      end

      % Read in the variables
      fprintf('  Reading: %s (%s)\n', InFile, InVname);

      VAR = squeeze(h5read(InFile, InVname));
      X   = squeeze(h5read(InFile, '/x_coords'));
      Y   = squeeze(h5read(InFile, '/y_coords'));
      Z   = squeeze(h5read(InFile, '/z_coords'));
      T   = squeeze(h5read(InFile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % If a hydrometeor number concentration, then convert #/kg to #/cm^2.
      % VAR is (x,z) so repeat RhoAir along the x dimension.
      if (~isempty(regexp(InVname, 'cloud_num')) || ...
          ~isempty(regexp(InVname, 'rain_num'))  || ...
          ~isempty(regexp(InVname, 'pris_num'))  || ...
          ~isempty(regexp(InVname, 'snow_num'))  || ...
          ~isempty(regexp(InVname, 'aggr_num'))  || ...
          ~isempty(regexp(InVname, 'graup_num')) || ...
          ~isempty(regexp(InVname, 'hail_num')))
        VAR      = VAR .* repmat(RhoAir', [ Nx 1 ]) .* 1e-6;
      end

      % If this is the first set, write out the coordinates into the output file
      % so that subsequent variables can have these coordinates attached.
      if (icount == 1)
        Xname = '/x_coords';
        Yname = '/y_coords';
        Zname = '/z_coords';
        Tname = '/t_coords';
  
        CreateDimensionsXyzt(OutFile, X, Y, Z, [ 1 ], Xname, Yname, Zname, Tname);
        % Add COARDS annotations
        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
      end

      % Determine the size and outdims spec for the output using size():
      %     [ 1 1 ]             --> Vsize = 1,         DimOrder = { }
      %     [ n 1 ] or [ 1 n ]  --> Vsize = n,         DimOrder = { 'x' }
      %     [ n m ]             --> Vsize = size(VAR), DimOrder = { 'x' 'z' };
      %
      % size() always returns at least two values. At this point, only
      % accommodating up to scalar, vector or 2D.
      Vsize = size(VAR);
      if ((Vsize(1) == 1) && (Vsize(2) == 1))
        % VAR is [ 1 1 ], ie. a scalar value
        Vsize = 1;
        DimOrder = { };
      elseif ((Vsize(1) == 1) || (Vsize(2) == 1))
        % VAR is [ 1 n ] or [ n 1 ], ie. a vector value
        % This implies 1D in radius which is the x dimension
        Vsize = Nx;
        DimOrder = { 'x' };
      else
        % Var is [ n m ]
        % Vsize does not need to be modified
        % This implies 2D in radius and height which are the x and z dimensions
        DimOrder = { 'x' 'z' };
      end

      % Collect output vars into one case specific file.
      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
      h5create(OutFile, OutVname, Vsize);
      h5write(OutFile, OutVname, VAR);

      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);

      fprintf('\n');
    end
  end
end
