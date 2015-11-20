function [ ] = GenStormXsections()

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

%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_ps_d1_mass' '/all_ps_d1_mass' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_ps_d2_mass' '/all_ps_d2_mass' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_s_d1_mass' '/all_s_d1_mass' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_s_d2_mass' '/all_s_d2_mass' }
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
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_cloud'  '/all_s_dust_cloud'  }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_rain'   '/all_s_dust_rain'   }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_pris'   '/all_s_dust_pris'   }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_snow'   '/all_s_dust_snow'   }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_aggr'   '/all_s_dust_aggr'   }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_graup'  '/all_s_dust_graup'  }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_hail'   '/all_s_dust_hail'   }

    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_ps_cloud_mass' '/all_ps_cloud_mass' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_s_cloud_mass'  '/all_s_cloud_mass' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_ps_cloud_num'  '/all_ps_cloud_num' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_s_cloud_num'   '/all_s_cloud_num' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_ps_cloud_diam' '/all_ps_cloud_diam' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_s_cloud_diam'  '/all_s_cloud_diam' }

    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_ps_rain_mass'      '/all_ps_rain_mass' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_s_rain_mass'       '/all_s_rain_mass' }

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
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/all_ps_lhf_heat' '/all_ps_lhf_heat' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/all_s_lhf_cool'  '/all_s_lhf_cool' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/all_s_lhf_heat'  '/all_s_lhf_heat' }

    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/all_ps_lhv_cool' '/all_ps_lhv_cool' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/all_ps_lhv_heat' '/all_ps_lhv_heat' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/all_s_lhv_cool'  '/all_s_lhv_cool' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/all_s_lhv_heat'  '/all_s_lhv_heat' }
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

      ControlInFile = regexprep(XsectionList{iset}{1}, '<CASE>', ControlCase);

      OutDiffVname = sprintf('%s_diff', OutVname);

      % Read in the variables
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      fprintf('  Reading: %s (%s)\n', ControlInFile, InVname);

      VAR = squeeze(h5read(InFile, InVname));
      X   = squeeze(h5read(InFile, '/x_coords'));
      Y   = squeeze(h5read(InFile, '/y_coords'));
      Z   = squeeze(h5read(InFile, '/z_coords'));
      T   = squeeze(h5read(InFile, '/t_coords'));

      CNTL_VAR = squeeze(h5read(ControlInFile, InVname));

      DIFF_VAR = VAR - CNTL_VAR;

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

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

      fprintf('  Writing: %s (%s)\n', OutFile, OutDiffVname);
      h5create(OutFile, OutDiffVname, Vsize);
      h5write(OutFile, OutDiffVname, DIFF_VAR);

      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, OutDiffVname, DimOrder, Xname, Yname, Zname, Tname);

      fprintf('\n');
    end
  end
end
