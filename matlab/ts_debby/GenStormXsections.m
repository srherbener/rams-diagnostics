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
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_avg_speed_t' '/ps_speed_t' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_avg_speed_r' '/ps_speed_r' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_avg_speed_t' '/s_speed_t' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_avg_speed_r' '/s_speed_r' }

    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_max_speed10m' '/ps_max_speed10m' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_max_speed10m' '/s_max_speed10m' }

    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_max_speed25m' '/ps_max_speed25m' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_max_speed25m' '/s_max_speed25m' }

    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_avg_updraft' '/ps_updraft' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_avg_dndraft' '/ps_dndraft' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_avg_updraft' '/s_updraft' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_avg_dndraft' '/s_dndraft' }

    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/ps_avg_theta_e' '/ps_theta_e' }
    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/s_avg_theta_e' '/s_theta_e' }

    { 'DIAGS/hist_meas_pcprate_<CASE>.h5' '/ps_avg_pcprate' '/ps_pcprate' }
    { 'DIAGS/hist_meas_pcprate_<CASE>.h5' '/s_avg_pcprate' '/s_pcprate' }

    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/ps_avg_vint_cond' '/ps_vint_cond' }
    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/s_avg_vint_cond' '/s_vint_cond' }

    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_avg_d1_mass' '/ps_d1_mass' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_avg_d1_num'  '/ps_d1_num'  }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_avg_d2_mass' '/ps_d2_mass' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/ps_avg_d2_num'  '/ps_d2_num'  }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_avg_d1_mass' '/s_d1_mass' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_avg_d1_num'  '/s_d1_num'  }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_avg_d2_mass' '/s_d2_mass' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/s_avg_d2_num'  '/s_d2_num'  }

    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/ps_avg_cloud'  '/ps_cloud' }
    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/s_avg_cloud'  '/s_cloud' }

    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_avg_rain'  '/ps_rain' }
    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_avg_rain'  '/s_rain' }

    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_avg_pris'  '/ps_pris' }
    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_avg_pris'  '/s_pris' }

    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_avg_snow'  '/ps_snow' }
    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_avg_snow'  '/s_snow' }

    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_avg_aggr'  '/ps_aggr' }
    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_avg_aggr'  '/s_aggr' }

    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_avg_graup'  '/ps_graup' }
    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_avg_graup'  '/s_graup' }

    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_avg_hail'  '/ps_hail' }
    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_avg_hail'  '/s_hail' }

    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/ps_avg_tcond'  '/ps_tcond' }
    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/s_avg_tcond'  '/s_tcond' }

    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/ps_avg_lhf_cool'  '/ps_lhf_cool' }
    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/ps_avg_lhf_heat'  '/ps_lhf_heat' }
    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/s_avg_lhf_cool'  '/s_lhf_cool' }
    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/s_avg_lhf_heat'  '/s_lhf_heat' }

    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/ps_avg_lhv_cool'  '/ps_lhv_cool' }
    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/ps_avg_lhv_heat'  '/ps_lhv_heat' }
    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/s_avg_lhv_cool'  '/s_lhv_cool' }
    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/s_avg_lhv_heat'  '/s_lhv_heat' }
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

    for iset = 1:Nsets
      InFile   = regexprep(XsectionList{iset}{1}, '<CASE>', Case);
      InVname  = XsectionList{iset}{2};
      OutVname = XsectionList{iset}{3};

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
