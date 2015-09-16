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

  % Two time periods:
  %    Pre SAL: T = 10 to 30 h (after RI, before encounter SAL)
  %        SAL: T = 40 to 60 h (during SAL)

  PreSalTstart = 10;
  PreSalTend   = 30;
  SalTstart    = 40;
  SalTend      = 60;

  % Description of cross sections
  XsectionList = {
    % in_file in_var pre_sal_out_var sal_out_var
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_t_fa' '/ps_speed_t_fa' '/s_speed_t_fa' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_t_wm' '/ps_speed_t_wm' '/s_speed_t_wm' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_r_fa' '/ps_speed_r_fa' '/s_speed_r_fa' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_r_wm' '/ps_speed_r_wm' '/s_speed_r_wm' }

    { 'DIAGS/hist_meas_w_<CASE>.h5' '/avg_updraft_fa' '/ps_updraft_fa' '/s_updraft_fa' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/avg_updraft_wm' '/ps_updraft_wm' '/s_updraft_wm' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/avg_dndraft_fa' '/ps_dndraft_fa' '/s_dndraft_fa' }
    { 'DIAGS/hist_meas_w_<CASE>.h5' '/avg_dndraft_wm' '/ps_dndraft_wm' '/s_dndraft_wm' }

    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/avg_theta_e_fa' '/ps_theta_e_fa' '/s_theta_e_fa' }
    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/avg_theta_e_wm' '/ps_theta_e_wm' '/s_theta_e_wm' }

    { 'DIAGS/hist_meas_pcprate_<CASE>.h5' '/avg_pcprate_fa' '/ps_pcprate_fa' '/s_pcprate_fa' }
    { 'DIAGS/hist_meas_pcprate_<CASE>.h5' '/avg_pcprate_wm' '/ps_pcprate_wm' '/s_pcprate_wm' }

    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/avg_vint_cond_fa' '/ps_vint_cond_fa' '/s_vint_cond_fa' }
    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/avg_vint_cond_wm' '/ps_vint_cond_wm' '/s_vint_cond_wm' }

    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d1_mass_fa' '/ps_d1_mass_fa' '/s_d1_mass_fa' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d1_mass_wm' '/ps_d1_mass_wm' '/s_d1_mass_wm' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d1_num_fa'  '/ps_d1_num_fa'  '/s_d1_num_fa' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d1_num_wm'  '/ps_d1_num_wm'  '/s_d1_num_wm' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d2_mass_fa' '/ps_d2_mass_fa' '/s_d2_mass_fa' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d2_mass_wm' '/ps_d2_mass_wm' '/s_d2_mass_wm' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d2_num_fa'  '/ps_d2_num_fa'  '/s_d2_num_fa' }
    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d2_num_wm'  '/ps_d2_num_wm'  '/s_d2_num_wm' }

    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/avg_cloud_fa'  '/ps_cloud_fa'  '/s_cloud_fa' }
    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/avg_cloud_wm'  '/ps_cloud_wm'  '/s_cloud_wm' }

    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/avg_rain_fa'  '/ps_rain_fa'  '/s_rain_fa' }
    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/avg_rain_wm'  '/ps_rain_wm'  '/s_rain_wm' }

    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/avg_pris_fa'  '/ps_pris_fa'  '/s_pris_fa' }
    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/avg_pris_wm'  '/ps_pris_wm'  '/s_pris_wm' }

    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/avg_snow_fa'  '/ps_snow_fa'  '/s_snow_fa' }
    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/avg_snow_wm'  '/ps_snow_wm'  '/s_snow_wm' }

    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/avg_aggr_fa'  '/ps_aggr_fa'  '/s_aggr_fa' }
    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/avg_aggr_wm'  '/ps_aggr_wm'  '/s_aggr_wm' }

    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/avg_graup_fa'  '/ps_graup_fa'  '/s_graup_fa' }
    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/avg_graup_wm'  '/ps_graup_wm'  '/s_graup_wm' }

    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/avg_hail_fa'  '/ps_hail_fa'  '/s_hail_fa' }
    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/avg_hail_wm'  '/ps_hail_wm'  '/s_hail_wm' }

    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/avg_lhf_cool_fa'  '/ps_lhf_cool_fa'  '/s_lhf_cool_fa' }
    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/avg_lhf_cool_wm'  '/ps_lhf_cool_wm'  '/s_lhf_cool_wm' }
    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/avg_lhf_heat_fa'  '/ps_lhf_heat_fa'  '/s_lhf_heat_fa' }
    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/avg_lhf_heat_wm'  '/ps_lhf_heat_wm'  '/s_lhf_heat_wm' }

    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/avg_lhv_cool_fa'  '/ps_lhv_cool_fa'  '/s_lhv_cool_fa' }
    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/avg_lhv_cool_wm'  '/ps_lhv_cool_wm'  '/s_lhv_cool_wm' }
    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/avg_lhv_heat_fa'  '/ps_lhv_heat_fa'  '/s_lhv_heat_fa' }
    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/avg_lhv_heat_wm'  '/ps_lhv_heat_wm'  '/s_lhv_heat_wm' }
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
      InFile      = regexprep(XsectionList{iset}{1}, '<CASE>', Case);
      InVname     = XsectionList{iset}{2};
      PreSalVname = XsectionList{iset}{3};
      SalVname    = XsectionList{iset}{4};

      ControlInFile = regexprep(XsectionList{iset}{1}, '<CASE>', ControlCase);

      PreSalDiffVname = sprintf('%s_diff', PreSalVname);
      SalDiffVname    = sprintf('%s_diff', SalVname);

      % Read in the variables, split into pre-SAL and SAL sections
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      fprintf('  Reading: %s (%s)\n', ControlInFile, InVname);

      VAR = squeeze(h5read(InFile, InVname));
      X   = squeeze(h5read(InFile, '/x_coords'));
      Y   = squeeze(h5read(InFile, '/y_coords'));
      Z   = squeeze(h5read(InFile, '/z_coords'));
      T   = squeeze(h5read(InFile, '/t_coords'));

      CNTL_VAR = squeeze(h5read(ControlInFile, InVname));

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

      % Determine the number of dimensions. When you have a vector (size is [ 1 n ] or [ n 1 ]),
      % Ndims will be set to 2 when you really want it to be 1.
      Ndims = ndims(VAR);
      if (Ndims == 2)
        if ((size(VAR,1) == 1) || (size(VAR,2) == 1))
          Ndims = 1
        end
      end
    
      % Form time indices for selecting the pre-SAL and SAL time periods
      SIM_TIME = (T ./ 3600) - 42;  % in hours, 0 through 60

      PT1 = find(SIM_TIME >= PreSalTstart, 1, 'first');
      PT2 = find(SIM_TIME <= PreSalTend,   1, 'last');

      ST1 = find(SIM_TIME >= SalTstart, 1, 'first');
      ST2 = find(SIM_TIME <= SalTend,   1, 'last');

      % Want to average across time, which is always the last dimension.
      % Reshape the output variables to 4D for Grads viewing.
      % Radius is the x dimension.
      switch(Ndims)
        case 1
          % var is 1D, (t)
          % results is a scalar
          P_VAR = squeeze(mean(VAR(PT1:PT2)));
          S_VAR = squeeze(mean(VAR(ST1:ST2)));

          P_CNTL_VAR = squeeze(mean(CNTL_VAR(PT1:PT2)));
          S_CNTL_VAR = squeeze(mean(CNTL_VAR(ST1:ST2)));

          Vsize = 1;
          DimOrder = { };
        case 2
          % var is 2D, (r,t)
          % result is (x)
          P_VAR = squeeze(mean(VAR(:,PT1:PT2),2));
          S_VAR = squeeze(mean(VAR(:,ST1:ST2),2));

          P_CNTL_VAR = squeeze(mean(CNTL_VAR(:,PT1:PT2),2));
          S_CNTL_VAR = squeeze(mean(CNTL_VAR(:,ST1:ST2),2));

          Vsize = Nx;
          DimOrder = { 'x' };
        case 3
          % var is 3D, (r,z,t)
          % result is (x,z)
          P_VAR = squeeze(mean(VAR(:,:,PT1:PT2),3));
          S_VAR = squeeze(mean(VAR(:,:,ST1:ST2),3));

          P_CNTL_VAR = squeeze(mean(CNTL_VAR(:,:,PT1:PT2),3));
          S_CNTL_VAR = squeeze(mean(CNTL_VAR(:,:,ST1:ST2),3));

          Vsize = size(P_VAR);
          DimOrder = { 'x' 'z' };
      end

      P_DIFF_VAR = P_VAR - P_CNTL_VAR;
      S_DIFF_VAR = S_VAR - S_CNTL_VAR;

      % Collect output vars into one case specific file.
      fprintf('  Writing: %s (%s, %s)\n', OutFile, PreSalVname, SalVname);
      h5create(OutFile, PreSalVname, Vsize);
      h5create(OutFile, SalVname,    Vsize);

      h5write(OutFile, PreSalVname, P_VAR);
      h5write(OutFile, SalVname,    S_VAR);

      fprintf('  Writing: %s (%s, %s)\n', OutFile, PreSalDiffVname, SalDiffVname);
      h5create(OutFile, PreSalDiffVname, Vsize);
      h5create(OutFile, SalDiffVname,    Vsize);

      h5write(OutFile, PreSalDiffVname, P_DIFF_VAR);
      h5write(OutFile, SalDiffVname,    S_DIFF_VAR);


      AttachDimensionsXyzt(OutFile, PreSalVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, SalVname,    DimOrder, Xname, Yname, Zname, Tname);

      AttachDimensionsXyzt(OutFile, PreSalDiffVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, SalDiffVname,    DimOrder, Xname, Yname, Zname, Tname);

      fprintf('\n');
    end
  end
end
