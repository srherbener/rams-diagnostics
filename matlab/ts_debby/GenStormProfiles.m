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

  % Two time periods:
  %    Pre SAL: T = 10 to 30 h (after RI, before encounter SAL)
  %        SAL: T = 40 to 60 h (during SAL)

  PreSalTstart = 10;
  PreSalTend   = 30;
  SalTstart    = 40;
  SalTend      = 60;

  % Two radial regions:
  %   CORE: 0 100 km
  %   RB:   100 250 km
  CoreRstart = 0;
  CoreRend   = 100;
  RbRstart   = 100;
  RbRend     = 250;

  % Description of profiles
  ProfileList = {
    % in_file in_var pre_sal_core_var pre_sal_rb_var sal_core_var sal_rb_var
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_t_fa' '/ps_speed_t_fa' '/s_speed_t_fa' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_t_wm' '/ps_speed_t_wm' '/s_speed_t_wm' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_r_fa' '/ps_speed_r_fa' '/s_speed_r_fa' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_r_wm' '/ps_speed_r_wm' '/s_speed_r_wm' }
%
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/avg_updraft_fa' '/ps_updraft_fa' '/s_updraft_fa' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/avg_updraft_wm' '/ps_updraft_wm' '/s_updraft_wm' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/avg_dndraft_fa' '/ps_dndraft_fa' '/s_dndraft_fa' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/avg_dndraft_wm' '/ps_dndraft_wm' '/s_dndraft_wm' }
%
%    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/avg_theta_e_fa' '/ps_theta_e_fa' '/s_theta_e_fa' }
%    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/avg_theta_e_wm' '/ps_theta_e_wm' '/s_theta_e_wm' }
%
%    { 'DIAGS/hist_meas_pcprate_<CASE>.h5' '/avg_pcprate_fa' '/ps_pcprate_fa' '/s_pcprate_fa' }
%    { 'DIAGS/hist_meas_pcprate_<CASE>.h5' '/avg_pcprate_wm' '/ps_pcprate_wm' '/s_pcprate_wm' }
%
%    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/avg_vint_cond_fa' '/ps_vint_cond_fa' '/s_vint_cond_fa' }
%    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/avg_vint_cond_wm' '/ps_vint_cond_wm' '/s_vint_cond_wm' }
%
%    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d1_mass_fa' '/ps_d1_mass_fa' '/s_d1_mass_fa' }
%    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d1_mass_wm' '/ps_d1_mass_wm' '/s_d1_mass_wm' }
%    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d1_num_fa'  '/ps_d1_num_fa'  '/s_d1_num_fa' }
%    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d1_num_wm'  '/ps_d1_num_wm'  '/s_d1_num_wm' }
%    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d2_mass_fa' '/ps_d2_mass_fa' '/s_d2_mass_fa' }
%    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d2_mass_wm' '/ps_d2_mass_wm' '/s_d2_mass_wm' }
%    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d2_num_fa'  '/ps_d2_num_fa'  '/s_d2_num_fa' }
%    { 'DIAGS/hist_meas_dust_<CASE>.h5' '/avg_d2_num_wm'  '/ps_d2_num_wm'  '/s_d2_num_wm' }

    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/avg_cloud_fa'  '/ps_core_cloud_fa' '/ps_rb_cloud_fa' '/s_core_cloud_fa' '/s_rb_cloud_fa' }
    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/avg_cloud_wm'  '/ps_core_cloud_wm' '/ps_rb_cloud_wm' '/s_core_cloud_wm' '/s_rb_cloud_wm' }

%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/avg_rain_fa'  '/ps_rain_fa'  '/s_rain_fa' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/avg_rain_wm'  '/ps_rain_wm'  '/s_rain_wm' }

%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/avg_pris_fa'  '/ps_pris_fa'  '/s_pris_fa' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/avg_pris_wm'  '/ps_pris_wm'  '/s_pris_wm' }

%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/avg_snow_fa'  '/ps_snow_fa'  '/s_snow_fa' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/avg_snow_wm'  '/ps_snow_wm'  '/s_snow_wm' }

%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/avg_aggr_fa'  '/ps_aggr_fa'  '/s_aggr_fa' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/avg_aggr_wm'  '/ps_aggr_wm'  '/s_aggr_wm' }

%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/avg_graup_fa'  '/ps_graup_fa'  '/s_graup_fa' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/avg_graup_wm'  '/ps_graup_wm'  '/s_graup_wm' }

%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/avg_hail_fa'  '/ps_hail_fa'  '/s_hail_fa' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/avg_hail_wm'  '/ps_hail_wm'  '/s_hail_wm' }

%    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/avg_lhf_cool_fa'  '/ps_lhf_cool_fa'  '/s_lhf_cool_fa' }
%    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/avg_lhf_cool_wm'  '/ps_lhf_cool_wm'  '/s_lhf_cool_wm' }
%    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/avg_lhf_heat_fa'  '/ps_lhf_heat_fa'  '/s_lhf_heat_fa' }
%    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/avg_lhf_heat_wm'  '/ps_lhf_heat_wm'  '/s_lhf_heat_wm' }
%
%    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/avg_lhv_cool_fa'  '/ps_lhv_cool_fa'  '/s_lhv_cool_fa' }
%    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/avg_lhv_cool_wm'  '/ps_lhv_cool_wm'  '/s_lhv_cool_wm' }
%    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/avg_lhv_heat_fa'  '/ps_lhv_heat_fa'  '/s_lhv_heat_fa' }
%    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/avg_lhv_heat_wm'  '/ps_lhv_heat_wm'  '/s_lhv_heat_wm' }
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
      InFile          = regexprep(ProfileList{iset}{1}, '<CASE>', Case);
      InVname         = ProfileList{iset}{2};
      PreSalCoreVname = ProfileList{iset}{3};
      PreSalRbVname   = ProfileList{iset}{4};
      SalCoreVname    = ProfileList{iset}{5};
      SalRbVname      = ProfileList{iset}{6};

      ControlInFile = regexprep(ProfileList{iset}{1}, '<CASE>', ControlCase);

      PreSalCoreDiffVname = sprintf('%s_diff', PreSalCoreVname);
      PreSalRbDiffVname   = sprintf('%s_diff', PreSalRbVname);
      SalCoreDiffVname    = sprintf('%s_diff', SalCoreVname);
      SalRbDiffVname      = sprintf('%s_diff', SalRbVname);

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
  
        CreateDimensionsXyzt(OutFile, [ 1 ], Y, Z, [ 1 ], Xname, Yname, Zname, Tname);
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

      % Form radius indices for selecting the core and rainband regions
      RADIUS = X ./ 1000;  % in km

      CR1 = find(RADIUS >= CoreRstart, 1, 'first');
      CR2 = find(RADIUS <= CoreRend,   1, 'last');

      RR1 = find(RADIUS >= RbRstart, 1, 'first');
      RR2 = find(RADIUS <= RbRend,   1, 'last');

      % Want to average across radius and time
      % Radius is the x dimension.
      switch(Ndims)
        case 1
          % var is 1D, (t)
          % result is a scalar
          PC_VAR = squeeze(mean(VAR(PT1:PT2)));
          SC_VAR = squeeze(mean(VAR(ST1:ST2)));
          PR_VAR = PC_VAR;
          SR_VAR = SC_VAR;

          PC_CNTL_VAR = squeeze(mean(CNTL_VAR(PT1:PT2)));
          SC_CNTL_VAR = squeeze(mean(CNTL_VAR(ST1:ST2)));
          PR_CNTL_VAR = PC_CNTL_VAR;
          SR_CNTL_VAR = SC_CNTL_VAR;

          Vsize = 1;
          DimOrder = { };
        case 2
          % var is 2D, (r,t)
          % result is scalar
          PC_VAR = squeeze(mean(mean(VAR(CR1:CR2,PT1:PT2),2),1));
          SC_VAR = squeeze(mean(mean(VAR(CR1:CR2,ST1:ST2),2),1));
          PR_VAR = squeeze(mean(mean(VAR(RR1:RR2,PT1:PT2),2),1));
          SR_VAR = squeeze(mean(mean(VAR(RR1:RR2,ST1:ST2),2),1));

          PC_CNTL_VAR = squeeze(mean(mean(CNTL_VAR(CR1:CR2,PT1:PT2),2),1));
          SC_CNTL_VAR = squeeze(mean(mean(CNTL_VAR(CR1:CR2,ST1:ST2),2),1));
          PR_CNTL_VAR = squeeze(mean(mean(CNTL_VAR(RR1:RR2,PT1:PT2),2),1));
          SR_CNTL_VAR = squeeze(mean(mean(CNTL_VAR(RR1:RR2,ST1:ST2),2),1));

          Vsize = 1;
          DimOrder = { };
        case 3
          % var is 3D, (r,z,t)
          % result is (z)
          PC_VAR = squeeze(mean(mean(VAR(CR1:CR2,:,PT1:PT2),3),1));
          SC_VAR = squeeze(mean(mean(VAR(CR1:CR2,:,ST1:ST2),3),1));
          PR_VAR = squeeze(mean(mean(VAR(RR1:RR2,:,PT1:PT2),3),1));
          SR_VAR = squeeze(mean(mean(VAR(RR1:RR2,:,ST1:ST2),3),1));

          PC_CNTL_VAR = squeeze(mean(mean(CNTL_VAR(CR1:CR2,:,PT1:PT2),3),1));
          SC_CNTL_VAR = squeeze(mean(mean(CNTL_VAR(CR1:CR2,:,ST1:ST2),3),1));
          PR_CNTL_VAR = squeeze(mean(mean(CNTL_VAR(RR1:RR2,:,PT1:PT2),3),1));
          SR_CNTL_VAR = squeeze(mean(mean(CNTL_VAR(RR1:RR2,:,ST1:ST2),3),1));

          Vsize = Nz;
          DimOrder = { 'z' };
      end

      PC_DIFF_VAR = PC_VAR - PC_CNTL_VAR;
      SC_DIFF_VAR = SC_VAR - SC_CNTL_VAR;
      PR_DIFF_VAR = PR_VAR - PR_CNTL_VAR;
      SR_DIFF_VAR = SR_VAR - SR_CNTL_VAR;

      % Collect output vars into one case specific file.
      fprintf('  Writing: %s (%s, %s, %s, %s)\n', OutFile, PreSalCoreVname, PreSalRbVname, SalCoreVname, SalRbVname);
      h5create(OutFile, PreSalCoreVname, Vsize);
      h5create(OutFile, SalCoreVname,    Vsize);
      h5create(OutFile, PreSalRbVname,   Vsize);
      h5create(OutFile, SalRbVname,      Vsize);

      h5write(OutFile, PreSalCoreVname, PC_VAR);
      h5write(OutFile, SalCoreVname,    SC_VAR);
      h5write(OutFile, PreSalRbVname,   PR_VAR);
      h5write(OutFile, SalRbVname,      SR_VAR);

      fprintf('  Writing: %s (%s, %s, %s, %s)\n', OutFile, PreSalCoreDiffVname, PreSalRbVname, SalCoreDiffVname, SalRbVname);
      h5create(OutFile, PreSalCoreDiffVname, Vsize);
      h5create(OutFile, SalCoreDiffVname,    Vsize);
      h5create(OutFile, PreSalRbDiffVname,   Vsize);
      h5create(OutFile, SalRbDiffVname,      Vsize);

      h5write(OutFile, PreSalCoreDiffVname, PC_DIFF_VAR);
      h5write(OutFile, SalCoreDiffVname,    SC_DIFF_VAR);
      h5write(OutFile, PreSalRbDiffVname,   PR_DIFF_VAR);
      h5write(OutFile, SalRbDiffVname,      SR_DIFF_VAR);


      AttachDimensionsXyzt(OutFile, PreSalCoreVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, SalCoreVname,    DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, PreSalRbVname,   DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, SalRbVname,      DimOrder, Xname, Yname, Zname, Tname);

      AttachDimensionsXyzt(OutFile, PreSalCoreDiffVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, SalCoreDiffVname,    DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, PreSalRbDiffVname,   DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, SalRbDiffVname,      DimOrder, Xname, Yname, Zname, Tname);

      fprintf('\n');
    end
  end
end
