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
    { 'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_ps_speed_t' '/all_ps_speed_t', 'z' }
    { 'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_ps_speed_r' '/all_ps_speed_r', 'z' }
    { 'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_s_speed_t' '/all_s_speed_t', 'z' }
    { 'DIAGS/hist_meas_az_speed_<CASE>.h5' '/all_s_speed_r' '/all_s_speed_r', 'z' }

    { 'DIAGS/hist_meas_az_p_speed_<CASE>.h5' '/all_ps_speed_t' '/all_p_ps_speed_t', 'p' }
    { 'DIAGS/hist_meas_az_p_speed_<CASE>.h5' '/all_ps_speed_r' '/all_p_ps_speed_r', 'p' }
    { 'DIAGS/hist_meas_az_p_speed_<CASE>.h5' '/all_s_speed_t' '/all_p_s_speed_t', 'p' }
    { 'DIAGS/hist_meas_az_p_speed_<CASE>.h5' '/all_s_speed_r' '/all_p_s_speed_r', 'p' }

    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_ps_updraft' '/all_ps_updraft', 'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_ps_dndraft' '/all_ps_dndraft', 'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_s_updraft' '/all_s_updraft', 'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_s_dndraft' '/all_s_dndraft', 'z' }

    { 'DIAGS/hist_meas_az_p_w_<CASE>.h5' '/all_ps_updraft' '/all_p_ps_updraft', 'p' }
    { 'DIAGS/hist_meas_az_p_w_<CASE>.h5' '/all_ps_dndraft' '/all_p_ps_dndraft', 'p' }
    { 'DIAGS/hist_meas_az_p_w_<CASE>.h5' '/all_s_updraft' '/all_p_s_updraft', 'p' }
    { 'DIAGS/hist_meas_az_p_w_<CASE>.h5' '/all_s_dndraft' '/all_p_s_dndraft', 'p' }

    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/all_ps_theta_e' '/all_ps_theta_e', 'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/all_s_theta_e' '/all_s_theta_e', 'z' }

%    { 'DIAGS/hist_meas_az_theta_v_<CASE>.h5' '/all_ps_theta_v' '/all_ps_theta_v', 'z' }
%    { 'DIAGS/hist_meas_az_theta_v_<CASE>.h5' '/all_s_theta_v' '/all_s_theta_v', 'z' }

%    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/all_ps_theta'  '/all_ps_theta', 'z' }
%    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/all_s_theta'  '/all_s_theta', 'z' }

%    { 'DIAGS/hist_meas_az_theta_rho_<CASE>.h5' '/all_ps_theta_rho' '/all_ps_theta_rho', 'z' }
%    { 'DIAGS/hist_meas_az_theta_rho_<CASE>.h5' '/all_s_theta_rho' '/all_s_theta_rho', 'z' }

    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/all_ps_tempc'  '/all_ps_tempc', 'z' }
    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/all_s_tempc'  '/all_s_tempc', 'z' }

%    { 'DIAGS/hist_meas_az_dewptc_<CASE>.h5' '/all_ps_dewptc'  '/all_ps_dewptc', 'z' }
%    { 'DIAGS/hist_meas_az_dewptc_<CASE>.h5' '/all_s_dewptc'  '/all_s_dewptc', 'z' }

%    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_core_ps_buoy_pos_acc_xsect'  '/all_core_ps_buoy_pos_acc', 'z' }
%    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_core_s_buoy_pos_acc_xsect'   '/all_core_s_buoy_pos_acc' , 'z' }
%    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_rb_ps_buoy_pos_acc_xsect'    '/all_rb_ps_buoy_pos_acc'  , 'z' }
%    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_rb_s_buoy_pos_acc_xsect'     '/all_rb_s_buoy_pos_acc'   , 'z' }

%    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_core_ps_buoy_neg_acc_xsect'  '/all_core_ps_buoy_neg_acc', 'z' }
%    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_core_s_buoy_neg_acc_xsect'   '/all_core_s_buoy_neg_acc' , 'z' }
%    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_rb_ps_buoy_neg_acc_xsect'    '/all_rb_ps_buoy_neg_acc'  , 'z' }
%    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_rb_s_buoy_neg_acc_xsect'     '/all_rb_s_buoy_neg_acc'   , 'z' }

    { 'DIAGS/hist_meas_az_relhum_<CASE>.h5' '/all_ps_relhum'  '/all_ps_relhum', 'z' }
    { 'DIAGS/hist_meas_az_relhum_<CASE>.h5' '/all_s_relhum'  '/all_s_relhum', 'z' }

    { 'DIAGS/hist_meas_az_pcprate_<CASE>.h5' '/all_ps_pcprate' '/all_ps_pcprate', 'z' }
    { 'DIAGS/hist_meas_az_pcprate_<CASE>.h5' '/all_s_pcprate' '/all_s_pcprate', 'z' }

%    { 'DIAGS/hist_meas_az_vint_cond_<CASE>.h5' '/all_ps_vint_cond' '/all_ps_vint_cond', 'z' }
%    { 'DIAGS/hist_meas_az_vint_cond_<CASE>.h5' '/all_s_vint_cond' '/all_s_vint_cond', 'z' }

    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/all_ps_vapor' '/all_ps_vapor', 'z' }
    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/all_s_vapor' '/all_s_vapor', 'z' }
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/lead_ps_vapor' '/lead_ps_vapor', 'z' }
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/lead_s_vapor' '/lead_s_vapor', 'z' }

%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_ps_d1_num'  '/all_ps_d1_num' , 'z' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_ps_d2_num'  '/all_ps_d2_num' , 'z' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_s_d1_num'  '/all_s_d1_num' , 'z' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_s_d2_num'  '/all_s_d2_num' , 'z' }

%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_cloud' '/all_ps_dust_cloud', 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_rain'  '/all_ps_dust_rain' , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_pris'  '/all_ps_dust_pris' , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_snow'  '/all_ps_dust_snow' , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_aggr'  '/all_ps_dust_aggr' , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_graup' '/all_ps_dust_graup', 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_hail'  '/all_ps_dust_hail' , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_ps_dust_hydro' '/all_ps_dust_hydro', 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_cloud'  '/all_s_dust_cloud' , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_rain'   '/all_s_dust_rain'  , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_pris'   '/all_s_dust_pris'  , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_snow'   '/all_s_dust_snow'  , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_aggr'   '/all_s_dust_aggr'  , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_graup'  '/all_s_dust_graup' , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_hail'   '/all_s_dust_hail'  , 'z' }
%    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_s_dust_hydro'  '/all_s_dust_hydro' , 'z' }

    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_ps_cloud_mass' '/all_ps_cloud_mass', 'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_s_cloud_mass'  '/all_s_cloud_mass', 'z' }
%    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_ps_cloud_num'  '/all_ps_cloud_num', 'z' }
%    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_s_cloud_num'   '/all_s_cloud_num', 'z' }
%    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_ps_cloud_diam' '/all_ps_cloud_diam', 'z' }
%    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_s_cloud_diam'  '/all_s_cloud_diam', 'z' }
    { 'DIAGS/hist_meas_az_p_cloud_<CASE>.h5' '/all_ps_cloud_mass' '/all_p_ps_cloud_mass', 'p' }
    { 'DIAGS/hist_meas_az_p_cloud_<CASE>.h5' '/all_s_cloud_mass'  '/all_p_s_cloud_mass', 'p' }

    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_ps_rain_mass' '/all_ps_rain_mass', 'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_s_rain_mass'  '/all_s_rain_mass', 'z' }
%    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_ps_rain_num'  '/all_ps_rain_num', 'z' }
%    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_s_rain_num'   '/all_s_rain_num', 'z' }
%    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_ps_rain_diam' '/all_ps_rain_diam', 'z' }
%    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_s_rain_diam'  '/all_s_rain_diam', 'z' }

%    { 'DIAGS/hist_meas_az_pris_<CASE>.h5' '/all_ps_pris_mass'      '/all_ps_pris_mass', 'z' }
%    { 'DIAGS/hist_meas_az_pris_<CASE>.h5' '/all_s_pris_mass'       '/all_s_pris_mass', 'z' }
%
%    { 'DIAGS/hist_meas_az_snow_<CASE>.h5' '/all_ps_snow_mass'      '/all_ps_snow_mass', 'z' }
%    { 'DIAGS/hist_meas_az_snow_<CASE>.h5' '/all_s_snow_mass'       '/all_s_snow_mass', 'z' }
%
%    { 'DIAGS/hist_meas_az_aggr_<CASE>.h5' '/all_ps_aggr_mass'      '/all_ps_aggr_mass', 'z' }
%    { 'DIAGS/hist_meas_az_aggr_<CASE>.h5' '/all_s_aggr_mass'       '/all_s_aggr_mass', 'z' }
%
%    { 'DIAGS/hist_meas_az_graup_<CASE>.h5' '/all_ps_graup_mass'      '/all_ps_graup_mass', 'z' }
%    { 'DIAGS/hist_meas_az_graup_<CASE>.h5' '/all_s_graup_mass'       '/all_s_graup_mass', 'z' }
%
%    { 'DIAGS/hist_meas_az_hail_<CASE>.h5' '/all_ps_hail_mass'      '/all_ps_hail_mass', 'z' }
%    { 'DIAGS/hist_meas_az_hail_<CASE>.h5' '/all_s_hail_mass'       '/all_s_hail_mass', 'z' }
%
%    { 'DIAGS/hist_meas_az_tcond_<CASE>.h5' '/all_ps_tcond_mass' '/all_ps_tcond_mass', 'z' }
%    { 'DIAGS/hist_meas_az_tcond_<CASE>.h5' '/all_s_tcond_mass'  '/all_s_tcond_mass', 'z' }

%    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/all_ps_lhf_cool' '/all_ps_lhf_cool', 'z' }
%    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/all_s_lhf_cool'  '/all_s_lhf_cool', 'z' }

%    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/all_ps_lhf_heat' '/all_ps_lhf_heat', 'z' }
%    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/all_s_lhf_heat'  '/all_s_lhf_heat', 'z' }

%    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/all_ps_lhv_cool' '/all_ps_lhv_cool', 'z' }
%    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/all_s_lhv_cool'  '/all_s_lhv_cool', 'z' }

%    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/all_ps_lhv_heat' '/all_ps_lhv_heat', 'z' }
%    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/all_s_lhv_heat'  '/all_s_lhv_heat', 'z' }

    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/all_ps_liq_evap' '/all_ps_liq_evap', 'z' }
    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/all_s_liq_evap'  '/all_s_liq_evap', 'z' }

    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/all_ps_liq_cond' '/all_ps_liq_cond', 'z' }
    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/all_s_liq_cond'  '/all_s_liq_cond', 'z' }

    { 'DIAGS/hist_meas_az_ice_dep_<CASE>.h5' '/all_ps_ice_dep' '/all_ps_ice_dep', 'z' }
    { 'DIAGS/hist_meas_az_ice_dep_<CASE>.h5' '/all_s_ice_dep'  '/all_s_ice_dep', 'z' }

    { 'DIAGS/hist_meas_az_ice_sub_<CASE>.h5' '/all_ps_ice_sub' '/all_ps_ice_sub', 'z' }
    { 'DIAGS/hist_meas_az_ice_sub_<CASE>.h5' '/all_s_ice_sub'  '/all_s_ice_sub', 'z' }

%    { 'DIAGS/hist_meas_az_ice_melt_<CASE>.h5' '/all_ps_ice_melt' '/all_ps_ice_melt', 'z' }
%    { 'DIAGS/hist_meas_az_ice_melt_<CASE>.h5' '/all_s_ice_melt'  '/all_s_ice_melt', 'z' }

%    { 'DIAGS/hist_meas_az_cloud_evap_<CASE>.h5' '/all_ps_cloud_evap' '/all_ps_cloud_evap', 'z' }
%    { 'DIAGS/hist_meas_az_cloud_evap_<CASE>.h5' '/all_s_cloud_evap'  '/all_s_cloud_evap', 'z' }
%
%    { 'DIAGS/hist_meas_az_cloud_cond_<CASE>.h5' '/all_ps_cloud_cond' '/all_ps_cloud_cond', 'z' }
%    { 'DIAGS/hist_meas_az_cloud_cond_<CASE>.h5' '/all_s_cloud_cond'  '/all_s_cloud_cond', 'z' }
%
%    { 'DIAGS/hist_meas_az_rain_evap_<CASE>.h5' '/all_ps_rain_evap' '/all_ps_rain_evap', 'z' }
%    { 'DIAGS/hist_meas_az_rain_evap_<CASE>.h5' '/all_s_rain_evap'  '/all_s_rain_evap', 'z' }
%
%    { 'DIAGS/hist_meas_az_rain_cond_<CASE>.h5' '/all_ps_rain_cond' '/all_ps_rain_cond', 'z' }
%    { 'DIAGS/hist_meas_az_rain_cond_<CASE>.h5' '/all_s_rain_cond'  '/all_s_rain_cond', 'z' }
%
%    { 'DIAGS/hist_meas_az_pris_sub_<CASE>.h5' '/all_ps_pris_sub' '/all_ps_pris_sub', 'z' }
%    { 'DIAGS/hist_meas_az_pris_sub_<CASE>.h5' '/all_s_pris_sub'  '/all_s_pris_sub', 'z' }
%
%    { 'DIAGS/hist_meas_az_pris_dep_<CASE>.h5' '/all_ps_pris_dep' '/all_ps_pris_dep', 'z' }
%    { 'DIAGS/hist_meas_az_pris_dep_<CASE>.h5' '/all_s_pris_dep'  '/all_s_pris_dep', 'z' }
%
%    { 'DIAGS/hist_meas_az_snow_sub_<CASE>.h5' '/all_ps_snow_sub' '/all_ps_snow_sub', 'z' }
%    { 'DIAGS/hist_meas_az_snow_sub_<CASE>.h5' '/all_s_snow_sub'  '/all_s_snow_sub', 'z' }
%
%    { 'DIAGS/hist_meas_az_snow_dep_<CASE>.h5' '/all_ps_snow_dep' '/all_ps_snow_dep', 'z' }
%    { 'DIAGS/hist_meas_az_snow_dep_<CASE>.h5' '/all_s_snow_dep'  '/all_s_snow_dep', 'z' }
%
%    { 'DIAGS/hist_meas_az_aggr_sub_<CASE>.h5' '/all_ps_aggr_sub' '/all_ps_aggr_sub', 'z' }
%    { 'DIAGS/hist_meas_az_aggr_sub_<CASE>.h5' '/all_s_aggr_sub'  '/all_s_aggr_sub', 'z' }
%
%    { 'DIAGS/hist_meas_az_aggr_dep_<CASE>.h5' '/all_ps_aggr_dep' '/all_ps_aggr_dep', 'z' }
%    { 'DIAGS/hist_meas_az_aggr_dep_<CASE>.h5' '/all_s_aggr_dep'  '/all_s_aggr_dep', 'z' }
%
%    { 'DIAGS/hist_meas_az_graup_sub_<CASE>.h5' '/all_ps_graup_sub' '/all_ps_graup_sub', 'z' }
%    { 'DIAGS/hist_meas_az_graup_sub_<CASE>.h5' '/all_s_graup_sub'  '/all_s_graup_sub', 'z' }
%
%    { 'DIAGS/hist_meas_az_graup_dep_<CASE>.h5' '/all_ps_graup_dep' '/all_ps_graup_dep', 'z' }
%    { 'DIAGS/hist_meas_az_graup_dep_<CASE>.h5' '/all_s_graup_dep'  '/all_s_graup_dep', 'z' }
%
%    { 'DIAGS/hist_meas_az_hail_sub_<CASE>.h5' '/all_ps_hail_sub' '/all_ps_hail_sub', 'z' }
%    { 'DIAGS/hist_meas_az_hail_sub_<CASE>.h5' '/all_s_hail_sub'  '/all_s_hail_sub', 'z' }
%
%    { 'DIAGS/hist_meas_az_hail_dep_<CASE>.h5' '/all_ps_hail_dep' '/all_ps_hail_dep', 'z' }
%    { 'DIAGS/hist_meas_az_hail_dep_<CASE>.h5' '/all_s_hail_dep'  '/all_s_hail_dep', 'z' }
%
%    { 'DIAGS/hist_meas_az_cloud_rime_<CASE>.h5' '/all_ps_cloud_rime' '/all_ps_cloud_rime', 'z' }
%    { 'DIAGS/hist_meas_az_cloud_rime_<CASE>.h5' '/all_s_cloud_rime'  '/all_s_cloud_rime', 'z' }
%
%    { 'DIAGS/hist_meas_az_rain2ice_<CASE>.h5' '/all_ps_rain2ice' '/all_ps_rain2ice', 'z' }
%    { 'DIAGS/hist_meas_az_rain2ice_<CASE>.h5' '/all_s_rain2ice'  '/all_s_rain2ice', 'z' }

    };
  Nsets = length(XsectionList);

  Xname = '/x_coords';
  Yname = '/y_coords';
  Zname = '/z_coords';
  Pname = '/p_coords';
  Tname = '/t_coords';

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

    ixcoord = 0;
    iycoord = 0;
    izcoord = 0;
    ipcoord = 0;
    itcoord = 0;
    for iset = 1:Nsets
      InFile     = regexprep(XsectionList{iset}{1}, '<CASE>', Case);
      InVname    = XsectionList{iset}{2};
      OutVname   = XsectionList{iset}{3};
      VcoordType = XsectionList{iset}{4};

      % skip this profile set if doing dust and on a NODUST case
      if ((~isempty(regexp(Case, 'NODUST'))) && ...
          ((~isempty(regexp(InVname, 'd[12]_num'))) || (~isempty(regexp(InVname, 'dust_'))) || (~isempty(regexp(InVname, 'dustifn_')))))
        continue
      end

      % Read in the variables
      fprintf('  Reading: %s (%s)\n', InFile, InVname);

      VAR = squeeze(h5read(InFile, InVname));
      X   = squeeze(h5read(InFile, '/x_coords'));
      Y   = squeeze(h5read(InFile, '/y_coords'));
      Z   = squeeze(h5read(InFile, '/z_coords'));
      T   = squeeze(h5read(InFile, '/t_coords'));

      ixcoord = ixcoord + 1;
      iycoord = iycoord + 1;
      if (VcoordType == 'z')
        izcoord = izcoord + 1;
      elseif (VcoordType == 'p')
        ipcoord = ipcoord + 1;
      end
      itcoord = itcoord + 1;

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
      if (ixcoord == 1)
        CreateDimension(OutFile, X, Xname, 'x');
        NotateDimension(OutFile, Xname, 'x');
      end
      if (iycoord == 1)
        CreateDimension(OutFile, Y, Yname, 'y');
        NotateDimension(OutFile, Yname, 'y');
      end
      if (izcoord == 1)
        CreateDimension(OutFile, Z, Zname, 'z');
        NotateDimension(OutFile, Zname, 'z');
      end
      if (ipcoord == 1)
        CreateDimension(OutFile, Z, Pname, 'p');
        NotateDimension(OutFile, Pname, 'p');
      end
      if (itcoord == 1)
        CreateDimension(OutFile, T, Tname, 't');
        NotateDimension(OutFile, Tname, 't');
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

      if (VcoordType == 'z')
        AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
      elseif (VcoordType == 'p')
        AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Pname, Tname);
      end

      fprintf('\n');
    end
  end
end
