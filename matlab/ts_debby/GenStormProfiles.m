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

  % Description of profiles
  ProfileList = {
    % in_file in_var out_var pre_sal_rb_var sal_core_var sal_rb_var

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/core_ps_d1_num'  '/core_ps_d1_num'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/core_s_d1_num'   '/core_s_d1_num'   'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/rb_ps_d1_num'    '/rb_ps_d1_num'    'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/rb_s_d1_num'     '/rb_s_d1_num'     'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/env_ps_d1_num'   '/env_ps_d1_num'   'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/env_s_d1_num'    '/env_s_d1_num'    'z' }

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/core_ps_d2_num'  '/core_ps_d2_num'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/core_s_d2_num'   '/core_s_d2_num'   'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/rb_ps_d2_num'    '/rb_ps_d2_num'    'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/rb_s_d2_num'     '/rb_s_d2_num'     'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/env_ps_d2_num'   '/env_ps_d2_num'   'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/env_s_d2_num'    '/env_s_d2_num'    'z' }

    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/core_ps_cloud_num'  '/core_ps_cloud_num'  'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/core_s_cloud_num'   '/core_s_cloud_num'   'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/rb_ps_cloud_num'    '/rb_ps_cloud_num'    'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/rb_s_cloud_num'     '/rb_s_cloud_num'     'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/env_ps_cloud_num'   '/env_ps_cloud_num'   'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/env_s_cloud_num'    '/env_s_cloud_num'    'z' }

    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/core_ps_cloud_diam' '/core_ps_cloud_diam'  'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/core_s_cloud_diam'  '/core_s_cloud_diam'   'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/rb_ps_cloud_diam'   '/rb_ps_cloud_diam'    'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/rb_s_cloud_diam'    '/rb_s_cloud_diam'     'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/env_ps_cloud_diam'  '/env_ps_cloud_diam'   'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/env_s_cloud_diam'   '/env_s_cloud_diam'    'z' }

    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/core_ps_lhf_cool' '/core_ps_lhf_cool' 'z' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/core_s_lhf_cool'  '/core_s_lhf_cool'  'z' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/rb_ps_lhf_cool'   '/rb_ps_lhf_cool'   'z' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/rb_s_lhf_cool'    '/rb_s_lhf_cool'    'z' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/env_ps_lhf_cool'  '/env_ps_lhf_cool'  'z' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/env_s_lhf_cool'   '/env_s_lhf_cool'   'z' }

    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/core_ps_lhf_heat' '/core_ps_lhf_heat' 'z' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/core_s_lhf_heat'  '/core_s_lhf_heat'  'z' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/rb_ps_lhf_heat'   '/rb_ps_lhf_heat'   'z' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/rb_s_lhf_heat'    '/rb_s_lhf_heat'    'z' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/env_ps_lhf_heat'  '/env_ps_lhf_heat'  'z' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/env_s_lhf_heat'   '/env_s_lhf_heat'   'z' }

    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/core_ps_lhv_cool' '/core_ps_lhv_cool' 'z' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/core_s_lhv_cool'  '/core_s_lhv_cool'  'z' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/rb_ps_lhv_cool'   '/rb_ps_lhv_cool'   'z' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/rb_s_lhv_cool'    '/rb_s_lhv_cool'    'z' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/env_ps_lhv_cool'  '/env_ps_lhv_cool'  'z' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/env_s_lhv_cool'   '/env_s_lhv_cool'   'z' }

    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/core_ps_lhv_heat' '/core_ps_lhv_heat' 'z' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/core_s_lhv_heat'  '/core_s_lhv_heat'  'z' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/rb_ps_lhv_heat'   '/rb_ps_lhv_heat'   'z' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/rb_s_lhv_heat'    '/rb_s_lhv_heat'    'z' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/env_ps_lhv_heat'  '/env_ps_lhv_heat'  'z' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/env_s_lhv_heat'   '/env_s_lhv_heat'   'z' }




%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_core_speed_t' '/ps_core_speed_t'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_rb_speed_t' '/ps_rb_speed_t'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_core_speed_t' '/s_core_speed_t'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_rb_speed_t' '/s_rb_speed_t'  'z' }
%
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_core_speed_r' '/ps_core_speed_r'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_rb_speed_r' '/ps_rb_speed_r'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_core_speed_r' '/s_core_speed_r'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_rb_speed_r' '/s_rb_speed_r'  'z' }
%
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_core_updraft' '/ps_core_updraft'  'z' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_rb_updraft' '/ps_rb_updraft'  'z' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_core_updraft' '/s_core_updraft'  'z' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_rb_updraft' '/s_rb_updraft'  'z' }
%
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_core_dndraft' '/ps_core_dndraft'  'z' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_rb_dndraft' '/ps_rb_dndraft'  'z' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_core_dndraft' '/s_core_dndraft'  'z' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_rb_dndraft' '/s_rb_dndraft'  'z' }
%
%    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/ps_core_theta_e' '/ps_core_theta_e'  'z' }
%    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/ps_rb_theta_e' '/ps_rb_theta_e'  'z' }
%    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/s_core_theta_e' '/s_core_theta_e'  'z' }
%    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/s_rb_theta_e' '/s_rb_theta_e'  'z' }
%
%    { 'DIAGS/hist_meas_theta_<CASE>.h5' '/ps_core_theta' '/ps_core_theta'  'z' }
%    { 'DIAGS/hist_meas_theta_<CASE>.h5' '/ps_rb_theta' '/ps_rb_theta'  'z' }
%    { 'DIAGS/hist_meas_theta_<CASE>.h5' '/s_core_theta' '/s_core_theta'  'z' }
%    { 'DIAGS/hist_meas_theta_<CASE>.h5' '/s_rb_theta' '/s_rb_theta'  'z' }
%
%    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/ps_core_tempc' '/ps_core_tempc'  'z' }
%    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/ps_rb_tempc' '/ps_rb_tempc'  'z' }
%    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/s_core_tempc' '/s_core_tempc'  'z' }
%    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/s_rb_tempc' '/s_rb_tempc'  'z' }
%
%    { 'DIAGS/hist_meas_ts_tempc_<CASE>.h5' '/spath_ps_tempc'  '/spath_ps_tempc'  'z' }
%    { 'DIAGS/hist_meas_ts_tempc_<CASE>.h5' '/smaxcp_ps_tempc' '/smaxcp_ps_tempc' 'z' }
%    { 'DIAGS/hist_meas_ts_tempc_<CASE>.h5' '/spath_s_tempc'   '/spath_s_tempc'   'z' }
%    { 'DIAGS/hist_meas_ts_tempc_<CASE>.h5' '/smaxcp_s_tempc'  '/smaxcp_s_tempc'  'z' }
%
%    { 'DIAGS/hist_meas_az_pcprate_<CASE>.h5' '/ps_pcprate' '/ps_pcprate'  'r' }
%    { 'DIAGS/hist_meas_az_pcprate_<CASE>.h5' '/s_pcprate' '/s_pcprate'  'r' }
%
%    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/ps_vint_cond' '/ps_vint_cond'  'r' }
%    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/s_vint_cond' '/s_vint_cond'  'r' }
%
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/ps_core_vapor'      '/ps_core_vapor'  'z' }
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/ps_rb_vapor'        '/ps_rb_vapor'  'z' }
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/s_core_vapor'       '/s_core_vapor'  'z' }
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/s_rb_vapor'         '/s_rb_vapor'  'z' }
%
%    { 'DIAGS/hist_meas_ts_vapor_<CASE>.h5' '/spath_ps_vapor'  '/spath_ps_vapor'  'z' }
%    { 'DIAGS/hist_meas_ts_vapor_<CASE>.h5' '/smaxcp_ps_vapor' '/smaxcp_ps_vapor' 'z' }
%    { 'DIAGS/hist_meas_ts_vapor_<CASE>.h5' '/spath_s_vapor'   '/spath_s_vapor'   'z' }
%    { 'DIAGS/hist_meas_ts_vapor_<CASE>.h5' '/smaxcp_s_vapor'  '/smaxcp_s_vapor'  'z' }
%
%
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d1_mass'  '/spath_ps_d1_mass'  'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_ps_d1_mass' '/smaxcp_ps_rb_d1_mass'  'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d1_mass'   '/spath_s_d1_mass'  'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_s_d1_mass'  '/smaxcp_s_rb_d1_mass'  'z' }
%
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d1_num'   '/spath_ps_d1_num'   'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_ps_d1_num'  '/smaxcp_ps_rb_d1_num'   'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d1_num'    '/spath_s_d1_num'   'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_s_d1_num'   '/smaxcp_s_rb_d1_num'   'z' }
%
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d2_mass'  '/spath_ps_d2_mass'  'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_ps_d2_mass' '/smaxcp_ps_rb_d2_mass'  'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d2_mass'   '/spath_s_d2_mass'  'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_s_d2_mass'  '/smaxcp_s_rb_d2_mass'  'z' }
%
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d2_num'   '/spath_ps_d2_num'   'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_ps_d2_num'  '/smaxcp_ps_rb_d2_num'   'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d2_num'    '/spath_s_d2_num'   'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_s_d2_num'   '/smaxcp_s_rb_d2_num'   'z' }
%
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_core_rain'      '/ps_core_rain'   'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_rb_rain'        '/ps_rb_rain'   'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_core_rain'       '/s_core_rain'   'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_rb_rain'         '/s_rb_rain'   'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_core_rain_num'  '/ps_core_rain_num'  'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_rb_rain_num'    '/ps_rb_rain_num'  'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_core_rain_num'   '/s_core_rain_num'  'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_rb_rain_num'     '/s_rb_rain_num'  'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_core_rain_diam' '/ps_core_rain_diam'  'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_rb_rain_diam'   '/ps_rb_rain_diam'  'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_core_rain_diam'  '/s_core_rain_diam'  'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_rb_rain_diam'    '/s_rb_rain_diam'  'z' }
%
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_core_pris'      '/ps_core_pris'   'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_rb_pris'        '/ps_rb_pris'   'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_core_pris'       '/s_core_pris'   'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_rb_pris'         '/s_rb_pris'   'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_core_pris_num'  '/ps_core_pris_num'  'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_rb_pris_num'    '/ps_rb_pris_num'  'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_core_pris_num'   '/s_core_pris_num'  'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_rb_pris_num'     '/s_rb_pris_num'  'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_core_pris_diam' '/ps_core_pris_diam'  'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_rb_pris_diam'   '/ps_rb_pris_diam'  'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_core_pris_diam'  '/s_core_pris_diam'  'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_rb_pris_diam'    '/s_rb_pris_diam'  'z' }
%
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_core_snow'      '/ps_core_snow'   'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_rb_snow'        '/ps_rb_snow'   'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_core_snow'       '/s_core_snow'   'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_rb_snow'         '/s_rb_snow'   'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_core_snow_num'  '/ps_core_snow_num'  'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_rb_snow_num'    '/ps_rb_snow_num'  'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_core_snow_num'   '/s_core_snow_num'  'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_rb_snow_num'     '/s_rb_snow_num'  'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_core_snow_diam' '/ps_core_snow_diam'  'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_rb_snow_diam'   '/ps_rb_snow_diam'  'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_core_snow_diam'  '/s_core_snow_diam'  'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_rb_snow_diam'    '/s_rb_snow_diam'  'z' }
%
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_core_aggr'      '/ps_core_aggr'   'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_rb_aggr'        '/ps_rb_aggr'   'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_core_aggr'       '/s_core_aggr'   'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_rb_aggr'         '/s_rb_aggr'   'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_core_aggr_num'  '/ps_core_aggr_num'  'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_rb_aggr_num'    '/ps_rb_aggr_num'  'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_core_aggr_num'   '/s_core_aggr_num'  'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_rb_aggr_num'     '/s_rb_aggr_num'  'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_core_aggr_diam' '/ps_core_aggr_diam'  'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_rb_aggr_diam'   '/ps_rb_aggr_diam'  'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_core_aggr_diam'  '/s_core_aggr_diam'  'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_rb_aggr_diam'    '/s_rb_aggr_diam'  'z' }
%
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_core_graup'      '/ps_core_graup'   'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_rb_graup'        '/ps_rb_graup'   'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_core_graup'       '/s_core_graup'   'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_rb_graup'         '/s_rb_graup'   'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_core_graup_num'  '/ps_core_graup_num'  'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_rb_graup_num'    '/ps_rb_graup_num'  'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_core_graup_num'   '/s_core_graup_num'  'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_rb_graup_num'     '/s_rb_graup_num'  'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_core_graup_diam' '/ps_core_graup_diam'  'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_rb_graup_diam'   '/ps_rb_graup_diam'  'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_core_graup_diam'  '/s_core_graup_diam'  'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_rb_graup_diam'    '/s_rb_graup_diam'  'z' }
%
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_core_hail'      '/ps_core_hail'   'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_rb_hail'        '/ps_rb_hail'   'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_core_hail'       '/s_core_hail'   'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_rb_hail'         '/s_rb_hail'   'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_core_hail_num'  '/ps_core_hail_num'  'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_rb_hail_num'    '/ps_rb_hail_num'  'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_core_hail_num'   '/s_core_hail_num'  'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_rb_hail_num'     '/s_rb_hail_num'  'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_core_hail_diam' '/ps_core_hail_diam'  'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_rb_hail_diam'   '/ps_rb_hail_diam'  'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_core_hail_diam'  '/s_core_hail_diam'  'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_rb_hail_diam'    '/s_rb_hail_diam'  'z' }
%
%    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/ps_core_tcond'  '/ps_core_tcond'   'z' }
%    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/ps_rb_tcond'  '/ps_rb_tcond'   'z' }
%    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/s_core_tcond'  '/s_core_tcond'   'z' }
%    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/s_rb_tcond'  '/s_rb_tcond'   'z' }
%
%    % Env region profiles
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_env_speed_t' '/ps_env_speed_t'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_env_speed_t' '/s_env_speed_t'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/ps_env_speed_r' '/ps_env_speed_r'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/s_env_speed_r' '/s_env_speed_r'  'z' }
%
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_env_updraft' '/ps_env_updraft'  'z' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_env_updraft' '/s_env_updraft'  'z' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/ps_env_dndraft' '/ps_env_dndraft'  'z' }
%    { 'DIAGS/hist_meas_w_<CASE>.h5' '/s_env_dndraft' '/s_env_dndraft'  'z' }
%
%    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/ps_env_theta_e' '/ps_env_theta_e'  'z' }
%    { 'DIAGS/hist_meas_theta_e_<CASE>.h5' '/s_env_theta_e' '/s_env_theta_e'  'z' }
%
%    { 'DIAGS/hist_meas_theta_<CASE>.h5' '/ps_env_theta' '/ps_env_theta'  'z' }
%    { 'DIAGS/hist_meas_theta_<CASE>.h5' '/s_env_theta' '/s_env_theta'  'z' }
%
%    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/ps_env_tempc' '/ps_env_tempc'  'z' }
%    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/s_env_tempc' '/s_env_tempc'  'z' }
%
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/ps_env_vapor'        '/ps_env_vapor'  'z' }
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/s_env_vapor'         '/s_env_vapor'  'z' }
%
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/ps_env_d1_mass' '/ps_env_d1_mass'  'z' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/s_env_d1_mass' '/s_env_d1_mass'  'z' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/ps_env_d1_num'  '/ps_env_d1_num'   'z' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/s_env_d1_num'  '/s_env_d1_num'   'z' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/ps_env_d2_mass' '/ps_env_d2_mass'  'z' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/s_env_d2_mass' '/s_env_d2_mass'  'z' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/ps_env_d2_num'  '/ps_env_d2_num'   'z' }
%    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/s_env_d2_num'  '/s_env_d2_num'   'z' }
%
%    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/ps_env_cloud'        '/ps_env_cloud'  'z' }
%    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/s_env_cloud'         '/s_env_cloud'  'z' }
%    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/ps_env_cloud_num'    '/ps_env_cloud_num'  'z' }
%    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/s_env_cloud_num'     '/s_env_cloud_num'  'z' }
%    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/ps_env_cloud_diam'   '/ps_env_cloud_diam'  'z' }
%    { 'DIAGS/hist_meas_cloud_<CASE>.h5' '/s_env_cloud_diam'    '/s_env_cloud_diam'  'z' }
%
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_env_rain'        '/ps_env_rain'   'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_env_rain'         '/s_env_rain'   'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_env_rain_num'    '/ps_env_rain_num'  'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_env_rain_num'     '/s_env_rain_num'  'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/ps_env_rain_diam'   '/ps_env_rain_diam'  'z' }
%    { 'DIAGS/hist_meas_rain_<CASE>.h5' '/s_env_rain_diam'    '/s_env_rain_diam'  'z' }
%
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_env_pris'        '/ps_env_pris'   'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_env_pris'         '/s_env_pris'   'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_env_pris_num'    '/ps_env_pris_num'  'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_env_pris_num'     '/s_env_pris_num'  'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/ps_env_pris_diam'   '/ps_env_pris_diam'  'z' }
%    { 'DIAGS/hist_meas_pris_<CASE>.h5' '/s_env_pris_diam'    '/s_env_pris_diam'  'z' }
%
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_env_snow'        '/ps_env_snow'   'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_env_snow'         '/s_env_snow'   'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_env_snow_num'    '/ps_env_snow_num'  'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_env_snow_num'     '/s_env_snow_num'  'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/ps_env_snow_diam'   '/ps_env_snow_diam'  'z' }
%    { 'DIAGS/hist_meas_snow_<CASE>.h5' '/s_env_snow_diam'    '/s_env_snow_diam'  'z' }
%
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_env_aggr'        '/ps_env_aggr'   'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_env_aggr'         '/s_env_aggr'   'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_env_aggr_num'    '/ps_env_aggr_num'  'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_env_aggr_num'     '/s_env_aggr_num'  'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/ps_env_aggr_diam'   '/ps_env_aggr_diam'  'z' }
%    { 'DIAGS/hist_meas_aggr_<CASE>.h5' '/s_env_aggr_diam'    '/s_env_aggr_diam'  'z' }
%
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_env_graup'        '/ps_env_graup'   'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_env_graup'         '/s_env_graup'   'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_env_graup_num'    '/ps_env_graup_num'  'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_env_graup_num'     '/s_env_graup_num'  'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/ps_env_graup_diam'   '/ps_env_graup_diam'  'z' }
%    { 'DIAGS/hist_meas_graup_<CASE>.h5' '/s_env_graup_diam'    '/s_env_graup_diam'  'z' }
%
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_env_hail'        '/ps_env_hail'   'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_env_hail'         '/s_env_hail'   'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_env_hail_num'    '/ps_env_hail_num'  'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_env_hail_num'     '/s_env_hail_num'  'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/ps_env_hail_diam'   '/ps_env_hail_diam'  'z' }
%    { 'DIAGS/hist_meas_hail_<CASE>.h5' '/s_env_hail_diam'    '/s_env_hail_diam'  'z' }
%
%    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/ps_env_tcond'  '/ps_env_tcond'   'z' }
%    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/s_env_tcond'  '/s_env_tcond'   'z' }
%
%    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/ps_env_lhf_cool'  '/ps_env_lhf_cool'   'z' }
%    { 'DIAGS/hist_meas_lhf_cool_<CASE>.h5' '/s_env_lhf_cool'  '/s_env_lhf_cool'   'z' }
%
%    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/ps_env_lhf_heat'  '/ps_env_lhf_heat'   'z' }
%    { 'DIAGS/hist_meas_lhf_heat_<CASE>.h5' '/s_env_lhf_heat'  '/s_env_lhf_heat'   'z' }
%
%    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/ps_env_lhv_cool'  '/ps_env_lhv_cool'   'z' }
%    { 'DIAGS/hist_meas_lhv_cool_<CASE>.h5' '/s_env_lhv_cool'  '/s_env_lhv_cool'   'z' }
%
%    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/ps_env_lhv_heat'  '/ps_env_lhv_heat'   'z' }
%    { 'DIAGS/hist_meas_lhv_heat_<CASE>.h5' '/s_env_lhv_heat'  '/s_env_lhv_heat'   'z' }

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

    icount = 0;
    for iset = 1:Nsets
      InFile   = regexprep(ProfileList{iset}{1}, '<CASE>', Case);
      InVname  = ProfileList{iset}{2};
      OutVname = ProfileList{iset}{3};
      ProfDim  = ProfileList{iset}{4};

      % skip this profile set if doing dust and on a NODUST case
      if ((~isempty(regexp(Case, 'NODUST'))) && (~isempty(regexp(InVname, 'd[12]_num'))))
        continue
      else
        icount = icount + 1;
      end

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

      DIFF_VAR = VAR - CNTL_VAR;

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % If a hydrometeor number concentration, then convert #/kg to #/cm^2.
      if (~isempty(regexp(InVname, 'cloud_num')) || ...
          ~isempty(regexp(InVname, 'rain_num'))  || ...
          ~isempty(regexp(InVname, 'pris_num'))  || ...
          ~isempty(regexp(InVname, 'snow_num'))  || ...
          ~isempty(regexp(InVname, 'aggr_num'))  || ...
          ~isempty(regexp(InVname, 'graup_num')) || ...
          ~isempty(regexp(InVname, 'hail_num')))
        VAR      = VAR .* RhoAir .* 1e-6;
        DIFF_VAR = DIFF_VAR .* RhoAir .* 1e-6;
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
        switch(ProfDim)
          case 'r'
            Vsize = Nx;
            DimOrder = { 'x' };

          case 'z'
            Vsize = Nz;
            DimOrder = { 'z' };
        end
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
