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

    % Region: All
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_core_ps_d1_num'  '/all_core_ps_d1_num'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_core_s_d1_num'   '/all_core_s_d1_num'   'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_rb_ps_d1_num'    '/all_rb_ps_d1_num'    'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_rb_s_d1_num'     '/all_rb_s_d1_num'     'z' }

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_core_ps_d2_num'  '/all_core_ps_d2_num'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_core_s_d2_num'   '/all_core_s_d2_num'   'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_rb_ps_d2_num'    '/all_rb_ps_d2_num'    'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_rb_s_d2_num'     '/all_rb_s_d2_num'     'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_ps_dust_cloud'  '/all_core_ps_dust_cloud'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_s_dust_cloud'   '/all_core_s_dust_cloud'   'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_ps_dust_cloud'    '/all_rb_ps_dust_cloud'    'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_s_dust_cloud'     '/all_rb_s_dust_cloud'     'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_ps_dust_rain'  '/all_core_ps_dust_rain'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_s_dust_rain'   '/all_core_s_dust_rain'   'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_ps_dust_rain'    '/all_rb_ps_dust_rain'    'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_s_dust_rain'     '/all_rb_s_dust_rain'     'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_ps_dust_pris'  '/all_core_ps_dust_pris'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_s_dust_pris'   '/all_core_s_dust_pris'   'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_ps_dust_pris'    '/all_rb_ps_dust_pris'    'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_s_dust_pris'     '/all_rb_s_dust_pris'     'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_ps_dust_snow'  '/all_core_ps_dust_snow'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_s_dust_snow'   '/all_core_s_dust_snow'   'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_ps_dust_snow'    '/all_rb_ps_dust_snow'    'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_s_dust_snow'     '/all_rb_s_dust_snow'     'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_ps_dust_aggr'  '/all_core_ps_dust_aggr'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_s_dust_aggr'   '/all_core_s_dust_aggr'   'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_ps_dust_aggr'    '/all_rb_ps_dust_aggr'    'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_s_dust_aggr'     '/all_rb_s_dust_aggr'     'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_ps_dust_graup'  '/all_core_ps_dust_graup'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_s_dust_graup'   '/all_core_s_dust_graup'   'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_ps_dust_graup'    '/all_rb_ps_dust_graup'    'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_s_dust_graup'     '/all_rb_s_dust_graup'     'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_ps_dust_hail'  '/all_core_ps_dust_hail'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_s_dust_hail'   '/all_core_s_dust_hail'   'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_ps_dust_hail'    '/all_rb_ps_dust_hail'    'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_s_dust_hail'     '/all_rb_s_dust_hail'     'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_ps_dust_hydro'  '/all_core_ps_dust_hydro'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_core_s_dust_hydro'   '/all_core_s_dust_hydro'   'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_ps_dust_hydro'    '/all_rb_ps_dust_hydro'    'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_rb_s_dust_hydro'     '/all_rb_s_dust_hydro'     'z' }

    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_ps_dustifn_cloud'  '/all_core_ps_dustifn_cloud'  'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_s_dustifn_cloud'   '/all_core_s_dustifn_cloud'   'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_ps_dustifn_cloud'    '/all_rb_ps_dustifn_cloud'    'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_s_dustifn_cloud'     '/all_rb_s_dustifn_cloud'     'z' }

    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_ps_dustifn_rain'  '/all_core_ps_dustifn_rain'  'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_s_dustifn_rain'   '/all_core_s_dustifn_rain'   'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_ps_dustifn_rain'    '/all_rb_ps_dustifn_rain'    'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_s_dustifn_rain'     '/all_rb_s_dustifn_rain'     'z' }

    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_ps_dustifn_pris'  '/all_core_ps_dustifn_pris'  'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_s_dustifn_pris'   '/all_core_s_dustifn_pris'   'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_ps_dustifn_pris'    '/all_rb_ps_dustifn_pris'    'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_s_dustifn_pris'     '/all_rb_s_dustifn_pris'     'z' }

    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_ps_dustifn_snow'  '/all_core_ps_dustifn_snow'  'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_s_dustifn_snow'   '/all_core_s_dustifn_snow'   'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_ps_dustifn_snow'    '/all_rb_ps_dustifn_snow'    'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_s_dustifn_snow'     '/all_rb_s_dustifn_snow'     'z' }

    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_ps_dustifn_aggr'  '/all_core_ps_dustifn_aggr'  'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_s_dustifn_aggr'   '/all_core_s_dustifn_aggr'   'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_ps_dustifn_aggr'    '/all_rb_ps_dustifn_aggr'    'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_s_dustifn_aggr'     '/all_rb_s_dustifn_aggr'     'z' }

    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_ps_dustifn_graup'  '/all_core_ps_dustifn_graup'  'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_s_dustifn_graup'   '/all_core_s_dustifn_graup'   'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_ps_dustifn_graup'    '/all_rb_ps_dustifn_graup'    'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_s_dustifn_graup'     '/all_rb_s_dustifn_graup'     'z' }

    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_ps_dustifn_hail'  '/all_core_ps_dustifn_hail'  'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_core_s_dustifn_hail'   '/all_core_s_dustifn_hail'   'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_ps_dustifn_hail'    '/all_rb_ps_dustifn_hail'    'z' }
    { 'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5' '/all_rb_s_dustifn_hail'     '/all_rb_s_dustifn_hail'     'z' }

    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_core_ps_cloud_mass'  '/all_core_ps_cloud_mass'  'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_core_s_cloud_mass'   '/all_core_s_cloud_mass'   'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_rb_ps_cloud_mass'    '/all_rb_ps_cloud_mass'    'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_rb_s_cloud_mass'     '/all_rb_s_cloud_mass'     'z' }

    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_core_ps_cloud_num'  '/all_core_ps_cloud_num'  'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_core_s_cloud_num'   '/all_core_s_cloud_num'   'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_rb_ps_cloud_num'    '/all_rb_ps_cloud_num'    'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_rb_s_cloud_num'     '/all_rb_s_cloud_num'     'z' }

    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_core_ps_cloud_diam' '/all_core_ps_cloud_diam'  'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_core_s_cloud_diam'  '/all_core_s_cloud_diam'   'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_rb_ps_cloud_diam'   '/all_rb_ps_cloud_diam'    'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/all_rb_s_cloud_diam'    '/all_rb_s_cloud_diam'     'z' }

    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_core_ps_rain_mass'  '/all_core_ps_rain_mass'  'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_core_s_rain_mass'   '/all_core_s_rain_mass'   'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_rb_ps_rain_mass'    '/all_rb_ps_rain_mass'    'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_rb_s_rain_mass'     '/all_rb_s_rain_mass'     'z' }

    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_core_ps_rain_num'  '/all_core_ps_rain_num'  'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_core_s_rain_num'   '/all_core_s_rain_num'   'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_rb_ps_rain_num'    '/all_rb_ps_rain_num'    'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_rb_s_rain_num'     '/all_rb_s_rain_num'     'z' }

    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_core_ps_rain_diam' '/all_core_ps_rain_diam'  'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_core_s_rain_diam'  '/all_core_s_rain_diam'   'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_rb_ps_rain_diam'   '/all_rb_ps_rain_diam'    'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/all_rb_s_rain_diam'    '/all_rb_s_rain_diam'     'z' }

    { 'DIAGS/hist_meas_az_pris_<CASE>.h5' '/all_core_ps_pris_mass'  '/all_core_ps_pris_mass'  'z' }
    { 'DIAGS/hist_meas_az_pris_<CASE>.h5' '/all_core_s_pris_mass'   '/all_core_s_pris_mass'   'z' }
    { 'DIAGS/hist_meas_az_pris_<CASE>.h5' '/all_rb_ps_pris_mass'    '/all_rb_ps_pris_mass'    'z' }
    { 'DIAGS/hist_meas_az_pris_<CASE>.h5' '/all_rb_s_pris_mass'     '/all_rb_s_pris_mass'     'z' }

    { 'DIAGS/hist_meas_az_snow_<CASE>.h5' '/all_core_ps_snow_mass'  '/all_core_ps_snow_mass'  'z' }
    { 'DIAGS/hist_meas_az_snow_<CASE>.h5' '/all_core_s_snow_mass'   '/all_core_s_snow_mass'   'z' }
    { 'DIAGS/hist_meas_az_snow_<CASE>.h5' '/all_rb_ps_snow_mass'    '/all_rb_ps_snow_mass'    'z' }
    { 'DIAGS/hist_meas_az_snow_<CASE>.h5' '/all_rb_s_snow_mass'     '/all_rb_s_snow_mass'     'z' }

    { 'DIAGS/hist_meas_az_aggr_<CASE>.h5' '/all_core_ps_aggr_mass'  '/all_core_ps_aggr_mass'  'z' }
    { 'DIAGS/hist_meas_az_aggr_<CASE>.h5' '/all_core_s_aggr_mass'   '/all_core_s_aggr_mass'   'z' }
    { 'DIAGS/hist_meas_az_aggr_<CASE>.h5' '/all_rb_ps_aggr_mass'    '/all_rb_ps_aggr_mass'    'z' }
    { 'DIAGS/hist_meas_az_aggr_<CASE>.h5' '/all_rb_s_aggr_mass'     '/all_rb_s_aggr_mass'     'z' }

    { 'DIAGS/hist_meas_az_graup_<CASE>.h5' '/all_core_ps_graup_mass'  '/all_core_ps_graup_mass'  'z' }
    { 'DIAGS/hist_meas_az_graup_<CASE>.h5' '/all_core_s_graup_mass'   '/all_core_s_graup_mass'   'z' }
    { 'DIAGS/hist_meas_az_graup_<CASE>.h5' '/all_rb_ps_graup_mass'    '/all_rb_ps_graup_mass'    'z' }
    { 'DIAGS/hist_meas_az_graup_<CASE>.h5' '/all_rb_s_graup_mass'     '/all_rb_s_graup_mass'     'z' }

    { 'DIAGS/hist_meas_az_hail_<CASE>.h5' '/all_core_ps_hail_mass'  '/all_core_ps_hail_mass'  'z' }
    { 'DIAGS/hist_meas_az_hail_<CASE>.h5' '/all_core_s_hail_mass'   '/all_core_s_hail_mass'   'z' }
    { 'DIAGS/hist_meas_az_hail_<CASE>.h5' '/all_rb_ps_hail_mass'    '/all_rb_ps_hail_mass'    'z' }
    { 'DIAGS/hist_meas_az_hail_<CASE>.h5' '/all_rb_s_hail_mass'     '/all_rb_s_hail_mass'     'z' }

    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/all_core_ps_lhf_cool' '/all_core_ps_lhf_cool' 'z' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/all_core_s_lhf_cool'  '/all_core_s_lhf_cool'  'z' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/all_rb_ps_lhf_cool'   '/all_rb_ps_lhf_cool'   'z' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/all_rb_s_lhf_cool'    '/all_rb_s_lhf_cool'    'z' }

    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/all_core_ps_lhf_heat' '/all_core_ps_lhf_heat' 'z' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/all_core_s_lhf_heat'  '/all_core_s_lhf_heat'  'z' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/all_rb_ps_lhf_heat'   '/all_rb_ps_lhf_heat'   'z' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/all_rb_s_lhf_heat'    '/all_rb_s_lhf_heat'    'z' }

    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/all_core_ps_lhv_cool' '/all_core_ps_lhv_cool' 'z' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/all_core_s_lhv_cool'  '/all_core_s_lhv_cool'  'z' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/all_rb_ps_lhv_cool'   '/all_rb_ps_lhv_cool'   'z' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/all_rb_s_lhv_cool'    '/all_rb_s_lhv_cool'    'z' }

    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/all_core_ps_lhv_heat' '/all_core_ps_lhv_heat' 'z' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/all_core_s_lhv_heat'  '/all_core_s_lhv_heat'  'z' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/all_rb_ps_lhv_heat'   '/all_rb_ps_lhv_heat'   'z' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/all_rb_s_lhv_heat'    '/all_rb_s_lhv_heat'    'z' }

    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/all_core_ps_liq_evap' '/all_core_ps_liq_evap' 'z' }
    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/all_core_s_liq_evap'  '/all_core_s_liq_evap'  'z' }
    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/all_rb_ps_liq_evap'   '/all_rb_ps_liq_evap'   'z' }
    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/all_rb_s_liq_evap'    '/all_rb_s_liq_evap'    'z' }

    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/all_core_ps_liq_cond' '/all_core_ps_liq_cond' 'z' }
    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/all_core_s_liq_cond'  '/all_core_s_liq_cond'  'z' }
    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/all_rb_ps_liq_cond'   '/all_rb_ps_liq_cond'   'z' }
    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/all_rb_s_liq_cond'    '/all_rb_s_liq_cond'    'z' }

    { 'DIAGS/hist_meas_az_cloud_evap_<CASE>.h5' '/all_core_ps_cloud_evap' '/all_core_ps_cloud_evap' 'z' }
    { 'DIAGS/hist_meas_az_cloud_evap_<CASE>.h5' '/all_core_s_cloud_evap'  '/all_core_s_cloud_evap'  'z' }
    { 'DIAGS/hist_meas_az_cloud_evap_<CASE>.h5' '/all_rb_ps_cloud_evap'   '/all_rb_ps_cloud_evap'   'z' }
    { 'DIAGS/hist_meas_az_cloud_evap_<CASE>.h5' '/all_rb_s_cloud_evap'    '/all_rb_s_cloud_evap'    'z' }

    { 'DIAGS/hist_meas_az_cloud_cond_<CASE>.h5' '/all_core_ps_cloud_cond' '/all_core_ps_cloud_cond' 'z' }
    { 'DIAGS/hist_meas_az_cloud_cond_<CASE>.h5' '/all_core_s_cloud_cond'  '/all_core_s_cloud_cond'  'z' }
    { 'DIAGS/hist_meas_az_cloud_cond_<CASE>.h5' '/all_rb_ps_cloud_cond'   '/all_rb_ps_cloud_cond'   'z' }
    { 'DIAGS/hist_meas_az_cloud_cond_<CASE>.h5' '/all_rb_s_cloud_cond'    '/all_rb_s_cloud_cond'    'z' }

    { 'DIAGS/hist_meas_az_rain_evap_<CASE>.h5' '/all_core_ps_rain_evap' '/all_core_ps_rain_evap' 'z' }
    { 'DIAGS/hist_meas_az_rain_evap_<CASE>.h5' '/all_core_s_rain_evap'  '/all_core_s_rain_evap'  'z' }
    { 'DIAGS/hist_meas_az_rain_evap_<CASE>.h5' '/all_rb_ps_rain_evap'   '/all_rb_ps_rain_evap'   'z' }
    { 'DIAGS/hist_meas_az_rain_evap_<CASE>.h5' '/all_rb_s_rain_evap'    '/all_rb_s_rain_evap'    'z' }

    { 'DIAGS/hist_meas_az_rain_cond_<CASE>.h5' '/all_core_ps_rain_cond' '/all_core_ps_rain_cond' 'z' }
    { 'DIAGS/hist_meas_az_rain_cond_<CASE>.h5' '/all_core_s_rain_cond'  '/all_core_s_rain_cond'  'z' }
    { 'DIAGS/hist_meas_az_rain_cond_<CASE>.h5' '/all_rb_ps_rain_cond'   '/all_rb_ps_rain_cond'   'z' }
    { 'DIAGS/hist_meas_az_rain_cond_<CASE>.h5' '/all_rb_s_rain_cond'    '/all_rb_s_rain_cond'    'z' }

    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_core_ps_updraft' '/all_core_ps_updraft' 'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_core_s_updraft'  '/all_core_s_updraft'  'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_rb_ps_updraft'   '/all_rb_ps_updraft'   'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_rb_s_updraft'    '/all_rb_s_updraft'    'z' }

    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_core_ps_dndraft' '/all_core_ps_dndraft'  'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_core_s_dndraft'  '/all_core_s_dndraft'  'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_rb_ps_dndraft'   '/all_rb_ps_dndraft'   'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/all_rb_s_dndraft'    '/all_rb_s_dndraft'    'z' }

    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/all_core_ps_theta_e' '/all_core_ps_theta_e' 'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/all_core_s_theta_e'  '/all_core_s_theta_e'  'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/all_rb_ps_theta_e'   '/all_rb_ps_theta_e'   'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/all_rb_s_theta_e'    '/all_rb_s_theta_e'    'z' }

    { 'DIAGS/hist_meas_az_theta_rho_<CASE>.h5' '/all_core_ps_theta_rho' '/all_core_ps_theta_rho' 'z' }
    { 'DIAGS/hist_meas_az_theta_rho_<CASE>.h5' '/all_core_s_theta_rho'  '/all_core_s_theta_rho'  'z' }
    { 'DIAGS/hist_meas_az_theta_rho_<CASE>.h5' '/all_rb_ps_theta_rho'   '/all_rb_ps_theta_rho'   'z' }
    { 'DIAGS/hist_meas_az_theta_rho_<CASE>.h5' '/all_rb_s_theta_rho'    '/all_rb_s_theta_rho'    'z' }

    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/all_core_ps_theta' '/all_core_ps_theta' 'z' }
    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/all_core_s_theta'  '/all_core_s_theta'  'z' }
    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/all_rb_ps_theta'   '/all_rb_ps_theta'   'z' }
    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/all_rb_s_theta'    '/all_rb_s_theta'    'z' }

    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/all_core_ps_tempc' '/all_core_ps_tempc' 'z' }
    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/all_core_s_tempc'  '/all_core_s_tempc'  'z' }
    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/all_rb_ps_tempc'   '/all_rb_ps_tempc'   'z' }
    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/all_rb_s_tempc'    '/all_rb_s_tempc'    'z' }

    % Region: Lead
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/lead_core_ps_d1_num'  '/lead_core_ps_d1_num'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/lead_core_s_d1_num'   '/lead_core_s_d1_num'   'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/lead_rb_ps_d1_num'    '/lead_rb_ps_d1_num'    'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/lead_rb_s_d1_num'     '/lead_rb_s_d1_num'     'z' }

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/lead_core_ps_d2_num'  '/lead_core_ps_d2_num'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/lead_core_s_d2_num'   '/lead_core_s_d2_num'   'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/lead_rb_ps_d2_num'    '/lead_rb_ps_d2_num'    'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/lead_rb_s_d2_num'     '/lead_rb_s_d2_num'     'z' }

    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_core_ps_cloud_mass'  '/lead_core_ps_cloud_mass'  'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_core_s_cloud_mass'   '/lead_core_s_cloud_mass'   'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_rb_ps_cloud_mass'    '/lead_rb_ps_cloud_mass'    'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_rb_s_cloud_mass'     '/lead_rb_s_cloud_mass'     'z' }

    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_core_ps_cloud_num'  '/lead_core_ps_cloud_num'  'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_core_s_cloud_num'   '/lead_core_s_cloud_num'   'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_rb_ps_cloud_num'    '/lead_rb_ps_cloud_num'    'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_rb_s_cloud_num'     '/lead_rb_s_cloud_num'     'z' }

    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_core_ps_cloud_diam' '/lead_core_ps_cloud_diam'  'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_core_s_cloud_diam'  '/lead_core_s_cloud_diam'   'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_rb_ps_cloud_diam'   '/lead_rb_ps_cloud_diam'    'z' }
    { 'DIAGS/hist_meas_az_cloud_<CASE>.h5' '/lead_rb_s_cloud_diam'    '/lead_rb_s_cloud_diam'     'z' }

    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_core_ps_rain_mass'  '/lead_core_ps_rain_mass'  'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_core_s_rain_mass'   '/lead_core_s_rain_mass'   'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_rb_ps_rain_mass'    '/lead_rb_ps_rain_mass'    'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_rb_s_rain_mass'     '/lead_rb_s_rain_mass'     'z' }

    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_core_ps_rain_num'  '/lead_core_ps_rain_num'  'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_core_s_rain_num'   '/lead_core_s_rain_num'   'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_rb_ps_rain_num'    '/lead_rb_ps_rain_num'    'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_rb_s_rain_num'     '/lead_rb_s_rain_num'     'z' }

    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_core_ps_rain_diam' '/lead_core_ps_rain_diam'  'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_core_s_rain_diam'  '/lead_core_s_rain_diam'   'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_rb_ps_rain_diam'   '/lead_rb_ps_rain_diam'    'z' }
    { 'DIAGS/hist_meas_az_rain_<CASE>.h5' '/lead_rb_s_rain_diam'    '/lead_rb_s_rain_diam'     'z' }

    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/lead_core_ps_lhf_cool' '/lead_core_ps_lhf_cool' 'z' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/lead_core_s_lhf_cool'  '/lead_core_s_lhf_cool'  'z' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/lead_rb_ps_lhf_cool'   '/lead_rb_ps_lhf_cool'   'z' }
    { 'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5' '/lead_rb_s_lhf_cool'    '/lead_rb_s_lhf_cool'    'z' }

    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/lead_core_ps_lhf_heat' '/lead_core_ps_lhf_heat' 'z' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/lead_core_s_lhf_heat'  '/lead_core_s_lhf_heat'  'z' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/lead_rb_ps_lhf_heat'   '/lead_rb_ps_lhf_heat'   'z' }
    { 'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5' '/lead_rb_s_lhf_heat'    '/lead_rb_s_lhf_heat'    'z' }

    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/lead_core_ps_lhv_cool' '/lead_core_ps_lhv_cool' 'z' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/lead_core_s_lhv_cool'  '/lead_core_s_lhv_cool'  'z' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/lead_rb_ps_lhv_cool'   '/lead_rb_ps_lhv_cool'   'z' }
    { 'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5' '/lead_rb_s_lhv_cool'    '/lead_rb_s_lhv_cool'    'z' }

    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/lead_core_ps_lhv_heat' '/lead_core_ps_lhv_heat' 'z' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/lead_core_s_lhv_heat'  '/lead_core_s_lhv_heat'  'z' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/lead_rb_ps_lhv_heat'   '/lead_rb_ps_lhv_heat'   'z' }
    { 'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5' '/lead_rb_s_lhv_heat'    '/lead_rb_s_lhv_heat'    'z' }

    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/lead_core_ps_liq_evap' '/lead_core_ps_liq_evap' 'z' }
    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/lead_core_s_liq_evap'  '/lead_core_s_liq_evap'  'z' }
    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/lead_rb_ps_liq_evap'   '/lead_rb_ps_liq_evap'   'z' }
    { 'DIAGS/hist_meas_az_liq_evap_<CASE>.h5' '/lead_rb_s_liq_evap'    '/lead_rb_s_liq_evap'    'z' }

    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/lead_core_ps_liq_cond' '/lead_core_ps_liq_cond' 'z' }
    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/lead_core_s_liq_cond'  '/lead_core_s_liq_cond'  'z' }
    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/lead_rb_ps_liq_cond'   '/lead_rb_ps_liq_cond'   'z' }
    { 'DIAGS/hist_meas_az_liq_cond_<CASE>.h5' '/lead_rb_s_liq_cond'    '/lead_rb_s_liq_cond'    'z' }

    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/lead_core_ps_updraft' '/lead_core_ps_updraft' 'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/lead_core_s_updraft'  '/lead_core_s_updraft'  'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/lead_rb_ps_updraft'   '/lead_rb_ps_updraft'   'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/lead_rb_s_updraft'    '/lead_rb_s_updraft'    'z' }

    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/lead_core_ps_dndraft' '/lead_core_ps_dndraft'  'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/lead_core_s_dndraft'  '/lead_core_s_dndraft'  'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/lead_rb_ps_dndraft'   '/lead_rb_ps_dndraft'   'z' }
    { 'DIAGS/hist_meas_az_w_<CASE>.h5' '/lead_rb_s_dndraft'    '/lead_rb_s_dndraft'    'z' }

    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/lead_core_ps_theta_e' '/lead_core_ps_theta_e' 'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/lead_core_s_theta_e'  '/lead_core_s_theta_e'  'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/lead_rb_ps_theta_e'   '/lead_rb_ps_theta_e'   'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/lead_rb_s_theta_e'    '/lead_rb_s_theta_e'    'z' }

    { 'DIAGS/hist_meas_az_theta_rho_<CASE>.h5' '/lead_core_ps_theta_rho' '/lead_core_ps_theta_rho' 'z' }
    { 'DIAGS/hist_meas_az_theta_rho_<CASE>.h5' '/lead_core_s_theta_rho'  '/lead_core_s_theta_rho'  'z' }
    { 'DIAGS/hist_meas_az_theta_rho_<CASE>.h5' '/lead_rb_ps_theta_rho'   '/lead_rb_ps_theta_rho'   'z' }
    { 'DIAGS/hist_meas_az_theta_rho_<CASE>.h5' '/lead_rb_s_theta_rho'    '/lead_rb_s_theta_rho'    'z' }

    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/lead_core_ps_theta' '/lead_core_ps_theta' 'z' }
    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/lead_core_s_theta'  '/lead_core_s_theta'  'z' }
    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/lead_rb_ps_theta'   '/lead_rb_ps_theta'   'z' }
    { 'DIAGS/hist_meas_az_theta_<CASE>.h5' '/lead_rb_s_theta'    '/lead_rb_s_theta'    'z' }

    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/lead_core_ps_tempc' '/lead_core_ps_tempc' 'z' }
    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/lead_core_s_tempc'  '/lead_core_s_tempc'  'z' }
    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/lead_rb_ps_tempc'   '/lead_rb_ps_tempc'   'z' }
    { 'DIAGS/hist_meas_az_tempc_<CASE>.h5' '/lead_rb_s_tempc'    '/lead_rb_s_tempc'    'z' }

    % Region: spath - area inside SAL region, that is in the path of the storm
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d1_mass'  '/spath_ps_d1_mass' 'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d1_mass'   '/spath_s_d1_mass'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d1_num'   '/spath_ps_d1_num'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d1_num'    '/spath_s_d1_num'   'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d2_mass'  '/spath_ps_d2_mass' 'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d2_mass'   '/spath_s_d2_mass'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d2_num'   '/spath_ps_d2_num'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d2_num'    '/spath_s_d2_num'   'z' }

    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_tracer1'  '/spath_ps_tracer1' 'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_tracer1'   '/spath_s_tracer1'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_tracer2'  '/spath_ps_tracer2' 'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_tracer2'   '/spath_s_tracer2'  'z' }

    % Region: All, selected by precip rate > 0.01 mm/h
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/precip_all_core_ps_theta_e' '/precip_all_core_ps_theta_e' 'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/precip_all_core_s_theta_e'  '/precip_all_core_s_theta_e'  'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/precip_all_rb_ps_theta_e'   '/precip_all_rb_ps_theta_e'   'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/precip_all_rb_s_theta_e'    '/precip_all_rb_s_theta_e'    'z' }

    % Region: Lead, selected by precip rate > 0.01 mm/h
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/precip_lead_core_ps_theta_e' '/precip_lead_core_ps_theta_e' 'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/precip_lead_core_s_theta_e'  '/precip_lead_core_s_theta_e'  'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/precip_lead_rb_ps_theta_e'   '/precip_lead_rb_ps_theta_e'   'z' }
    { 'DIAGS/hist_meas_az_theta_e_<CASE>.h5' '/precip_lead_rb_s_theta_e'    '/precip_lead_rb_s_theta_e'    'z' }




%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/core_ps_speed_t' '/pcore_s_speed_t'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/rb_ps_speed_t' '/rb_ps_speed_t'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/core_s_speed_t' '/core_s_speed_t'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/rb_s_speed_t' '/rb_s_speed_t'  'z' }
%
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/core_ps_speed_r' '/pcore_s_speed_r'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/rb_ps_speed_r' '/rb_ps_speed_r'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/core_s_speed_r' '/core_s_speed_r'  'z' }
%    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/rb_s_speed_r' '/rb_s_speed_r'  'z' }
%
%    { 'DIAGS/hist_meas_az_pcprate_<CASE>.h5' '/ps_pcprate' '/ps_pcprate'  'r' }
%    { 'DIAGS/hist_meas_az_pcprate_<CASE>.h5' '/s_pcprate' '/s_pcprate'  'r' }
%
%    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/ps_vint_cond' '/ps_vint_cond'  'r' }
%    { 'DIAGS/hist_meas_vint_cond_<CASE>.h5' '/s_vint_cond' '/s_vint_cond'  'r' }
%
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/core_ps_vapor'      '/pcore_s_vapor'  'z' }
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/rb_ps_vapor'        '/rb_ps_vapor'  'z' }
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/core_s_vapor'       '/core_s_vapor'  'z' }
%    { 'DIAGS/hist_meas_az_vapor_<CASE>.h5' '/rb_s_vapor'         '/rb_s_vapor'  'z' }
%
%    { 'DIAGS/hist_meas_ts_vapor_<CASE>.h5' '/spath_ps_vapor'  '/spath_ps_vapor'  'z' }
%    { 'DIAGS/hist_meas_ts_vapor_<CASE>.h5' '/smaxcp_ps_vapor' '/smaxcp_ps_vapor' 'z' }
%    { 'DIAGS/hist_meas_ts_vapor_<CASE>.h5' '/spath_s_vapor'   '/spath_s_vapor'   'z' }
%    { 'DIAGS/hist_meas_ts_vapor_<CASE>.h5' '/smaxcp_s_vapor'  '/smaxcp_s_vapor'  'z' }
%
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_ps_d1_mass' '/smaxcp_rb_ps_d1_mass'  'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_s_d1_mass'  '/smaxcp_rb_s_d1_mass'  'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_ps_d1_num'  '/smaxcp_rb_ps_d1_num'   'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_s_d1_num'   '/smaxcp_rb_s_d1_num'   'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_ps_d2_mass' '/smaxcp_rb_ps_d2_mass'  'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_s_d2_mass'  '/smaxcp_rb_s_d2_mass'  'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_ps_d2_num'  '/smaxcp_rb_ps_d2_num'   'z' }
%    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/smaxcp_s_d2_num'   '/smaxcp_rb_s_d2_num'   'z' }
%
%    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/core_ps_tcond'  '/pcore_s_tcond'   'z' }
%    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/rb_ps_tcond'  '/rb_ps_tcond'   'z' }
%    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/core_s_tcond'  '/core_s_tcond'   'z' }
%    { 'DIAGS/hist_meas_tcond_<CASE>.h5' '/rb_s_tcond'  '/rb_s_tcond'   'z' }
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
      if ((~isempty(regexp(Case, 'NODUST'))) && ...
          ((~isempty(regexp(InVname, 'd[12]_num'))) || ...
           (~isempty(regexp(InVname, 'd[12]_mass'))) || ...
           (~isempty(regexp(InVname, '_dust_'))) || ...
           (~isempty(regexp(InVname, '_dustifn_'))) || ...
           (~isempty(regexp(InVname, '_tracer[12]')))))
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
