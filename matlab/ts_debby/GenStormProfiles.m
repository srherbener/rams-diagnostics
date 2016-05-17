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
    { 'DIAGS/hist_meas_az_ccn_<CASE>.h5' '/all_whole_ps_ccn_mass' '/all_whole_ps_ccn_mass' 'z' }
    { 'DIAGS/hist_meas_az_ccn_<CASE>.h5' '/all_whole_s_ccn_mass'  '/all_whole_s_ccn_mass'  'z' }
    { 'DIAGS/hist_meas_az_ccn_<CASE>.h5' '/all_whole_i_ccn_mass'  '/all_whole_i_ccn_mass'  'z' }
    { 'DIAGS/hist_meas_az_ccn_<CASE>.h5' '/all_whole_m_ccn_mass'  '/all_whole_m_ccn_mass'  'z' }
    { 'DIAGS/hist_meas_az_ccn_<CASE>.h5' '/all_whole_f_ccn_mass'  '/all_whole_f_ccn_mass'  'z' }

    { 'DIAGS/hist_meas_az_ra_<CASE>.h5' '/all_whole_ps_ra_mass' '/all_whole_ps_ra_mass' 'z' }
    { 'DIAGS/hist_meas_az_ra_<CASE>.h5' '/all_whole_s_ra_mass'  '/all_whole_s_ra_mass'  'z' }
    { 'DIAGS/hist_meas_az_ra_<CASE>.h5' '/all_whole_i_ra_mass'  '/all_whole_i_ra_mass'  'z' }
    { 'DIAGS/hist_meas_az_ra_<CASE>.h5' '/all_whole_m_ra_mass'  '/all_whole_m_ra_mass'  'z' }
    { 'DIAGS/hist_meas_az_ra_<CASE>.h5' '/all_whole_f_ra_mass'  '/all_whole_f_ra_mass'  'z' }

    { 'DIAGS/hist_meas_az_aero_<CASE>.h5' '/all_whole_ps_aero_mass' '/all_whole_ps_aero_mass' 'z' }
    { 'DIAGS/hist_meas_az_aero_<CASE>.h5' '/all_whole_s_aero_mass'  '/all_whole_s_aero_mass'  'z' }
    { 'DIAGS/hist_meas_az_aero_<CASE>.h5' '/all_whole_i_aero_mass'  '/all_whole_i_aero_mass'  'z' }
    { 'DIAGS/hist_meas_az_aero_<CASE>.h5' '/all_whole_m_aero_mass'  '/all_whole_m_aero_mass'  'z' }
    { 'DIAGS/hist_meas_az_aero_<CASE>.h5' '/all_whole_f_aero_mass'  '/all_whole_f_aero_mass'  'z' }

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_ps_d1_num' '/all_whole_ps_d1_num' 'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_s_d1_num'  '/all_whole_s_d1_num'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_core_ps_d1_num'  '/all_core_ps_d1_num'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_core_s_d1_num'   '/all_core_s_d1_num'   'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_rb_ps_d1_num'    '/all_rb_ps_d1_num'    'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_rb_s_d1_num'     '/all_rb_s_d1_num'     'z' }

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_ps_d2_num' '/all_whole_ps_d2_num' 'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_s_d2_num'  '/all_whole_s_d2_num'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_core_ps_d2_num'  '/all_core_ps_d2_num'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_core_s_d2_num'   '/all_core_s_d2_num'   'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_rb_ps_d2_num'    '/all_rb_ps_d2_num'    'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_rb_s_d2_num'     '/all_rb_s_d2_num'     'z' }

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_ps_d1_mass' '/all_whole_ps_d1_mass' 'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_s_d1_mass'  '/all_whole_s_d1_mass'  'z' }

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_ps_d2_mass' '/all_whole_ps_d2_mass' 'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_s_d2_mass'  '/all_whole_s_d2_mass'  'z' }

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_ps_trdust1_diff' '/all_whole_ps_trdust1_diff' 'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_s_trdust1_diff'  '/all_whole_s_trdust1_diff'  'z' }

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_ps_trdust2_diff' '/all_whole_ps_trdust2_diff' 'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_s_trdust2_diff'  '/all_whole_s_trdust2_diff'  'z' }

    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_ps_dust_mass' '/all_whole_ps_dust_mass' 'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_s_dust_mass'  '/all_whole_s_dust_mass'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_i_dust_mass'  '/all_whole_i_dust_mass'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_m_dust_mass'  '/all_whole_m_dust_mass'  'z' }
    { 'DIAGS/hist_meas_az_dust_<CASE>.h5' '/all_whole_f_dust_mass'  '/all_whole_f_dust_mass'  'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_ps_dust_cloud' '/all_whole_ps_dust_cloud' 'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_s_dust_cloud'  '/all_whole_s_dust_cloud'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_i_dust_cloud'  '/all_whole_i_dust_cloud'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_m_dust_cloud'  '/all_whole_m_dust_cloud'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_f_dust_cloud'  '/all_whole_f_dust_cloud'  'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_ps_dust_rain' '/all_whole_ps_dust_rain' 'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_s_dust_rain'  '/all_whole_s_dust_rain'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_i_dust_rain'  '/all_whole_i_dust_rain'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_m_dust_rain'  '/all_whole_m_dust_rain'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_f_dust_rain'  '/all_whole_f_dust_rain'  'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_ps_dust_pris' '/all_whole_ps_dust_pris' 'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_s_dust_pris'  '/all_whole_s_dust_pris'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_i_dust_pris'  '/all_whole_i_dust_pris'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_m_dust_pris'  '/all_whole_m_dust_pris'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_f_dust_pris'  '/all_whole_f_dust_pris'  'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_ps_dust_snow' '/all_whole_ps_dust_snow' 'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_s_dust_snow'  '/all_whole_s_dust_snow'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_i_dust_snow'  '/all_whole_i_dust_snow'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_m_dust_snow'  '/all_whole_m_dust_snow'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_f_dust_snow'  '/all_whole_f_dust_snow'  'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_ps_dust_aggr' '/all_whole_ps_dust_aggr' 'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_s_dust_aggr'  '/all_whole_s_dust_aggr'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_i_dust_aggr'  '/all_whole_i_dust_aggr'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_m_dust_aggr'  '/all_whole_m_dust_aggr'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_f_dust_aggr'  '/all_whole_f_dust_aggr'  'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_ps_dust_graup' '/all_whole_ps_dust_graup' 'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_s_dust_graup'  '/all_whole_s_dust_graup'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_i_dust_graup'  '/all_whole_i_dust_graup'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_m_dust_graup'  '/all_whole_m_dust_graup'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_f_dust_graup'  '/all_whole_f_dust_graup'  'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_ps_dust_hail' '/all_whole_ps_dust_hail' 'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_s_dust_hail'  '/all_whole_s_dust_hail'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_i_dust_hail'  '/all_whole_i_dust_hail'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_m_dust_hail'  '/all_whole_m_dust_hail'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_f_dust_hail'  '/all_whole_f_dust_hail'  'z' }

    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_ps_dust_hydro' '/all_whole_ps_dust_hydro' 'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_s_dust_hydro'  '/all_whole_s_dust_hydro'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_i_dust_hydro'  '/all_whole_i_dust_hydro'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_m_dust_hydro'  '/all_whole_m_dust_hydro'  'z' }
    { 'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5' '/all_whole_f_dust_hydro'  '/all_whole_f_dust_hydro'  'z' }

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

    { 'DIAGS/hist_meas_az_pris_sub_<CASE>.h5' '/all_core_ps_pris_sub' '/all_core_ps_pris_sub' 'z' }
    { 'DIAGS/hist_meas_az_pris_sub_<CASE>.h5' '/all_core_s_pris_sub'  '/all_core_s_pris_sub'  'z' }
    { 'DIAGS/hist_meas_az_pris_sub_<CASE>.h5' '/all_rb_ps_pris_sub'   '/all_rb_ps_pris_sub'   'z' }
    { 'DIAGS/hist_meas_az_pris_sub_<CASE>.h5' '/all_rb_s_pris_sub'    '/all_rb_s_pris_sub'    'z' }

    { 'DIAGS/hist_meas_az_pris_dep_<CASE>.h5' '/all_core_ps_pris_dep' '/all_core_ps_pris_dep' 'z' }
    { 'DIAGS/hist_meas_az_pris_dep_<CASE>.h5' '/all_core_s_pris_dep'  '/all_core_s_pris_dep'  'z' }
    { 'DIAGS/hist_meas_az_pris_dep_<CASE>.h5' '/all_rb_ps_pris_dep'   '/all_rb_ps_pris_dep'   'z' }
    { 'DIAGS/hist_meas_az_pris_dep_<CASE>.h5' '/all_rb_s_pris_dep'    '/all_rb_s_pris_dep'    'z' }

    { 'DIAGS/hist_meas_az_snow_sub_<CASE>.h5' '/all_core_ps_snow_sub' '/all_core_ps_snow_sub' 'z' }
    { 'DIAGS/hist_meas_az_snow_sub_<CASE>.h5' '/all_core_s_snow_sub'  '/all_core_s_snow_sub'  'z' }
    { 'DIAGS/hist_meas_az_snow_sub_<CASE>.h5' '/all_rb_ps_snow_sub'   '/all_rb_ps_snow_sub'   'z' }
    { 'DIAGS/hist_meas_az_snow_sub_<CASE>.h5' '/all_rb_s_snow_sub'    '/all_rb_s_snow_sub'    'z' }

    { 'DIAGS/hist_meas_az_snow_dep_<CASE>.h5' '/all_core_ps_snow_dep' '/all_core_ps_snow_dep' 'z' }
    { 'DIAGS/hist_meas_az_snow_dep_<CASE>.h5' '/all_core_s_snow_dep'  '/all_core_s_snow_dep'  'z' }
    { 'DIAGS/hist_meas_az_snow_dep_<CASE>.h5' '/all_rb_ps_snow_dep'   '/all_rb_ps_snow_dep'   'z' }
    { 'DIAGS/hist_meas_az_snow_dep_<CASE>.h5' '/all_rb_s_snow_dep'    '/all_rb_s_snow_dep'    'z' }

    { 'DIAGS/hist_meas_az_aggr_sub_<CASE>.h5' '/all_core_ps_aggr_sub' '/all_core_ps_aggr_sub' 'z' }
    { 'DIAGS/hist_meas_az_aggr_sub_<CASE>.h5' '/all_core_s_aggr_sub'  '/all_core_s_aggr_sub'  'z' }
    { 'DIAGS/hist_meas_az_aggr_sub_<CASE>.h5' '/all_rb_ps_aggr_sub'   '/all_rb_ps_aggr_sub'   'z' }
    { 'DIAGS/hist_meas_az_aggr_sub_<CASE>.h5' '/all_rb_s_aggr_sub'    '/all_rb_s_aggr_sub'    'z' }

    { 'DIAGS/hist_meas_az_aggr_dep_<CASE>.h5' '/all_core_ps_aggr_dep' '/all_core_ps_aggr_dep' 'z' }
    { 'DIAGS/hist_meas_az_aggr_dep_<CASE>.h5' '/all_core_s_aggr_dep'  '/all_core_s_aggr_dep'  'z' }
    { 'DIAGS/hist_meas_az_aggr_dep_<CASE>.h5' '/all_rb_ps_aggr_dep'   '/all_rb_ps_aggr_dep'   'z' }
    { 'DIAGS/hist_meas_az_aggr_dep_<CASE>.h5' '/all_rb_s_aggr_dep'    '/all_rb_s_aggr_dep'    'z' }

    { 'DIAGS/hist_meas_az_graup_sub_<CASE>.h5' '/all_core_ps_graup_sub' '/all_core_ps_graup_sub' 'z' }
    { 'DIAGS/hist_meas_az_graup_sub_<CASE>.h5' '/all_core_s_graup_sub'  '/all_core_s_graup_sub'  'z' }
    { 'DIAGS/hist_meas_az_graup_sub_<CASE>.h5' '/all_rb_ps_graup_sub'   '/all_rb_ps_graup_sub'   'z' }
    { 'DIAGS/hist_meas_az_graup_sub_<CASE>.h5' '/all_rb_s_graup_sub'    '/all_rb_s_graup_sub'    'z' }

    { 'DIAGS/hist_meas_az_graup_dep_<CASE>.h5' '/all_core_ps_graup_dep' '/all_core_ps_graup_dep' 'z' }
    { 'DIAGS/hist_meas_az_graup_dep_<CASE>.h5' '/all_core_s_graup_dep'  '/all_core_s_graup_dep'  'z' }
    { 'DIAGS/hist_meas_az_graup_dep_<CASE>.h5' '/all_rb_ps_graup_dep'   '/all_rb_ps_graup_dep'   'z' }
    { 'DIAGS/hist_meas_az_graup_dep_<CASE>.h5' '/all_rb_s_graup_dep'    '/all_rb_s_graup_dep'    'z' }

    { 'DIAGS/hist_meas_az_hail_sub_<CASE>.h5' '/all_core_ps_hail_sub' '/all_core_ps_hail_sub' 'z' }
    { 'DIAGS/hist_meas_az_hail_sub_<CASE>.h5' '/all_core_s_hail_sub'  '/all_core_s_hail_sub'  'z' }
    { 'DIAGS/hist_meas_az_hail_sub_<CASE>.h5' '/all_rb_ps_hail_sub'   '/all_rb_ps_hail_sub'   'z' }
    { 'DIAGS/hist_meas_az_hail_sub_<CASE>.h5' '/all_rb_s_hail_sub'    '/all_rb_s_hail_sub'    'z' }

    { 'DIAGS/hist_meas_az_hail_dep_<CASE>.h5' '/all_core_ps_hail_dep' '/all_core_ps_hail_dep' 'z' }
    { 'DIAGS/hist_meas_az_hail_dep_<CASE>.h5' '/all_core_s_hail_dep'  '/all_core_s_hail_dep'  'z' }
    { 'DIAGS/hist_meas_az_hail_dep_<CASE>.h5' '/all_rb_ps_hail_dep'   '/all_rb_ps_hail_dep'   'z' }
    { 'DIAGS/hist_meas_az_hail_dep_<CASE>.h5' '/all_rb_s_hail_dep'    '/all_rb_s_hail_dep'    'z' }

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

    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_core_ps_buoy_pos_acc' '/all_core_ps_buoy_pos_acc' 'z' }
    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_core_s_buoy_pos_acc'  '/all_core_s_buoy_pos_acc'  'z' }
    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_rb_ps_buoy_pos_acc'   '/all_rb_ps_buoy_pos_acc'   'z' }
    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_rb_s_buoy_pos_acc'    '/all_rb_s_buoy_pos_acc'    'z' }
    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_core_ps_buoy_neg_acc' '/all_core_ps_buoy_neg_acc' 'z' }
    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_core_s_buoy_neg_acc'  '/all_core_s_buoy_neg_acc'  'z' }
    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_rb_ps_buoy_neg_acc'   '/all_rb_ps_buoy_neg_acc'   'z' }
    { 'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5' '/all_rb_s_buoy_neg_acc'    '/all_rb_s_buoy_neg_acc'    'z' }

    % Region: spath - area inside SAL region, that is in the path of the storm
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_i_d1_mass'   '/spath_i_d1_mass'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d1_mass'  '/spath_ps_d1_mass' 'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d1_mass'   '/spath_s_d1_mass'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_i_d1_num'    '/spath_i_d1_num'   'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d1_num'   '/spath_ps_d1_num'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d1_num'    '/spath_s_d1_num'   'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_i_d2_mass'   '/spath_i_d2_mass'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d2_mass'  '/spath_ps_d2_mass' 'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d2_mass'   '/spath_s_d2_mass'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_i_d2_num'    '/spath_i_d2_num'   'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_d2_num'   '/spath_ps_d2_num'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_d2_num'    '/spath_s_d2_num'   'z' }

    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_i_tracer1'   '/spath_i_tracer1'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_tracer1'  '/spath_ps_tracer1' 'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_tracer1'   '/spath_s_tracer1'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_i_tracer2'   '/spath_i_tracer2'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_ps_tracer2'  '/spath_ps_tracer2' 'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_s_tracer2'   '/spath_s_tracer2'  'z' }

    { 'DIAGS/hist_meas_ts_ccn_<CASE>.h5' '/spath_i_ccn_mass'   '/spath_i_ccn_mass'  'z' }
    { 'DIAGS/hist_meas_ts_ccn_<CASE>.h5' '/spath_m_ccn_mass'   '/spath_m_ccn_mass'  'z' }
    { 'DIAGS/hist_meas_ts_ccn_<CASE>.h5' '/spath_f_ccn_mass'   '/spath_f_ccn_mass'  'z' }

    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_i_dust_mass'   '/spath_i_dust_mass'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_m_dust_mass'   '/spath_m_dust_mass'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_f_dust_mass'   '/spath_f_dust_mass'  'z' }

    { 'DIAGS/hist_meas_ts_ra_<CASE>.h5' '/spath_i_ra_mass'   '/spath_i_ra_mass'  'z' }
    { 'DIAGS/hist_meas_ts_ra_<CASE>.h5' '/spath_m_ra_mass'   '/spath_m_ra_mass'  'z' }
    { 'DIAGS/hist_meas_ts_ra_<CASE>.h5' '/spath_f_ra_mass'   '/spath_f_ra_mass'  'z' }

    { 'DIAGS/hist_meas_ts_aero_<CASE>.h5' '/spath_i_aero_mass'   '/spath_i_aero_mass'  'z' }
    { 'DIAGS/hist_meas_ts_aero_<CASE>.h5' '/spath_m_aero_mass'   '/spath_m_aero_mass'  'z' }
    { 'DIAGS/hist_meas_ts_aero_<CASE>.h5' '/spath_f_aero_mass'   '/spath_f_aero_mass'  'z' }

    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_i_tracer_mass'   '/spath_i_tracer_mass'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_m_tracer_mass'   '/spath_m_tracer_mass'  'z' }
    { 'DIAGS/hist_meas_ts_dust_<CASE>.h5' '/spath_f_tracer_mass'   '/spath_f_tracer_mass'  'z' }

    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/spath_i_dust_hydro'   '/spath_i_dust_hydro'  'z' }
    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/spath_m_dust_hydro'   '/spath_m_dust_hydro'  'z' }
    { 'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5' '/spath_f_dust_hydro'   '/spath_f_dust_hydro'  'z' }

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

      OnFinalProf = ~isempty(regexp(InVname, '_f_'));

      % skip this profile if input file or dataset is missing (eg., NODUST cases don't
      % have d1_num, d1_mass, etc. extractions)
      try
        HINFO = h5info(InFile, InVname);
      catch
        fprintf('  Input file/dataset does not exist, skipping: %s (%s)\n', InFile, InVname);
        continue
      end
      icount = icount + 1;
      
      ControlInFile = regexprep(ProfileList{iset}{1}, '<CASE>', ControlCase);

      OutDiffVname  = sprintf('%s_diff',  OutVname);
      OutDeltaVname = sprintf('%s_delta', OutVname);

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

      % If we are on a "final" profile, create the difference between initial and final
      % call this quatity "delta"
      if (OnFinalProf)
        DeltaVname = regexprep(InVname, '_f_', '_i_');

        fprintf('  Reading: %s (%s)\n', InFile, DeltaVname);
        INIT_VAR = squeeze(h5read(InFile, DeltaVname));
        
        INIT_VAR(isnan(INIT_VAR)) = 0;
        DELTA_VAR = VAR - INIT_VAR;
      end

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
        if (OnFinalProf)
          DELTA_VAR = DELTA_VAR .* RhoAir .* 1e-6;
        end
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

      if (OnFinalProf)
        fprintf('  Writing: %s (%s)\n', OutFile, OutDeltaVname);
        h5create(OutFile, OutDeltaVname, Vsize);
        h5write(OutFile, OutDeltaVname, DELTA_VAR);
        AttachDimensionsXyzt(OutFile, OutDeltaVname, DimOrder, Xname, Yname, Zname, Tname);
      end

      fprintf('\n');
    end
  end
end
