function [ ] = GenTsdHistMeas()

  % make sure output directory exists
  Ddir = 'DIAGS';  % coordinate this with output file names below
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % list of simulation cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_DUST'
    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  % Description of measurements
  %   <measure_name> <measure_list> <out_file>
  %
  %   where <measure_list> is one or more of:
  %     <in_file> <var_name> <reduction_method> <arg_for_reduction_method> <out_var_name> <rt_reduction> <r_start> <r_end> <t_start> <t_end> <bin_select_op> <bin_select_val>
  %         <rt_reduction>:
  %             'all': only reduce bins, result is (r,z,t)
  %             'xsection': reduce time and bins, result is (r,z)
  %             'profile': reduce time, radius and bins, result is (z)
  %         <r_start, r_end> radius in km -> full range is 0 450
  %         <t_start, t_end> sim time in h  -> full range is 0 60
  MeasSets = {
    % storm speed measurements
    {
      'Speed'
      {
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/avg_speed'         'all'        0 450   0 60 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_avg_speed'      'xsection'   0 450  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_avg_speed'       'xsection'   0 450  40 60 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_core_avg_speed' 'profile'    0 100  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_core_avg_speed'  'profile'    0 100  40 60 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_rb_avg_speed'   'profile'  100 250  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_rb_avg_speed'    'profile'  100 250  40 60 'ge' 0   }

        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/avg_speed_t'         'all'        0 450   0 60 'ge' -100   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_avg_speed_t'      'xsection'   0 450  10 30 'ge' -100   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_avg_speed_t'       'xsection'   0 450  40 60 'ge' -100   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_core_avg_speed_t' 'profile'    0 100  10 30 'ge' -100   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_core_avg_speed_t'  'profile'    0 100  40 60 'ge' -100   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_rb_avg_speed_t'   'profile'  100 250  10 30 'ge' -100   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_rb_avg_speed_t'    'profile'  100 250  40 60 'ge' -100   }

        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/avg_speed_r'         'all'        0 450   0 60 'ge' -100   }
        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_avg_speed_r'      'xsection'   0 450  10 30 'ge' -100   }
        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_avg_speed_r'       'xsection'   0 450  40 60 'ge' -100   }
        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_core_avg_speed_r' 'profile'    0 100  10 30 'ge' -100   }
        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_core_avg_speed_r'  'profile'    0 100  40 60 'ge' -100   }
        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_rb_avg_speed_r'   'profile'  100 250  10 30 'ge' -100   }
        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_rb_avg_speed_r'    'profile'  100 250  40 60 'ge' -100   }

        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'farea'   0.99 '/max_speed'         'all'        0 450   0 60 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'farea'   0.99 '/ps_max_speed'      'xsection'   0 450  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'farea'   0.99 '/s_max_speed'       'xsection'   0 450  40 60 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'farea'   0.99 '/ps_core_max_speed' 'profile'    0 100  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'farea'   0.99 '/s_core_max_speed'  'profile'    0 100  40 60 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'farea'   0.99 '/ps_rb_max_speed'   'profile'  100 250  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'farea'   0.99 '/s_rb_max_speed'    'profile'  100 250  40 60 'ge' 0   }

        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'farea'   0.99 '/max_speed_t'         'all'        0 450   0 60 'ge' 0   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'farea'   0.99 '/ps_max_speed_t'      'xsection'   0 450  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'farea'   0.99 '/s_max_speed_t'       'xsection'   0 450  40 60 'ge' 0   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'farea'   0.99 '/ps_core_max_speed_t' 'profile'    0 100  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'farea'   0.99 '/s_core_max_speed_t'  'profile'    0 100  40 60 'ge' 0   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'farea'   0.99 '/ps_rb_max_speed_t'   'profile'  100 250  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'farea'   0.99 '/s_rb_max_speed_t'    'profile'  100 250  40 60 'ge' 0   }

        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'farea'   0.99 '/max_speed10m'         'all'        0 450   0 60 'ge' 0   }
        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'farea'   0.99 '/ps_max_speed10m'      'xsection'   0 450  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'farea'   0.99 '/s_max_speed10m'       'xsection'   0 450  40 60 'ge' 0   }
        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'farea'   0.99 '/ps_core_max_speed10m' 'profile'    0 100  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'farea'   0.99 '/s_core_max_speed10m'  'profile'    0 100  40 60 'ge' 0   }
        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'farea'   0.99 '/ps_rb_max_speed10m'   'profile'  100 250  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'farea'   0.99 '/s_rb_max_speed10m'    'profile'  100 250  40 60 'ge' 0   }

        { 'AzAveragedData/hist_speed25m_<CASE>.h5' '/speed25m' 'farea'   0.99 '/max_speed25m'         'all'        0 450   0 60 'ge' 0   }
        { 'AzAveragedData/hist_speed25m_<CASE>.h5' '/speed25m' 'farea'   0.99 '/ps_max_speed25m'      'xsection'   0 450  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed25m_<CASE>.h5' '/speed25m' 'farea'   0.99 '/s_max_speed25m'       'xsection'   0 450  40 60 'ge' 0   }
        { 'AzAveragedData/hist_speed25m_<CASE>.h5' '/speed25m' 'farea'   0.99 '/ps_core_max_speed25m' 'profile'    0 100  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed25m_<CASE>.h5' '/speed25m' 'farea'   0.99 '/s_core_max_speed25m'  'profile'    0 100  40 60 'ge' 0   }
        { 'AzAveragedData/hist_speed25m_<CASE>.h5' '/speed25m' 'farea'   0.99 '/ps_rb_max_speed25m'   'profile'  100 250  10 30 'ge' 0   }
        { 'AzAveragedData/hist_speed25m_<CASE>.h5' '/speed25m' 'farea'   0.99 '/s_rb_max_speed25m'    'profile'  100 250  40 60 'ge' 0   }
      }
      'DIAGS/hist_meas_speed_<CASE>.h5'
    }

%    {
%      'Pressure'
%      {
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/avg_press'         'all'        0 450   0 60 'ge' 0   }
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/ps_avg_press'      'xsection'   0 450  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/s_avg_press'       'xsection'   0 450  40 60 'ge' 0   }
%
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'farea'   0.01 '/min_press'         'all'        0 450   0 60 'ge' 0   }
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'farea'   0.01 '/ps_min_press'      'xsection'   0 450  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'farea'   0.01 '/s_min_press'       'xsection'   0 450  40 60 'ge' 0   }
%      }
%      'DIAGS/hist_meas_press_<CASE>.h5'
%    }
%
%    % precip rate measurements
%    {
%      'Precip Rate'
%      {
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/avg_pcprate'         'all'        0 450   0 60 'ge' 0   }
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/ps_avg_pcprate'      'xsection'   0 450  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/s_avg_pcprate'       'xsection'   0 450  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/ps_core_avg_pcprate' 'profile'    0 100  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/s_core_avg_pcprate'  'profile'    0 100  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/ps_rb_avg_pcprate'   'profile'  100 250  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/s_rb_avg_pcprate'    'profile'  100 250  40 60 'ge' 0   }
%      }
%      'DIAGS/hist_meas_pcprate_<CASE>.h5'
%    }
%
%    % vertially integrated condensate measurements
%    {
%      'Vert Cond'
%      {
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/avg_vint_cond'         'all'        0 450   0 60 'ge' 0   }
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/ps_avg_vint_cond'      'xsection'   0 450  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/s_avg_vint_cond'       'xsection'   0 450  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/ps_core_avg_vint_cond' 'profile'    0 100  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/s_core_avg_vint_cond'  'profile'    0 100  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/ps_rb_avg_vint_cond'   'profile'  100 250  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/s_rb_avg_vint_cond'    'profile'  100 250  40 60 'ge' 0   }
%      }
%      'DIAGS/hist_meas_vint_cond_<CASE>.h5'
%    }
%
%    % vertical velocity measurements
%    {
%      'Vertical Velocity'
%      {
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/avg_updraft'         'all'        0 450   0 60 'ge' 0.01   }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_avg_updraft'      'xsection'   0 450  10 30 'ge' 0.01   }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_avg_updraft'       'xsection'   0 450  40 60 'ge' 0.01   }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_core_avg_updraft' 'profile'    0 100  10 30 'ge' 0.01   }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_core_avg_updraft'  'profile'    0 100  40 60 'ge' 0.01   }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_rb_avg_updraft'   'profile'  100 250  10 30 'ge' 0.01   }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_rb_avg_updraft'    'profile'  100 250  40 60 'ge' 0.01   }
%
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.99 '/max_updraft'         'all'        0 450   0 60 'ge' 0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.99 '/ps_max_updraft'      'xsection'   0 450  10 30 'ge' 0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.99 '/s_max_updraft'       'xsection'   0 450  40 60 'ge' 0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.99 '/ps_core_max_updraft' 'profile'    0 100  10 30 'ge' 0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.99 '/s_core_max_updraft'  'profile'    0 100  40 60 'ge' 0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.99 '/ps_rb_max_updraft'   'profile'  100 250  10 30 'ge' 0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.99 '/s_rb_max_updraft'    'profile'  100 250  40 60 'ge' 0.05   }
%
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/avg_dndraft'         'all'        0 450   0 60 'le' -0.01   }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_avg_dndraft'      'xsection'   0 450  10 30 'le' -0.01   }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_avg_dndraft'       'xsection'   0 450  40 60 'le' -0.01   }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_core_avg_dndraft' 'profile'    0 100  10 30 'le' -0.01   }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_core_avg_dndraft'  'profile'    0 100  40 60 'le' -0.01   }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_rb_avg_dndraft'   'profile'  100 250  10 30 'le' -0.01   }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_rb_avg_dndraft'    'profile'  100 250  40 60 'le' -0.01   }
%
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.01 '/max_dndraft'         'all'        0 450   0 60 'le' -0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.01 '/ps_max_dndraft'      'xsection'   0 450  10 30 'le' -0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.01 '/s_max_dndraft'       'xsection'   0 450  40 60 'le' -0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.01 '/ps_core_max_dndraft' 'profile'    0 100  10 30 'le' -0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.01 '/s_core_max_dndraft'  'profile'    0 100  40 60 'le' -0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.01 '/ps_rb_max_dndraft'   'profile'  100 250  10 30 'le' -0.05   }
%        { 'AzAveragedData/hist_w_<CASE>.h5'  '/w' 'farea'   0.01 '/s_rb_max_dndraft'    'profile'  100 250  40 60 'le' -0.05   }
%      }
%      'DIAGS/hist_meas_w_<CASE>.h5'
%    }
%
%    % theta_e measurements
%    {
%      'Theta-E'
%      {
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/avg_theta_e'         'all'        0 450   0 60 'ge' 0   }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/ps_avg_theta_e'      'xsection'   0 450  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/s_avg_theta_e'       'xsection'   0 450  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/ps_core_avg_theta_e' 'profile'    0 100  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/s_core_avg_theta_e'  'profile'    0 100  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/ps_rb_avg_theta_e'   'profile'  100 250  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/s_rb_avg_theta_e'    'profile'  100 250  40 60 'ge' 0   }
%      }
%      'DIAGS/hist_meas_theta_e_<CASE>.h5'
%    }
%
%    % dust measurements
%    {
%      'Dust'
%      {
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/avg_d1_mass'         'all'        0 450   0 60 'ge' 0   }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/ps_avg_d1_mass'      'xsection'   0 450  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/s_avg_d1_mass'       'xsection'   0 450  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/ps_core_avg_d1_mass' 'profile'    0 100  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/s_core_avg_d1_mass'  'profile'    0 100  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/ps_rb_avg_d1_mass'   'profile'  100 250  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/s_rb_avg_d1_mass'    'profile'  100 250  40 60 'ge' 0   }
%
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/avg_d1_num'         'all'        0 450   0 60 'ge' 0   }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/ps_avg_d1_num'      'xsection'   0 450  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/s_avg_d1_num'       'xsection'   0 450  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/ps_core_avg_d1_num' 'profile'    0 100  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/s_core_avg_d1_num'  'profile'    0 100  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/ps_rb_avg_d1_num'   'profile'  100 250  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/s_rb_avg_d1_num'    'profile'  100 250  40 60 'ge' 0   }
%
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/avg_d2_mass'         'all'        0 450   0 60 'ge' 0   }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/ps_avg_d2_mass'      'xsection'   0 450  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/s_avg_d2_mass'       'xsection'   0 450  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/ps_core_avg_d2_mass' 'profile'    0 100  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/s_core_avg_d2_mass'  'profile'    0 100  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/ps_rb_avg_d2_mass'   'profile'  100 250  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/s_rb_avg_d2_mass'    'profile'  100 250  40 60 'ge' 0   }
%
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/avg_d2_num'         'all'        0 450   0 60 'ge' 0   }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/ps_avg_d2_num'      'xsection'   0 450  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/s_avg_d2_num'       'xsection'   0 450  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/ps_core_avg_d2_num' 'profile'    0 100  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/s_core_avg_d2_num'  'profile'    0 100  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/ps_rb_avg_d2_num'   'profile'  100 250  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/s_rb_avg_d2_num'    'profile'  100 250  40 60 'ge' 0   }
%      }
%      'DIAGS/hist_meas_dust_<CASE>.h5'
%    }
%
%    % cloud
%    {
%      'Cloud'
%     {
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/avg_cloud'         'all'        0 450   0 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/ps_avg_cloud'      'xsection'   0 450  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/s_avg_cloud'       'xsection'   0 450  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/ps_core_avg_cloud' 'profile'    0 100  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/s_core_avg_cloud'  'profile'    0 100  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/ps_rb_avg_cloud'   'profile'  100 250  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/s_rb_avg_cloud'    'profile'  100 250  40 60 'ge' 0.01 }
%
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/avg_cloud_num'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/ps_avg_cloud_num'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/s_avg_cloud_num'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/ps_core_avg_cloud_num' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/s_core_avg_cloud_num'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/ps_rb_avg_cloud_num'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/s_rb_avg_cloud_num'    'profile'  100 250  40 60 'ge' 0 }
%
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/avg_cloud_diam'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/ps_avg_cloud_diam'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/s_avg_cloud_diam'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/ps_core_avg_cloud_diam' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/s_core_avg_cloud_diam'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/ps_rb_avg_cloud_diam'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/s_rb_avg_cloud_diam'    'profile'  100 250  40 60 'ge' 0 }
%      }
%      'DIAGS/hist_meas_cloud_<CASE>.h5'
%    }
%
%    % rain
%    {
%      'Rain'
%     {
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/avg_rain'         'all'        0 450   0 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/ps_avg_rain'      'xsection'   0 450  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/s_avg_rain'       'xsection'   0 450  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/ps_core_avg_rain' 'profile'    0 100  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/s_core_avg_rain'  'profile'    0 100  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/ps_rb_avg_rain'   'profile'  100 250  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/s_rb_avg_rain'    'profile'  100 250  40 60 'ge' 0.01 }
%
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/avg_rain_num'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/ps_avg_rain_num'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/s_avg_rain_num'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/ps_core_avg_rain_num' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/s_core_avg_rain_num'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/ps_rb_avg_rain_num'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/s_rb_avg_rain_num'    'profile'  100 250  40 60 'ge' 0 }
%
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/avg_rain_diam'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/ps_avg_rain_diam'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/s_avg_rain_diam'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/ps_core_avg_rain_diam' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/s_core_avg_rain_diam'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/ps_rb_avg_rain_diam'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/s_rb_avg_rain_diam'    'profile'  100 250  40 60 'ge' 0 }
%      }
%      'DIAGS/hist_meas_rain_<CASE>.h5'
%    }
%
%    % pristine ice
%    {
%      'Pristine'
%     {
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/avg_pris'         'all'        0 450   0 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/ps_avg_pris'      'xsection'   0 450  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/s_avg_pris'       'xsection'   0 450  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/ps_core_avg_pris' 'profile'    0 100  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/s_core_avg_pris'  'profile'    0 100  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/ps_rb_avg_pris'   'profile'  100 250  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/s_rb_avg_pris'    'profile'  100 250  40 60 'ge' 0.01 }
%
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/avg_pris_num'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/ps_avg_pris_num'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/s_avg_pris_num'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/ps_core_avg_pris_num' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/s_core_avg_pris_num'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/ps_rb_avg_pris_num'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/s_rb_avg_pris_num'    'profile'  100 250  40 60 'ge' 0 }
%
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/avg_pris_diam'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/ps_avg_pris_diam'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/s_avg_pris_diam'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/ps_core_avg_pris_diam' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/s_core_avg_pris_diam'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/ps_rb_avg_pris_diam'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/s_rb_avg_pris_diam'    'profile'  100 250  40 60 'ge' 0 }
%      }
%      'DIAGS/hist_meas_pris_<CASE>.h5'
%    }
%
%    % aggregates
%    {
%      'Aggregates'
%     {
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/avg_aggr'         'all'        0 450   0 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/ps_avg_aggr'      'xsection'   0 450  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/s_avg_aggr'       'xsection'   0 450  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/ps_core_avg_aggr' 'profile'    0 100  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/s_core_avg_aggr'  'profile'    0 100  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/ps_rb_avg_aggr'   'profile'  100 250  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/s_rb_avg_aggr'    'profile'  100 250  40 60 'ge' 0.01 }
%
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/avg_aggr_num'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/ps_avg_aggr_num'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/s_avg_aggr_num'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/ps_core_avg_aggr_num' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/s_core_avg_aggr_num'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/ps_rb_avg_aggr_num'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/s_rb_avg_aggr_num'    'profile'  100 250  40 60 'ge' 0 }
%
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/avg_aggr_diam'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/ps_avg_aggr_diam'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/s_avg_aggr_diam'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/ps_core_avg_aggr_diam' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/s_core_avg_aggr_diam'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/ps_rb_avg_aggr_diam'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/s_rb_avg_aggr_diam'    'profile'  100 250  40 60 'ge' 0 }
%      }
%      'DIAGS/hist_meas_aggr_<CASE>.h5'
%    }
%
%    % snow
%    {
%      'Snow'
%     {
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/avg_snow'         'all'        0 450   0 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/ps_avg_snow'      'xsection'   0 450  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/s_avg_snow'       'xsection'   0 450  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/ps_core_avg_snow' 'profile'    0 100  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/s_core_avg_snow'  'profile'    0 100  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/ps_rb_avg_snow'   'profile'  100 250  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/s_rb_avg_snow'    'profile'  100 250  40 60 'ge' 0.01 }
%
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/avg_snow_num'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/ps_avg_snow_num'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/s_avg_snow_num'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/ps_core_avg_snow_num' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/s_core_avg_snow_num'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/ps_rb_avg_snow_num'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/s_rb_avg_snow_num'    'profile'  100 250  40 60 'ge' 0 }
%
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/avg_snow_diam'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/ps_avg_snow_diam'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/s_avg_snow_diam'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/ps_core_avg_snow_diam' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/s_core_avg_snow_diam'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/ps_rb_avg_snow_diam'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/s_rb_avg_snow_diam'    'profile'  100 250  40 60 'ge' 0 }
%      }
%      'DIAGS/hist_meas_snow_<CASE>.h5'
%    }
%
%    % graupel
%    {
%      'Graupel'
%     {
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/avg_graup'         'all'        0 450   0 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/ps_avg_graup'      'xsection'   0 450  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/s_avg_graup'       'xsection'   0 450  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/ps_core_avg_graup' 'profile'    0 100  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/s_core_avg_graup'  'profile'    0 100  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/ps_rb_avg_graup'   'profile'  100 250  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/s_rb_avg_graup'    'profile'  100 250  40 60 'ge' 0.01 }
%
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/avg_graup_num'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/ps_avg_graup_num'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/s_avg_graup_num'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/ps_core_avg_graup_num' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/s_core_avg_graup_num'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/ps_rb_avg_graup_num'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/s_rb_avg_graup_num'    'profile'  100 250  40 60 'ge' 0 }
%
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/avg_graup_diam'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/ps_avg_graup_diam'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/s_avg_graup_diam'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/ps_core_avg_graup_diam' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/s_core_avg_graup_diam'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/ps_rb_avg_graup_diam'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/s_rb_avg_graup_diam'    'profile'  100 250  40 60 'ge' 0 }
%      }
%      'DIAGS/hist_meas_graup_<CASE>.h5'
%    }
%
%    % hail
%    {
%      'Hail'
%     {
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/avg_hail'         'all'        0 450   0 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/ps_avg_hail'      'xsection'   0 450  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/s_avg_hail'       'xsection'   0 450  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/ps_core_avg_hail' 'profile'    0 100  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/s_core_avg_hail'  'profile'    0 100  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/ps_rb_avg_hail'   'profile'  100 250  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/s_rb_avg_hail'    'profile'  100 250  40 60 'ge' 0.01 }
%
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/avg_hail_num'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/ps_avg_hail_num'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/s_avg_hail_num'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/ps_core_avg_hail_num' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/s_core_avg_hail_num'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/ps_rb_avg_hail_num'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/s_rb_avg_hail_num'    'profile'  100 250  40 60 'ge' 0 }
%
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/avg_hail_diam'         'all'        0 450   0 60 'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/ps_avg_hail_diam'      'xsection'   0 450  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/s_avg_hail_diam'       'xsection'   0 450  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/ps_core_avg_hail_diam' 'profile'    0 100  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/s_core_avg_hail_diam'  'profile'    0 100  40 60 'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/ps_rb_avg_hail_diam'   'profile'  100 250  10 30 'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/s_rb_avg_hail_diam'    'profile'  100 250  40 60 'ge' 0 }
%      }
%      'DIAGS/hist_meas_hail_<CASE>.h5'
%    }
%
%    % total condensate
%    {
%      'Total Condensate'
%     {
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/avg_tcond'         'all'        0 450   0 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/ps_avg_tcond'      'xsection'   0 450  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/s_avg_tcond'       'xsection'   0 450  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/ps_core_avg_tcond' 'profile'    0 100  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/s_core_avg_tcond'  'profile'    0 100  40 60 'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/ps_rb_avg_tcond'   'profile'  100 250  10 30 'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/s_rb_avg_tcond'    'profile'  100 250  40 60 'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_tcond_<CASE>.h5'
%    }
%
%    % cooling via latent heat of freezing
%    {
%      'LHF Cooling'
%     {
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/avg_lhf_cool'         'all'        0 450   0 60 'le' 0   }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/ps_avg_lhf_cool'      'xsection'   0 450  10 30 'le' 0   }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/s_avg_lhf_cool'       'xsection'   0 450  40 60 'le' 0   }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/ps_core_avg_lhf_cool' 'profile'    0 100  10 30 'le' 0   }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/s_core_avg_lhf_cool'  'profile'    0 100  40 60 'le' 0   }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/ps_rb_avg_lhf_cool'   'profile'  100 250  10 30 'le' 0   }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/s_rb_avg_lhf_cool'    'profile'  100 250  40 60 'le' 0   }
%      }
%      'DIAGS/hist_meas_lhf_cool_<CASE>.h5'
%    }
%
%    % heating via latent heat of freezing
%    {
%      'LHF Heating'
%     {
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/avg_lhf_heat'         'all'        0 450   0 60 'ge' 0   }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/ps_avg_lhf_heat'      'xsection'   0 450  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/s_avg_lhf_heat'       'xsection'   0 450  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/ps_core_avg_lhf_heat' 'profile'    0 100  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/s_core_avg_lhf_heat'  'profile'    0 100  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/ps_rb_avg_lhf_heat'   'profile'  100 250  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/s_rb_avg_lhf_heat'    'profile'  100 250  40 60 'ge' 0   }
%      }
%      'DIAGS/hist_meas_lhf_heat_<CASE>.h5'
%    }
%
%    % cooling via latent heat of vaporization
%    {
%      'LHV Cooling'
%     {
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/avg_lhv_cool'         'all'        0 450   0 60 'le' 0   }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/ps_avg_lhv_cool'      'xsection'   0 450  10 30 'le' 0   }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/s_avg_lhv_cool'       'xsection'   0 450  40 60 'le' 0   }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/ps_core_avg_lhv_cool' 'profile'    0 100  10 30 'le' 0   }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/s_core_avg_lhv_cool'  'profile'    0 100  40 60 'le' 0   }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/ps_rb_avg_lhv_cool'   'profile'  100 250  10 30 'le' 0   }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/s_rb_avg_lhv_cool'    'profile'  100 250  40 60 'le' 0   }
%      }
%      'DIAGS/hist_meas_lhv_cool_<CASE>.h5'
%    }
%
%    % heating via latent heat of vaporization
%    {
%      'LHV Heating'
%     {
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/avg_lhv_heat'         'all'        0 450   0 60 'ge' 0   }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/ps_avg_lhv_heat'      'xsection'   0 450  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/s_avg_lhv_heat'       'xsection'   0 450  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/ps_core_avg_lhv_heat' 'profile'    0 100  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/s_core_avg_lhv_heat'  'profile'    0 100  40 60 'ge' 0   }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/ps_rb_avg_lhv_heat'   'profile'  100 250  10 30 'ge' 0   }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/s_rb_avg_lhv_heat'    'profile'  100 250  40 60 'ge' 0   }
%      }
%      'DIAGS/hist_meas_lhv_heat_<CASE>.h5'
%    }

    };

  Nsets = length(MeasSets);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating histogram measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');


    for iset = 1:Nsets
      MeasName     = MeasSets{iset}{1};
      MeasList     = MeasSets{iset}{2};
      OutFtemplate = MeasSets{iset}{3};

      fprintf('    Measurement set: %s\n', MeasName);
      fprintf('\n');

      % Put all measurements into one file per case
      % If the file exists, remove it so that the HDF5 commands
      % can effectively re-create datasets.
      OutFile = regexprep(OutFtemplate, '<CASE>', Case);
      if (exist(OutFile, 'file') == 2)
        delete(OutFile);
      end
  
      Nmeas = length(MeasList);
      for imeas = 1:Nmeas
        Ftemplate = MeasList{imeas}{1};
        Vname     = MeasList{imeas}{2};
        Rmethod   = MeasList{imeas}{3};
        Param     = MeasList{imeas}{4};
        OutVname  = MeasList{imeas}{5};
        Rspec     = MeasList{imeas}{6};
        Rstart    = MeasList{imeas}{7};
        Rend      = MeasList{imeas}{8};
        Tstart    = MeasList{imeas}{9};
        Tend      = MeasList{imeas}{10};
        SelectOp  = MeasList{imeas}{11};
        SelectVal = MeasList{imeas}{12};
  
        InFile = regexprep(Ftemplate, '<CASE>', Case);
  
        fprintf('      Reading: %s (%s)\n', InFile, Vname);
        fprintf('        Reduction spec: %s\n', Rspec);
        fprintf('        Radius reduction range: %.2f to %.2f\n', Rstart, Rend);
        fprintf('        Time reduction range: %.2f to %.2f\n', Tstart, Tend);
        fprintf('        Bin reduction method: %s (%.2f)\n', Rmethod, Param);
        fprintf('        Bin selection: %s %.2f\n', SelectOp, SelectVal);
  
        % Read in data which will be 4D -> (x,y,z,t)
        %
        %     x --> radial bands
        %     y --> histogram bins
        %     z --> height
        %     t --> time
        %
        HDATA = squeeze(h5read(InFile, Vname));
        BINS  = squeeze(h5read(InFile, '/y_coords'));
  
        % Assume same r,z,t values for all measurements
        X    = squeeze(h5read(InFile, '/x_coords'));
        Y    = 1; % dummy dimension for output (y represents bins which always get reduced)
        Z    = squeeze(h5read(InFile, '/z_coords'));
        T    = squeeze(h5read(InFile, '/t_coords'));
  
        Nx = length(X);
        Ny = length(Y);
        Nz = length(Z);
        Nt = length(T);

        % If first measurement, then write out coordinates for later
        % use in attaching vars to them.
        if (imeas == 1)
          Xname = '/x_coords';
          Yname = '/y_coords';
          Zname = '/z_coords';
          Tname = '/t_coords';
  
          CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
          % Add COARDS annotations
          NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
        end

        % Calculate indices for selecting by radius, bins and time.
        % Do selection on bin values (ge or le)
        % Select all bin values by default
        R = X ./ 1000;
        R1 = find(R >= Rstart, 1, 'first');
        R2 = find(R <= Rend  , 1, 'last');

        B1 = 1;         
        B2 = length(BINS);
        if (strcmp(SelectOp, 'ge'))
          B1 = find(BINS >= SelectVal, 1, 'first');
        end
        if (strcmp(SelectOp, 'le'))
          B2 = find(BINS <= SelectVal, 1, 'last');
        end

        SIM_TIME = (T ./ 3600) - 42;
        T1 = find(SIM_TIME >= Tstart, 1, 'first');
        T2 = find(SIM_TIME <= Tend  , 1, 'last');

        % Reduce the histograms via three different styles:
        %   all: reduce only bins, result is (r,z,t)
        %   xsection: reduce bins and time, result is (r,z)
        %   profile: reduce bins, radius and time, result is (z)
        switch(Rspec)
          case 'all'
            if (Nz == 1)
              MEAS = squeeze(ReduceHists(HDATA(:,B1:B2,:), 2, BINS(B1:B2), Rmethod, Param));
              OutSize = [ Nx Nt ];
              DimOrder = { 'x' 't' };
            else
              MEAS = squeeze(ReduceHists(HDATA(:,B1:B2,:,:), 2, BINS(B1:B2), Rmethod, Param));
              OutSize = [ Nx Nz Nt ];
              DimOrder = { 'x' 'z' 't' };
            end
            
          case 'xsection'
            if (Nz == 1)
              HDATA = squeeze(sum(HDATA(:,:,T1:T2), 3)); % sum up counts in time dimension
              MEAS = squeeze(ReduceHists(HDATA(:,B1:B2), 2, BINS(B1:B2), Rmethod, Param));
              OutSize = Nx;
              DimOrder = { 'x' };
            else
              HDATA = squeeze(sum(HDATA(:,:,:,T1:T2), 4)); % sum up counts in time dimension
              MEAS = squeeze(ReduceHists(HDATA(:,B1:B2,:), 2, BINS(B1:B2), Rmethod, Param));
              OutSize = [ Nx Nz ];
              DimOrder = { 'x' 'z' };
            end

          case 'profile'
            if (Nz == 1)
              HDATA = squeeze(sum(HDATA(:,:,T1:T2), 3)); % sum up counts in time dimension
              HDATA = squeeze(sum(HDATA(R1:R2,:), 1)); % sum up counts in radius dimension
              MEAS = squeeze(ReduceHist1D(HDATA(B1:B2)', BINS(B1:B2), Rmethod, Param));
              OutSize = 1;
              DimOrder = { };
            else
              HDATA = squeeze(sum(HDATA(:,:,:,T1:T2), 4)); % sum up counts in time dimension
              HDATA = squeeze(sum(HDATA(R1:R2,:,:), 1)); % sum up counts in radius dimension
              MEAS = squeeze(ReduceHists(HDATA(B1:B2,:), 1, BINS(B1:B2), Rmethod, Param));
              OutSize = Nz;
              DimOrder = { 'z' };
            end

          case default
            fprintf('\n');
            fprintf('ERROR: unrecognized reduction spec (%s), skipping this measurement\n', Rspec);
            continue;
        end

        % Write out measurement
        fprintf('      Writing: %s (%s)\n', OutFile, OutVname)
        fprintf('\n');
  
        % Write out measurement (x,z,t)
        % attach code below.
        h5create(OutFile, OutVname, OutSize);
        h5write(OutFile, OutVname, MEAS);

        % Attach dimensions
        AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
      end % measurements
    end % sets
  end % cases
end % function
