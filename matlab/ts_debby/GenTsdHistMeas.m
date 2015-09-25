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

  % Temporal ranges
  TstartPreSal = 10;
  TendPreSal   = 30;

  TstartSal = 40;
  TendSal   = 60;

  % Radial ranges
  RstartCore =   0;
  RendCore   = 100;

  RstartRband = 100;
  RendRband   = 250;

  % Height ranges
  Zsfc = 0;
  Ztop = 1.5; % km

  % Input is a histogram based upon either 2D or 3D field
  %   2D - input will be of the form (r,b,t)
  %   3D - input will be of the form (r,b,z,t)
  %
  %   where: r - radius
  %          b - bins (histogram counts)
  %          z - height
  %          t - time

  % Output can be specified as one of the following reduction types:
  %
  %     type      2D result       3D result
  %   xsection      N/A            (r,z)
  %   z_profile     N/A             (z)
  %   r_profile     (r)             (r)
  %   
  %   xsection_ts   N/A            (r,z,t)
  %   z_profile_ts  N/A             (z,t)
  %   r_profile_ts  (r,t)           (r,t)
  %   
  % The idea for reduction is to combine (sum) all of the histogram counts first, then
  % reduce the histogram counts once (as opposed to reducing the histogram counts first
  % followed by averaging to combine along the other dimensions).
  %
  % Another notion is to avoid combining counts along the z axis, i.e., preserve
  % levels. This means that to get an r_profile from a 3D field, a single level
  % must be selected. 
  %
  % The temporal, radial and height ranges above are used for reduction. Temporal and
  % radial ranges are used for summing histogram counts, and height is used for selecting
  % a single level.
  %
  % Range specs for summing counts:
  %
  %    Radial range
  %      core            RstartCore, RendCore
  %      rband           RstartRband, RendRband
  %
  %    Temporal range
  %       pre_sal        TstartPreSal, TendPreSal
  %       sal            TstartSal, TendSal
  %
  % Range specs for selecting a level:
  %
  %    Height range
  %       sfc            Zsfc
  %       maxlev         Go through histogram reduction steps, then select level that
  %                      contains the maximum value AND lies between Zsfc and Ztop.
  %       minlev         Same as maxlev, except select the level with the minumum value.

  % Description of measurements
  %   <measure_name> <measure_list> <out_file>
  %
  %   where <measure_list> is one or more of:
  %     <in_file> <var_name> <reduction_method> <arg_for_reduction_method> <out_var_name> <reduction> <r_range> <t_range> <z_range> <bin_select_op> <bin_select_val>
  %     ranges are defined above

  MeasSets = {
%    % storm speed measurements
%    {
%      'Speed'
%      {
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/speed_ts'           'xsection_ts' ''      ''        ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_speed'           'xsection'    ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_speed'            'xsection'    ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_core_speed'      'zprofile'    'core'  'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_core_speed'       'zprofile'    'core'  'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_rb_speed'        'zprofile'    'rband' 'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_rb_speed'         'zprofile'    'rband' 'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_speed_sfc'       'rprofile'    ''      'pre_sal' 'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_speed_sfc'        'rprofile'    ''      'sal'     'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_speed_maxlev'    'rprofile'    ''      'pre_sal' 'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_speed_maxlev'     'rprofile'    ''      'sal'     'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/speed_sfc_ts'       'rprofile_ts' ''      ''        'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/speed_maxlev_ts'    'rprofile_ts' ''      ''        'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_speed_sfc_ts'    'rprofile_ts' ''      'pre_sal' 'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_speed_maxlev_ts' 'rprofile_ts' ''      'pre_sal' 'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_speed_sfc_ts'     'rprofile_ts' ''      'sal'     'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_speed_maxlev_ts'  'rprofile_ts' ''      'sal'     'maxlev' 'ge' 0   }
%
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/speed_t_ts'           'xsection_ts' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_speed_t'           'xsection'    ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_speed_t'            'xsection'    ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_core_speed_t'      'zprofile'    'core'  'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_core_speed_t'       'zprofile'    'core'  'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_rb_speed_t'        'zprofile'    'rband' 'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_rb_speed_t'         'zprofile'    'rband' 'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_speed_t_sfc'       'rprofile'    ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_speed_t_sfc'        'rprofile'    ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_speed_t_maxlev'    'rprofile'    ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_speed_t_maxlev'     'rprofile'    ''      'sal'     'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/speed_t_sfc_ts'       'rprofile_ts' ''      ''        'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/speed_t_maxlev_ts'    'rprofile_ts' ''      ''        'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_speed_t_sfc_ts'    'rprofile_ts' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_speed_t_maxlev_ts' 'rprofile_ts' ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_speed_t_sfc_ts'     'rprofile_ts' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_speed_t_maxlev_ts'  'rprofile_ts' ''      'sal'     'maxlev' 'ge' -100 }
%
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/speed_r_ts'           'xsection_ts' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_speed_r'           'xsection'    ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_speed_r'            'xsection'    ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_core_speed_r'      'zprofile'    'core'  'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_core_speed_r'       'zprofile'    'core'  'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_rb_speed_r'        'zprofile'    'rband' 'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_rb_speed_r'         'zprofile'    'rband' 'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_speed_r_sfc'       'rprofile'    ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_speed_r_sfc'        'rprofile'    ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_speed_r_minlev'    'rprofile'    ''      'pre_sal' 'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_speed_r_minlev'     'rprofile'    ''      'sal'     'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/speed_r_sfc_ts'       'rprofile_ts' ''      ''        'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/speed_r_minlev_ts'    'rprofile_ts' ''      ''        'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_speed_r_sfc_ts'    'rprofile_ts' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_speed_r_maxlev_ts' 'rprofile_ts' ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_speed_r_sfc_ts'     'rprofile_ts' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_speed_r_maxlev_ts'  'rprofile_ts' ''      'sal'     'maxlev' 'ge' -100 }
%
%        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/ps_speed10m'      'rprofile'    ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/s_speed10m'       'rprofile'    ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/speed10m_ts'      'rprofile_ts' ''      ''        ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/ps_speed10m_ts'   'rprofile_ts' ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/s_speed10m_ts'    'rprofile_ts' ''      'sal'     ''       'ge' 0   }
%      }
%      'DIAGS/hist_meas_speed_<CASE>.h5'
%    }

    {
      'Pressure'
      {
        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/press_ts'         'xsection_ts' ''      ''        ''       'ge' 0   }
        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/ps_press'         'xsection'    ''      'pre_sal' ''       'ge' 0   }
        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/s_press'          'xsection'    ''      'sal'     ''       'ge' 0   }
        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/ps_press_sfc'     'rprofile'    ''      'pre_sal' 'sfc'    'ge' 0   }
        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/s_press_sfc'      'rprofile'    ''      'sal'     'sfc'    'ge' 0   }
        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/press_sfc_ts'     'rprofile_ts' ''      ''        'sfc'    'ge' 0   }
        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/ps_press_sfc_ts'  'rprofile_ts' ''      'pre_sal' 'sfc'    'ge' 0   }
        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/s_press_sfc_ts'   'rprofile_ts' ''      'sal'     'sfc'    'ge' 0   }
      }
      'DIAGS/hist_meas_press_<CASE>.h5'
    }

    % precip rate measurements
    {
      'Precip Rate'
      {
        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/ps_pcprate'      'rprofile'    ''      'pre_sal' ''       'ge' 0   }
        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/s_pcprate'       'rprofile'    ''      'sal'     ''       'ge' 0   }
        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/pcprate_ts'      'rprofile_ts' ''      ''        ''       'ge' 0   }
        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/ps_pcprate_ts'   'rprofile_ts' ''      'pre_sal' ''       'ge' 0   }
        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/s_pcprate_ts'    'rprofile_ts' ''      'sal'     ''       'ge' 0   }
      }
      'DIAGS/hist_meas_pcprate_<CASE>.h5'
    }

%    % vertially integrated condensate measurements
%    {
%      'Vert Cond'
%      {
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/ps_vint_cond'      'rprofile'    ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/s_vint_cond'       'rprofile'    ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/vint_cond_ts'      'rprofile_ts' ''      ''        ''       'ge' 0   }
%      }
%      'DIAGS/hist_meas_vint_cond_<CASE>.h5'
%    }
%
%    % vertical velocity measurements
%    {
%      'Vertical Velocity'
%      {
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/updraft_ts'         'xsection_ts' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_updraft'         'xsection'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_updraft'          'xsection'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/core_updraft_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/rb_updraft_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_core_updraft'    'zprofile'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_core_updraft'     'zprofile'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_rb_updraft'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_rb_updraft'       'zprofile'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/dndraft_ts'         'xsection_ts' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_dndraft'         'xsection'    ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_dndraft'          'xsection'    ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/core_dndraft_ts'    'zprofile_ts' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/rb_dndraft_ts'      'zprofile_ts' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_core_dndraft'    'zprofile'    'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_core_dndraft'     'zprofile'    'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_rb_dndraft'      'zprofile'    'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_rb_dndraft'       'zprofile'    'rband' 'sal'     ''       'le' -0.01 }
%      }
%      'DIAGS/hist_meas_w_<CASE>.h5'
%    }
%
%    % theta_e measurements
%    {
%      'Theta-E'
%      {
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/theta_e_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/ps_theta_e'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/s_theta_e'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/core_theta_e_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/rb_theta_e_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/ps_core_theta_e'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/s_core_theta_e'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/ps_rb_theta_e'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/s_rb_theta_e'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_theta_e_<CASE>.h5'
%    }
%
%    % dust measurements
%    {
%      'Dust'
%      {
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/d1_mass_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/ps_d1_mass'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/s_d1_mass'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/core_d1_mass_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/rb_d1_mass_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/ps_core_d1_mass'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/s_core_d1_mass'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/ps_rb_d1_mass'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/s_rb_d1_mass'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/d1_num_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/ps_d1_num'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/s_d1_num'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/core_d1_num_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/rb_d1_num_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/ps_core_d1_num'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/s_core_d1_num'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/ps_rb_d1_num'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/s_rb_d1_num'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/d2_mass_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/ps_d2_mass'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/s_d2_mass'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/core_d2_mass_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/rb_d2_mass_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/ps_core_d2_mass'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/s_core_d2_mass'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/ps_rb_d2_mass'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/s_rb_d2_mass'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/d2_num_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/ps_d2_num'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/s_d2_num'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/core_d2_num_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/rb_d2_num_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/ps_core_d2_num'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/s_core_d2_num'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/ps_rb_d2_num'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/s_rb_d2_num'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_dust_<CASE>.h5'
%    }
%
%    % cloud
%    {
%      'Cloud'
%      {
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/cloud_ts'         'xsection_ts' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/ps_cloud'         'xsection'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/s_cloud'          'xsection'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/core_cloud_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/rb_cloud_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/ps_core_cloud'    'zprofile'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/s_core_cloud'     'zprofile'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/ps_rb_cloud'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/s_rb_cloud'       'zprofile'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/cloud_num_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/ps_cloud_num'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/s_cloud_num'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/core_cloud_num_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/rb_cloud_num_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/ps_core_cloud_num'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/s_core_cloud_num'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/ps_rb_cloud_num'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/s_rb_cloud_num'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/cloud_diam_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/ps_cloud_diam'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/s_cloud_diam'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/core_cloud_diam_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/rb_cloud_diam_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/ps_core_cloud_diam'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/s_core_cloud_diam'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/ps_rb_cloud_diam'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/s_rb_cloud_diam'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_cloud_<CASE>.h5'
%    }
%
%    % rain
%    {
%      'Rain'
%      {
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/rain_ts'         'xsection_ts' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/ps_rain'         'xsection'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/s_rain'          'xsection'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/core_rain_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/rb_rain_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/ps_core_rain'    'zprofile'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/s_core_rain'     'zprofile'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/ps_rb_rain'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/s_rb_rain'       'zprofile'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/rain_num_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/ps_rain_num'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/s_rain_num'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/core_rain_num_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/rb_rain_num_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/ps_core_rain_num'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/s_core_rain_num'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/ps_rb_rain_num'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/s_rb_rain_num'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/rain_diam_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/ps_rain_diam'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/s_rain_diam'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/core_rain_diam_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/rb_rain_diam_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/ps_core_rain_diam'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/s_core_rain_diam'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/ps_rb_rain_diam'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/s_rb_rain_diam'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_rain_<CASE>.h5'
%    }
%
%    % pristine ice
%    {
%      'Pristine'
%      {
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/pris_ts'         'xsection_ts' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/ps_pris'         'xsection'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/s_pris'          'xsection'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/core_pris_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/rb_pris_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/ps_core_pris'    'zprofile'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/s_core_pris'     'zprofile'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/ps_rb_pris'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/s_rb_pris'       'zprofile'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/pris_num_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/ps_pris_num'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/s_pris_num'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/core_pris_num_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/rb_pris_num_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/ps_core_pris_num'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/s_core_pris_num'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/ps_rb_pris_num'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/s_rb_pris_num'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/pris_diam_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/ps_pris_diam'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/s_pris_diam'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/core_pris_diam_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/rb_pris_diam_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/ps_core_pris_diam'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/s_core_pris_diam'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/ps_rb_pris_diam'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/s_rb_pris_diam'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_pris_<CASE>.h5'
%    }
%
%    % aggregates
%    {
%      'Aggregates'
%      {
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/aggr_ts'         'xsection_ts' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/ps_aggr'         'xsection'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/s_aggr'          'xsection'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/core_aggr_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/rb_aggr_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/ps_core_aggr'    'zprofile'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/s_core_aggr'     'zprofile'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/ps_rb_aggr'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/s_rb_aggr'       'zprofile'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/aggr_num_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/ps_aggr_num'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/s_aggr_num'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/core_aggr_num_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/rb_aggr_num_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/ps_core_aggr_num'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/s_core_aggr_num'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/ps_rb_aggr_num'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/s_rb_aggr_num'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/aggr_diam_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/ps_aggr_diam'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/s_aggr_diam'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/core_aggr_diam_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/rb_aggr_diam_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/ps_core_aggr_diam'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/s_core_aggr_diam'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/ps_rb_aggr_diam'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/s_rb_aggr_diam'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_aggr_<CASE>.h5'
%    }
%
%    % snow
%    {
%      'Snow'
%      {
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/snow_ts'         'xsection_ts' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/ps_snow'         'xsection'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/s_snow'          'xsection'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/core_snow_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/rb_snow_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/ps_core_snow'    'zprofile'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/s_core_snow'     'zprofile'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/ps_rb_snow'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/s_rb_snow'       'zprofile'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/snow_num_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/ps_snow_num'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/s_snow_num'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/core_snow_num_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/rb_snow_num_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/ps_core_snow_num'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/s_core_snow_num'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/ps_rb_snow_num'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/s_rb_snow_num'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/snow_diam_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/ps_snow_diam'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/s_snow_diam'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/core_snow_diam_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/rb_snow_diam_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/ps_core_snow_diam'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/s_core_snow_diam'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/ps_rb_snow_diam'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/s_rb_snow_diam'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_snow_<CASE>.h5'
%    }
%
%    % graupel
%    {
%      'Graupel'
%      {
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/graup_ts'         'xsection_ts' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/ps_graup'         'xsection'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/s_graup'          'xsection'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/core_graup_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/rb_graup_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/ps_core_graup'    'zprofile'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/s_core_graup'     'zprofile'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/ps_rb_graup'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/s_rb_graup'       'zprofile'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/graup_num_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/ps_graup_num'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/s_graup_num'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/core_graup_num_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/rb_graup_num_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/ps_core_graup_num'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/s_core_graup_num'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/ps_rb_graup_num'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/s_rb_graup_num'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/graup_diam_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/ps_graup_diam'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/s_graup_diam'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/core_graup_diam_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/rb_graup_diam_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/ps_core_graup_diam'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/s_core_graup_diam'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/ps_rb_graup_diam'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/s_rb_graup_diam'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_graup_<CASE>.h5'
%    }
%
%    % hail
%    {
%      'Hail'
%      {
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/hail_ts'         'xsection_ts' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/ps_hail'         'xsection'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/s_hail'          'xsection'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/core_hail_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/rb_hail_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/ps_core_hail'    'zprofile'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/s_core_hail'     'zprofile'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/ps_rb_hail'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/s_rb_hail'       'zprofile'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/hail_num_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/ps_hail_num'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/s_hail_num'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/core_hail_num_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/rb_hail_num_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/ps_core_hail_num'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/s_core_hail_num'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/ps_rb_hail_num'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/s_rb_hail_num'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/hail_diam_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/ps_hail_diam'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/s_hail_diam'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/core_hail_diam_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/rb_hail_diam_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/ps_core_hail_diam'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/s_core_hail_diam'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/ps_rb_hail_diam'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/s_rb_hail_diam'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_hail_<CASE>.h5'
%    }
%
%    % total condensate
%    {
%      'Total Condensate'
%     {
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/tcond_ts'         'xsection_ts' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/ps_tcond'         'xsection'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/s_tcond'          'xsection'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/core_tcond_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/rb_tcond_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/ps_core_tcond'    'zprofile'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/s_core_tcond'     'zprofile'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/ps_rb_tcond'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/s_rb_tcond'       'zprofile'    'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_tcond_<CASE>.h5'
%    }
%
%    % cooling via latent heat of freezing
%    {
%      'LHF Cooling'
%     {
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/lhf_cool_ts'         'xsection_ts' ''      ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/ps_lhf_cool'         'xsection'    ''      'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/s_lhf_cool'          'xsection'    ''      'sal'     ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/core_lhf_cool_ts'    'zprofile_ts' 'core'  ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/rb_lhf_cool_ts'      'zprofile_ts' 'rband' ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/ps_core_lhf_cool'    'zprofile'    'core'  'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/s_core_lhf_cool'     'zprofile'    'core'  'sal'     ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/ps_rb_lhf_cool'      'zprofile'    'rband' 'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/s_rb_lhf_cool'       'zprofile'    'rband' 'sal'     ''       'le' 0 }
%      }
%      'DIAGS/hist_meas_lhf_cool_<CASE>.h5'
%    }
%
%    % heating via latent heat of freezing
%    {
%      'LHF Heating'
%     {
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/lhf_heat_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/ps_lhf_heat'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/s_lhf_heat'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/core_lhf_heat_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/rb_lhf_heat_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/ps_core_lhf_heat'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/s_core_lhf_heat'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/ps_rb_lhf_heat'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/s_rb_lhf_heat'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_lhf_heat_<CASE>.h5'
%    }
%
%    % cooling via latent heat of vaporization
%    {
%      'LHV Cooling'
%     {
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/lhv_cool_ts'         'xsection_ts' ''      ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/ps_lhv_cool'         'xsection'    ''      'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/s_lhv_cool'          'xsection'    ''      'sal'     ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/core_lhv_cool_ts'    'zprofile_ts' 'core'  ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/rb_lhv_cool_ts'      'zprofile_ts' 'rband' ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/ps_core_lhv_cool'    'zprofile'    'core'  'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/s_core_lhv_cool'     'zprofile'    'core'  'sal'     ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/ps_rb_lhv_cool'      'zprofile'    'rband' 'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/s_rb_lhv_cool'       'zprofile'    'rband' 'sal'     ''       'le' 0 }
%      }
%      'DIAGS/hist_meas_lhv_cool_<CASE>.h5'
%    }
%
%    % heating via latent heat of vaporization
%    {
%      'LHV Heating'
%     {
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/lhv_heat_ts'         'xsection_ts' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/ps_lhv_heat'         'xsection'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/s_lhv_heat'          'xsection'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/core_lhv_heat_ts'    'zprofile_ts' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/rb_lhv_heat_ts'      'zprofile_ts' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/ps_core_lhv_heat'    'zprofile'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/s_core_lhv_heat'     'zprofile'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/ps_rb_lhv_heat'      'zprofile'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/s_rb_lhv_heat'       'zprofile'    'rband' 'sal'     ''       'ge' 0 }
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
        Rrange    = MeasList{imeas}{7};
        Trange    = MeasList{imeas}{8};
        Zrange    = MeasList{imeas}{9};
        SelectOp  = MeasList{imeas}{10};
        SelectVal = MeasList{imeas}{11};
  
        InFile = regexprep(Ftemplate, '<CASE>', Case);

        fprintf('      Reading: %s (%s)\n', InFile, Vname);
        fprintf('        Reduction spec: %s\n', Rspec);
        fprintf('        Bin reduction method: %s (%.2f)\n', Rmethod, Param);
        fprintf('        Bin selection: %s %.2f\n', SelectOp, SelectVal);
  
        % Read in data which will be of two forms: (x,y,z,t) or (x,y,t)
        %
        %     x --> radial bands, r
        %     y --> histogram bins, b
        %     z --> height
        %     t --> time
        %
        HDATA = squeeze(h5read(InFile, Vname));
        X     = squeeze(h5read(InFile, '/x_coords'));
        Y     = squeeze(h5read(InFile, '/y_coords'));
        Z     = squeeze(h5read(InFile, '/z_coords'));
        T     = squeeze(h5read(InFile, '/t_coords'));
  
        % Determine dimensionality of input
        Ndims = ndims(HDATA) - 1;  % number of dimensions minus the time dimension
                                   %  Ndims == 2 -> input is (r,b,t)
                                   %  Ndims == 3 -> input is (r,b,z,t)

        % Determine indices corresponding to r,b,z,t range specs
        % Default is entire range of dimension sizes
        % RADIUS
        R = X ./ 1000;
        R1 = 1;
        R2 = length(R);
        switch(Rrange)
          case 'core'
            R1 = find(R >= RstartCore, 1, 'first');
            R2 = find(R <= RendCore,   1, 'last');

          case 'rband'
            R1 = find(R >= RstartRband, 1, 'first');
            R2 = find(R <= RendRband,   1, 'last');
        end
        fprintf('        Radius selection range: %.2f to %.2f (%d:%d)\n', R(R1), R(R2), R1, R2);

        % BINS
        B = Y;
        B1 = 1;         
        B2 = length(B);
        switch(SelectOp)
          case 'ge'
            B1 = find(B >= SelectVal, 1, 'first');

          case 'le'
            B2 = find(B <= SelectVal, 1, 'last');
        end
        fprintf('        Bin selection range: %.2f to %.2f (%d:%d)\n', B(B1), B(B2), B1, B2);

        % HEIGHT
        H = Z ./ 1000;
        Z1 = 1;
        Z2 = length(H);
        if (Ndims == 3)
          switch(Zrange)
            case { 'maxlev' 'minlev' 'sfc' }
              Z1 = find(H >= Zsfc, 1, 'first');
              Z2 = find(H <= Ztop, 1, 'last');
          end
        end
        fprintf('        Height selection range: %.2f to %.2f (%d:%d)\n', H(Z1), H(Z2), Z1, Z2);

        % SIM TIME
        ST = (T ./ 3600) - 42;
        T1 = 1;
        T2 = length(ST);
        switch(Trange)
          case 'pre_sal'
            T1 = find(ST >= TstartPreSal, 1, 'first');
            T2 = find(ST <= TendPreSal,   1, 'last');

          case 'sal'
            T1 = find(ST >= TstartSal, 1, 'first');
            T2 = find(ST <= TendSal,   1, 'last');
        end
        fprintf('        Time selection range: %.2f to %.2f (%d:%d)\n', ST(T1), ST(T2), T1, T2);

        % Trim down the coordinates according to the selection indices
        R  = R(R1:R2);
        B  = B(B1:B2);
        H  = H(Z1:Z2);
        ST = ST(T1:T2);

        Nr = length(R);
        Nz = length(H);
        Nt = length(ST);

        % Trim down the input data according to the selection indices
        if (Ndims == 2)
          HDATA = HDATA(R1:R2,B1:B2,T1:T2);
        else
          HDATA = HDATA(R1:R2,B1:B2,Z1:Z2,T1:T2);
        end
  
        % If first measurement, then write out coordinates for later
        % use in attaching vars to them. Write out the original coordinates,
        % except for Y (histogram bin counts) which will always be reduced. The
        % original coordinates with their original lengths will always be the
        % appropriate vectors to use for attaching to variables.
        if (imeas == 1)
          Xname = '/x_coords';
          Yname = '/y_coords';
          Zname = '/z_coords';
          Tname = '/t_coords';
  
          CreateDimensionsXyzt(OutFile, X, [ 1 ], Z, T, Xname, Yname, Zname, Tname);
          % Add COARDS annotations
          NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
        end

        % Reduce the histograms according to reduction spec
        if (Ndims == 2)
          % 2D field input, HDATA is (r,b,t)
          switch(Rspec)
            case { 'xsection' 'xsection_ts' 'zprofile' 'zprofile_ts' }
              fprintf('\n');
              fprintf('ERROR: Cannot use "%s" reduction with 2D field input, skipping this measurement\n', Rspec);
              continue;

            case 'rprofile'
              % (r,b,t) -> (r)
              HDATA = squeeze(sum(HDATA, 3)); % sum up counts in time dimension
              MEAS = squeeze(ReduceHists(HDATA, 2, B, Rmethod, Param));
              OutSize = Nr;
              DimOrder = { 'x' };

            case 'rprofile_ts'
              % (r,b,t) -> (r,t)
              MEAS = squeeze(ReduceHists(HDATA, 2, B, Rmethod, Param));
              OutSize = [ Nr Nt ];
              DimOrder = { 'x' 't' };

            case default
            fprintf('\n');
            fprintf('ERROR: unrecognized reduction spec (%s), skipping this measurement\n', Rspec);
            continue;
          end
        else
          % 3D field input, HDATA is (r,b,z,t)
          switch(Rspec)
            case 'xsection'
              % (r,b,z,t) -> (r,z)
              HDATA = squeeze(sum(HDATA, 4)); % sum up counts in time dimension
              MEAS = squeeze(ReduceHists(HDATA, 2, B, Rmethod, Param));
              OutSize = [ Nr Nz ];
              DimOrder = { 'x' 'z' };

            case 'zprofile'
              % (r,b,z,t) -> (z)
              HDATA = squeeze(sum(HDATA, 4)); % sum up counts in time dimension
              HDATA = squeeze(sum(HDATA, 1)); % sum up counts in radius dimension
              MEAS = squeeze(ReduceHists(HDATA, 1, B, Rmethod, Param));
              OutSize = Nz;
              DimOrder = { 'z' };

            case 'rprofile'
              % (r,b,z,t) -> (r)
              HDATA = squeeze(sum(HDATA, 4)); % sum up counts in time dimension
              MEAS = squeeze(ReduceHists(HDATA, 2, B, Rmethod, Param));

              % MEAS is (r,z) at this point, with z trimmed down to the range we
              % are restricting the search for min or max values. Also z == 1
              % represents the surface. Determine which level needs to be selected.
              switch(Zrange)
                case 'sfc'
                  ZSEL = 1;

                case 'minlev'
                  [ SVAL ZSEL ] = min(min(MEAS,[],1));

                case 'maxlev'
                  [ SVAL ZSEL ] = max(max(MEAS,[],1));
              end
              fprintf('        Z level selected for rprofile: %.2f (%d)\n', Z(ZSEL), ZSEL);
              MEAS = squeeze(MEAS(:,ZSEL));
              OutSize = Nr;
              DimOrder = { 'x' };

            case 'xsection_ts'
              % (r,b,z,t) -> (r,z,t)
              MEAS = squeeze(ReduceHists(HDATA, 2, B, Rmethod, Param));
              OutSize = [ Nr Nz Nt ];
              DimOrder = { 'x' 'z' 't' };

            case 'zprofile_ts'
              % (r,b,z,t) -> (z,t)
              HDATA = squeeze(sum(HDATA, 1)); % sum up counts in radius dimension
              MEAS = squeeze(ReduceHists(HDATA, 1, B, Rmethod, Param));
              OutSize = [ Nz Nt ];
              DimOrder = { 'z' 't' };

            case 'rprofile_ts'
              % (r,b,z,t) -> (r,t)
              MEAS = squeeze(ReduceHists(HDATA, 2, B, Rmethod, Param));

              % MEAS is (r,z,t) at this point, with z trimmed down to the range we
              % are restricting the search for min or max values. Also z == 1
              % represents the surface. Determine which level needs to be selected.
              switch(Zrange)
                case 'sfc'
                  ZSEL = 1;

                case 'minlev'
                  [ SVAL ZSEL ] = min(min(min(MEAS,[],3),[],1));

                case 'maxlev'
                  [ SVAL ZSEL ] = max(max(max(MEAS,[],3),[],1));
              end
              fprintf('        Z level selected for rprofile: %.2f (%d)\n', Z(ZSEL), ZSEL);
              MEAS = squeeze(MEAS(:,ZSEL,:));
              OutSize = [ Nr Nt ];
              DimOrder = { 'x' 't' };

            case default
              fprintf('\n');
              fprintf('ERROR: unrecognized reduction spec (%s), skipping this measurement\n', Rspec);
              continue;
          end
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
