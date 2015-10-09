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
  RstartCore =  50;
  RendCore   = 150;

  RstartRband = 200;
  RendRband   = 300;

  RstartEnv = 350;
  RendEnv   = 450;

  RstartLead = 300;
  RendLead   = 500;

  % Height ranges
  Zsfc = 0;
  Ztop = 1.5; % km

  % Input is a histogram based upon either 2D or 3D field
  %
  % Azimuthal averaged data (azavg output):
  %   2D - input will be of the form (r,b,t)
  %   3D - input will be of the form (r,b,z,t)
  %
  % Time series averaged data (tsavg output):
  %   2D - input will be of the form (b,t)
  %   3D - input will be of the form (b,z,t)
  %
  %   where: r - radius
  %          b - bins (histogram counts)
  %          z - height
  %          t - time

  % The input form and reduction is specified by a 4 letter string. The string follows
  % the convention of azavg and tsavg output always writing a 4 dimension variable of
  % the form (x,y,z,t). The reduction spec uses 4 letters that correspond to the
  % xyzt dimension order of the input. If a lower case string is used, then that means
  % to reduce that dimension. If an underscore is used, that means that this input
  % dimension (one of x,y,z,t) is not used. For example, to read in a 3D field from azavg, and reduce
  % to a height profile, use 'rbZt'. This means:
  %
  %   input dim     quatity     reduce
  %       x          radius       yes
  %       y          bins         yes
  %       z          height       no
  %       t          time         yes
  %
  % To produce a time series of height-radius cross sections from azavg 3D output,
  % use 'RbZT'. This means:
  %
  %   input dim     quatity     reduce
  %       x          radius       no
  %       y          bins         yes
  %       z          height       no
  %       t          time         no
  %
  % To produce a radial profile from a 2D azavg output, use 'Rb_t'. This means:
  %
  %
  %   input dim     quatity     reduce
  %       x          radius       no
  %       y          bins         yes
  %       z          unused       ---
  %       t          time         yes
  %
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
  %      core            RstartCore,  RendCore
  %      rband           RstartRband, RendRband
  %      env             RstartEnv,   RendEnv
  %
  %    Temporal range
  %       pre_sal        TstartPreSal, TendPreSal
  %       sal            TstartSal,    TendSal
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
  %     <in_file> <var_name> <reduction_method> <arg_for_reduction_method> <out_var_name> <reduce_spec> <r_range> <t_range> <z_range> <bin_select_op> <bin_select_val>
  %       ranges are defined above
  %       <reduce_spec> describes input dimensions (and order) with upper case letters on the dims that are to be reduced.
  %           <reduce_spec> = 4-character string that corresponds to 'xyzt' that describes which dimensions are radius, bins, height and time,
  %                           as well as sepcifiying which dimensions are to be reduced (lower case).
  %                For example, 'rb_T' means that the input has three dimensions (after squeezing), radius is x, bins are y, and t (z will disappaer
  %                after squeezing), making the input data (r,b,t). The lower case 'r' and 'b' say to reduce the r and b dimensions, leaving
  %                the output a vector in dimension t.
  %
  %                'RbZT' means radius is x, bins are y, height is z and time is t, reduce only the bins, leaving the output (r,z,t).

  MeasSets = {
    % storm speed measurements
%    {
%      'Speed'
%      {
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/speed_ts'           'RbZT' ''      ''        ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_speed'           'RbZt'    ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_speed'            'RbZt'    ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_core_speed'      'rbZt'    'core'  'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_core_speed'       'rbZt'    'core'  'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_rb_speed'        'rbZt'    'rband' 'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_rb_speed'         'rbZt'    'rband' 'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_env_speed'       'rbZt'    'env'   'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_env_speed'        'rbZt'    'env'   'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_speed_sfc'       'Rbzt'    ''      'pre_sal' 'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_speed_sfc'        'Rbzt'    ''      'sal'     'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_speed_maxlev'    'Rbzt'    ''      'pre_sal' 'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_speed_maxlev'     'Rbzt'    ''      'sal'     'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/speed_sfc_ts'       'RbzT' ''      ''        'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/speed_maxlev_ts'    'RbzT' ''      ''        'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_speed_sfc_ts'    'RbzT' ''      'pre_sal' 'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/ps_speed_maxlev_ts' 'RbzT' ''      'pre_sal' 'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_speed_sfc_ts'     'RbzT' ''      'sal'     'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/s_speed_maxlev_ts'  'RbzT' ''      'sal'     'maxlev' 'ge' 0   }
%
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/speed_t_ts'           'RbZT' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_speed_t'           'RbZt'    ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_speed_t'            'RbZt'    ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_core_speed_t'      'rbZt'    'core'  'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_core_speed_t'       'rbZt'    'core'  'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_rb_speed_t'        'rbZt'    'rband' 'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_rb_speed_t'         'rbZt'    'rband' 'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_env_speed_t'       'rbZt'    'env'   'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_env_speed_t'        'rbZt'    'env'   'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_speed_t_sfc'       'Rbzt'    ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_speed_t_sfc'        'Rbzt'    ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_speed_t_maxlev'    'Rbzt'    ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_speed_t_maxlev'     'Rbzt'    ''      'sal'     'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/speed_t_sfc_ts'       'RbzT' ''      ''        'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/speed_t_maxlev_ts'    'RbzT' ''      ''        'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_speed_t_sfc_ts'    'RbzT' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/ps_speed_t_maxlev_ts' 'RbzT' ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_speed_t_sfc_ts'     'RbzT' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/s_speed_t_maxlev_ts'  'RbzT' ''      'sal'     'maxlev' 'ge' -100 }
%
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/speed_r_ts'           'RbZT' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_speed_r'           'RbZt'    ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_speed_r'            'RbZt'    ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_core_speed_r'      'rbZt'    'core'  'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_core_speed_r'       'rbZt'    'core'  'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_rb_speed_r'        'rbZt'    'rband' 'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_rb_speed_r'         'rbZt'    'rband' 'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_env_speed_r'       'rbZt'    'env'   'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_env_speed_r'        'rbZt'    'env'   'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_speed_r_sfc'       'Rbzt'    ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_speed_r_sfc'        'Rbzt'    ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_speed_r_minlev'    'Rbzt'    ''      'pre_sal' 'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_speed_r_minlev'     'Rbzt'    ''      'sal'     'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/speed_r_sfc_ts'       'RbzT' ''      ''        'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/speed_r_minlev_ts'    'RbzT' ''      ''        'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_speed_r_sfc_ts'    'RbzT' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/ps_speed_r_maxlev_ts' 'RbzT' ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_speed_r_sfc_ts'     'RbzT' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/s_speed_r_maxlev_ts'  'RbzT' ''      'sal'     'maxlev' 'ge' -100 }
%
%        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/ps_speed10m'      'Rb_t'    ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/s_speed10m'       'Rb_t'    ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/speed10m_ts'      'Rb_T' ''      ''        ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/ps_speed10m_ts'   'Rb_T' ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/s_speed10m_ts'    'Rb_T' ''      'sal'     ''       'ge' 0   }
%      }
%      'DIAGS/hist_meas_speed_<CASE>.h5'
%    }
%
%    {
%      'Pressure'
%      {
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/press_ts'         'RbZT' ''      ''        ''       'ge' 0   }
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/ps_press'         'RbZt'    ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/s_press'          'RbZt'    ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/ps_press_sfc'     'Rbzt'    ''      'pre_sal' 'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/s_press_sfc'      'Rbzt'    ''      'sal'     'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/press_sfc_ts'     'RbzT' ''      ''        'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/ps_press_sfc_ts'  'RbzT' ''      'pre_sal' 'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/s_press_sfc_ts'   'RbzT' ''      'sal'     'sfc'    'ge' 0   }
%      }
%      'DIAGS/hist_meas_press_<CASE>.h5'
%    }
%
%    % precip rate measurements
%    {
%      'Precip Rate Azavg'
%      {
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/ps_pcprate'      'Rb_t'  ''      'pre_sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/s_pcprate'       'Rb_t'  ''      'sal'     ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/pcprate_ts'      'Rb_T'  ''      ''        ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/ps_pcprate_ts'   'Rb_T'  ''      'pre_sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/s_pcprate_ts'    'Rb_T'  ''      'sal'     ''       'ge' 0.001 }
%      }
%      'DIAGS/hist_meas_az_pcprate_<CASE>.h5'
%    }
%
%    {
%      'Precip Rate Tsavg'
%      {
%        % sample region
%        { 'TsAveragedData/hist_spath_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/spath_pcprate_ts'      'b__T'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/spath_ps_pcprate_ts'   'b__T'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/spath_s_pcprate_ts'    'b__T'   ''      'sal'     ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_spath_pcprate'    'B__t'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_spath_ps_pcprate' 'B__t'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_spath_s_pcprate'  'B__t'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_smaxcp_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/smaxcp_pcprate_ts'      'b__T'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_smaxcp_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/smaxcp_ps_pcprate_ts'   'b__T'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_smaxcp_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/smaxcp_s_pcprate_ts'    'b__T'   ''      'sal'     ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_smaxcp_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_smaxcp_pcprate'    'B__t'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_smaxcp_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_smaxcp_ps_pcprate' 'B__t'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_smaxcp_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_smaxcp_s_pcprate'  'B__t'   ''      'sal'     ''       'ge' 0.001 }
%
%        % lead region
%        { 'TsAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/lead_pcprate_ts'      'b__T'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/lead_ps_pcprate_ts'   'b__T'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/lead_s_pcprate_ts'    'b__T'   ''      'sal'     ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_lead_pcprate'    'B__t'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_lead_ps_pcprate' 'B__t'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_lead_s_pcprate'  'B__t'   ''      'sal'     ''       'ge' 0.001 }
%      }
%      'DIAGS/hist_meas_ts_pcprate_<CASE>.h5'
%    }
%
%    % vertially integrated condensate measurements
%    {
%      'Vert Cond'
%      {
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/ps_vint_cond'      'Rb_t'    ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/s_vint_cond'       'Rb_t'    ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/vint_cond_ts'      'Rb_T' ''      ''        ''       'ge' 0   }
%      }
%      'DIAGS/hist_meas_vint_cond_<CASE>.h5'
%    }
%
%    % vertical velocity measurements
%    {
%      'Vertical Velocity'
%      {
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/updraft_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_updraft'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_updraft'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/core_updraft_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/rb_updraft_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/env_updraft_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_core_updraft'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_core_updraft'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_rb_updraft'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_rb_updraft'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_env_updraft'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_env_updraft'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/dndraft_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_dndraft'         'RbZt'    ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_dndraft'          'RbZt'    ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/core_dndraft_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/rb_dndraft_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/env_dndraft_ts'     'rbZT' 'env'   ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_core_dndraft'    'rbZt'    'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_core_dndraft'     'rbZt'    'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_rb_dndraft'      'rbZt'    'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_rb_dndraft'       'rbZt'    'rband' 'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/ps_env_dndraft'     'rbZt'    'env'   'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/s_env_dndraft'      'rbZt'    'env'   'sal'     ''       'le' -0.01 }
%      }
%      'DIAGS/hist_meas_w_<CASE>.h5'
%    }
%
%    % theta_e measurements
%    {
%      'Theta-E'
%      {
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/theta_e_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/ps_theta_e'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/s_theta_e'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/core_theta_e_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/rb_theta_e_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/env_theta_e_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/ps_core_theta_e'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/s_core_theta_e'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/ps_rb_theta_e'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/s_rb_theta_e'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/ps_env_theta_e'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/s_env_theta_e'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_theta_e_<CASE>.h5'
%    }
%
%    % dust measurements
%    {
%      'Dust Azavg'
%      {
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/d1_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/ps_d1_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/s_d1_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/core_d1_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/rb_d1_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/env_d1_mass_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/ps_core_d1_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/s_core_d1_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/ps_rb_d1_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/s_rb_d1_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/ps_env_d1_mass'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/s_env_d1_mass'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/d1_num_ts'         'RbZT' ''      ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/ps_d1_num'         'RbZt'    ''      'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/s_d1_num'          'RbZt'    ''      'sal'     ''       'ge' 1 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/core_d1_num_ts'    'rbZT' 'core'  ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/rb_d1_num_ts'      'rbZT' 'rband' ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/env_d1_num_ts'     'rbZT' 'env'   ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/ps_core_d1_num'    'rbZt'    'core'  'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/s_core_d1_num'     'rbZt'    'core'  'sal'     ''       'ge' 1 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/ps_rb_d1_num'      'rbZt'    'rband' 'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/s_rb_d1_num'       'rbZt'    'rband' 'sal'     ''       'ge' 1 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/ps_env_d1_num'     'rbZt'    'env'   'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/s_env_d1_num'      'rbZt'    'env'   'sal'     ''       'ge' 1 }
%
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/d2_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/ps_d2_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/s_d2_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/core_d2_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/rb_d2_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/env_d2_mass_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/ps_core_d2_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/s_core_d2_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/ps_rb_d2_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/s_rb_d2_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/ps_env_d2_mass'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/s_env_d2_mass'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/d2_num_ts'         'RbZT' ''      ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/ps_d2_num'         'RbZt'    ''      'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/s_d2_num'          'RbZt'    ''      'sal'     ''       'ge' 1 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/core_d2_num_ts'    'rbZT' 'core'  ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/rb_d2_num_ts'      'rbZT' 'rband' ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/env_d2_num_ts'     'rbZT' 'env'   ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/ps_core_d2_num'    'rbZt'    'core'  'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/s_core_d2_num'     'rbZt'    'core'  'sal'     ''       'ge' 1 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/ps_rb_d2_num'      'rbZt'    'rband' 'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/s_rb_d2_num'       'rbZt'    'rband' 'sal'     ''       'ge' 1 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/ps_env_d2_num'     'rbZt'    'env'   'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/s_env_d2_num'      'rbZt'    'env'   'sal'     ''       'ge' 1 }
%
%        % lead region
%        { 'AzAveragedData/hist_lead_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/l_d1_mass_ts'      'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/l_ps_d1_mass'      'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/l_s_d1_mass'       'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/lead_d1_mass_ts'   'rbZT' 'lead'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/ps_lead_d1_mass'   'rbZt'    'lead'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/s_lead_d1_mass'    'rbZt'    'lead'  'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/l_d1_num_ts'       'RbZT' ''      ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_lead_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/l_ps_d1_num'       'RbZt'    ''      'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_lead_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/l_s_d1_num'        'RbZt'    ''      'sal'     ''       'ge' 1 }
%        { 'AzAveragedData/hist_lead_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/lead_d1_num_ts'    'rbZT' 'lead'  ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_lead_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/ps_lead_d1_num'    'rbZt'    'lead'  'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_lead_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/s_lead_d1_num'     'rbZt'    'lead'  'sal'     ''       'ge' 1 }
%
%        { 'AzAveragedData/hist_lead_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/l_d2_mass_ts'      'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/l_ps_d2_mass'      'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/l_s_d2_mass'       'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/lead_d2_mass_ts'   'rbZT' 'lead'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/ps_lead_d2_mass'   'rbZt'    'lead'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/s_lead_d2_mass'    'rbZt'    'lead'  'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_d2_num_<CASE>.h5'  '/d2_num'  'wtmean'  0.0  '/l_d2_num_ts'       'RbZT' ''      ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_lead_d2_num_<CASE>.h5'  '/d2_num'  'wtmean'  0.0  '/l_ps_d2_num'       'RbZt'    ''      'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_lead_d2_num_<CASE>.h5'  '/d2_num'  'wtmean'  0.0  '/l_s_d2_num'        'RbZt'    ''      'sal'     ''       'ge' 1 }
%        { 'AzAveragedData/hist_lead_d2_num_<CASE>.h5'  '/d2_num'  'wtmean'  0.0  '/lead_d2_num_ts'    'rbZT' 'lead'  ''        ''       'ge' 1 }
%        { 'AzAveragedData/hist_lead_d2_num_<CASE>.h5'  '/d2_num'  'wtmean'  0.0  '/ps_lead_d2_num'    'rbZt'    'lead'  'pre_sal' ''       'ge' 1 }
%        { 'AzAveragedData/hist_lead_d2_num_<CASE>.h5'  '/d2_num'  'wtmean'  0.0  '/s_lead_d2_num'     'rbZt'    'lead'  'sal'     ''       'ge' 1 }
%      }
%      'DIAGS/hist_meas_az_dust_<CASE>.h5'
%    }
%
%    {
%      'Dust Tsavg'
%      {
%        % sample regions in SAL
%        { 'TsAveragedData/hist_spath_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/spath_d1_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/spath_ps_d1_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/spath_s_d1_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_smaxcp_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/smaxcp_d1_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_smaxcp_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/smaxcp_ps_d1_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_smaxcp_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/smaxcp_s_d1_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/spath_d1_num_ts'  'b_ZT'   ''      ''        ''       'ge' 1 }
%        { 'TsAveragedData/hist_spath_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/spath_ps_d1_num'  'b_Zt'   ''      'pre_sal' ''       'ge' 1 }
%        { 'TsAveragedData/hist_spath_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/spath_s_d1_num'   'b_Zt'   ''      'sal'     ''       'ge' 1 }
%
%        { 'TsAveragedData/hist_smaxcp_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/smaxcp_d1_num_ts'  'b_ZT'   ''      ''        ''       'ge' 1 }
%        { 'TsAveragedData/hist_smaxcp_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/smaxcp_ps_d1_num'  'b_Zt'   ''      'pre_sal' ''       'ge' 1 }
%        { 'TsAveragedData/hist_smaxcp_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/smaxcp_s_d1_num'   'b_Zt'   ''      'sal'     ''       'ge' 1 }
%
%        { 'TsAveragedData/hist_spath_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/spath_d2_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/spath_ps_d2_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/spath_s_d2_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_smaxcp_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/smaxcp_d2_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_smaxcp_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/smaxcp_ps_d2_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_smaxcp_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/smaxcp_s_d2_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/spath_d2_num_ts'  'b_ZT'   ''      ''        ''       'ge' 1 }
%        { 'TsAveragedData/hist_spath_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/spath_ps_d2_num'  'b_Zt'   ''      'pre_sal' ''       'ge' 1 }
%        { 'TsAveragedData/hist_spath_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/spath_s_d2_num'   'b_Zt'   ''      'sal'     ''       'ge' 1 }
%
%        { 'TsAveragedData/hist_smaxcp_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/smaxcp_d2_num_ts'  'b_ZT'   ''      ''        ''       'ge' 1 }
%        { 'TsAveragedData/hist_smaxcp_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/smaxcp_ps_d2_num'  'b_Zt'   ''      'pre_sal' ''       'ge' 1 }
%        { 'TsAveragedData/hist_smaxcp_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/smaxcp_s_d2_num'   'b_Zt'   ''      'sal'     ''       'ge' 1 }
%      }
%      'DIAGS/hist_meas_ts_dust_<CASE>.h5'
%    }
%
%    % cloud
%    {
%      'Cloud'
%      {
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/cloud_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/ps_cloud'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/s_cloud'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/core_cloud_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/rb_cloud_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/env_cloud_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/ps_core_cloud'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/s_core_cloud'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/ps_rb_cloud'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/s_rb_cloud'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/ps_env_cloud'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/s_env_cloud'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/cloud_num_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/ps_cloud_num'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/s_cloud_num'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/core_cloud_num_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/rb_cloud_num_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/env_cloud_num_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/ps_core_cloud_num'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/s_core_cloud_num'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/ps_rb_cloud_num'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/s_rb_cloud_num'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/ps_env_cloud_num'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/s_env_cloud_num'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/cloud_diam_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/ps_cloud_diam'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/s_cloud_diam'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/core_cloud_diam_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/rb_cloud_diam_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/env_cloud_diam_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/ps_core_cloud_diam'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/s_core_cloud_diam'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/ps_rb_cloud_diam'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/s_rb_cloud_diam'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/ps_env_cloud_diam'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/s_env_cloud_diam'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_cloud_<CASE>.h5'
%    }
%
%    % rain
%    {
%      'Rain'
%      {
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/rain_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/ps_rain'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/s_rain'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/core_rain_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/rb_rain_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/env_rain_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/ps_core_rain'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/s_core_rain'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/ps_rb_rain'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/s_rb_rain'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/ps_env_rain'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/s_env_rain'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/rain_num_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/ps_rain_num'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/s_rain_num'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/core_rain_num_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/rb_rain_num_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/env_rain_num_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/ps_core_rain_num'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/s_core_rain_num'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/ps_rb_rain_num'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/s_rb_rain_num'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/ps_env_rain_num'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/s_env_rain_num'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/rain_diam_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/ps_rain_diam'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/s_rain_diam'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/core_rain_diam_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/rb_rain_diam_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/env_rain_diam_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/ps_core_rain_diam'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/s_core_rain_diam'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/ps_rb_rain_diam'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/s_rb_rain_diam'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/ps_env_rain_diam'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/s_env_rain_diam'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_rain_<CASE>.h5'
%    }
%
%    % pristine ice
%    {
%      'Pristine'
%      {
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/pris_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/ps_pris'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/s_pris'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/core_pris_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/rb_pris_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/env_pris_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/ps_core_pris'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/s_core_pris'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/ps_rb_pris'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/s_rb_pris'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/ps_env_pris'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/s_env_pris'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/pris_num_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/ps_pris_num'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/s_pris_num'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/core_pris_num_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/rb_pris_num_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/env_pris_num_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/ps_core_pris_num'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/s_core_pris_num'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/ps_rb_pris_num'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/s_rb_pris_num'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/ps_env_pris_num'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_num_<CASE>.h5' '/pris_num' 'wtmean'  0.0  '/s_env_pris_num'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/pris_diam_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/ps_pris_diam'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/s_pris_diam'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/core_pris_diam_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/rb_pris_diam_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/env_pris_diam_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/ps_core_pris_diam'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/s_core_pris_diam'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/ps_rb_pris_diam'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/s_rb_pris_diam'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/ps_env_pris_diam'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_pris_diam_<CASE>.h5' '/pris_diam' 'wtmean'  0.0  '/s_env_pris_diam'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_pris_<CASE>.h5'
%    }
%
%    % aggregates
%    {
%      'Aggregates'
%      {
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/aggr_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/ps_aggr'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/s_aggr'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/core_aggr_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/rb_aggr_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/env_aggr_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/ps_core_aggr'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/s_core_aggr'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/ps_rb_aggr'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/s_rb_aggr'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/ps_env_aggr'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/s_env_aggr'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/aggr_num_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/ps_aggr_num'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/s_aggr_num'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/core_aggr_num_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/rb_aggr_num_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/env_aggr_num_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/ps_core_aggr_num'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/s_core_aggr_num'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/ps_rb_aggr_num'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/s_rb_aggr_num'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/ps_env_aggr_num'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_num_<CASE>.h5' '/aggr_num' 'wtmean'  0.0  '/s_env_aggr_num'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/aggr_diam_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/ps_aggr_diam'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/s_aggr_diam'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/core_aggr_diam_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/rb_aggr_diam_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/env_aggr_diam_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/ps_core_aggr_diam'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/s_core_aggr_diam'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/ps_rb_aggr_diam'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/s_rb_aggr_diam'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/ps_env_aggr_diam'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_aggr_diam_<CASE>.h5' '/aggr_diam' 'wtmean'  0.0  '/s_env_aggr_diam'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_aggr_<CASE>.h5'
%    }
%
%    % snow
%    {
%      'Snow'
%      {
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/snow_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/ps_snow'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/s_snow'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/core_snow_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/rb_snow_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/env_snow_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/ps_core_snow'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/s_core_snow'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/ps_rb_snow'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/s_rb_snow'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/ps_env_snow'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/s_env_snow'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/snow_num_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/ps_snow_num'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/s_snow_num'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/core_snow_num_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/rb_snow_num_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/env_snow_num_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/ps_core_snow_num'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/s_core_snow_num'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/ps_rb_snow_num'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/s_rb_snow_num'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/ps_env_snow_num'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_num_<CASE>.h5' '/snow_num' 'wtmean'  0.0  '/s_env_snow_num'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/snow_diam_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/ps_snow_diam'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/s_snow_diam'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/core_snow_diam_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/rb_snow_diam_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/env_snow_diam_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/ps_core_snow_diam'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/s_core_snow_diam'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/ps_rb_snow_diam'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/s_rb_snow_diam'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/ps_env_snow_diam'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_snow_diam_<CASE>.h5' '/snow_diam' 'wtmean'  0.0  '/s_env_snow_diam'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_snow_<CASE>.h5'
%    }
%
%    % graupel
%    {
%      'Graupel'
%      {
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/graup_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/ps_graup'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/s_graup'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/core_graup_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/rb_graup_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/env_graup_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/ps_core_graup'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/s_core_graup'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/ps_rb_graup'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/s_rb_graup'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/ps_env_graup'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/s_env_graup'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/graup_num_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/ps_graup_num'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/s_graup_num'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/core_graup_num_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/rb_graup_num_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/env_graup_num_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/ps_core_graup_num'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/s_core_graup_num'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/ps_rb_graup_num'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/s_rb_graup_num'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/ps_env_graup_num'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_num_<CASE>.h5' '/graup_num' 'wtmean'  0.0  '/s_env_graup_num'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/graup_diam_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/ps_graup_diam'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/s_graup_diam'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/core_graup_diam_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/rb_graup_diam_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/env_graup_diam_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/ps_core_graup_diam'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/s_core_graup_diam'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/ps_rb_graup_diam'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/s_rb_graup_diam'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/ps_env_graup_diam'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_graup_diam_<CASE>.h5' '/graup_diam' 'wtmean'  0.0  '/s_env_graup_diam'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_graup_<CASE>.h5'
%    }
%
%    % hail
%    {
%      'Hail'
%      {
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/hail_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/ps_hail'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/s_hail'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/core_hail_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/rb_hail_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/env_hail_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/ps_core_hail'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/s_core_hail'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/ps_rb_hail'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/s_rb_hail'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/ps_env_hail'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/s_env_hail'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/hail_num_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/ps_hail_num'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/s_hail_num'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/core_hail_num_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/rb_hail_num_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/env_hail_num_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/ps_core_hail_num'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/s_core_hail_num'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/ps_rb_hail_num'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/s_rb_hail_num'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/ps_env_hail_num'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_num_<CASE>.h5' '/hail_num' 'wtmean'  0.0  '/s_env_hail_num'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/hail_diam_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/ps_hail_diam'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/s_hail_diam'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/core_hail_diam_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/rb_hail_diam_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/env_hail_diam_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/ps_core_hail_diam'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/s_core_hail_diam'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/ps_rb_hail_diam'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/s_rb_hail_diam'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/ps_env_hail_diam'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_hail_diam_<CASE>.h5' '/hail_diam' 'wtmean'  0.0  '/s_env_hail_diam'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_hail_<CASE>.h5'
%    }
%
%    % total condensate
%    {
%      'Total Condensate'
%     {
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/tcond_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/ps_tcond'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/s_tcond'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/core_tcond_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/rb_tcond_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/env_tcond_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/ps_core_tcond'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/s_core_tcond'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/ps_rb_tcond'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/s_rb_tcond'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/ps_env_tcond'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/s_env_tcond'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_tcond_<CASE>.h5'
%    }
%
%    % cooling via latent heat of freezing
%    {
%      'LHF Cooling'
%     {
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/lhf_cool_ts'         'RbZT' ''      ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/ps_lhf_cool'         'RbZt'    ''      'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/s_lhf_cool'          'RbZt'    ''      'sal'     ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/core_lhf_cool_ts'    'rbZT' 'core'  ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/rb_lhf_cool_ts'      'rbZT' 'rband' ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/env_lhf_cool_ts'     'rbZT' 'env'   ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/ps_core_lhf_cool'    'rbZt'    'core'  'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/s_core_lhf_cool'     'rbZt'    'core'  'sal'     ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/ps_rb_lhf_cool'      'rbZt'    'rband' 'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/s_rb_lhf_cool'       'rbZt'    'rband' 'sal'     ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/ps_env_lhf_cool'     'rbZt'    'env'   'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/s_env_lhf_cool'      'rbZt'    'env'   'sal'     ''       'le' 0 }
%      }
%      'DIAGS/hist_meas_lhf_cool_<CASE>.h5'
%    }
%
%    % heating via latent heat of freezing
%    {
%      'LHF Heating'
%     {
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/lhf_heat_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/ps_lhf_heat'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/s_lhf_heat'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/core_lhf_heat_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/rb_lhf_heat_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/env_lhf_heat_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/ps_core_lhf_heat'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/s_core_lhf_heat'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/ps_rb_lhf_heat'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/s_rb_lhf_heat'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/ps_env_lhf_heat'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/s_env_lhf_heat'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_lhf_heat_<CASE>.h5'
%    }
%
%    % cooling via latent heat of vaporization
%    {
%      'LHV Cooling'
%     {
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/lhv_cool_ts'         'RbZT' ''      ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/ps_lhv_cool'         'RbZt'    ''      'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/s_lhv_cool'          'RbZt'    ''      'sal'     ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/core_lhv_cool_ts'    'rbZT' 'core'  ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/rb_lhv_cool_ts'      'rbZT' 'rband' ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/env_lhv_cool_ts'     'rbZT' 'env'   ''        ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/ps_core_lhv_cool'    'rbZt'    'core'  'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/s_core_lhv_cool'     'rbZt'    'core'  'sal'     ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/ps_rb_lhv_cool'      'rbZt'    'rband' 'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/s_rb_lhv_cool'       'rbZt'    'rband' 'sal'     ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/ps_env_lhv_cool'     'rbZt'    'env'   'pre_sal' ''       'le' 0 }
%        { 'AzAveragedData/hist_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/s_env_lhv_cool'      'rbZt'    'env'   'sal'     ''       'le' 0 }
%      }
%      'DIAGS/hist_meas_lhv_cool_<CASE>.h5'
%    }
%
%    % heating via latent heat of vaporization
%    {
%      'LHV Heating'
%     {
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/lhv_heat_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/ps_lhv_heat'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/s_lhv_heat'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/core_lhv_heat_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/rb_lhv_heat_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/env_lhv_heat_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/ps_core_lhv_heat'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/s_core_lhv_heat'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/ps_rb_lhv_heat'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/s_rb_lhv_heat'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/ps_env_lhv_heat'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/s_env_lhv_heat'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_lhv_heat_<CASE>.h5'
%    }
%
%    % liquid condensation
%    {
%      'Liquid Condensation'
%     {
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/liq_cond_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/ps_liq_cond'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/s_liq_cond'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/core_liq_cond_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/rb_liq_cond_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/env_liq_cond_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/ps_core_liq_cond'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/s_core_liq_cond'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/ps_rb_liq_cond'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/s_rb_liq_cond'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/ps_env_liq_cond'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/s_env_liq_cond'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_liq_cond_<CASE>.h5'
%    }
%
%    % liquid evaporation
%    {
%      'Liquid Evaporation'
%     {
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/liq_evap_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/ps_liq_evap'         'RbZt'    ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/s_liq_evap'          'RbZt'    ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/core_liq_evap_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/rb_liq_evap_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/env_liq_evap_ts'     'rbZT' 'env'   ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/ps_core_liq_evap'    'rbZt'    'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/s_core_liq_evap'     'rbZt'    'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/ps_rb_liq_evap'      'rbZt'    'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/s_rb_liq_evap'       'rbZt'    'rband' 'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/ps_env_liq_evap'     'rbZt'    'env'   'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/s_env_liq_evap'      'rbZt'    'env'   'sal'     ''       'le' -0.01 }
%      }
%      'DIAGS/hist_meas_liq_evap_<CASE>.h5'
%    }
%
%    % ice deposition
%    {
%      'Ice Deposition'
%     {
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/ice_dep_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/ps_ice_dep'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/s_ice_dep'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/core_ice_dep_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/rb_ice_dep_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/env_ice_dep_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/ps_core_ice_dep'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/s_core_ice_dep'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/ps_rb_ice_dep'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/s_rb_ice_dep'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/ps_env_ice_dep'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/s_env_ice_dep'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_ice_dep_<CASE>.h5'
%    }
%
%    % ice sublimation
%    {
%      'Ice Sublimation'
%     {
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/ice_sub_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/ps_ice_sub'         'RbZt'    ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/s_ice_sub'          'RbZt'    ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/core_ice_sub_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/rb_ice_sub_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/env_ice_sub_ts'     'rbZT' 'env'   ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/ps_core_ice_sub'    'rbZt'    'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/s_core_ice_sub'     'rbZt'    'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/ps_rb_ice_sub'      'rbZt'    'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/s_rb_ice_sub'       'rbZt'    'rband' 'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/ps_env_ice_sub'     'rbZt'    'env'   'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/s_env_ice_sub'      'rbZt'    'env'   'sal'     ''       'le' -0.01 }
%      }
%      'DIAGS/hist_meas_ice_sub_<CASE>.h5'
%    }
%
%    % ice melting
%    {
%      'Ice Melting'
%     {
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/ice_melt_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/ps_ice_melt'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/s_ice_melt'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/core_ice_melt_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/rb_ice_melt_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/env_ice_melt_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/ps_core_ice_melt'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/s_core_ice_melt'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/ps_rb_ice_melt'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/s_rb_ice_melt'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/ps_env_ice_melt'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/s_env_ice_melt'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_ice_melt_<CASE>.h5'
%    }
%
%    % riming, cloud to ice
%    {
%      'Cloud Riming'
%     {
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/cloud_rime_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/ps_cloud_rime'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/s_cloud_rime'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/core_cloud_rime_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/rb_cloud_rime_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/env_cloud_rime_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/ps_core_cloud_rime'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/s_core_cloud_rime'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/ps_rb_cloud_rime'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/s_rb_cloud_rime'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/ps_env_cloud_rime'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/s_env_cloud_rime'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_cloud_rime_<CASE>.h5'
%    }
%
%    % rain to ice
%    {
%      'Rain to Ice'
%     {
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/rain2ice_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/ps_rain2ice'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/s_rain2ice'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/core_rain2ice_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/rb_rain2ice_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/env_rain2ice_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/ps_core_rain2ice'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/s_core_rain2ice'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/ps_rb_rain2ice'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/s_rb_rain2ice'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/ps_env_rain2ice'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/s_env_rain2ice'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_rain2ice_<CASE>.h5'
%    }
%
%    % vapor mixing ratio
%    {
%      'Vapor'
%      {
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/vapor_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/ps_vapor'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/s_vapor'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/core_vapor_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/rb_vapor_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/env_vapor_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/ps_core_vapor'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/s_core_vapor'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/ps_rb_vapor'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/s_rb_vapor'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/ps_env_vapor'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/s_env_vapor'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_vapor_<CASE>.h5'
%    }
%
%    % theta measurements
%    {
%      'Theta-E'
%      {
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/theta_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/ps_theta'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/s_theta'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/core_theta_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/rb_theta_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/env_theta_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/ps_core_theta'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/s_core_theta'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/ps_rb_theta'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/s_rb_theta'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/ps_env_theta'     'rbZt'    'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/s_env_theta'      'rbZt'    'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_theta_<CASE>.h5'
%    }

    % tempc measurements
    {
      'Temperature Azavg'
      {
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/tempc_ts'         'RbZT'   ''      ''        ''       'ge' -100 }
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/ps_tempc'         'RbZt'   ''      'pre_sal' ''       'ge' -100 }
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/s_tempc'          'RbZt'   ''      'sal'     ''       'ge' -100 }
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/core_tempc_ts'    'rbZT'   'core'  ''        ''       'ge' -100 }
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/rb_tempc_ts'      'rbZT'   'rband' ''        ''       'ge' -100 }
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/env_tempc_ts'     'rbZT'   'env'   ''        ''       'ge' -100 }
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/ps_core_tempc'    'rbZt'   'core'  'pre_sal' ''       'ge' -100 }
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/s_core_tempc'     'rbZt'   'core'  'sal'     ''       'ge' -100 }
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/ps_rb_tempc'      'rbZt'   'rband' 'pre_sal' ''       'ge' -100 }
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/s_rb_tempc'       'rbZt'   'rband' 'sal'     ''       'ge' -100 }
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/ps_env_tempc'     'rbZt'   'env'   'pre_sal' ''       'ge' -100 }
        { 'AzAveragedData/hist_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/s_env_tempc'      'rbZt'   'env'   'sal'     ''       'ge' -100 }

        % lead region
        { 'AzAveragedData/hist_lead_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/l_tempc_ts'      'RbZT'   ''      ''        ''       'ge' -100 }
        { 'AzAveragedData/hist_lead_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/l_ps_tempc'      'RbZt'   ''      'pre_sal' ''       'ge' -100 }
        { 'AzAveragedData/hist_lead_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/l_s_tempc'       'RbZt'   ''      'sal'     ''       'ge' -100 }
        { 'AzAveragedData/hist_lead_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/lead_tempc_ts'   'rbZT'   'lead'  ''        ''       'ge' -100 }
        { 'AzAveragedData/hist_lead_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/ps_lead_tempc'   'rbZt'   'lead'  'pre_sal' ''       'ge' -100 }
        { 'AzAveragedData/hist_lead_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/s_lead_tempc'    'rbZt'   'lead'  'sal'     ''       'ge' -100 }
      }
      'DIAGS/hist_meas_az_tempc_<CASE>.h5'
    }

    {
      'Temperature Tsavg'
      {
        % sample regions in SAL
        { 'TsAveragedData/hist_spath_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/spath_tempc_ts'  'b_ZT'   ''      ''        ''       'ge' -100 }
        { 'TsAveragedData/hist_spath_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/spath_ps_tempc'  'b_Zt'   ''      'pre_sal' ''       'ge' -100 }
        { 'TsAveragedData/hist_spath_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/spath_s_tempc'   'b_Zt'   ''      'sal'     ''       'ge' -100 }

        { 'TsAveragedData/hist_smaxcp_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/smaxcp_tempc_ts'  'b_ZT'   ''      ''        ''       'ge' -100 }
        { 'TsAveragedData/hist_smaxcp_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/smaxcp_ps_tempc'  'b_Zt'   ''      'pre_sal' ''       'ge' -100 }
        { 'TsAveragedData/hist_smaxcp_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/smaxcp_s_tempc'   'b_Zt'   ''      'sal'     ''       'ge' -100 }
      }
      'DIAGS/hist_meas_ts_tempc_<CASE>.h5'
    }

    {
      'Cold Pools Tsavg'
      {
        % lead region
        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/lead_cpools_ts'      'b__T'   ''      ''        ''       'ge' -100 }
        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/lead_ps_cpools_ts'   'b__T'   ''      'pre_sal' ''       'ge' -100 }
        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/lead_s_cpools_ts'    'b__T'   ''      'sal'     ''       'ge' -100 }
        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_lead_cpools'    'B__t'   ''      ''        ''       'ge' -100 }
        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_lead_ps_cpools' 'B__t'   ''      'pre_sal' ''       'ge' -100 }
        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_lead_s_cpools'  'B__t'   ''      'sal'     ''       'ge' -100 }

        % sample regions in SAL
        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/spath_cpools_ts'      'b__T'   ''      ''        ''       'ge' -100 }
        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/spath_ps_cpools_ts'   'b__T'   ''      'pre_sal' ''       'ge' -100 }
        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/spath_s_cpools_ts'    'b__T'   ''      'sal'     ''       'ge' -100 }
        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_spath_cpools'    'B__t'   ''      ''        ''       'ge' -100 }
        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_spath_ps_cpools' 'B__t'   ''      'pre_sal' ''       'ge' -100 }
        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_spath_s_cpools'  'B__t'   ''      'sal'     ''       'ge' -100 }

        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/smaxcp_cpools_ts'      'b__T'   ''      ''        ''       'ge' -100 }
        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/smaxcp_ps_cpools_ts'   'b__T'   ''      'pre_sal' ''       'ge' -100 }
        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/smaxcp_s_cpools_ts'    'b__T'   ''      'sal'     ''       'ge' -100 }
        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_smaxcp_cpools'    'B__t'   ''      ''        ''       'ge' -100 }
        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_smaxcp_ps_cpools' 'B__t'   ''      'pre_sal' ''       'ge' -100 }
        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_smaxcp_s_cpools'  'B__t'   ''      'sal'     ''       'ge' -100 }

      }
      'DIAGS/hist_meas_ts_cpools_<CASE>.h5'
    }

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
        Rspec     = num2cell(MeasList{imeas}{6});
        Rrange    = MeasList{imeas}{7};
        Trange    = MeasList{imeas}{8};
        Zrange    = MeasList{imeas}{9};
        SelectOp  = MeasList{imeas}{10};
        SelectVal = MeasList{imeas}{11};

        InFile = regexprep(Ftemplate, '<CASE>', Case);

        % Parse the reduction spec. Always has 4 characters that correspond to x, y, z, t.
        Rindex = 0;  % if index remains zero, then corresponding quantity is not in the input
        Bindex = 0;
        Zindex = 0;
        Tindex = 0;

        Rreduce = false;
        Breduce = false;
        Zreduce = false;
        Treduce = false;

        InForm = {};
        Ndims = 0;
        InDims = { 'x' 'y' 'z' 't' }';
        for i = 1:length(Rspec)
          switch(Rspec{i})
            case { 'r' 'R' }
              RinVar = InDims{i};
              Ndims = Ndims + 1;
              Rindex = Ndims;
              InForm{Ndims} = 'r';
              if (strcmp(Rspec{i}, 'r'))
                Rreduce = true;
              end

            case { 'b' 'B' }
              BinVar = InDims{i};
              Ndims = Ndims + 1;
              Bindex = Ndims;
              InForm{Ndims} = 'b';
              if (strcmp(Rspec{i}, 'b'))
                Breduce = true;
              end

            case { 'z' 'Z' }
              ZinVar = InDims{i};
              Ndims = Ndims + 1;
              Zindex = Ndims;
              InForm{Ndims} = 'z';
              if (strcmp(Rspec{i}, 'z'))
                Zreduce = true;
              end

            case { 't' 'T' }
              TinVar = InDims{i};
              Ndims = Ndims + 1;
              Tindex = Ndims;
              InForm{Ndims} = 't';
              if (strcmp(Rspec{i}, 't'))
                Treduce = true;
              end
          end
        end

        fprintf('      Reading: %s (%s)\n', InFile, Vname);
        fprintf('        Input form: (');
        fprintf('%s', InForm{1});
        for i = 2:Ndims
          fprintf(',%s', InForm{i});
        end
        fprintf(')\n');
        fprintf('\n');

        HDATA = squeeze(h5read(InFile, Vname));
        X     = squeeze(h5read(InFile, '/x_coords'));
        Y     = squeeze(h5read(InFile, '/y_coords'));
        Z     = squeeze(h5read(InFile, '/z_coords'));
        T     = squeeze(h5read(InFile, '/t_coords'));
  
        % Determine indices corresponding to r,b,z,t range specs
        % Default is entire range of dimension sizes
        clear R1;
        clear R2;
        clear B1;
        clear B2;
        clear Z1
        clear Z2
        clear T1;
        clear T2;

        fprintf('        Selection:\n');
        % RADIUS
        if (Rindex > 0)
          % convert to km
          switch(RinVar)
            case 'x'
              R = X ./ 1000;
            case 'y'
              R = Y ./ 1000;
            case 'z'
              R = Z ./ 1000;
            case 't'
              R = T ./ 1000;
          end
          R1 = 1;
          R2 = length(R);
          switch(Rrange)
            case 'core'
              R1 = find(R >= RstartCore, 1, 'first');
              R2 = find(R <= RendCore,   1, 'last');

            case 'rband'
              R1 = find(R >= RstartRband, 1, 'first');
              R2 = find(R <= RendRband,   1, 'last');
  
            case 'env'
              R1 = find(R >= RstartEnv, 1, 'first');
              R2 = find(R <= RendEnv,   1, 'last');
          end

          fprintf('          Radius: "%s" -> %.2f to %.2f (%d:%d)\n', Rrange, R(R1), R(R2), R1, R2);

          % trim the radius coordinates
          R  = R(R1:R2);
          Nr = length(R);

          % add to selection spec
          if (Rindex == 1)
            SelectSpec = sprintf('%d:%d', R1, R2);
          else
            SelectSpec = sprintf('%s,%d:%d', SelectSpec, R1, R2);
          end
        end

        % BINS
        if (Bindex > 0)
          switch(BinVar)
            case 'x'
              B = X;
            case 'y'
              B = Y;
            case 'z'
              B = Z;
            case 't'
              B = T;
          end
          B1 = 1;         
          B2 = length(B);
          switch(SelectOp)
            case 'ge'
              B1 = find(B >= SelectVal, 1, 'first');
  
            case 'le'
              B2 = find(B <= SelectVal, 1, 'last');
          end

          fprintf('          Bin: "%s %.2f" -> %.2f to %.2f (%d:%d)\n', SelectOp, SelectVal, B(B1), B(B2), B1, B2);

          % trim the bins coordinates
          B  = B(B1:B2);
          Nb = length(B);

          % add to selection spec
          if (Bindex == 1)
            SelectSpec = sprintf('%d:%d', B1, B2);
          else
            SelectSpec = sprintf('%s,%d:%d', SelectSpec, B1, B2);
          end
        end

        % HEIGHT
        if (Zindex > 0)
          % convert to km
          switch(ZinVar)
            case 'x'
              H = X ./ 1000;
            case 'y'
              H = Y ./ 1000;
            case 'z'
              H = Z ./ 1000;
            case 't'
              H = T ./ 1000;
          end
          Z1 = 1;
          Z2 = length(H);
          switch(Zrange)
            case { 'maxlev' 'minlev' 'sfc' }
              Z1 = find(H >= Zsfc, 1, 'first');
              Z2 = find(H <= Ztop, 1, 'last');
          end

          fprintf('          Height: "%s" -> %.2f to %.2f (%d:%d)\n', Zrange, H(Z1), H(Z2), Z1, Z2);

          % trim the height coordinates
          H  = H(Z1:Z2);
          Nz = length(H);

          % add to selection spec
          if (Zindex == 1)
            SelectSpec = sprintf('%d:%d', Z1, Z2);
          else
            SelectSpec = sprintf('%s,%d:%d', SelectSpec, Z1, Z2);
          end
        end

        % SIM TIME
        if (Tindex > 0)
          % convert to hours starting with zero
          switch(TinVar)
            case 'x'
              ST = (X ./ 3600) - 42;
            case 'y'
              ST = (Y ./ 3600) - 42;
            case 'z'
              ST = (Z ./ 3600) - 42;
            case 't'
              ST = (T ./ 3600) - 42;
          end
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

          fprintf('          Time: "%s" -> %.2f to %.2f (%d:%d)\n', Trange, ST(T1), ST(T2), T1, T2);

          % trim the sim time coordinates
          ST  = ST(T1:T2);
          Nt = length(ST);

          % add to selection spec
          if (Tindex == 1)
            SelectSpec = sprintf('%d:%d', T1, T2);
          else
            SelectSpec = sprintf('%s,%d:%d', SelectSpec, T1, T2);
          end
        end

        % Trim down the input data according to the selection indices
        SelectCmd = sprintf('HDATA(%s)', SelectSpec);
        MEAS = eval(SelectCmd);
        fprintf('\n');
  
        % If first measurement, then write out coordinates for later
        % use in attaching vars to them. Write out the original coordinates.
        % The original coordinates with their original lengths will always be the
        % appropriate vectors to use for attaching to variables.
        if (imeas == 1)
          Xname = '/x_coords';
          Yname = '/y_coords';
          Zname = '/z_coords';
          Tname = '/t_coords';
  
          CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
          % Add COARDS annotations
          NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
        end

        % Reduce the input data. Go in this order:
        %    1. Sum histogram counts (dims r and t)
        %    2. Reduce histogram counts (dim b)
        %    3. Select z level (dim z)
        %
        % Keep original dimensions intact until finished so that [RBZT]index vars
        % can remain constant.
        fprintf('        Reduction:\n');

        if (Rreduce)
          fprintf('          Summing radial counts\n');
          MEAS = squeeze(sum(MEAS, Rindex));

          % Adjust indices, number of dimensions
          if (Bindex > Rindex)
            Bindex = Bindex - 1;
          end
          if (Zindex > Rindex)
            Zindex = Zindex - 1;
          end
          if (Tindex > Rindex)
            Tindex = Tindex - 1;
          end
          Rindex = 0;
          Ndims = Ndims - 1;
        end

        if (Treduce)
          fprintf('          Summing temporal counts\n');
          MEAS = squeeze(sum(MEAS, Tindex));

          % Adjust indices, number of dimensions
          if (Rindex > Tindex)
            Rindex = Rindex - 1;
          end
          if (Bindex > Tindex)
            Bindex = Bindex - 1;
          end
          if (Zindex > Tindex)
            Zindex = Zindex - 1;
          end
          Tindex = 0;
          Ndims = Ndims - 1;
        end

        if (Breduce)
          fprintf('          Reducing bin counts\n');
          MEAS = squeeze(ReduceHists(MEAS, Bindex, B, Rmethod, Param));

          % Adjust indices, number of dimensions
          if (Rindex > Bindex)
            Rindex = Rindex - 1;
          end
          if (Zindex > Bindex)
            Zindex = Zindex - 1;
          end
          if (Tindex > Bindex)
            Tindex = Tindex - 1;
          end
          Bindex = 0;
          Ndims = Ndims - 1;
        end

        if (Zreduce)
          fprintf('          Reducing height\n');
          % Z has been trimmed down to the range we are restricting the search
          % for min or max values, therefore z == 1 represents the surface.
          % Determine which level needs to be selected.
          switch(Zrange)
            case 'sfc'
              ZSEL = 1;

            case 'minlev'
              % 1. grab linear index into MEAS where min occurs
              % 2. translate linear index back to 4D indices, (r,b,z,t)
              %    doesn't matter if MEAS is less than 4D, the extra dims get 1's
              % 3. pick off the index that corresponds to the z dimension
              [ MVAL MLOC ] = min(MEAS(:));
              [ I1 I2 I3 I4 ] = ind2sub(size(MEAS),MLOC);
              MIND = [ I1 I2 I3 I4 ];
              ZSEL = MIND(Zindex);

            case 'maxlev'
              % 1. grab linear index into MEAS where max occurs
              % 2. translate linear index back to 4D indices, (r,b,z,t)
              %    doesn't matter if MEAS is less than 4D, the extra dims get 1's
              % 3. pick off the index that corresponds to the z dimension
              [ MVAL MLOC ] = max(MEAS(:));
              [ I1 I2 I3 I4 ] = ind2sub(size(MEAS),MLOC);
              MIND = [ I1 I2 I3 I4 ];
              ZSEL = MIND(Zindex);
          end
          fprintf('            Z level selected for reduction %s: %.2f (%d)\n', Zrange, H(ZSEL), ZSEL);

          % select the single z level, thus reducing the z dimension
          for i = 1:Ndims
            if (i == Zindex)
              DimSpecs{i} = sprintf('%d', ZSEL);
            else
              DimSpecs{i} = ':';
            end
          end
          SelectCmd = sprintf('squeeze(MEAS(%s))', strjoin(DimSpecs, ','));
          MEAS = eval(SelectCmd);

          % Adjust indices, number of dimensions
          if (Rindex > Zindex)
            Rindex = Rindex - 1;
          end
          if (Bindex > Zindex)
            Bindex = Bindex - 1;
          end
          if (Tindex > Zindex)
            Tindex = Tindex - 1;
          end
          Zindex = 0;
          Ndims = Ndims - 1;
        end
        fprintf('\n');

        % Write out measurement
        fprintf('      Writing: %s (%s)\n', OutFile, OutVname)
        fprintf('\n');
  
        % Figure out output dimensions, sizes
        clear OutSize;
        clear DimOrder;

        Oindex = 0;

        for i = 1:length(Rspec)
          switch(Rspec{i})
            case('R')
              Oindex = Oindex + 1;
              OutSize(Oindex) = Nr;
              DimOrder{Oindex} = RinVar;

            case('B')
              Oindex = Oindex + 1;
              OutSize(Oindex) = Nb;
              DimOrder{Oindex} = BinVar;

            case('Z')
              Oindex = Oindex + 1;
              OutSize(Oindex) = Nz;
              DimOrder{Oindex} = ZinVar;

            case('T')
              Oindex = Oindex + 1;
              OutSize(Oindex) = Nt;
              DimOrder{Oindex} = TinVar;

          end
        end

        h5create(OutFile, OutVname, OutSize);
        h5write(OutFile, OutVname, MEAS);

        % Attach dimensions
        AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
      end % measurements
    end % sets
  end % cases
end % function
