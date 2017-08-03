function [ ] = GenTsdHistMeas()

  % make sure output directory exists
  Ddir = 'DIAGS';  % coordinate this with output file names below
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % list of simulation cases
  CaseList = {
     'TSD_SD_1KM'
%     'TSD_SD'
%     'TSD_SD_1G'
%    'TSD_SAL_DUST'
%    'TSD_SAL_NODUST'
%    'TSD_NONSAL_DUST'
%    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  % Temporal ranges
  TstartInit = 0;
  TendInit = 2;

  TstartMid = 34;
  TendMid = 36;
  
  TstartFinal = 58;
  TendFinal = 60;


  TstartPreSal = 10;
  TendPreSal   = 30;

  TstartSal = 40;
  TendSal   = 60;


  % Radial ranges
  RstartCore =   0;
  RendCore   = 200;

  RstartRband = 250;
  RendRband   = 450;

  RstartEnv = 450;
  RendEnv   = 600;

  % Height ranges
  Zsfc = 0;
  Ztop = 2; % km

  Psfc = 1010; % mb
  Ptop = 800;

  % For converting input time values.
  TimeScale = 1 / 3600;  % seconds to hours

  %TimeOffset = 42;      % time (h) at start of sim
  %TimeOffset = 12;
  TimeOffset = 0;

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
  %       init           TstartInit,   TendInit
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Azimuthally averaged data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%% Vertical coords -> pressure %%%%%%

    %%%%%%% lite data %%%%%%%

%    % theta
%    {
%      'Theta Azavg'
%      {
%        { 'AzAveragedData/hist_all_p_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_theta_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_ps_theta'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_s_theta'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_core_theta_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_core_ps_theta'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_core_s_theta'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_rb_theta_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_rb_ps_theta'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_rb_s_theta'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_p_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_theta_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_ps_theta'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_s_theta'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_core_theta_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_core_ps_theta'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_core_s_theta'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_rb_theta_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_rb_ps_theta'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_rb_s_theta'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_p_theta_lite_<CASE>.h5'
%    }
%
%    % theta_e
%    {
%      'Theta-E Azavg'
%      {
%        { 'AzAveragedData/hist_all_p_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_theta_e_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_ps_theta_e'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_s_theta_e'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_theta_e_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_ps_theta_e'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_s_theta_e'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_theta_e_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_ps_theta_e'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_s_theta_e'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_p_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_theta_es_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_ps_theta_es'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_s_theta_es'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_core_theta_es_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_core_ps_theta_es'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_core_s_theta_es'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_rb_theta_es_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_rb_ps_theta_es'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_rb_s_theta_es'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_p_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_theta_e_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_ps_theta_e'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_s_theta_e'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_core_theta_e_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_core_ps_theta_e'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_core_s_theta_e'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_rb_theta_e_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_rb_ps_theta_e'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_rb_s_theta_e'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_p_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_theta_es_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_ps_theta_es'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_s_theta_es'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_core_theta_es_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_core_ps_theta_es'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_core_s_theta_es'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_rb_theta_es_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_rb_ps_theta_es'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_rb_s_theta_es'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_p_theta_e_lite_<CASE>.h5'
%    }
%
%    % entropy
%    {
%      'Entropy Azavg'
%      {
%        { 'AzAveragedData/hist_all_p_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_entropy_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_ps_entropy'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_s_entropy'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_core_entropy_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_core_ps_entropy'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_core_s_entropy'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_rb_entropy_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_rb_ps_entropy'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_rb_s_entropy'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_p_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_entropy_s_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_ps_entropy_s'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_s_entropy_s'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_core_entropy_s_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_core_ps_entropy_s'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_core_s_entropy_s'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_rb_entropy_s_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_rb_ps_entropy_s'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_rb_s_entropy_s'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_p_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_entropy_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_ps_entropy'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_s_entropy'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_core_entropy_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_core_ps_entropy'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_core_s_entropy'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_rb_entropy_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_rb_ps_entropy'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_rb_s_entropy'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_p_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_entropy_s_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_ps_entropy_s'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_s_entropy_s'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_core_entropy_s_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_core_ps_entropy_s'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_core_s_entropy_s'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_rb_entropy_s_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_rb_ps_entropy_s'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_rb_s_entropy_s'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_p_entropy_lite_<CASE>.h5'
%    }
%
%    %%%%%%% full resolution data %%%%%%%
%
%    % storm speed measurements
%    {
%      'Speed'
%      {
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_speed_t_ts'           'RbZT' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t_ts'        'RbZT' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t_ts'         'RbZT' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t'           'RbZt' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t'            'RbZt' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_whole_speed_t'        'rbZT' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_whole_ps_speed_t'     'rbZt' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_whole_s_speed_t'      'rbZt' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_core_ps_speed_t'      'rbZt' 'core'  'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_core_s_speed_t'       'rbZt' 'core'  'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_rb_ps_speed_t'        'rbZt' 'rband' 'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_rb_s_speed_t'         'rbZt' 'rband' 'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t_sfc'       'Rbzt' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t_sfc'        'Rbzt' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t_maxlev'    'Rbzt' ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t_maxlev'     'Rbzt' ''      'sal'     'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_speed_t_sfc_ts'       'RbzT' ''      ''        'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_speed_t_maxlev_ts'    'RbzT' ''      ''        'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t_sfc_ts'    'RbzT' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t_maxlev_ts' 'RbzT' ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t_sfc_ts'     'RbzT' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t_maxlev_ts'  'RbzT' ''      'sal'     'maxlev' 'ge' -100 }
%
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_speed_r_ts'           'RbZT' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r_ts'        'RbZT' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r_ts'         'RbZT' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r'           'RbZt' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r'            'RbZt' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_whole_speed_r'        'rbZT' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_whole_ps_speed_r'     'rbZt' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_whole_s_speed_r'      'rbZt' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_core_ps_speed_r'      'rbZt' 'core'  'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_core_s_speed_r'       'rbZt' 'core'  'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_rb_ps_speed_r'        'rbZt' 'rband' 'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_rb_s_speed_r'         'rbZt' 'rband' 'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r_sfc'       'Rbzt' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r_sfc'        'Rbzt' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r_minlev'    'Rbzt' ''      'pre_sal' 'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r_minlev'     'Rbzt' ''      'sal'     'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_speed_r_sfc_ts'       'RbzT' ''      ''        'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_speed_r_minlev_ts'    'RbzT' ''      ''        'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r_sfc_ts'    'RbzT' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r_maxlev_ts' 'RbzT' ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r_sfc_ts'     'RbzT' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_p_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r_maxlev_ts'  'RbzT' ''      'sal'     'maxlev' 'ge' -100 }
%
%      }
%      'DIAGS/hist_meas_az_p_speed_<CASE>.h5'
%    }
%
%    % vertical velocity measurements
%    {
%      'Vertical Velocity'
%      {
%        { 'AzAveragedData/hist_all_p_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_updraft_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_ps_updraft'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_s_updraft'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_updraft_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_ps_updraft'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_s_updraft'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_updraft_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_ps_updraft'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_s_updraft'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_p_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_dndraft_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_p_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_ps_dndraft'         'RbZt' ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_p_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_s_dndraft'          'RbZt' ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_p_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_dndraft_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_p_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_ps_dndraft'    'rbZt' 'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_p_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_s_dndraft'     'rbZt' 'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_p_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_dndraft_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_p_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_ps_dndraft'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_p_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_s_dndraft'       'rbZt' 'rband' 'sal'     ''       'le' -0.01 }
%
%      }
%      'DIAGS/hist_meas_az_p_w_<CASE>.h5'
%    }
%
%    % vapor mixing ratio
%    {
%      'Vapor Azavg'
%      {
%        { 'AzAveragedData/hist_all_p_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_vapor_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_ps_vapor'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_s_vapor'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_core_vapor_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_rb_vapor_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_core_ps_vapor'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_core_s_vapor'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_rb_ps_vapor'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_rb_s_vapor'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%
%      }
%      'DIAGS/hist_meas_az_p_vapor_<CASE>.h5'
%    }
%
%    % cloud
%    {
%      'Cloud'
%      {
%        { 'AzAveragedData/hist_all_p_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_cloud_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_ps_cloud_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_s_cloud_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_core_cloud_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_core_ps_cloud_mass'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_core_s_cloud_mass'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_rb_cloud_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_rb_ps_cloud_mass'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_rb_s_cloud_mass'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%      }
%      'DIAGS/hist_meas_az_p_cloud_<CASE>.h5'
%    }
%
%    % rain
%    {
%      'Rain'
%      {
%        { 'AzAveragedData/hist_all_p_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_rain_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_ps_rain_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_s_rain_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_core_rain_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_rb_rain_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_core_ps_rain_mass'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_core_s_rain_mass'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_rb_ps_rain_mass'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_p_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_rb_s_rain_mass'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%      }
%      'DIAGS/hist_meas_az_p_rain_<CASE>.h5'
%    }
%
%    % precip rate measurements
%    {
%      'Precip Rate Azavg'
%      {
%        { 'AzAveragedData/hist_all_p_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_pcprate_ts'      'Rb_T'  ''      ''        ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_all_p_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_ps_pcprate_ts'   'Rb_T'  ''      'pre_sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_all_p_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_s_pcprate_ts'    'Rb_T'  ''      'sal'     ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_all_p_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_ps_pcprate'      'Rb_t'  ''      'pre_sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_all_p_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_s_pcprate'       'Rb_t'  ''      'sal'     ''       'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_p_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_all_core_ps_pcprate' 'rB_t'  'core' 'pre_sal' ''   'ge' 0.001 }
%        { 'AzAveragedData/hist_all_p_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_all_core_s_pcprate'  'rB_t'  'core' 'sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_all_p_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_all_rb_ps_pcprate'   'rB_t'  'rb'   'pre_sal' ''   'ge' 0.001 }
%        { 'AzAveragedData/hist_all_p_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_all_rb_s_pcprate'    'rB_t'  'rb'   'sal' ''       'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_p_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_whole_pcprate_ts'      'rb_T'  ''      ''        ''       'ge' 0.001 }
%      }
%      'DIAGS/hist_meas_az_p_pcprate_<CASE>.h5'
%    }
%
%
%    % theta_e measurements
%    {
%      'Theta-E Azavg'
%      {
%        { 'AzAveragedData/hist_all_p_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_theta_e_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_ps_theta_e'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_s_theta_e'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_theta_e_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_ps_theta_e'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_s_theta_e'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_theta_e_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_ps_theta_e'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_s_theta_e'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%      }
%      'DIAGS/hist_meas_az_p_theta_e_<CASE>.h5'
%    }
%
%    % tempc measurements
%    {
%      'Temperature Azavg'
%      {
%        { 'AzAveragedData/hist_all_p_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_tempc_ts'         'RbZT'   ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_ps_tempc'         'RbZt'   ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_s_tempc'          'RbZt'   ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_core_tempc_ts'    'rbZT'   'core'  ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_rb_tempc_ts'      'rbZT'   'rband' ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_core_ps_tempc'    'rbZt'   'core'  'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_core_s_tempc'     'rbZt'   'core'  'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_rb_ps_tempc'      'rbZt'   'rband' 'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_p_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_rb_s_tempc'       'rbZt'   'rband' 'sal'     ''       'ge' -100 }
%      }
%      'DIAGS/hist_meas_az_p_tempc_<CASE>.h5'
%    }
%
%    % relhum measurements
%    {
%      'RH Azavg'
%      {
%        { 'AzAveragedData/hist_all_p_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_relhum_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_ps_relhum'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_s_relhum'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_core_relhum_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_rb_relhum_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_core_ps_relhum'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_core_s_relhum'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_rb_ps_relhum'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_p_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_rb_s_relhum'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_p_relhum_<CASE>.h5'
%    }
%
%    % liquid evaporation
%    {
%      'Liquid Evaporation'
%     {
%        { 'AzAveragedData/hist_all_p_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_liq_evap_ts'         'RbZT' ''      ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_ps_liq_evap'         'RbZt' ''      'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_s_liq_evap'          'RbZt' ''      'sal'     ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_core_liq_evap_ts'    'rbZT' 'core'  ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_core_ps_liq_evap'    'rbZt' 'core'  'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_core_s_liq_evap'     'rbZt' 'core'  'sal'     ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_rb_liq_evap_ts'      'rbZT' 'rband' ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_rb_ps_liq_evap'      'rbZt' 'rband' 'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_rb_s_liq_evap'       'rbZt' 'rband' 'sal'     ''       'lt' 0 }
%
%      }
%      'DIAGS/hist_meas_az_p_liq_evap_<CASE>.h5'
%    }
%
%    % Liquid condensation
%    {
%      'Liquid Condensation'
%     {
%        { 'AzAveragedData/hist_all_p_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_liq_cond_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_ps_liq_cond'         'RbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_s_liq_cond'          'RbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_core_liq_cond_ts'    'rbZT' 'core'  ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_core_ps_liq_cond'    'rbZt' 'core'  'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_core_s_liq_cond'     'rbZt' 'core'  'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_rb_liq_cond_ts'      'rbZT' 'rband' ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_rb_ps_liq_cond'      'rbZt' 'rband' 'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_rb_s_liq_cond'       'rbZt' 'rband' 'sal'     ''       'gt' 0 }
%
%      }
%      'DIAGS/hist_meas_az_p_liq_cond_<CASE>.h5'
%    }
%
%    % ice deposition
%    {
%      'Ice Deposition'
%     {
%        { 'AzAveragedData/hist_all_p_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_ice_dep_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_ps_ice_dep'         'RbZt'    ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_s_ice_dep'          'RbZt'    ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_core_ice_dep_ts'    'rbZT' 'core'  ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_rb_ice_dep_ts'      'rbZT' 'rband' ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_core_ps_ice_dep'    'rbZt'    'core'  'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_core_s_ice_dep'     'rbZt'    'core'  'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_rb_ps_ice_dep'      'rbZt'    'rband' 'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_rb_s_ice_dep'       'rbZt'    'rband' 'sal'     ''       'gt' 0 }
%
%      }
%      'DIAGS/hist_meas_az_p_ice_dep_<CASE>.h5'
%    }
%
%    % ice sublimation
%    {
%      'Ice Sublimation'
%     {
%        { 'AzAveragedData/hist_all_p_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_ice_sub_ts'         'RbZT' ''      ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_ps_ice_sub'         'RbZt'    ''      'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_s_ice_sub'          'RbZt'    ''      'sal'     ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_core_ice_sub_ts'    'rbZT' 'core'  ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_rb_ice_sub_ts'      'rbZT' 'rband' ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_core_ps_ice_sub'    'rbZt'    'core'  'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_core_s_ice_sub'     'rbZt'    'core'  'sal'     ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_rb_ps_ice_sub'      'rbZt'    'rband' 'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_p_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_rb_s_ice_sub'       'rbZt'    'rband' 'sal'     ''       'lt' 0 }
%
%      }
%      'DIAGS/hist_meas_az_p_ice_sub_<CASE>.h5'
%    }
%
%
%
%
%
%    %%%%%% Vertical coords -> height %%%%%%
%
%    %%%%%% lite data %%%%%%%%
%
%    % theta
%    {
%      'Theta Azavg'
%      {
%        { 'AzAveragedData/hist_all_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_theta_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_ps_theta'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_s_theta'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_core_theta_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_core_ps_theta'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_core_s_theta'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_rb_theta_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_rb_ps_theta'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_rb_s_theta'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_theta_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_ps_theta'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_s_theta'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_core_theta_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_core_ps_theta'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_core_s_theta'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_rb_theta_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_rb_ps_theta'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_nv_lite_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_nv_rb_s_theta'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_theta_lite_<CASE>.h5'
%    }
%
%    % theta_e
%    {
%      'Theta-E Azavg'
%      {
%        { 'AzAveragedData/hist_all_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_theta_e_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_ps_theta_e'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_s_theta_e'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_theta_e_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_ps_theta_e'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_s_theta_e'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_theta_e_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_ps_theta_e'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_s_theta_e'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_theta_es_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_ps_theta_es'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_s_theta_es'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_core_theta_es_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_core_ps_theta_es'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_core_s_theta_es'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_rb_theta_es_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_rb_ps_theta_es'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_rb_s_theta_es'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_theta_e_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_ps_theta_e'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_s_theta_e'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_core_theta_e_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_core_ps_theta_e'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_core_s_theta_e'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_rb_theta_e_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_rb_ps_theta_e'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_nv_lite_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_nv_rb_s_theta_e'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_theta_es_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_ps_theta_es'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_s_theta_es'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_core_theta_es_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_core_ps_theta_es'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_core_s_theta_es'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_rb_theta_es_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_rb_ps_theta_es'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_es_nv_lite_<CASE>.h5' '/theta_es' 'wtmean'  0.0  '/all_nv_rb_s_theta_es'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_theta_e_lite_<CASE>.h5'
%    }
%
%    % entropy
%    {
%      'Entropy Azavg'
%      {
%        { 'AzAveragedData/hist_all_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_entropy_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_ps_entropy'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_s_entropy'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_core_entropy_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_core_ps_entropy'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_core_s_entropy'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_rb_entropy_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_rb_ps_entropy'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_rb_s_entropy'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_entropy_s_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_ps_entropy_s'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_s_entropy_s'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_core_entropy_s_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_core_ps_entropy_s'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_core_s_entropy_s'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_rb_entropy_s_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_rb_ps_entropy_s'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_rb_s_entropy_s'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_entropy_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_ps_entropy'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_s_entropy'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_core_entropy_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_core_ps_entropy'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_core_s_entropy'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_rb_entropy_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_rb_ps_entropy'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_nv_lite_<CASE>.h5' '/entropy' 'wtmean'  0.0  '/all_nv_rb_s_entropy'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_all_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_entropy_s_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_ps_entropy_s'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_s_entropy_s'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_core_entropy_s_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_core_ps_entropy_s'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_core_s_entropy_s'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_rb_entropy_s_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_rb_ps_entropy_s'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_entropy_s_nv_lite_<CASE>.h5' '/entropy_s' 'wtmean'  0.0  '/all_nv_rb_s_entropy_s'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_entropy_lite_<CASE>.h5'
%    }
%
%    %%%%%% full resolution data %%%%%%%%
%    % vertical velocity measurements
%    {
%      'Vertical Velocity'
%      {
%        { 'AzAveragedData/hist_all_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_updraft_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_ps_updraft'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_s_updraft'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_updraft_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_ps_updraft'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_s_updraft'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_updraft_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_ps_updraft'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_s_updraft'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_dndraft_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_ps_dndraft'         'RbZt' ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_s_dndraft'          'RbZt' ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_dndraft_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_ps_dndraft'    'rbZt' 'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_core_s_dndraft'     'rbZt' 'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_dndraft_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_ps_dndraft'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/all_rb_s_dndraft'       'rbZt' 'rband' 'sal'     ''       'le' -0.01 }
%
%        { 'AzAveragedData/hist_lead_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_updraft_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_ps_updraft'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_s_updraft'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_core_updraft_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_core_ps_updraft'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_core_s_updraft'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_rb_updraft_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_rb_ps_updraft'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_up_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_rb_s_updraft'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_dndraft_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_ps_dndraft'         'RbZt' ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_s_dndraft'          'RbZt' ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_core_dndraft_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_core_ps_dndraft'    'rbZt' 'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_core_s_dndraft'     'rbZt' 'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_rb_dndraft_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_rb_ps_dndraft'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_dn_<CASE>.h5' '/w' 'wtmean'  0.0  '/lead_rb_s_dndraft'       'rbZt' 'rband' 'sal'     ''       'le' -0.01 }
%      }
%      'DIAGS/hist_meas_az_w_<CASE>.h5'
%    }

    % dust measurements
    {
      'Dust Azavg'
      {
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_d1_num_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_ps_d1_num'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_s_d1_num'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_whole_d1_num_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_whole_ps_d1_num'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_whole_s_d1_num'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_whole_i_d1_num'    'rbZt' ''      'init'    ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_whole_m_d1_num'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_whole_f_d1_num'    'rbZt' ''      'final'   ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_core_ps_d1_num'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_core_s_d1_num'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_rb_ps_d1_num'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_num_<CASE>.h5' '/d1_num' 'wtmean'  0.0  '/all_rb_s_d1_num'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/all_d1_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/all_ps_d1_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/all_s_d1_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/all_whole_d1_mass_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/all_whole_ps_d1_mass'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/all_whole_s_d1_mass'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/all_whole_i_d1_mass'    'rbZt' ''      'init'    ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/all_whole_m_d1_mass'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/all_whole_f_d1_mass'    'rbZt' ''      'final'   ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/all_trdust1_diff_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/all_ps_trdust1_diff'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/all_s_trdust1_diff'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/all_whole_trdust1_diff_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/all_whole_ps_trdust1_diff'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/all_whole_s_trdust1_diff'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/all_whole_i_trdust1_diff'    'rbZt' ''      'init'    ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/all_whole_m_trdust1_diff'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/all_whole_f_trdust1_diff'    'rbZt' ''      'final'   ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_d2_num_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_ps_d2_num'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_s_d2_num'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_whole_d2_num_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_whole_ps_d2_num'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_whole_s_d2_num'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_whole_i_d2_num'    'rbZt' ''      'init'    ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_whole_m_d2_num'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_whole_f_d2_num'    'rbZt' ''      'final'   ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_core_ps_d2_num'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_core_s_d2_num'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_rb_ps_d2_num'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_num_<CASE>.h5' '/d2_num' 'wtmean'  0.0  '/all_rb_s_d2_num'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/all_d2_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/all_ps_d2_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/all_s_d2_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/all_whole_d2_mass_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/all_whole_ps_d2_mass'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/all_whole_s_d2_mass'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/all_whole_i_d2_mass'    'rbZt' ''      'init'    ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/all_whole_m_d2_mass'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/all_whole_f_d2_mass'    'rbZt' ''      'final'   ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/all_trdust2_diff_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/all_ps_trdust2_diff'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/all_s_trdust2_diff'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/all_whole_trdust2_diff_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/all_whole_ps_trdust2_diff'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/all_whole_s_trdust2_diff'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/all_whole_i_trdust2_diff'    'rbZt' ''      'init'    ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/all_whole_m_trdust2_diff'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/all_whole_f_trdust2_diff'    'rbZt' ''      'final'   ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_dust_sfc_<CASE>.h5' '/dust_sfc' 'wtmean'  0.0  '/all_dust_sfc_ts'         'Rb_T' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_sfc_<CASE>.h5' '/dust_sfc' 'wtmean'  0.0  '/all_ps_dust_sfc'         'Rb_t' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_sfc_<CASE>.h5' '/dust_sfc' 'wtmean'  0.0  '/all_s_dust_sfc'          'Rb_t' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_sfc_<CASE>.h5' '/dust_sfc' 'wtmean'  0.0  '/all_whole_dust_sfc_ts'   'rb_T' ''      ''        ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/all_dust_adv_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/all_ps_dust_adv'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/all_s_dust_adv'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/all_whole_dust_adv_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/all_whole_ps_dust_adv'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/all_whole_s_dust_adv'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/all_whole_i_dust_adv'    'rbZt' ''      'init'    ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/all_whole_m_dust_adv'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/all_whole_f_dust_adv'    'rbZt' ''      'final'   ''       'ge' 0.01 }
%
        { 'AzAveragedData/hist_all_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/all_dust_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/all_ps_dust_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/all_s_dust_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/all_whole_ps_dust_mass'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/all_whole_s_dust_mass'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/all_whole_i_dust_mass'    'rbZt' ''      'init'    ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/all_whole_m_dust_mass'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/all_whole_f_dust_mass'    'rbZt' ''      'final'   ''       'ge' 0.01 }

%        { 'AzAveragedData/hist_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/tcond_dust_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/tcond_ps_dust_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/tcond_s_dust_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/tcond_whole_dust_mass_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/tcond_whole_ps_dust_mass'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/tcond_whole_s_dust_mass'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/tcond_whole_i_dust_mass'    'rbZt' ''      'init'    ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/tcond_whole_m_dust_mass'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/tcond_whole_f_dust_mass'    'rbZt' ''      'final'   ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/clear_dust_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/clear_ps_dust_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/clear_s_dust_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/clear_whole_dust_mass_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/clear_whole_ps_dust_mass'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/clear_whole_s_dust_mass'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/clear_whole_i_dust_mass'    'rbZt' ''      'init'    ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/clear_whole_m_dust_mass'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/clear_whole_f_dust_mass'    'rbZt' ''      'final'   ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_dust_mass_<CASE>.h5'   '/dust_mass'   'wtmean'  0.0  '/all_whole_dust_mass_ts'     'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_tracer_mass_<CASE>.h5' '/tracer_mass' 'wtmean'  0.0  '/all_whole_tracer_mass_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
      }
      'DIAGS/hist_meas_az_dust_<CASE>.h5'
    }
%
%    % ccn measurements
%    {
%      'CCN Azavg'
%      {
%        { 'AzAveragedData/hist_all_ccn_num_<CASE>.h5' '/ccn_num' 'wtmean'  0.0  '/all_ccn_num_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_num_<CASE>.h5' '/ccn_num' 'wtmean'  0.0  '/all_ps_ccn_num'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_num_<CASE>.h5' '/ccn_num' 'wtmean'  0.0  '/all_s_ccn_num'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_num_<CASE>.h5' '/ccn_num' 'wtmean'  0.0  '/all_whole_ccn_num_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_num_<CASE>.h5' '/ccn_num' 'wtmean'  0.0  '/all_whole_ps_ccn_num'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_num_<CASE>.h5' '/ccn_num' 'wtmean'  0.0  '/all_whole_s_ccn_num'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_num_<CASE>.h5' '/ccn_num' 'wtmean'  0.0  '/all_whole_i_ccn_num'    'rbZt' ''      'init'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_num_<CASE>.h5' '/ccn_num' 'wtmean'  0.0  '/all_whole_m_ccn_num'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_num_<CASE>.h5' '/ccn_num' 'wtmean'  0.0  '/all_whole_f_ccn_num'    'rbZt' ''      'final'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/all_ccn_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/all_ps_ccn_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/all_s_ccn_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/all_whole_ps_ccn_mass'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/all_whole_s_ccn_mass'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/all_whole_i_ccn_mass'    'rbZt' ''      'init'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/all_whole_m_ccn_mass'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/all_whole_f_ccn_mass'    'rbZt' ''      'final'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/tcond_ccn_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/tcond_ps_ccn_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/tcond_s_ccn_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/tcond_whole_ccn_mass_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/tcond_whole_ps_ccn_mass'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/tcond_whole_s_ccn_mass'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/tcond_whole_i_ccn_mass'    'rbZt' ''      'init'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/tcond_whole_m_ccn_mass'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/tcond_whole_f_ccn_mass'    'rbZt' ''      'final'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/clear_ccn_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/clear_ps_ccn_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/clear_s_ccn_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/clear_whole_ccn_mass_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/clear_whole_ps_ccn_mass'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/clear_whole_s_ccn_mass'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/clear_whole_i_ccn_mass'    'rbZt' ''      'init'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/clear_whole_m_ccn_mass'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/clear_whole_f_ccn_mass'    'rbZt' ''      'final'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/all_whole_ccn_mass_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_ccn_<CASE>.h5'
%    }

    % Regenerated aerosols measurements
    {
      'Regen Azavg'
      {
        { 'AzAveragedData/hist_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/all_ra_mass_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/all_ps_ra_mass'         'RbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/all_s_ra_mass'          'RbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/all_whole_ps_ra_mass'   'rbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/all_whole_s_ra_mass'    'rbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/all_whole_i_ra_mass'    'rbZt' ''      'init'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/all_whole_m_ra_mass'    'rbZt' ''      'mid'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/all_whole_f_ra_mass'    'rbZt' ''      'final'     ''       'gt' 0 }
%
%        { 'AzAveragedData/hist_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/tcond_ra_mass_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/tcond_ps_ra_mass'         'RbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/tcond_s_ra_mass'          'RbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/tcond_whole_ra_mass_ts'   'rbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/tcond_whole_ps_ra_mass'   'rbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/tcond_whole_s_ra_mass'    'rbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/tcond_whole_i_ra_mass'    'rbZt' ''      'init'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/tcond_whole_m_ra_mass'    'rbZt' ''      'mid'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/tcond_whole_f_ra_mass'    'rbZt' ''      'final'     ''       'gt' 0 }
%
%        { 'AzAveragedData/hist_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/clear_ra_mass_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/clear_ps_ra_mass'         'RbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/clear_s_ra_mass'          'RbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/clear_whole_ra_mass_ts'   'rbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/clear_whole_ps_ra_mass'   'rbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/clear_whole_s_ra_mass'    'rbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/clear_whole_i_ra_mass'    'rbZt' ''      'init'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/clear_whole_m_ra_mass'    'rbZt' ''      'mid'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/clear_whole_f_ra_mass'    'rbZt' ''      'final'     ''       'gt' 0 }
%
%        { 'AzAveragedData/hist_all_ra1_num_<CASE>.h5' '/ra1_num' 'wtmean'  0.0  '/all_ra1_num_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_num_<CASE>.h5' '/ra1_num' 'wtmean'  0.0  '/all_ps_ra1_num'         'RbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_num_<CASE>.h5' '/ra1_num' 'wtmean'  0.0  '/all_s_ra1_num'          'RbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_num_<CASE>.h5' '/ra1_num' 'wtmean'  0.0  '/all_whole_ra1_num_ts'   'rbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_num_<CASE>.h5' '/ra1_num' 'wtmean'  0.0  '/all_whole_ps_ra1_num'   'rbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_num_<CASE>.h5' '/ra1_num' 'wtmean'  0.0  '/all_whole_s_ra1_num'    'rbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_num_<CASE>.h5' '/ra1_num' 'wtmean'  0.0  '/all_whole_i_ra1_num'    'rbZt' ''      'init'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_num_<CASE>.h5' '/ra1_num' 'wtmean'  0.0  '/all_whole_m_ra1_num'    'rbZt' ''      'mid'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_num_<CASE>.h5' '/ra1_num' 'wtmean'  0.0  '/all_whole_f_ra1_num'    'rbZt' ''      'final'     ''       'gt' 0 }
%
%        { 'AzAveragedData/hist_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/all_ra1_mass_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/all_ps_ra1_mass'         'RbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/all_s_ra1_mass'          'RbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/all_whole_ra1_mass_ts'   'rbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/all_whole_ps_ra1_mass'   'rbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/all_whole_s_ra1_mass'    'rbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/all_whole_i_ra1_mass'    'rbZt' ''      'init'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/all_whole_m_ra1_mass'    'rbZt' ''      'mid'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/all_whole_f_ra1_mass'    'rbZt' ''      'final'     ''       'gt' 0 }
%
%        { 'AzAveragedData/hist_all_ra2_num_<CASE>.h5' '/ra2_num' 'wtmean'  0.0  '/all_ra2_num_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_num_<CASE>.h5' '/ra2_num' 'wtmean'  0.0  '/all_ps_ra2_num'         'RbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_num_<CASE>.h5' '/ra2_num' 'wtmean'  0.0  '/all_s_ra2_num'          'RbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_num_<CASE>.h5' '/ra2_num' 'wtmean'  0.0  '/all_whole_ra2_num_ts'   'rbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_num_<CASE>.h5' '/ra2_num' 'wtmean'  0.0  '/all_whole_ps_ra2_num'   'rbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_num_<CASE>.h5' '/ra2_num' 'wtmean'  0.0  '/all_whole_s_ra2_num'    'rbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_num_<CASE>.h5' '/ra2_num' 'wtmean'  0.0  '/all_whole_i_ra2_num'    'rbZt' ''      'init'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_num_<CASE>.h5' '/ra2_num' 'wtmean'  0.0  '/all_whole_m_ra2_num'    'rbZt' ''      'mid'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_num_<CASE>.h5' '/ra2_num' 'wtmean'  0.0  '/all_whole_f_ra2_num'    'rbZt' ''      'fianl'     ''       'gt' 0 }
%
%        { 'AzAveragedData/hist_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/all_ra2_mass_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/all_ps_ra2_mass'         'RbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/all_s_ra2_mass'          'RbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/all_whole_ra2_mass_ts'   'rbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/all_whole_ps_ra2_mass'   'rbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/all_whole_s_ra2_mass'    'rbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/all_whole_i_ra2_mass'    'rbZt' ''      'init'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/all_whole_m_ra2_mass'    'rbZt' ''      'mid'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/all_whole_f_ra2_mass'    'rbZt' ''      'final'     ''       'gt' 0 }
%
%        { 'AzAveragedData/hist_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/all_whole_ra_mass_ts'   'rbZT' ''      ''        ''       'gt' 0 }
%
%        { 'AzAveragedData/hist_cloud_all_ra_mass_<CASE>.h5'     '/ra_mass' 'wtmean'  0.0  '/cloud_all_ra_mass_ts'           'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ra_mass_<CASE>.h5'     '/ra_mass' 'wtmean'  0.0  '/cloud_all_whole_ra_mass_ts'     'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ra_mass_<CASE>.h5'     '/ra_mass' 'wtmean'  0.0  '/cloud_all_whole_ps_ra_mass'     'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ra_mass_<CASE>.h5'     '/ra_mass' 'wtmean'  0.0  '/cloud_all_whole_s_ra_mass'      'rbZt' '' 'sal'     '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/not_cloud_all_ra_mass_ts'       'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/not_cloud_all_whole_ra_mass_ts' 'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/not_cloud_all_whole_ps_ra_mass' 'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/not_cloud_all_whole_s_ra_mass'  'rbZt' '' 'sal'     '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra_mass_<CASE>.h5'     '/ra_mass' 'wtmean'  0.0  '/tcond_all_ra_mass_ts'           'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra_mass_<CASE>.h5'     '/ra_mass' 'wtmean'  0.0  '/tcond_all_whole_ra_mass_ts'     'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra_mass_<CASE>.h5'     '/ra_mass' 'wtmean'  0.0  '/tcond_all_whole_ps_ra_mass'     'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra_mass_<CASE>.h5'     '/ra_mass' 'wtmean'  0.0  '/tcond_all_whole_s_ra_mass'      'rbZt' '' 'sal'     '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/not_tcond_all_ra_mass_ts'       'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/not_tcond_all_whole_ra_mass_ts' 'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/not_tcond_all_whole_ps_ra_mass' 'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/not_tcond_all_whole_s_ra_mass'  'rbZt' '' 'sal'     '' '' 'gt' 0 }
%
%        { 'AzAveragedData/hist_cloud_all_ra1_mass_<CASE>.h5'     '/ra1_mass' 'wtmean'  0.0  '/cloud_all_ra1_mass_ts'           'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ra1_mass_<CASE>.h5'     '/ra1_mass' 'wtmean'  0.0  '/cloud_all_whole_ra1_mass_ts'     'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ra1_mass_<CASE>.h5'     '/ra1_mass' 'wtmean'  0.0  '/cloud_all_whole_ps_ra1_mass'     'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ra1_mass_<CASE>.h5'     '/ra1_mass' 'wtmean'  0.0  '/cloud_all_whole_s_ra1_mass'      'rbZt' '' 'sal'     '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/not_cloud_all_ra1_mass_ts'       'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/not_cloud_all_whole_ra1_mass_ts' 'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/not_cloud_all_whole_ps_ra1_mass' 'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/not_cloud_all_whole_s_ra1_mass'  'rbZt' '' 'sal'     '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra1_mass_<CASE>.h5'     '/ra1_mass' 'wtmean'  0.0  '/tcond_all_ra1_mass_ts'           'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra1_mass_<CASE>.h5'     '/ra1_mass' 'wtmean'  0.0  '/tcond_all_whole_ra1_mass_ts'     'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra1_mass_<CASE>.h5'     '/ra1_mass' 'wtmean'  0.0  '/tcond_all_whole_ps_ra1_mass'     'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra1_mass_<CASE>.h5'     '/ra1_mass' 'wtmean'  0.0  '/tcond_all_whole_s_ra1_mass'      'rbZt' '' 'sal'     '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/not_tcond_all_ra1_mass_ts'       'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/not_tcond_all_whole_ra1_mass_ts' 'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/not_tcond_all_whole_ps_ra1_mass' 'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/not_tcond_all_whole_s_ra1_mass'  'rbZt' '' 'sal'     '' '' 'gt' 0 }
%
%        { 'AzAveragedData/hist_cloud_all_ra2_mass_<CASE>.h5'     '/ra2_mass' 'wtmean'  0.0  '/cloud_all_ra2_mass_ts'           'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ra2_mass_<CASE>.h5'     '/ra2_mass' 'wtmean'  0.0  '/cloud_all_whole_ra2_mass_ts'     'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ra2_mass_<CASE>.h5'     '/ra2_mass' 'wtmean'  0.0  '/cloud_all_whole_ps_ra2_mass'     'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ra2_mass_<CASE>.h5'     '/ra2_mass' 'wtmean'  0.0  '/cloud_all_whole_s_ra2_mass'      'rbZt' '' 'sal'     '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/not_cloud_all_ra2_mass_ts'       'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/not_cloud_all_whole_ra2_mass_ts' 'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/not_cloud_all_whole_ps_ra2_mass' 'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/not_cloud_all_whole_s_ra2_mass'  'rbZt' '' 'sal'     '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra2_mass_<CASE>.h5'     '/ra2_mass' 'wtmean'  0.0  '/tcond_all_ra2_mass_ts'           'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra2_mass_<CASE>.h5'     '/ra2_mass' 'wtmean'  0.0  '/tcond_all_whole_ra2_mass_ts'     'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra2_mass_<CASE>.h5'     '/ra2_mass' 'wtmean'  0.0  '/tcond_all_whole_ps_ra2_mass'     'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ra2_mass_<CASE>.h5'     '/ra2_mass' 'wtmean'  0.0  '/tcond_all_whole_s_ra2_mass'      'rbZt' '' 'sal'     '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/not_tcond_all_ra2_mass_ts'       'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/not_tcond_all_whole_ra2_mass_ts' 'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/not_tcond_all_whole_ps_ra2_mass' 'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/not_tcond_all_whole_s_ra2_mass'  'rbZt' '' 'sal'     '' 'gt' 0 }
      }
      'DIAGS/hist_meas_az_ra_<CASE>.h5'
    }

%    % "aero" (dust + regen) measurements
%    {
%      'AERO Azavg'
%      {
%        { 'AzAveragedData/hist_all_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/all_aero_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/all_ps_aero_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/all_s_aero_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/all_whole_ps_aero_mass'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/all_whole_s_aero_mass'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/all_whole_i_aero_mass'    'rbZt' ''      'init'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/all_whole_m_aero_mass'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/all_whole_f_aero_mass'    'rbZt' ''      'final'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/tcond_aero_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/tcond_ps_aero_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/tcond_s_aero_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/tcond_whole_aero_mass_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/tcond_whole_ps_aero_mass'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/tcond_whole_s_aero_mass'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/tcond_whole_i_aero_mass'    'rbZt' ''      'init'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/tcond_whole_m_aero_mass'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/tcond_whole_f_aero_mass'    'rbZt' ''      'final'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/clear_aero_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/clear_ps_aero_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/clear_s_aero_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/clear_whole_aero_mass_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/clear_whole_ps_aero_mass'   'rbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/clear_whole_s_aero_mass'    'rbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/clear_whole_i_aero_mass'    'rbZt' ''      'init'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/clear_whole_m_aero_mass'    'rbZt' ''      'mid'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/clear_whole_f_aero_mass'    'rbZt' ''      'final'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/all_whole_aero_mass_ts'   'rbZT' ''      ''        ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_aero_<CASE>.h5'
%    }
%
%    % cloud
%    {
%      'Cloud'
%      {
%        { 'AzAveragedData/hist_all_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_cloud_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_ps_cloud_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_s_cloud_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_core_cloud_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_core_ps_cloud_mass'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_core_s_cloud_mass'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_rb_cloud_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_rb_ps_cloud_mass'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/all_rb_s_cloud_mass'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/all_cloud_num_ts'         'RbZT' ''      ''        ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_all_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/all_ps_cloud_num'         'RbZt' ''      'pre_sal' ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_all_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/all_s_cloud_num'          'RbZt' ''      'sal'     ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_all_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/all_core_cloud_num_ts'    'rbZT' 'core'  ''        ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_all_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/all_core_ps_cloud_num'    'rbZt' 'core'  'pre_sal' ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_all_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/all_core_s_cloud_num'     'rbZt' 'core'  'sal'     ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_all_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/all_rb_cloud_num_ts'      'rbZT' 'rband' ''        ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_all_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/all_rb_ps_cloud_num'      'rbZt' 'rband' 'pre_sal' ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_all_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/all_rb_s_cloud_num'       'rbZt' 'rband' 'sal'     ''       'ge' 1e6 }
%
%        { 'AzAveragedData/hist_all_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/all_cloud_diam_ts'         'RbZT' ''      ''        ''       'ge' 20 }
%        { 'AzAveragedData/hist_all_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/all_ps_cloud_diam'         'RbZt' ''      'pre_sal' ''       'ge' 20 }
%        { 'AzAveragedData/hist_all_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/all_s_cloud_diam'          'RbZt' ''      'sal'     ''       'ge' 20 }
%        { 'AzAveragedData/hist_all_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/all_core_cloud_diam_ts'    'rbZT' 'core'  ''        ''       'ge' 20 }
%        { 'AzAveragedData/hist_all_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/all_core_ps_cloud_diam'    'rbZt' 'core'  'pre_sal' ''       'ge' 20 }
%        { 'AzAveragedData/hist_all_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/all_core_s_cloud_diam'     'rbZt' 'core'  'sal'     ''       'ge' 20 }
%        { 'AzAveragedData/hist_all_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/all_rb_cloud_diam_ts'      'rbZT' 'rband' ''        ''       'ge' 20 }
%        { 'AzAveragedData/hist_all_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/all_rb_ps_cloud_diam'      'rbZt' 'rband' 'pre_sal' ''       'ge' 20 }
%        { 'AzAveragedData/hist_all_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/all_rb_s_cloud_diam'       'rbZt' 'rband' 'sal'     ''       'ge' 20 }
%
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_cloud_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_ps_cloud_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_s_cloud_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_core_cloud_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_core_ps_cloud_mass'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_core_s_cloud_mass'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_rb_cloud_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_rb_ps_cloud_mass'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_rb_s_cloud_mass'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_env_cloud_mass_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_env_ps_cloud_mass'     'rbZt' 'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_cloud_<CASE>.h5' '/cloud' 'wtmean'  0.0  '/lead_env_s_cloud_mass'      'rbZt' 'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/lead_cloud_num_ts'         'RbZT' ''      ''        ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_lead_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/lead_ps_cloud_num'         'RbZt' ''      'pre_sal' ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_lead_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/lead_s_cloud_num'          'RbZt' ''      'sal'     ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_lead_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/lead_core_cloud_num_ts'    'rbZT' 'core'  ''        ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_lead_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/lead_core_ps_cloud_num'    'rbZt' 'core'  'pre_sal' ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_lead_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/lead_core_s_cloud_num'     'rbZt' 'core'  'sal'     ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_lead_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/lead_rb_cloud_num_ts'      'rbZT' 'rband' ''        ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_lead_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/lead_rb_ps_cloud_num'      'rbZt' 'rband' 'pre_sal' ''       'ge' 1e6 }
%        { 'AzAveragedData/hist_lead_cloud_num_<CASE>.h5' '/cloud_num' 'wtmean'  0.0  '/lead_rb_s_cloud_num'       'rbZt' 'rband' 'sal'     ''       'ge' 1e6 }
%
%        { 'AzAveragedData/hist_lead_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/lead_cloud_diam_ts'         'RbZT' ''      ''        ''       'ge' 20 }
%        { 'AzAveragedData/hist_lead_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/lead_ps_cloud_diam'         'RbZt' ''      'pre_sal' ''       'ge' 20 }
%        { 'AzAveragedData/hist_lead_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/lead_s_cloud_diam'          'RbZt' ''      'sal'     ''       'ge' 20 }
%        { 'AzAveragedData/hist_lead_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/lead_core_cloud_diam_ts'    'rbZT' 'core'  ''        ''       'ge' 20 }
%        { 'AzAveragedData/hist_lead_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/lead_core_ps_cloud_diam'    'rbZt' 'core'  'pre_sal' ''       'ge' 20 }
%        { 'AzAveragedData/hist_lead_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/lead_core_s_cloud_diam'     'rbZt' 'core'  'sal'     ''       'ge' 20 }
%        { 'AzAveragedData/hist_lead_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/lead_rb_cloud_diam_ts'      'rbZT' 'rband' ''        ''       'ge' 20 }
%        { 'AzAveragedData/hist_lead_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/lead_rb_ps_cloud_diam'      'rbZt' 'rband' 'pre_sal' ''       'ge' 20 }
%        { 'AzAveragedData/hist_lead_cloud_diam_<CASE>.h5' '/cloud_diam' 'wtmean'  0.0  '/lead_rb_s_cloud_diam'       'rbZt' 'rband' 'sal'     ''       'ge' 20 }
%      }
%      'DIAGS/hist_meas_az_cloud_<CASE>.h5'
%    }
%
%    % theta_v measurements
%    {
%      'Theta-V Azavg'
%      {
%        { 'AzAveragedData/hist_all_theta_v_<CASE>.h5' '/theta_v' 'wtmean'  0.0  '/all_theta_v_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_v_<CASE>.h5' '/theta_v' 'wtmean'  0.0  '/all_ps_theta_v'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_v_<CASE>.h5' '/theta_v' 'wtmean'  0.0  '/all_s_theta_v'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_v_<CASE>.h5' '/theta_v' 'wtmean'  0.0  '/all_core_theta_v_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_v_<CASE>.h5' '/theta_v' 'wtmean'  0.0  '/all_core_ps_theta_v'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_v_<CASE>.h5' '/theta_v' 'wtmean'  0.0  '/all_core_s_theta_v'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_v_<CASE>.h5' '/theta_v' 'wtmean'  0.0  '/all_rb_theta_v_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_v_<CASE>.h5' '/theta_v' 'wtmean'  0.0  '/all_rb_ps_theta_v'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_v_<CASE>.h5' '/theta_v' 'wtmean'  0.0  '/all_rb_s_theta_v'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_theta_v_<CASE>.h5'
%    }
%
%
%    % theta_e measurements
%    {
%      'Theta-E Azavg'
%      {
%        { 'AzAveragedData/hist_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_theta_e_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_ps_theta_e'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_s_theta_e'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_theta_e_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_ps_theta_e'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_core_s_theta_e'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_theta_e_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_ps_theta_e'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/all_rb_s_theta_e'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_theta_e_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_ps_theta_e'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_s_theta_e'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_core_theta_e_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_core_ps_theta_e'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_core_s_theta_e'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_rb_theta_e_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_rb_ps_theta_e'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_rb_s_theta_e'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_env_theta_e_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_env_ps_theta_e'     'rbZt' 'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/lead_env_s_theta_e'      'rbZt' 'env'   'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_precip_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_all_theta_e_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_all_ps_theta_e'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_all_s_theta_e'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_all_core_theta_e_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_all_core_ps_theta_e'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_all_core_s_theta_e'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_all_rb_theta_e_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_all_rb_ps_theta_e'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_all_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_all_rb_s_theta_e'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_precip_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_lead_theta_e_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_lead_ps_theta_e'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_lead_s_theta_e'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_lead_core_theta_e_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_lead_core_ps_theta_e'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_lead_core_s_theta_e'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_lead_rb_theta_e_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_lead_rb_ps_theta_e'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_precip_lead_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/precip_lead_rb_s_theta_e'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_theta_e_<CASE>.h5'
%    }
%
%    % theta_rho measurements
%    {
%      'Theta-Rho Azavg'
%      {
%        { 'AzAveragedData/hist_all_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/all_theta_rho_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/all_ps_theta_rho'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/all_s_theta_rho'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/all_core_theta_rho_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/all_core_ps_theta_rho'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/all_core_s_theta_rho'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/all_rb_theta_rho_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/all_rb_ps_theta_rho'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/all_rb_s_theta_rho'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_theta_rho_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_ps_theta_rho'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_s_theta_rho'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_core_theta_rho_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_core_ps_theta_rho'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_core_s_theta_rho'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_rb_theta_rho_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_rb_ps_theta_rho'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_rb_s_theta_rho'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_env_theta_rho_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_env_ps_theta_rho'     'rbZt' 'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_rho_<CASE>.h5' '/theta_rho' 'wtmean'  0.0  '/lead_env_s_theta_rho'      'rbZt' 'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_theta_rho_<CASE>.h5'
%    }
%
%    % theta measurements
%    {
%      'Theta Azavg'
%      {
%        { 'AzAveragedData/hist_all_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_theta_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_ps_theta'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_s_theta'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_core_theta_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_core_ps_theta'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_core_s_theta'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_rb_theta_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_rb_ps_theta'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/all_rb_s_theta'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_theta_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_ps_theta'         'RbZt' ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_s_theta'          'RbZt' ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_core_theta_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_core_ps_theta'    'rbZt' 'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_core_s_theta'     'rbZt' 'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_rb_theta_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_rb_ps_theta'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_rb_s_theta'       'rbZt' 'rband' 'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_env_theta_ts'     'rbZT' 'env'   ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_env_ps_theta'     'rbZt' 'env'   'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/lead_env_s_theta'      'rbZt' 'env'   'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_theta_<CASE>.h5'
%    }
%
%
%    % precip rate measurements
%    {
%      'Precip Rate Azavg'
%      {
%        { 'AzAveragedData/hist_all_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_pcprate_ts'      'Rb_T'  ''      ''        ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_all_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_ps_pcprate_ts'   'Rb_T'  ''      'pre_sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_all_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_s_pcprate_ts'    'Rb_T'  ''      'sal'     ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_all_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_ps_pcprate'      'Rb_t'  ''      'pre_sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_all_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_s_pcprate'       'Rb_t'  ''      'sal'     ''       'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_all_core_ps_pcprate' 'rB_t'  'core' 'pre_sal' ''   'ge' 0.001 }
%        { 'AzAveragedData/hist_all_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_all_core_s_pcprate'  'rB_t'  'core' 'sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_all_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_all_rb_ps_pcprate'   'rB_t'  'rb'   'pre_sal' ''   'ge' 0.001 }
%        { 'AzAveragedData/hist_all_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_all_rb_s_pcprate'    'rB_t'  'rb'   'sal' ''       'ge' 0.001 }
%
%        { 'AzAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/lead_pcprate_ts'      'Rb_T'  ''      ''        ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/lead_ps_pcprate_ts'   'Rb_T'  ''      'pre_sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/lead_s_pcprate_ts'    'Rb_T'  ''      'sal'     ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/lead_ps_pcprate'      'Rb_t'  ''      'pre_sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/lead_s_pcprate'       'Rb_t'  ''      'sal'     ''       'ge' 0.001 }
%
%        { 'AzAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_lead_core_ps_pcprate' 'rB_t'  'core' 'pre_sal' ''   'ge' 0.001 }
%        { 'AzAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_lead_core_s_pcprate'  'rB_t'  'core' 'sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_lead_rb_ps_pcprate'   'rB_t'  'rb'   'pre_sal' ''   'ge' 0.001 }
%        { 'AzAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_lead_rb_s_pcprate'    'rB_t'  'rb'   'sal' ''       'ge' 0.001 }
%        { 'AzAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_lead_env_ps_pcprate'  'rB_t'  'env'  'pre_sal' ''   'ge' 0.001 }
%        { 'AzAveragedData/hist_lead_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/hist_lead_env_s_pcprate'   'rB_t'  'env'  'sal' ''       'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/all_whole_pcprate_ts'      'rb_T'  ''      ''        ''       'ge' 0.001 }
%      }
%      'DIAGS/hist_meas_az_pcprate_<CASE>.h5'
%    }
%
%    % storm speed measurements
%    {
%      'Speed'
%      {
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_speed_ts'           'RbZT' ''      ''        ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_ps_speed_ts'        'RbZT' ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_s_speed_ts'         'RbZT' ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_ps_speed'           'RbZt' ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_s_speed'            'RbZt' ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_core_ps_speed'      'rbZt' 'core'  'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_core_s_speed'       'rbZt' 'core'  'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_rb_ps_speed'        'rbZt' 'rband' 'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_rb_s_speed'         'rbZt' 'rband' 'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_ps_speed_sfc'       'Rbzt' ''      'pre_sal' 'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_s_speed_sfc'        'Rbzt' ''      'sal'     'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_ps_speed_maxlev'    'Rbzt' ''      'pre_sal' 'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_s_speed_maxlev'     'Rbzt' ''      'sal'     'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_speed_sfc_ts'       'RbzT' ''      ''        'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_speed_maxlev_ts'    'RbzT' ''      ''        'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_ps_speed_sfc_ts'    'RbzT' ''      'pre_sal' 'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_ps_speed_maxlev_ts' 'RbzT' ''      'pre_sal' 'maxlev' 'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_s_speed_sfc_ts'     'RbzT' ''      'sal'     'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_all_speed_<CASE>.h5' '/speed' 'wtmean'  0.0  '/all_s_speed_maxlev_ts'  'RbzT' ''      'sal'     'maxlev' 'ge' 0   }
%
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_speed_t_ts'           'RbZT' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t_ts'        'RbZT' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t_ts'         'RbZT' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t'           'RbZt' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t'            'RbZt' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_whole_speed_t'        'rbZT' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_whole_ps_speed_t'     'rbZt' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_whole_s_speed_t'      'rbZt' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_core_ps_speed_t'      'rbZt' 'core'  'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_core_s_speed_t'       'rbZt' 'core'  'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_rb_ps_speed_t'        'rbZt' 'rband' 'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_rb_s_speed_t'         'rbZt' 'rband' 'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t_sfc'       'Rbzt' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t_sfc'        'Rbzt' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t_maxlev'    'Rbzt' ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t_maxlev'     'Rbzt' ''      'sal'     'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_speed_t_sfc_ts'       'RbzT' ''      ''        'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_speed_t_maxlev_ts'    'RbzT' ''      ''        'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t_sfc_ts'    'RbzT' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_ps_speed_t_maxlev_ts' 'RbzT' ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t_sfc_ts'     'RbzT' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_t_<CASE>.h5' '/speed_t' 'wtmean'  0.0  '/all_s_speed_t_maxlev_ts'  'RbzT' ''      'sal'     'maxlev' 'ge' -100 }
%
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_speed_r_ts'           'RbZT' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r_ts'        'RbZT' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r_ts'         'RbZT' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r'           'RbZt' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r'            'RbZt' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_whole_speed_r'        'rbZT' ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_whole_ps_speed_r'     'rbZt' ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_whole_s_speed_r'      'rbZt' ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_core_ps_speed_r'      'rbZt' 'core'  'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_core_s_speed_r'       'rbZt' 'core'  'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_rb_ps_speed_r'        'rbZt' 'rband' 'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_rb_s_speed_r'         'rbZt' 'rband' 'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r_sfc'       'Rbzt' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r_sfc'        'Rbzt' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r_minlev'    'Rbzt' ''      'pre_sal' 'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r_minlev'     'Rbzt' ''      'sal'     'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_speed_r_sfc_ts'       'RbzT' ''      ''        'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_speed_r_minlev_ts'    'RbzT' ''      ''        'minlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r_sfc_ts'    'RbzT' ''      'pre_sal' 'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_ps_speed_r_maxlev_ts' 'RbzT' ''      'pre_sal' 'maxlev' 'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r_sfc_ts'     'RbzT' ''      'sal'     'sfc'    'ge' -100 }
%        { 'AzAveragedData/hist_all_speed_r_<CASE>.h5' '/speed_r' 'wtmean'  0.0  '/all_s_speed_r_maxlev_ts'  'RbzT' ''      'sal'     'maxlev' 'ge' -100 }
%
%        { 'AzAveragedData/hist_all_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/all_ps_speed10m'      'Rb_t' ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/all_s_speed10m'       'Rb_t' ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/all_speed10m_ts'      'Rb_T' ''      ''        ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/all_ps_speed10m_ts'   'Rb_T' ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_speed10m_<CASE>.h5' '/speed10m' 'wtmean'  0.0  '/all_s_speed10m_ts'    'Rb_T' ''      'sal'     ''       'ge' 0   }
%      }
%      'DIAGS/hist_meas_az_speed_<CASE>.h5'
%    }
%
%    {
%      'Pressure'
%      {
%        { 'AzAveragedData/hist_all_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/all_press_ts'         'RbZT' ''      ''        ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/all_ps_press'         'RbZt' ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/all_s_press'          'RbZt' ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/all_ps_press_sfc'     'Rbzt' ''      'pre_sal' 'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_all_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/all_s_press_sfc'      'Rbzt' ''      'sal'     'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_all_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/all_press_sfc_ts'     'RbzT' ''      ''        'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_all_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/all_ps_press_sfc_ts'  'RbzT' ''      'pre_sal' 'sfc'    'ge' 0   }
%        { 'AzAveragedData/hist_all_press_<CASE>.h5' '/press' 'wtmean'  0.0  '/all_s_press_sfc_ts'   'RbzT' ''      'sal'     'sfc'    'ge' 0   }
%      }
%      'DIAGS/hist_meas_az_press_<CASE>.h5'
%    }
%
%    % rain
%    {
%      'Rain'
%      {
%        { 'AzAveragedData/hist_all_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_rain_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_ps_rain_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_s_rain_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_core_rain_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_rb_rain_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_core_ps_rain_mass'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_core_s_rain_mass'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_rb_ps_rain_mass'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/all_rb_s_rain_mass'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/all_rain_num_ts'         'RbZT' ''      ''        ''       'ge' 1000 }
%        { 'AzAveragedData/hist_all_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/all_ps_rain_num'         'RbZt' ''      'pre_sal' ''       'ge' 1000 }
%        { 'AzAveragedData/hist_all_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/all_s_rain_num'          'RbZt' ''      'sal'     ''       'ge' 1000 }
%        { 'AzAveragedData/hist_all_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/all_core_rain_num_ts'    'rbZT' 'core'  ''        ''       'ge' 1000 }
%        { 'AzAveragedData/hist_all_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/all_rb_rain_num_ts'      'rbZT' 'rband' ''        ''       'ge' 1000 }
%        { 'AzAveragedData/hist_all_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/all_core_ps_rain_num'    'rbZt' 'core'  'pre_sal' ''       'ge' 1000 }
%        { 'AzAveragedData/hist_all_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/all_core_s_rain_num'     'rbZt' 'core'  'sal'     ''       'ge' 1000 }
%        { 'AzAveragedData/hist_all_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/all_rb_ps_rain_num'      'rbZt' 'rband' 'pre_sal' ''       'ge' 1000 }
%        { 'AzAveragedData/hist_all_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/all_rb_s_rain_num'       'rbZt' 'rband' 'sal'     ''       'ge' 1000 }
%
%        { 'AzAveragedData/hist_all_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/all_rain_diam_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/all_ps_rain_diam'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/all_s_rain_diam'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/all_core_rain_diam_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/all_rb_rain_diam_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/all_core_ps_rain_diam'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/all_core_s_rain_diam'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/all_rb_ps_rain_diam'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/all_rb_s_rain_diam'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_rain_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_ps_rain_mass'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_s_rain_mass'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_core_rain_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_rb_rain_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_env_rain_mass_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_core_ps_rain_mass'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_core_s_rain_mass'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_rb_ps_rain_mass'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_rb_s_rain_mass'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_env_ps_rain_mass'     'rbZt' 'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_<CASE>.h5' '/rain' 'wtmean'  0.0  '/lead_env_s_rain_mass'      'rbZt' 'env'   'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/lead_rain_num_ts'         'RbZT' ''      ''        ''       'ge' 1000 }
%        { 'AzAveragedData/hist_lead_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/lead_ps_rain_num'         'RbZt' ''      'pre_sal' ''       'ge' 1000 }
%        { 'AzAveragedData/hist_lead_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/lead_s_rain_num'          'RbZt' ''      'sal'     ''       'ge' 1000 }
%        { 'AzAveragedData/hist_lead_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/lead_core_rain_num_ts'    'rbZT' 'core'  ''        ''       'ge' 1000 }
%        { 'AzAveragedData/hist_lead_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/lead_rb_rain_num_ts'      'rbZT' 'rband' ''        ''       'ge' 1000 }
%        { 'AzAveragedData/hist_lead_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/lead_core_ps_rain_num'    'rbZt' 'core'  'pre_sal' ''       'ge' 1000 }
%        { 'AzAveragedData/hist_lead_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/lead_core_s_rain_num'     'rbZt' 'core'  'sal'     ''       'ge' 1000 }
%        { 'AzAveragedData/hist_lead_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/lead_rb_ps_rain_num'      'rbZt' 'rband' 'pre_sal' ''       'ge' 1000 }
%        { 'AzAveragedData/hist_lead_rain_num_<CASE>.h5' '/rain_num' 'wtmean'  0.0  '/lead_rb_s_rain_num'       'rbZt' 'rband' 'sal'     ''       'ge' 1000 }
%
%        { 'AzAveragedData/hist_lead_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/lead_rain_diam_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/lead_ps_rain_diam'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/lead_s_rain_diam'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/lead_core_rain_diam_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/lead_rb_rain_diam_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/lead_core_ps_rain_diam'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/lead_core_s_rain_diam'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/lead_rb_ps_rain_diam'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_diam_<CASE>.h5' '/rain_diam' 'wtmean'  0.0  '/lead_rb_s_rain_diam'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_rain_<CASE>.h5'
%    }
%
%    % pristine ice
%    {
%      'Pristine'
%      {
%        { 'AzAveragedData/hist_all_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/all_pris_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/all_ps_pris_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/all_s_pris_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/all_core_pris_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/all_rb_pris_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/all_core_ps_pris_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/all_core_s_pris_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/all_rb_ps_pris_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/all_rb_s_pris_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/lead_pris_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/lead_ps_pris_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/lead_s_pris_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/lead_core_pris_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/lead_rb_pris_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/lead_core_ps_pris_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/lead_core_s_pris_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/lead_rb_ps_pris_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_pris_<CASE>.h5' '/pris' 'wtmean'  0.0  '/lead_rb_s_pris_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_pris_<CASE>.h5'
%    }
%
%    % aggregates
%    {
%      'Aggregates'
%      {
%        { 'AzAveragedData/hist_all_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/all_aggr_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/all_ps_aggr_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/all_s_aggr_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/all_core_aggr_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/all_rb_aggr_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/all_core_ps_aggr_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/all_core_s_aggr_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/all_rb_ps_aggr_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/all_rb_s_aggr_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/lead_aggr_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/lead_ps_aggr_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/lead_s_aggr_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/lead_core_aggr_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/lead_rb_aggr_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/lead_core_ps_aggr_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/lead_core_s_aggr_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/lead_rb_ps_aggr_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_aggr_<CASE>.h5' '/aggr' 'wtmean'  0.0  '/lead_rb_s_aggr_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_aggr_<CASE>.h5'
%    }
%
%    % snow
%    {
%      'Snow'
%      {
%        { 'AzAveragedData/hist_all_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/all_snow_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/all_ps_snow_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/all_s_snow_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/all_core_snow_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/all_rb_snow_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/all_core_ps_snow_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/all_core_s_snow_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/all_rb_ps_snow_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/all_rb_s_snow_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/lead_snow_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/lead_ps_snow_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/lead_s_snow_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/lead_core_snow_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/lead_rb_snow_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/lead_core_ps_snow_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/lead_core_s_snow_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/lead_rb_ps_snow_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_snow_<CASE>.h5' '/snow' 'wtmean'  0.0  '/lead_rb_s_snow_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_snow_<CASE>.h5'
%    }
%
%    % graupel
%    {
%      'Graupel'
%      {
%        { 'AzAveragedData/hist_all_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/all_graup_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/all_ps_graup_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/all_s_graup_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/all_core_graup_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/all_rb_graup_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/all_core_ps_graup_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/all_core_s_graup_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/all_rb_ps_graup_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/all_rb_s_graup_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/lead_graup_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/lead_ps_graup_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/lead_s_graup_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/lead_core_graup_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/lead_rb_graup_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/lead_core_ps_graup_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/lead_core_s_graup_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/lead_rb_ps_graup_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_graup_<CASE>.h5' '/graup' 'wtmean'  0.0  '/lead_rb_s_graup_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_graup_<CASE>.h5'
%    }
%
%    % hail
%    {
%      'Hail'
%      {
%        { 'AzAveragedData/hist_all_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/all_hail_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/all_ps_hail_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/all_s_hail_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/all_core_hail_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/all_rb_hail_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/all_core_ps_hail_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/all_core_s_hail_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/all_rb_ps_hail_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/all_rb_s_hail_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/lead_hail_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/lead_ps_hail_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/lead_s_hail_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/lead_core_hail_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/lead_rb_hail_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/lead_core_ps_hail_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/lead_core_s_hail_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/lead_rb_ps_hail_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_hail_<CASE>.h5' '/hail' 'wtmean'  0.0  '/lead_rb_s_hail_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_hail_<CASE>.h5'
%    }
%
%    % total condensate
%    {
%      'Total Condensate'
%      {
%        { 'AzAveragedData/hist_all_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/all_tcond_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/all_ps_tcond_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/all_s_tcond_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/all_core_tcond_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/all_rb_tcond_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/all_core_ps_tcond_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/all_core_s_tcond_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/all_rb_ps_tcond_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/all_rb_s_tcond_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/lead_tcond_mass_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/lead_ps_tcond_mass'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/lead_s_tcond_mass'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/lead_core_tcond_mass_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/lead_rb_tcond_mass_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/lead_core_ps_tcond_mass'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/lead_core_s_tcond_mass'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/lead_rb_ps_tcond_mass'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_tcond_<CASE>.h5' '/tcond' 'wtmean'  0.0  '/lead_rb_s_tcond_mass'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_tcond_<CASE>.h5'
%    }
%
%    % vapor mixing ratio
%    {
%      'Vapor Azavg'
%      {
%        { 'AzAveragedData/hist_all_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_vapor_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_ps_vapor'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_s_vapor'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_core_vapor_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_rb_vapor_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_core_ps_vapor'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_core_s_vapor'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_rb_ps_vapor'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/all_rb_s_vapor'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_vapor_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_ps_vapor'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_s_vapor'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_core_vapor_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_rb_vapor_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_env_vapor_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_core_ps_vapor'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_core_s_vapor'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_rb_ps_vapor'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_rb_s_vapor'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_env_ps_vapor'     'rbZt'    'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_env_s_vapor'      'rbZt'    'env'   'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_vapor_<CASE>.h5'
%    }
%
%    % tempc measurements
%    {
%      'Temperature Azavg'
%      {
%        { 'AzAveragedData/hist_all_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_tempc_ts'         'RbZT'   ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_ps_tempc'         'RbZt'   ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_s_tempc'          'RbZt'   ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_core_tempc_ts'    'rbZT'   'core'  ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_rb_tempc_ts'      'rbZT'   'rband' ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_core_ps_tempc'    'rbZt'   'core'  'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_core_s_tempc'     'rbZt'   'core'  'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_rb_ps_tempc'      'rbZt'   'rband' 'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/all_rb_s_tempc'       'rbZt'   'rband' 'sal'     ''       'ge' -100 }
%      }
%      'DIAGS/hist_meas_az_tempc_<CASE>.h5'
%    }
%
%    % dewptc measurements
%    {
%      'Dew Point Azavg'
%      {
%        { 'AzAveragedData/hist_all_dewptc_<CASE>.h5' '/dewptc' 'wtmean'  0.0  '/all_dewptc_ts'         'RbZT'   ''      ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_dewptc_<CASE>.h5' '/dewptc' 'wtmean'  0.0  '/all_ps_dewptc'         'RbZt'   ''      'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_dewptc_<CASE>.h5' '/dewptc' 'wtmean'  0.0  '/all_s_dewptc'          'RbZt'   ''      'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_dewptc_<CASE>.h5' '/dewptc' 'wtmean'  0.0  '/all_core_dewptc_ts'    'rbZT'   'core'  ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_dewptc_<CASE>.h5' '/dewptc' 'wtmean'  0.0  '/all_rb_dewptc_ts'      'rbZT'   'rband' ''        ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_dewptc_<CASE>.h5' '/dewptc' 'wtmean'  0.0  '/all_core_ps_dewptc'    'rbZt'   'core'  'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_dewptc_<CASE>.h5' '/dewptc' 'wtmean'  0.0  '/all_core_s_dewptc'     'rbZt'   'core'  'sal'     ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_dewptc_<CASE>.h5' '/dewptc' 'wtmean'  0.0  '/all_rb_ps_dewptc'      'rbZt'   'rband' 'pre_sal' ''       'ge' -100 }
%        { 'AzAveragedData/hist_all_dewptc_<CASE>.h5' '/dewptc' 'wtmean'  0.0  '/all_rb_s_dewptc'       'rbZt'   'rband' 'sal'     ''       'ge' -100 }
%      }
%      'DIAGS/hist_meas_az_dewptc_<CASE>.h5'
%    }
%
%    % vertially integrated liquid measurements
%    {
%      'Vert Liquid'
%      {
%        { 'AzAveragedData/hist_all_vint_liq_<CASE>.h5' '/vint_liq' 'wtmean'  0.0  '/all_ps_vint_liq'      'Rb_t'    ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_vint_liq_<CASE>.h5' '/vint_liq' 'wtmean'  0.0  '/all_s_vint_liq'       'Rb_t'    ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_vint_liq_<CASE>.h5' '/vint_liq' 'wtmean'  0.0  '/all_vint_liq_ts'      'Rb_T' ''      ''        ''       'ge' 0   }
%      }
%      'DIAGS/hist_meas_az_vint_liq_<CASE>.h5'
%    }
%
%    % vertially integrated ice measurements
%    {
%      'Vert Liquid'
%      {
%        { 'AzAveragedData/hist_all_vint_ice_<CASE>.h5' '/vint_ice' 'wtmean'  0.0  '/all_ps_vint_ice'      'Rb_t'    ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_vint_ice_<CASE>.h5' '/vint_ice' 'wtmean'  0.0  '/all_s_vint_ice'       'Rb_t'    ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_vint_ice_<CASE>.h5' '/vint_ice' 'wtmean'  0.0  '/all_vint_ice_ts'      'Rb_T' ''      ''        ''       'ge' 0   }
%      }
%      'DIAGS/hist_meas_az_vint_ice_<CASE>.h5'
%    }
%
%    % vertially integrated condensate measurements
%    {
%      'Vert Cond'
%      {
%        { 'AzAveragedData/hist_all_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/all_ps_vint_cond'      'Rb_t'    ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/all_s_vint_cond'       'Rb_t'    ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_vint_cond_<CASE>.h5' '/vint_cond' 'wtmean'  0.0  '/all_vint_cond_ts'      'Rb_T' ''      ''        ''       'ge' 0   }
%      }
%      'DIAGS/hist_meas_az_vint_cond_<CASE>.h5'
%    }
%
%    % vertially integrated vapor measurements
%    {
%      'Vert Vapor'
%      {
%        { 'AzAveragedData/hist_all_vint_vapor_<CASE>.h5' '/vint_vapor' 'wtmean'  0.0  '/all_ps_vint_vapor'      'Rb_t'    ''      'pre_sal' ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_vint_vapor_<CASE>.h5' '/vint_vapor' 'wtmean'  0.0  '/all_s_vint_vapor'       'Rb_t'    ''      'sal'     ''       'ge' 0   }
%        { 'AzAveragedData/hist_all_vint_vapor_<CASE>.h5' '/vint_vapor' 'wtmean'  0.0  '/all_vint_vapor_ts'      'Rb_T' ''      ''        ''       'ge' 0   }
%      }
%      'DIAGS/hist_meas_az_vint_vapor_<CASE>.h5'
%    }
%
%    % dust measurements
%    {
%      'Dust In Hydrometeors Azavg'
%      {
%        % Region: lead
%        { 'AzAveragedData/hist_all_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/all_dust_cloud_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/all_ps_dust_cloud'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/all_s_dust_cloud'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/all_whole_dust_cloud_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/all_whole_ps_dust_cloud'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/all_whole_s_dust_cloud'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/all_whole_i_dust_cloud'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/all_whole_m_dust_cloud'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/all_whole_f_dust_cloud'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/all_dust_rain_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/all_ps_dust_rain'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/all_s_dust_rain'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/all_whole_dust_rain_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/all_whole_ps_dust_rain'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/all_whole_s_dust_rain'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/all_whole_i_dust_rain'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/all_whole_m_dust_rain'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/all_whole_f_dust_rain'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/all_dust_pris_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/all_ps_dust_pris'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/all_s_dust_pris'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/all_whole_dust_pris_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/all_whole_ps_dust_pris'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/all_whole_s_dust_pris'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/all_whole_i_dust_pris'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/all_whole_m_dust_pris'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/all_whole_f_dust_pris'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/all_dust_snow_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/all_ps_dust_snow'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/all_s_dust_snow'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/all_whole_dust_snow_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/all_whole_ps_dust_snow'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/all_whole_s_dust_snow'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/all_whole_i_dust_snow'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/all_whole_m_dust_snow'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/all_whole_f_dust_snow'    'rbZt' ''      'fina'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/all_dust_aggr_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/all_ps_dust_aggr'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/all_s_dust_aggr'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/all_whole_dust_aggr_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/all_whole_ps_dust_aggr'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/all_whole_s_dust_aggr'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/all_whole_i_dust_aggr'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/all_whole_m_dust_aggr'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/all_whole_f_dust_aggr'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/all_dust_graup_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/all_ps_dust_graup'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/all_s_dust_graup'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/all_whole_dust_graup_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/all_whole_ps_dust_graup'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/all_whole_s_dust_graup'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/all_whole_i_dust_graup'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/all_whole_m_dust_graup'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/all_whole_f_dust_graup'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/all_dust_hail_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/all_ps_dust_hail'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/all_s_dust_hail'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/all_whole_dust_hail_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/all_whole_ps_dust_hail'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/all_whole_s_dust_hail'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/all_whole_i_dust_hail'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/all_whole_m_dust_hail'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/all_whole_f_dust_hail'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/all_dust_hydro_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/all_ps_dust_hydro'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/all_s_dust_hydro'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/all_whole_ps_dust_hydro'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/all_whole_s_dust_hydro'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/all_whole_i_dust_hydro'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/all_whole_m_dust_hydro'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/all_whole_f_dust_hydro'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/tcond_dust_cloud_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/tcond_ps_dust_cloud'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/tcond_s_dust_cloud'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/tcond_whole_dust_cloud_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/tcond_whole_ps_dust_cloud'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/tcond_whole_s_dust_cloud'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/tcond_whole_i_dust_cloud'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/tcond_whole_m_dust_cloud'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/tcond_whole_f_dust_cloud'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/tcond_dust_rain_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/tcond_ps_dust_rain'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/tcond_s_dust_rain'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/tcond_whole_dust_rain_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/tcond_whole_ps_dust_rain'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/tcond_whole_s_dust_rain'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/tcond_whole_i_dust_rain'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/tcond_whole_m_dust_rain'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/tcond_whole_f_dust_rain'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/tcond_dust_pris_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/tcond_ps_dust_pris'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/tcond_s_dust_pris'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/tcond_whole_dust_pris_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/tcond_whole_ps_dust_pris'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/tcond_whole_s_dust_pris'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/tcond_whole_i_dust_pris'    'rbZt' ''      'inti'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/tcond_whole_m_dust_pris'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/tcond_whole_f_dust_pris'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/tcond_dust_snow_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/tcond_ps_dust_snow'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/tcond_s_dust_snow'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/tcond_whole_dust_snow_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/tcond_whole_ps_dust_snow'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/tcond_whole_s_dust_snow'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/tcond_whole_i_dust_snow'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/tcond_whole_m_dust_snow'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/tcond_whole_f_dust_snow'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/tcond_dust_aggr_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/tcond_ps_dust_aggr'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/tcond_s_dust_aggr'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/tcond_whole_dust_aggr_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/tcond_whole_ps_dust_aggr'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/tcond_whole_s_dust_aggr'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/tcond_whole_i_dust_aggr'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/tcond_whole_m_dust_aggr'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/tcond_whole_f_dust_aggr'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/tcond_dust_graup_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/tcond_ps_dust_graup'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/tcond_s_dust_graup'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/tcond_whole_dust_graup_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/tcond_whole_ps_dust_graup'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/tcond_whole_s_dust_graup'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/tcond_whole_i_dust_graup'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/tcond_whole_m_dust_graup'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/tcond_whole_f_dust_graup'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/tcond_dust_hail_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/tcond_ps_dust_hail'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/tcond_s_dust_hail'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/tcond_whole_dust_hail_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/tcond_whole_ps_dust_hail'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/tcond_whole_s_dust_hail'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/tcond_whole_i_dust_hail'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/tcond_whole_m_dust_hail'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/tcond_whole_f_dust_hail'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/tcond_dust_hydro_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/tcond_ps_dust_hydro'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/tcond_s_dust_hydro'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/tcond_whole_dust_hydro_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/tcond_whole_ps_dust_hydro'   'rbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/tcond_whole_s_dust_hydro'    'rbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/tcond_whole_i_dust_hydro'    'rbZt' ''      'init'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/tcond_whole_m_dust_hydro'    'rbZt' ''      'mid'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/tcond_whole_f_dust_hydro'    'rbZt' ''      'final'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/all_whole_dust_hydro_ts'   'rbZT' ''      ''        ''    'ge' 0.001 }
%      }
%      'DIAGS/hist_meas_az_dust_hydro_<CASE>.h5'
%    }
%
%    % dustifn measurements
%    {
%      'Dust As IFN In Hydrometeors Azavg'
%      {
%        % Region: all
%        { 'AzAveragedData/hist_all_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/all_dustifn_cloud_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/all_ps_dustifn_cloud'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/all_s_dustifn_cloud'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/all_core_dustifn_cloud_ts'    'rbZT' 'core'  ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/all_rb_dustifn_cloud_ts'      'rbZT' 'rband' ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/all_core_ps_dustifn_cloud'    'rbZt' 'core'  'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/all_core_s_dustifn_cloud'     'rbZt' 'core'  'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/all_rb_ps_dustifn_cloud'      'rbZt' 'rband' 'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/all_rb_s_dustifn_cloud'       'rbZt' 'rband' 'sal'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/all_dustifn_rain_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/all_ps_dustifn_rain'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/all_s_dustifn_rain'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/all_core_dustifn_rain_ts'    'rbZT' 'core'  ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/all_rb_dustifn_rain_ts'      'rbZT' 'rband' ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/all_core_ps_dustifn_rain'    'rbZt' 'core'  'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/all_core_s_dustifn_rain'     'rbZt' 'core'  'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/all_rb_ps_dustifn_rain'      'rbZt' 'rband' 'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/all_rb_s_dustifn_rain'       'rbZt' 'rband' 'sal'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/all_dustifn_pris_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/all_ps_dustifn_pris'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/all_s_dustifn_pris'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/all_core_dustifn_pris_ts'    'rbZT' 'core'  ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/all_rb_dustifn_pris_ts'      'rbZT' 'rband' ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/all_core_ps_dustifn_pris'    'rbZt' 'core'  'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/all_core_s_dustifn_pris'     'rbZt' 'core'  'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/all_rb_ps_dustifn_pris'      'rbZt' 'rband' 'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/all_rb_s_dustifn_pris'       'rbZt' 'rband' 'sal'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/all_dustifn_snow_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/all_ps_dustifn_snow'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/all_s_dustifn_snow'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/all_core_dustifn_snow_ts'    'rbZT' 'core'  ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/all_rb_dustifn_snow_ts'      'rbZT' 'rband' ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/all_core_ps_dustifn_snow'    'rbZt' 'core'  'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/all_core_s_dustifn_snow'     'rbZt' 'core'  'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/all_rb_ps_dustifn_snow'      'rbZt' 'rband' 'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/all_rb_s_dustifn_snow'       'rbZt' 'rband' 'sal'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/all_dustifn_aggr_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/all_ps_dustifn_aggr'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/all_s_dustifn_aggr'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/all_core_dustifn_aggr_ts'    'rbZT' 'core'  ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/all_rb_dustifn_aggr_ts'      'rbZT' 'rband' ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/all_core_ps_dustifn_aggr'    'rbZt' 'core'  'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/all_core_s_dustifn_aggr'     'rbZt' 'core'  'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/all_rb_ps_dustifn_aggr'      'rbZt' 'rband' 'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/all_rb_s_dustifn_aggr'       'rbZt' 'rband' 'sal'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/all_dustifn_graup_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/all_ps_dustifn_graup'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/all_s_dustifn_graup'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/all_core_dustifn_graup_ts'    'rbZT' 'core'  ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/all_rb_dustifn_graup_ts'      'rbZT' 'rband' ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/all_core_ps_dustifn_graup'    'rbZt' 'core'  'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/all_core_s_dustifn_graup'     'rbZt' 'core'  'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/all_rb_ps_dustifn_graup'      'rbZt' 'rband' 'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/all_rb_s_dustifn_graup'       'rbZt' 'rband' 'sal'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/all_dustifn_hail_ts'         'RbZT' ''      ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/all_ps_dustifn_hail'         'RbZt' ''      'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/all_s_dustifn_hail'          'RbZt' ''      'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/all_core_dustifn_hail_ts'    'rbZT' 'core'  ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/all_rb_dustifn_hail_ts'      'rbZT' 'rband' ''        ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/all_core_ps_dustifn_hail'    'rbZt' 'core'  'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/all_core_s_dustifn_hail'     'rbZt' 'core'  'sal'     ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/all_rb_ps_dustifn_hail'      'rbZt' 'rband' 'pre_sal' ''    'ge' 0.001 }
%        { 'AzAveragedData/hist_all_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/all_rb_s_dustifn_hail'       'rbZt' 'rband' 'sal'     ''    'ge' 0.001 }
%
%        { 'AzAveragedData/hist_all_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/all_dustifn_hydro_ts'         'RbZT' ''      ''        ''    'ge' 1e-6 }
%        { 'AzAveragedData/hist_all_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/all_ps_dustifn_hydro'         'RbZt' ''      'pre_sal' ''    'ge' 1e-6 }
%        { 'AzAveragedData/hist_all_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/all_s_dustifn_hydro'          'RbZt' ''      'sal'     ''    'ge' 1e-6 }
%        { 'AzAveragedData/hist_all_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/all_core_dustifn_hydro_ts'    'rbZT' 'core'  ''        ''    'ge' 1e-6 }
%        { 'AzAveragedData/hist_all_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/all_rb_dustifn_hydro_ts'      'rbZT' 'rband' ''        ''    'ge' 1e-6 }
%        { 'AzAveragedData/hist_all_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/all_core_ps_dustifn_hydro'    'rbZt' 'core'  'pre_sal' ''    'ge' 1e-6 }
%        { 'AzAveragedData/hist_all_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/all_core_s_dustifn_hydro'     'rbZt' 'core'  'sal'     ''    'ge' 1e-6 }
%        { 'AzAveragedData/hist_all_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/all_rb_ps_dustifn_hydro'      'rbZt' 'rband' 'pre_sal' ''    'ge' 1e-6 }
%        { 'AzAveragedData/hist_all_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/all_rb_s_dustifn_hydro'       'rbZt' 'rband' 'sal'     ''    'ge' 1e-6 }
%      }
%      'DIAGS/hist_meas_az_dustifn_hydro_<CASE>.h5'
%    }
%
%    % cloud condensation
%    {
%      'Cloud Condensation'
%     {
%        { 'AzAveragedData/hist_all_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/all_cloud_cond_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/all_ps_cloud_cond'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/all_s_cloud_cond'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/all_core_cloud_cond_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/all_rb_cloud_cond_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/all_core_ps_cloud_cond'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/all_core_s_cloud_cond'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/all_rb_ps_cloud_cond'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/all_rb_s_cloud_cond'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_cloud_cond_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_ps_cloud_cond'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_s_cloud_cond'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_core_cloud_cond_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_rb_cloud_cond_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_env_cloud_cond_ts'     'rbZT' 'env'   ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_core_ps_cloud_cond'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_core_s_cloud_cond'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_rb_ps_cloud_cond'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_rb_s_cloud_cond'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_env_ps_cloud_cond'     'rbZt' 'env'   'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_cloud_cond_<CASE>.h5' '/cloud_cond' 'wtmean'  0.0  '/lead_env_s_cloud_cond'      'rbZt' 'env'   'sal'     ''       'ge' 0.1 }
%
%        { 'AzAveragedData/hist_cloud_all_cld_cond_<CASE>.h5'     '/cld_cond' 'wtmean'  0.0  '/cloud_all_cld_cond_ts'           'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_cld_cond_<CASE>.h5'     '/cld_cond' 'wtmean'  0.0  '/cloud_all_whole_cld_cond_ts'     'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_cld_cond_<CASE>.h5'     '/cld_cond' 'wtmean'  0.0  '/cloud_all_whole_ps_cld_cond'     'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_cld_cond_<CASE>.h5'     '/cld_cond' 'wtmean'  0.0  '/cloud_all_whole_s_cld_cond'      'rbZt' '' 'sal'     '' 'gt' 0 }
%
%        { 'AzAveragedData/hist_not_cloud_all_cld_cond_<CASE>.h5' '/cld_cond' 'wtmean'  0.0  '/not_cloud_all_cld_cond_ts'       'RbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_cld_cond_<CASE>.h5' '/cld_cond' 'wtmean'  0.0  '/not_cloud_all_whole_cld_cond_ts' 'rbZT' '' ''        '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_cld_cond_<CASE>.h5' '/cld_cond' 'wtmean'  0.0  '/not_cloud_all_whole_ps_cld_cond' 'rbZt' '' 'pre_sal' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_cld_cond_<CASE>.h5' '/cld_cond' 'wtmean'  0.0  '/not_cloud_all_whole_s_cld_cond'  'rbZt' '' 'sal'     '' 'gt' 0 }
%      }
%      'DIAGS/hist_meas_az_cloud_cond_<CASE>.h5'
%    }
%
%    % cloud evaporation
%    {
%      'Cloud Evaporation'
%     {
%        { 'AzAveragedData/hist_all_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/all_cloud_evap_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/all_ps_cloud_evap'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/all_s_cloud_evap'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/all_core_cloud_evap_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/all_rb_cloud_evap_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/all_core_ps_cloud_evap'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/all_core_s_cloud_evap'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/all_rb_ps_cloud_evap'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/all_rb_s_cloud_evap'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_cloud_evap_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_ps_cloud_evap'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_s_cloud_evap'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_core_cloud_evap_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_rb_cloud_evap_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_env_cloud_evap_ts'     'rbZT' 'env'   ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_core_ps_cloud_evap'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_core_s_cloud_evap'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_rb_ps_cloud_evap'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_rb_s_cloud_evap'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_env_ps_cloud_evap'     'rbZt' 'env'   'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_cloud_evap_<CASE>.h5' '/cloud_evap' 'wtmean'  0.0  '/lead_env_s_cloud_evap'      'rbZt' 'env'   'sal'     ''       'le' -0.1 }
%
%        { 'AzAveragedData/hist_cloud_all_cld_evap_<CASE>.h5'     '/cld_evap' 'wtmean'  0.0  '/cloud_all_cld_evap_ts'           'RbZT' '' ''        '' 'lt' 0 }
%        { 'AzAveragedData/hist_cloud_all_cld_evap_<CASE>.h5'     '/cld_evap' 'wtmean'  0.0  '/cloud_all_whole_cld_evap_ts'     'rbZT' '' ''        '' 'lt' 0 }
%        { 'AzAveragedData/hist_cloud_all_cld_evap_<CASE>.h5'     '/cld_evap' 'wtmean'  0.0  '/cloud_all_whole_ps_cld_evap'     'rbZt' '' 'pre_sal' '' 'lt' 0 }
%        { 'AzAveragedData/hist_cloud_all_cld_evap_<CASE>.h5'     '/cld_evap' 'wtmean'  0.0  '/cloud_all_whole_s_cld_evap'      'rbZt' '' 'sal'     '' 'lt' 0 }
%
%        { 'AzAveragedData/hist_not_cloud_all_cld_evap_<CASE>.h5' '/cld_evap' 'wtmean'  0.0  '/not_cloud_all_cld_evap_ts'       'RbZT' '' ''        '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_cld_evap_<CASE>.h5' '/cld_evap' 'wtmean'  0.0  '/not_cloud_all_whole_cld_evap_ts' 'rbZT' '' ''        '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_cld_evap_<CASE>.h5' '/cld_evap' 'wtmean'  0.0  '/not_cloud_all_whole_ps_cld_evap' 'rbZt' '' 'pre_sal' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_cld_evap_<CASE>.h5' '/cld_evap' 'wtmean'  0.0  '/not_cloud_all_whole_s_cld_evap'  'rbZt' '' 'sal'     '' 'lt' 0 }
%      }
%      'DIAGS/hist_meas_az_cloud_evap_<CASE>.h5'
%    }
%
%    % rain condensation
%    {
%      'Rain Condensation'
%     {
%        { 'AzAveragedData/hist_all_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/all_rain_cond_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/all_ps_rain_cond'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/all_s_rain_cond'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/all_core_rain_cond_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/all_rb_rain_cond_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/all_core_ps_rain_cond'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/all_core_s_rain_cond'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/all_rb_ps_rain_cond'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/all_rb_s_rain_cond'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_rain_cond_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_ps_rain_cond'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_s_rain_cond'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_core_rain_cond_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_rb_rain_cond_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_env_rain_cond_ts'     'rbZT' 'env'   ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_core_ps_rain_cond'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_core_s_rain_cond'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_rb_ps_rain_cond'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_rb_s_rain_cond'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_env_ps_rain_cond'     'rbZt' 'env'   'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rain_cond_<CASE>.h5' '/rain_cond' 'wtmean'  0.0  '/lead_env_s_rain_cond'      'rbZt' 'env'   'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_rain_cond_<CASE>.h5'
%    }
%
%    % rain evaporation
%    {
%      'Rain Evaporation'
%     {
%        { 'AzAveragedData/hist_all_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/all_rain_evap_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/all_ps_rain_evap'         'RbZt' ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/all_s_rain_evap'          'RbZt' ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/all_core_rain_evap_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/all_rb_rain_evap_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/all_core_ps_rain_evap'    'rbZt' 'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/all_core_s_rain_evap'     'rbZt' 'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/all_rb_ps_rain_evap'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/all_rb_s_rain_evap'       'rbZt' 'rband' 'sal'     ''       'le' -0.01 }
%
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_rain_evap_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_ps_rain_evap'         'RbZt' ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_s_rain_evap'          'RbZt' ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_core_rain_evap_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_rb_rain_evap_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_env_rain_evap_ts'     'rbZT' 'env'   ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_core_ps_rain_evap'    'rbZt' 'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_core_s_rain_evap'     'rbZt' 'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_rb_ps_rain_evap'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_rb_s_rain_evap'       'rbZt' 'rband' 'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_env_ps_rain_evap'     'rbZt' 'env'   'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rain_evap_<CASE>.h5' '/rain_evap' 'wtmean'  0.0  '/lead_env_s_rain_evap'      'rbZt' 'env'   'sal'     ''       'le' -0.01 }
%      }
%      'DIAGS/hist_meas_az_rain_evap_<CASE>.h5'
%    }
%
%    % pristine deposition
%    {
%      'Pristine Deposition'
%     {
%        { 'AzAveragedData/hist_all_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/all_pris_dep_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/all_ps_pris_dep'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/all_s_pris_dep'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/all_core_pris_dep_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/all_rb_pris_dep_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/all_core_ps_pris_dep'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/all_core_s_pris_dep'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/all_rb_ps_pris_dep'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/all_rb_s_pris_dep'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%
%        { 'AzAveragedData/hist_lead_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/lead_pris_dep_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/lead_ps_pris_dep'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/lead_s_pris_dep'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/lead_core_pris_dep_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/lead_rb_pris_dep_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/lead_core_ps_pris_dep'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/lead_core_s_pris_dep'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/lead_rb_ps_pris_dep'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_pris_dep_<CASE>.h5' '/pris_dep' 'wtmean'  0.0  '/lead_rb_s_pris_dep'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%      }
%      'DIAGS/hist_meas_az_pris_dep_<CASE>.h5'
%    }
%
%    % pristine sublimation
%    {
%      'Pristine Sublimation'
%     {
%        { 'AzAveragedData/hist_all_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/all_pris_sub_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/all_ps_pris_sub'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/all_s_pris_sub'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/all_core_pris_sub_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/all_rb_pris_sub_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/all_core_ps_pris_sub'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/all_core_s_pris_sub'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/all_rb_ps_pris_sub'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/all_rb_s_pris_sub'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%
%        { 'AzAveragedData/hist_lead_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/lead_pris_sub_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/lead_ps_pris_sub'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/lead_s_pris_sub'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/lead_core_pris_sub_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/lead_rb_pris_sub_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/lead_core_ps_pris_sub'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/lead_core_s_pris_sub'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/lead_rb_ps_pris_sub'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_pris_sub_<CASE>.h5' '/pris_sub' 'wtmean'  0.0  '/lead_rb_s_pris_sub'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%      }
%      'DIAGS/hist_meas_az_pris_sub_<CASE>.h5'
%    }
%
%
%    % snow deposition
%    {
%      'Snow Deposition'
%     {
%        { 'AzAveragedData/hist_all_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/all_snow_dep_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/all_ps_snow_dep'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/all_s_snow_dep'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/all_core_snow_dep_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/all_rb_snow_dep_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/all_core_ps_snow_dep'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/all_core_s_snow_dep'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/all_rb_ps_snow_dep'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/all_rb_s_snow_dep'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%
%        { 'AzAveragedData/hist_lead_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/lead_snow_dep_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/lead_ps_snow_dep'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/lead_s_snow_dep'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/lead_core_snow_dep_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/lead_rb_snow_dep_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/lead_core_ps_snow_dep'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/lead_core_s_snow_dep'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/lead_rb_ps_snow_dep'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_snow_dep_<CASE>.h5' '/snow_dep' 'wtmean'  0.0  '/lead_rb_s_snow_dep'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%      }
%      'DIAGS/hist_meas_az_snow_dep_<CASE>.h5'
%    }
%
%    % snow sublimation
%    {
%      'Snow Sublimation'
%     {
%        { 'AzAveragedData/hist_all_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/all_snow_sub_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/all_ps_snow_sub'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/all_s_snow_sub'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/all_core_snow_sub_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/all_rb_snow_sub_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/all_core_ps_snow_sub'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/all_core_s_snow_sub'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/all_rb_ps_snow_sub'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/all_rb_s_snow_sub'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%
%        { 'AzAveragedData/hist_lead_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/lead_snow_sub_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/lead_ps_snow_sub'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/lead_s_snow_sub'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/lead_core_snow_sub_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/lead_rb_snow_sub_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/lead_core_ps_snow_sub'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/lead_core_s_snow_sub'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/lead_rb_ps_snow_sub'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_snow_sub_<CASE>.h5' '/snow_sub' 'wtmean'  0.0  '/lead_rb_s_snow_sub'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%      }
%      'DIAGS/hist_meas_az_snow_sub_<CASE>.h5'
%    }
%
%
%    % aggregate deposition
%    {
%      'Aggregate Deposition'
%     {
%        { 'AzAveragedData/hist_all_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/all_aggr_dep_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/all_ps_aggr_dep'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/all_s_aggr_dep'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/all_core_aggr_dep_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/all_rb_aggr_dep_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/all_core_ps_aggr_dep'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/all_core_s_aggr_dep'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/all_rb_ps_aggr_dep'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/all_rb_s_aggr_dep'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%
%        { 'AzAveragedData/hist_lead_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/lead_aggr_dep_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/lead_ps_aggr_dep'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/lead_s_aggr_dep'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/lead_core_aggr_dep_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/lead_rb_aggr_dep_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/lead_core_ps_aggr_dep'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/lead_core_s_aggr_dep'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/lead_rb_ps_aggr_dep'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_aggr_dep_<CASE>.h5' '/aggr_dep' 'wtmean'  0.0  '/lead_rb_s_aggr_dep'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%      }
%      'DIAGS/hist_meas_az_aggr_dep_<CASE>.h5'
%    }
%
%    % aggregate sublimation
%    {
%      'Aggregate Sublimation'
%     {
%        { 'AzAveragedData/hist_all_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/all_aggr_sub_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/all_ps_aggr_sub'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/all_s_aggr_sub'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/all_core_aggr_sub_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/all_rb_aggr_sub_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/all_core_ps_aggr_sub'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/all_core_s_aggr_sub'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/all_rb_ps_aggr_sub'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/all_rb_s_aggr_sub'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%
%        { 'AzAveragedData/hist_lead_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/lead_aggr_sub_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/lead_ps_aggr_sub'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/lead_s_aggr_sub'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/lead_core_aggr_sub_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/lead_rb_aggr_sub_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/lead_core_ps_aggr_sub'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/lead_core_s_aggr_sub'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/lead_rb_ps_aggr_sub'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_aggr_sub_<CASE>.h5' '/aggr_sub' 'wtmean'  0.0  '/lead_rb_s_aggr_sub'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%      }
%      'DIAGS/hist_meas_az_aggr_sub_<CASE>.h5'
%    }
%
%
%    % graupel deposition
%    {
%      'Graupel Deposition'
%     {
%        { 'AzAveragedData/hist_all_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/all_graup_dep_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/all_ps_graup_dep'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/all_s_graup_dep'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/all_core_graup_dep_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/all_rb_graup_dep_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/all_core_ps_graup_dep'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/all_core_s_graup_dep'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/all_rb_ps_graup_dep'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/all_rb_s_graup_dep'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%
%        { 'AzAveragedData/hist_lead_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/lead_graup_dep_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/lead_ps_graup_dep'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/lead_s_graup_dep'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/lead_core_graup_dep_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/lead_rb_graup_dep_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/lead_core_ps_graup_dep'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/lead_core_s_graup_dep'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/lead_rb_ps_graup_dep'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_graup_dep_<CASE>.h5' '/graup_dep' 'wtmean'  0.0  '/lead_rb_s_graup_dep'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%      }
%      'DIAGS/hist_meas_az_graup_dep_<CASE>.h5'
%    }
%
%    % graupel sublimation
%    {
%      'Graupel Sublimation'
%     {
%        { 'AzAveragedData/hist_all_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/all_graup_sub_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/all_ps_graup_sub'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/all_s_graup_sub'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/all_core_graup_sub_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/all_rb_graup_sub_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/all_core_ps_graup_sub'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/all_core_s_graup_sub'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/all_rb_ps_graup_sub'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/all_rb_s_graup_sub'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%
%        { 'AzAveragedData/hist_lead_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/lead_graup_sub_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/lead_ps_graup_sub'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/lead_s_graup_sub'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/lead_core_graup_sub_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/lead_rb_graup_sub_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/lead_core_ps_graup_sub'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/lead_core_s_graup_sub'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/lead_rb_ps_graup_sub'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_graup_sub_<CASE>.h5' '/graup_sub' 'wtmean'  0.0  '/lead_rb_s_graup_sub'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%      }
%      'DIAGS/hist_meas_az_graup_sub_<CASE>.h5'
%    }
%
%
%    % hail deposition
%    {
%      'Hail Deposition'
%     {
%        { 'AzAveragedData/hist_all_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/all_hail_dep_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/all_ps_hail_dep'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/all_s_hail_dep'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/all_core_hail_dep_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/all_rb_hail_dep_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/all_core_ps_hail_dep'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/all_core_s_hail_dep'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/all_rb_ps_hail_dep'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_all_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/all_rb_s_hail_dep'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%
%        { 'AzAveragedData/hist_lead_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/lead_hail_dep_ts'         'RbZT' ''      ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/lead_ps_hail_dep'         'RbZt' ''      'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/lead_s_hail_dep'          'RbZt' ''      'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/lead_core_hail_dep_ts'    'rbZT' 'core'  ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/lead_rb_hail_dep_ts'      'rbZT' 'rband' ''        ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/lead_core_ps_hail_dep'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/lead_core_s_hail_dep'     'rbZt' 'core'  'sal'     ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/lead_rb_ps_hail_dep'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.1 }
%        { 'AzAveragedData/hist_lead_hail_dep_<CASE>.h5' '/hail_dep' 'wtmean'  0.0  '/lead_rb_s_hail_dep'       'rbZt' 'rband' 'sal'     ''       'ge' 0.1 }
%      }
%      'DIAGS/hist_meas_az_hail_dep_<CASE>.h5'
%    }
%
%    % hail sublimation
%    {
%      'Hail Sublimation'
%     {
%        { 'AzAveragedData/hist_all_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/all_hail_sub_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/all_ps_hail_sub'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/all_s_hail_sub'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/all_core_hail_sub_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/all_rb_hail_sub_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/all_core_ps_hail_sub'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/all_core_s_hail_sub'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/all_rb_ps_hail_sub'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_all_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/all_rb_s_hail_sub'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%
%        { 'AzAveragedData/hist_lead_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/lead_hail_sub_ts'         'RbZT' ''      ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/lead_ps_hail_sub'         'RbZt' ''      'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/lead_s_hail_sub'          'RbZt' ''      'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/lead_core_hail_sub_ts'    'rbZT' 'core'  ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/lead_rb_hail_sub_ts'      'rbZT' 'rband' ''        ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/lead_core_ps_hail_sub'    'rbZt' 'core'  'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/lead_core_s_hail_sub'     'rbZt' 'core'  'sal'     ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/lead_rb_ps_hail_sub'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.1 }
%        { 'AzAveragedData/hist_lead_hail_sub_<CASE>.h5' '/hail_sub' 'wtmean'  0.0  '/lead_rb_s_hail_sub'       'rbZt' 'rband' 'sal'     ''       'le' -0.1 }
%      }
%      'DIAGS/hist_meas_az_hail_sub_<CASE>.h5'
%    }
%
%
%    % cooling via latent heat of freezing
%    {
%      'LHF Cooling'
%     {
%        { 'AzAveragedData/hist_all_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/all_lhf_cool_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/all_ps_lhf_cool'         'RbZt' ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/all_s_lhf_cool'          'RbZt' ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/all_core_lhf_cool_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/all_core_ps_lhf_cool'    'rbZt' 'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/all_core_s_lhf_cool'     'rbZt' 'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/all_rb_lhf_cool_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/all_rb_ps_lhf_cool'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/all_rb_s_lhf_cool'       'rbZt' 'rband' 'sal'     ''       'le' -0.01 }
%
%        { 'AzAveragedData/hist_lead_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/lead_lhf_cool_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/lead_ps_lhf_cool'         'RbZt' ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/lead_s_lhf_cool'          'RbZt' ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/lead_core_lhf_cool_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/lead_core_ps_lhf_cool'    'rbZt' 'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/lead_core_s_lhf_cool'     'rbZt' 'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/lead_rb_lhf_cool_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/lead_rb_ps_lhf_cool'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhf_cool_<CASE>.h5' '/lhf_cool' 'wtmean'  0.0  '/lead_rb_s_lhf_cool'       'rbZt' 'rband' 'sal'     ''       'le' -0.01 }
%      }
%      'DIAGS/hist_meas_az_lhf_cool_<CASE>.h5'
%    }
%
%    % heating via latent heat of freezing
%    {
%      'LHF Heating'
%     {
%        { 'AzAveragedData/hist_all_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/all_lhf_heat_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/all_ps_lhf_heat'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/all_s_lhf_heat'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/all_core_lhf_heat_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/all_core_ps_lhf_heat'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/all_core_s_lhf_heat'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/all_rb_lhf_heat_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/all_rb_ps_lhf_heat'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/all_rb_s_lhf_heat'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/lead_lhf_heat_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/lead_ps_lhf_heat'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/lead_s_lhf_heat'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/lead_core_lhf_heat_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/lead_core_ps_lhf_heat'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/lead_core_s_lhf_heat'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/lead_rb_lhf_heat_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/lead_rb_ps_lhf_heat'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhf_heat_<CASE>.h5' '/lhf_heat' 'wtmean'  0.0  '/lead_rb_s_lhf_heat'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_lhf_heat_<CASE>.h5'
%    }
%
%    % cooling via latent heat of vaporization
%    {
%      'LHV Cooling'
%     {
%        { 'AzAveragedData/hist_all_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/all_lhv_cool_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/all_ps_lhv_cool'         'RbZt' ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/all_s_lhv_cool'          'RbZt' ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/all_core_lhv_cool_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/all_core_ps_lhv_cool'    'rbZt' 'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/all_core_s_lhv_cool'     'rbZt' 'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/all_rb_lhv_cool_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/all_rb_ps_lhv_cool'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_all_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/all_rb_s_lhv_cool'       'rbZt' 'rband' 'sal'     ''       'le' -0.01 }
%
%        { 'AzAveragedData/hist_lead_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/lead_lhv_cool_ts'         'RbZT' ''      ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/lead_ps_lhv_cool'         'RbZt' ''      'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/lead_s_lhv_cool'          'RbZt' ''      'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/lead_core_lhv_cool_ts'    'rbZT' 'core'  ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/lead_core_ps_lhv_cool'    'rbZt' 'core'  'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/lead_core_s_lhv_cool'     'rbZt' 'core'  'sal'     ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/lead_rb_lhv_cool_ts'      'rbZT' 'rband' ''        ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/lead_rb_ps_lhv_cool'      'rbZt' 'rband' 'pre_sal' ''       'le' -0.01 }
%        { 'AzAveragedData/hist_lead_lhv_cool_<CASE>.h5' '/lhv_cool' 'wtmean'  0.0  '/lead_rb_s_lhv_cool'       'rbZt' 'rband' 'sal'     ''       'le' -0.01 }
%      }
%      'DIAGS/hist_meas_az_lhv_cool_<CASE>.h5'
%    }
%
%    % heating via latent heat of vaporization
%    {
%      'LHV Heating'
%     {
%        { 'AzAveragedData/hist_all_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/all_lhv_heat_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/all_ps_lhv_heat'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/all_s_lhv_heat'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/all_core_lhv_heat_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/all_core_ps_lhv_heat'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/all_core_s_lhv_heat'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/all_rb_lhv_heat_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/all_rb_ps_lhv_heat'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/all_rb_s_lhv_heat'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/lead_lhv_heat_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/lead_ps_lhv_heat'         'RbZt' ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/lead_s_lhv_heat'          'RbZt' ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/lead_core_lhv_heat_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/lead_core_ps_lhv_heat'    'rbZt' 'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/lead_core_s_lhv_heat'     'rbZt' 'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/lead_rb_lhv_heat_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/lead_rb_ps_lhv_heat'      'rbZt' 'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_lhv_heat_<CASE>.h5' '/lhv_heat' 'wtmean'  0.0  '/lead_rb_s_lhv_heat'       'rbZt' 'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_lhv_heat_<CASE>.h5'
%    }
%
%    % liquid evaporation
%    {
%      'Liquid Evaporation'
%     {
%        { 'AzAveragedData/hist_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_liq_evap_ts'         'RbZT' ''      ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_ps_liq_evap'         'RbZt' ''      'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_s_liq_evap'          'RbZt' ''      'sal'     ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_core_liq_evap_ts'    'rbZT' 'core'  ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_core_ps_liq_evap'    'rbZt' 'core'  'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_core_s_liq_evap'     'rbZt' 'core'  'sal'     ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_rb_liq_evap_ts'      'rbZT' 'rband' ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_rb_ps_liq_evap'      'rbZt' 'rband' 'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/all_rb_s_liq_evap'       'rbZt' 'rband' 'sal'     ''       'lt' 0 }
%
%        { 'AzAveragedData/hist_lead_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/lead_liq_evap_ts'         'RbZT' ''      ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_lead_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/lead_ps_liq_evap'         'RbZt' ''      'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_lead_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/lead_s_liq_evap'          'RbZt' ''      'sal'     ''       'lt' 0 }
%        { 'AzAveragedData/hist_lead_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/lead_core_liq_evap_ts'    'rbZT' 'core'  ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_lead_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/lead_core_ps_liq_evap'    'rbZt' 'core'  'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_lead_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/lead_core_s_liq_evap'     'rbZt' 'core'  'sal'     ''       'lt' 0 }
%        { 'AzAveragedData/hist_lead_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/lead_rb_liq_evap_ts'      'rbZT' 'rband' ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_lead_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/lead_rb_ps_liq_evap'      'rbZt' 'rband' 'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_lead_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/lead_rb_s_liq_evap'       'rbZt' 'rband' 'sal'     ''       'lt' 0 }
%
%        { 'AzAveragedData/hist_cloud_all_liq_evap_<CASE>.h5'     '/liq_evap' 'wtmean'  0.0  '/cloud_all_liq_evap_ts'           'RbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/not_cloud_all_liq_evap_ts'       'RbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_tcond_all_liq_evap_<CASE>.h5'     '/liq_evap' 'wtmean'  0.0  '/tcond_all_liq_evap_ts'           'RbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/not_tcond_all_liq_evap_ts'       'RbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_cloud_all_liq_evap_<CASE>.h5'     '/liq_evap' 'wtmean'  0.0  '/cloud_all_whole_liq_evap_ts'     'rbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/not_cloud_all_whole_liq_evap_ts' 'rbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_tcond_all_liq_evap_<CASE>.h5'     '/liq_evap' 'wtmean'  0.0  '/tcond_all_whole_liq_evap_ts'     'rbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/not_tcond_all_whole_liq_evap_ts' 'rbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_cloud_all_liq_evap_<CASE>.h5'     '/liq_evap' 'wtmean'  0.0  '/cloud_all_whole_ps_liq_evap'     'rbZt' 'pre_sal' '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/not_cloud_all_whole_ps_liq_evap' 'rbZt' 'pre_sal' '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_tcond_all_liq_evap_<CASE>.h5'     '/liq_evap' 'wtmean'  0.0  '/tcond_all_whole_ps_liq_evap'     'rbZt' 'pre_sal' '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/not_tcond_all_whole_ps_liq_evap' 'rbZt' 'pre_sal' '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_cloud_all_liq_evap_<CASE>.h5'     '/liq_evap' 'wtmean'  0.0  '/cloud_all_whole_s_liq_evap'      'rbZt' 'sal'     '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/not_cloud_all_whole_s_liq_evap'  'rbZt' 'sal'     '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_tcond_all_liq_evap_<CASE>.h5'     '/liq_evap' 'wtmean'  0.0  '/tcond_all_whole_s_liq_evap'      'rbZt' 'sal'     '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_liq_evap_<CASE>.h5' '/liq_evap' 'wtmean'  0.0  '/not_tcond_all_whole_s_liq_evap'  'rbZt' 'sal'     '' '' 'lt' 0 }
%      }
%      'DIAGS/hist_meas_az_liq_evap_<CASE>.h5'
%    }
%
%    % Liquid condensation
%    {
%      'Liquid Condensation'
%     {
%        { 'AzAveragedData/hist_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_liq_cond_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_ps_liq_cond'         'RbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_s_liq_cond'          'RbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_core_liq_cond_ts'    'rbZT' 'core'  ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_core_ps_liq_cond'    'rbZt' 'core'  'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_core_s_liq_cond'     'rbZt' 'core'  'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_rb_liq_cond_ts'      'rbZT' 'rband' ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_rb_ps_liq_cond'      'rbZt' 'rband' 'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/all_rb_s_liq_cond'       'rbZt' 'rband' 'sal'     ''       'gt' 0 }
%
%        { 'AzAveragedData/hist_lead_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/lead_liq_cond_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_lead_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/lead_ps_liq_cond'         'RbZt' ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_lead_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/lead_s_liq_cond'          'RbZt' ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_lead_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/lead_core_liq_cond_ts'    'rbZT' 'core'  ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_lead_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/lead_core_ps_liq_cond'    'rbZt' 'core'  'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_lead_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/lead_core_s_liq_cond'     'rbZt' 'core'  'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_lead_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/lead_rb_liq_cond_ts'      'rbZT' 'rband' ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_lead_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/lead_rb_ps_liq_cond'      'rbZt' 'rband' 'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_lead_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/lead_rb_s_liq_cond'       'rbZt' 'rband' 'sal'     ''       'gt' 0 }
%
%        { 'AzAveragedData/hist_cloud_all_liq_cond_<CASE>.h5'     '/liq_cond' 'wtmean'  0.0  '/cloud_all_liq_cond_ts'           'RbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/not_cloud_all_liq_cond_ts'       'RbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_liq_cond_<CASE>.h5'     '/liq_cond' 'wtmean'  0.0  '/tcond_all_liq_cond_ts'           'RbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/not_tcond_all_liq_cond_ts'       'RbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_liq_cond_<CASE>.h5'     '/liq_cond' 'wtmean'  0.0  '/cloud_all_whole_liq_cond_ts'     'rbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/not_cloud_all_whole_liq_cond_ts' 'rbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_liq_cond_<CASE>.h5'     '/liq_cond' 'wtmean'  0.0  '/tcond_all_liq_whole_cond_ts'     'rbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/not_tcond_all_whole_liq_cond_ts' 'rbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_liq_cond_<CASE>.h5'     '/liq_cond' 'wtmean'  0.0  '/cloud_all_whole_ps_liq_cond'     'rbZt' 'pre_sal' '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/not_cloud_all_whole_ps_liq_cond' 'rbZt' 'pre_sal' '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_liq_cond_<CASE>.h5'     '/liq_cond' 'wtmean'  0.0  '/tcond_all_whole_ps_liq_cond'     'rbZt' 'pre_sal' '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/not_tcond_all_whole_ps_liq_cond' 'rbZt' 'pre_sal' '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_liq_cond_<CASE>.h5'     '/liq_cond' 'wtmean'  0.0  '/cloud_all_whole_s_liq_cond'      'rbZt' 'sal'     '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/not_cloud_all_whole_s_liq_cond'  'rbZt' 'sal'     '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_liq_cond_<CASE>.h5'     '/liq_cond' 'wtmean'  0.0  '/tcond_all_whole_s_liq_cond'      'rbZt' 'sal'     '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_liq_cond_<CASE>.h5' '/liq_cond' 'wtmean'  0.0  '/not_tcond_all_whole_s_liq_cond'  'rbZt' 'sal'     '' '' 'gt' 0 }
%      }
%      'DIAGS/hist_meas_az_liq_cond_<CASE>.h5'
%    }
%
%    % ice deposition
%    {
%      'Ice Deposition'
%     {
%        { 'AzAveragedData/hist_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_ice_dep_ts'         'RbZT' ''      ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_ps_ice_dep'         'RbZt'    ''      'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_s_ice_dep'          'RbZt'    ''      'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_core_ice_dep_ts'    'rbZT' 'core'  ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_rb_ice_dep_ts'      'rbZT' 'rband' ''        ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_core_ps_ice_dep'    'rbZt'    'core'  'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_core_s_ice_dep'     'rbZt'    'core'  'sal'     ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_rb_ps_ice_dep'      'rbZt'    'rband' 'pre_sal' ''       'gt' 0 }
%        { 'AzAveragedData/hist_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/all_rb_s_ice_dep'       'rbZt'    'rband' 'sal'     ''       'gt' 0 }
%
%        { 'AzAveragedData/hist_cloud_all_ice_dep_<CASE>.h5'     '/ice_dep' 'wtmean'  0.0  '/cloud_all_ice_dep_ts'           'RbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/not_cloud_all_ice_dep_ts'       'RbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ice_dep_<CASE>.h5'     '/ice_dep' 'wtmean'  0.0  '/tcond_all_ice_dep_ts'           'RbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/not_tcond_all_ice_dep_ts'       'RbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ice_dep_<CASE>.h5'     '/ice_dep' 'wtmean'  0.0  '/cloud_all_whole_ice_dep_ts'     'rbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/not_cloud_whole_all_ice_dep_ts' 'rbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ice_dep_<CASE>.h5'     '/ice_dep' 'wtmean'  0.0  '/tcond_all_whole_ice_dep_ts'     'rbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/not_tcond_all_whole_ice_dep_ts' 'rbZT' ''        '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ice_dep_<CASE>.h5'     '/ice_dep' 'wtmean'  0.0  '/cloud_all_whole_ps_ice_dep'     'rbZt' 'pre_sal' '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/not_cloud_all_whole_ps_ice_dep' 'rbZt' 'pre_sal' '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ice_dep_<CASE>.h5'     '/ice_dep' 'wtmean'  0.0  '/tcond_all_whole_ps_ice_dep'     'rbZt' 'pre_sal' '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/not_tcond_all_whole_ps_ice_dep' 'rbZt' 'pre_sal' '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ice_dep_<CASE>.h5'     '/ice_dep' 'wtmean'  0.0  '/cloud_all_whole_s_ice_dep'      'rbZt' 'sal'     '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/not_cloud_all_whole_s_ice_dep'  'rbZt' 'sal'     '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ice_dep_<CASE>.h5'     '/ice_dep' 'wtmean'  0.0  '/tcond_all_whole_s_ice_dep'      'rbZt' 'sal'     '' '' 'gt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ice_dep_<CASE>.h5' '/ice_dep' 'wtmean'  0.0  '/not_tcond_all_whole_s_ice_dep'  'rbZt' 'sal'     '' '' 'gt' 0 }
%      }
%      'DIAGS/hist_meas_az_ice_dep_<CASE>.h5'
%    }
%
%    % ice sublimation
%    {
%      'Ice Sublimation'
%     {
%        { 'AzAveragedData/hist_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_ice_sub_ts'         'RbZT' ''      ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_ps_ice_sub'         'RbZt'    ''      'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_s_ice_sub'          'RbZt'    ''      'sal'     ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_core_ice_sub_ts'    'rbZT' 'core'  ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_rb_ice_sub_ts'      'rbZT' 'rband' ''        ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_core_ps_ice_sub'    'rbZt'    'core'  'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_core_s_ice_sub'     'rbZt'    'core'  'sal'     ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_rb_ps_ice_sub'      'rbZt'    'rband' 'pre_sal' ''       'lt' 0 }
%        { 'AzAveragedData/hist_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/all_rb_s_ice_sub'       'rbZt'    'rband' 'sal'     ''       'lt' 0 }
%
%        { 'AzAveragedData/hist_cloud_all_ice_sub_<CASE>.h5'     '/ice_sub' 'wtmean'  0.0  '/cloud_all_ice_sub_ts'           'RbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/not_cloud_all_ice_sub_ts'       'RbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ice_sub_<CASE>.h5'     '/ice_sub' 'wtmean'  0.0  '/tcond_all_ice_sub_ts'           'RbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/not_tcond_all_ice_sub_ts'       'RbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ice_sub_<CASE>.h5'     '/ice_sub' 'wtmean'  0.0  '/cloud_all_whole_ice_sub_ts'     'rbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/not_cloud_all_whole_ice_sub_ts' 'rbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ice_sub_<CASE>.h5'     '/ice_sub' 'wtmean'  0.0  '/tcond_all_whole_ice_sub_ts'     'rbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/not_tcond_all_whole_ice_sub_ts' 'rbZT' ''        '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ice_sub_<CASE>.h5'     '/ice_sub' 'wtmean'  0.0  '/cloud_all_whole_ps_ice_sub'     'rbZt' 'pre_sal' '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/not_cloud_all_whole_ps_ice_sub' 'rbZt' 'pre_sal' '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ice_sub_<CASE>.h5'     '/ice_sub' 'wtmean'  0.0  '/tcond_all_whole_ps_ice_sub'     'rbZt' 'pre_sal' '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/not_tcond_all_whole_ps_ice_sub' 'rbZt' 'pre_sal' '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_cloud_all_ice_sub_<CASE>.h5'     '/ice_sub' 'wtmean'  0.0  '/cloud_all_whole_s_ice_sub'      'rbZt' 'sal'     '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_cloud_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/not_cloud_all_whole_s_ice_sub'  'rbZt' 'sal'     '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_tcond_all_ice_sub_<CASE>.h5'     '/ice_sub' 'wtmean'  0.0  '/tcond_all_whole_s_ice_sub'      'rbZt' 'sal'     '' '' 'lt' 0 }
%        { 'AzAveragedData/hist_not_tcond_all_ice_sub_<CASE>.h5' '/ice_sub' 'wtmean'  0.0  '/not_tcond_all_whole_s_ice_sub'  'rbZt' 'sal'     '' '' 'lt' 0 }
%      }
%      'DIAGS/hist_meas_az_ice_sub_<CASE>.h5'
%    }
%
%    % ice melting
%    {
%      'Ice Melting'
%     {
%        { 'AzAveragedData/hist_all_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/all_ice_melt_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/all_ps_ice_melt'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/all_s_ice_melt'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/all_core_ice_melt_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/all_rb_ice_melt_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/all_core_ps_ice_melt'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/all_core_s_ice_melt'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/all_rb_ps_ice_melt'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_ice_melt_<CASE>.h5' '/ice_melt' 'wtmean'  0.0  '/all_rb_s_ice_melt'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_ice_melt_<CASE>.h5'
%    }
%
%    % riming, cloud to ice
%    {
%      'Cloud Riming'
%     {
%        { 'AzAveragedData/hist_all_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/all_cloud_rime_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/all_ps_cloud_rime'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/all_s_cloud_rime'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/all_core_cloud_rime_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/all_rb_cloud_rime_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/all_core_ps_cloud_rime'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/all_core_s_cloud_rime'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/all_rb_ps_cloud_rime'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_cloud_rime_<CASE>.h5' '/cloud_rime' 'wtmean'  0.0  '/all_rb_s_cloud_rime'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_cloud_rime_<CASE>.h5'
%    }
%
%    % rain to ice
%    {
%      'Rain to Ice'
%     {
%        { 'AzAveragedData/hist_all_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/all_rain2ice_ts'         'RbZT' ''      ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/all_ps_rain2ice'         'RbZt'    ''      'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/all_s_rain2ice'          'RbZt'    ''      'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/all_core_rain2ice_ts'    'rbZT' 'core'  ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/all_rb_rain2ice_ts'      'rbZT' 'rband' ''        ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/all_core_ps_rain2ice'    'rbZt'    'core'  'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/all_core_s_rain2ice'     'rbZt'    'core'  'sal'     ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/all_rb_ps_rain2ice'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rain2ice_<CASE>.h5' '/rain2ice' 'wtmean'  0.0  '/all_rb_s_rain2ice'       'rbZt'    'rband' 'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_az_rain2ice_<CASE>.h5'
%    }
%
%    % relhum measurements
%    {
%      'RH Azavg'
%      {
%        { 'AzAveragedData/hist_all_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_relhum_ts'         'RbZT' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_ps_relhum'         'RbZt'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_s_relhum'          'RbZt'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_core_relhum_ts'    'rbZT' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_rb_relhum_ts'      'rbZT' 'rband' ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_core_ps_relhum'    'rbZt'    'core'  'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_core_s_relhum'     'rbZt'    'core'  'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_rb_ps_relhum'      'rbZt'    'rband' 'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/all_rb_s_relhum'       'rbZt'    'rband' 'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_relhum_<CASE>.h5'
%    }
%
%    % acceleration due to buoyancy
%    {
%      'Buoyancy acceleration'
%      {
%        { 'AzAveragedData/hist_all_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/all_core_buoy_pos_acc_ts'       'rbZT' 'core'  ''        ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_all_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/all_core_ps_buoy_pos_acc_xsect' 'RbZt' 'core'  'pre_sal' ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_all_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/all_core_s_buoy_pos_acc_xsect'  'RbZt' 'core'  'sal'     ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_all_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/all_core_ps_buoy_pos_acc'       'rbZt' 'core'  'pre_sal' ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_all_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/all_core_s_buoy_pos_acc'        'rbZt' 'core'  'sal'     ''  'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/all_core_buoy_neg_acc_ts'       'rbZT' 'core'  ''        ''  'le' -0.01 }
%        { 'AzAveragedData/hist_all_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/all_core_ps_buoy_neg_acc_xsect' 'RbZt' 'core'  'pre_sal' ''  'le' -0.01 }
%        { 'AzAveragedData/hist_all_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/all_core_s_buoy_neg_acc_xsect'  'RbZt' 'core'  'sal'     ''  'le' -0.01 }
%        { 'AzAveragedData/hist_all_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/all_core_ps_buoy_neg_acc'       'rbZt' 'core'  'pre_sal' ''  'le' -0.01 }
%        { 'AzAveragedData/hist_all_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/all_core_s_buoy_neg_acc'        'rbZt' 'core'  'sal'     ''  'le' -0.01 }
%
%        { 'AzAveragedData/hist_all_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/all_rb_buoy_pos_acc_ts'       'rbZT' 'rband'  ''        ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/all_rb_ps_buoy_pos_acc_xsect' 'RbZt' 'rband'  'pre_sal' ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/all_rb_s_buoy_pos_acc_xsect'  'RbZt' 'rband'  'sal'     ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/all_rb_ps_buoy_pos_acc'       'rbZt' 'rband'  'pre_sal' ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_all_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/all_rb_s_buoy_pos_acc'        'rbZt' 'rband'  'sal'     ''  'ge' 0.01 }
%
%        { 'AzAveragedData/hist_all_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/all_rb_buoy_neg_acc_ts'       'rbZT' 'rband'  ''        ''  'le' -0.01 }
%        { 'AzAveragedData/hist_all_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/all_rb_ps_buoy_neg_acc_xsect' 'RbZt' 'rband'  'pre_sal' ''  'le' -0.01 }
%        { 'AzAveragedData/hist_all_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/all_rb_s_buoy_neg_acc_xsect'  'RbZt' 'rband'  'sal'     ''  'le' -0.01 }
%        { 'AzAveragedData/hist_all_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/all_rb_ps_buoy_neg_acc'       'rbZt' 'rband'  'pre_sal' ''  'le' -0.01 }
%        { 'AzAveragedData/hist_all_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/all_rb_s_buoy_neg_acc'        'rbZt' 'rband'  'sal'     ''  'le' -0.01 }
%
%        { 'AzAveragedData/hist_lead_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/lead_core_buoy_pos_acc_ts'       'rbZT' 'core'  ''        ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/lead_core_ps_buoy_pos_acc_xsect' 'RbZt' 'core'  'pre_sal' ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/lead_core_s_buoy_pos_acc_xsect'  'RbZt' 'core'  'sal'     ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/lead_core_ps_buoy_pos_acc'       'rbZt' 'core'  'pre_sal' ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/lead_core_s_buoy_pos_acc'        'rbZt' 'core'  'sal'     ''  'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/lead_core_buoy_neg_acc_ts'       'rbZT' 'core'  ''        ''  'le' -0.01 }
%        { 'AzAveragedData/hist_lead_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/lead_core_ps_buoy_neg_acc_xsect' 'RbZt' 'core'  'pre_sal' ''  'le' -0.01 }
%        { 'AzAveragedData/hist_lead_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/lead_core_s_buoy_neg_acc_xsect'  'RbZt' 'core'  'sal'     ''  'le' -0.01 }
%        { 'AzAveragedData/hist_lead_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/lead_core_ps_buoy_neg_acc'       'rbZt' 'core'  'pre_sal' ''  'le' -0.01 }
%        { 'AzAveragedData/hist_lead_core_buoy_acc_<CASE>.h5' '/core_buoy_acc' 'wtmean'  0.0  '/lead_core_s_buoy_neg_acc'        'rbZt' 'core'  'sal'     ''  'le' -0.01 }
%
%        { 'AzAveragedData/hist_lead_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/lead_rb_buoy_pos_acc_ts'       'rbZT' 'rband'  ''        ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/lead_rb_ps_buoy_pos_acc_xsect' 'RbZt' 'rband'  'pre_sal' ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/lead_rb_s_buoy_pos_acc_xsect'  'RbZt' 'rband'  'sal'     ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/lead_rb_ps_buoy_pos_acc'       'rbZt' 'rband'  'pre_sal' ''  'ge' 0.01 }
%        { 'AzAveragedData/hist_lead_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/lead_rb_s_buoy_pos_acc'        'rbZt' 'rband'  'sal'     ''  'ge' 0.01 }
%
%        { 'AzAveragedData/hist_lead_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/lead_rb_buoy_neg_acc_ts'       'rbZT' 'rband'  ''        ''  'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/lead_rb_ps_buoy_neg_acc_xsect' 'RbZt' 'rband'  'pre_sal' ''  'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/lead_rb_s_buoy_neg_acc_xsect'  'RbZt' 'rband'  'sal'     ''  'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/lead_rb_ps_buoy_neg_acc'       'rbZt' 'rband'  'pre_sal' ''  'le' -0.01 }
%        { 'AzAveragedData/hist_lead_rband_buoy_acc_<CASE>.h5' '/rband_buoy_acc' 'wtmean'  0.0  '/lead_rb_s_buoy_neg_acc'        'rbZt' 'rband'  'sal'     ''  'le' -0.01 }
%      }
%      'DIAGS/hist_meas_az_buoy_acc_<CASE>.h5'
%    }
%
%    % zonal vertical wind shear
%    {
%      'Zonal vertical shear'
%      {
%        { 'AzAveragedData/hist_all_zonal_vshear_<CASE>.h5' '/vshear' 'wtmean'  0.0  '/all_vshear_ts'         'Rb_T' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_zonal_vshear_<CASE>.h5' '/vshear' 'wtmean'  0.0  '/all_ps_vshear'         'Rb_t'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_zonal_vshear_<CASE>.h5' '/vshear' 'wtmean'  0.0  '/all_s_vshear'          'Rb_t'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_zonal_vshear_<CASE>.h5' '/vshear' 'wtmean'  0.0  '/all_core_vshear_ts'    'rb_T' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_all_zonal_vshear_<CASE>.h5' '/vshear' 'wtmean'  0.0  '/all_rb_vshear_ts'      'rb_T' 'rband' ''        ''       'ge' 0 }
%
%        { 'AzAveragedData/hist_lead_zonal_vshear_<CASE>.h5' '/vshear' 'wtmean'  0.0  '/lead_vshear_ts'         'Rb_T' ''      ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_zonal_vshear_<CASE>.h5' '/vshear' 'wtmean'  0.0  '/lead_ps_vshear'         'Rb_t'    ''      'pre_sal' ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_zonal_vshear_<CASE>.h5' '/vshear' 'wtmean'  0.0  '/lead_s_vshear'          'Rb_t'    ''      'sal'     ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_zonal_vshear_<CASE>.h5' '/vshear' 'wtmean'  0.0  '/lead_core_vshear_ts'    'rb_T' 'core'  ''        ''       'ge' 0 }
%        { 'AzAveragedData/hist_lead_zonal_vshear_<CASE>.h5' '/vshear' 'wtmean'  0.0  '/lead_rb_vshear_ts'      'rb_T' 'rband' ''        ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_az_zonal_vshear_<CASE>.h5'
%    }
%
%
%
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    % Spatially averaged data
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    % Dust
%    {
%      'Dust Tsavg'
%      {
%        % sample regions in SAL
%        { 'TsAveragedData/hist_spath_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/spath_d1_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/spath_i_d1_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/spath_m_d1_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/spath_f_d1_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/spath_ps_d1_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d1_mass_<CASE>.h5' '/d1_mass' 'wtmean'  0.0  '/spath_s_d1_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/spath_d1_num_ts'   'b_ZT'   ''      ''        ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/spath_i_d1_num'    'b_Zt'   ''      'init'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/spath_m_d1_num'    'b_Zt'   ''      'mid'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/spath_f_d1_num'    'b_Zt'   ''      'final'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/spath_ps_d1_num'   'b_Zt'   ''      'pre_sal' ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_d1_num_<CASE>.h5'  '/d1_num'  'wtmean'  0.0  '/spath_s_d1_num'    'b_Zt'   ''      'sal'     ''       'ge' 10   }
%
%        { 'TsAveragedData/hist_spath_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/spath_d2_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/spath_i_d2_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/spath_m_d2_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/spath_f_d2_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/spath_ps_d2_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d2_mass_<CASE>.h5' '/d2_mass' 'wtmean'  0.0  '/spath_s_d2_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_d2_num_<CASE>.h5'  '/d2_num' 'wtmean'   0.0  '/spath_d2_num_ts'   'b_ZT'   ''      ''        ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_d2_num_<CASE>.h5'  '/d2_num' 'wtmean'   0.0  '/spath_i_d2_num'    'b_Zt'   ''      'init'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_d2_num_<CASE>.h5'  '/d2_num' 'wtmean'   0.0  '/spath_m_d2_num'    'b_Zt'   ''      'mid'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_d2_num_<CASE>.h5'  '/d2_num' 'wtmean'   0.0  '/spath_f_d2_num'    'b_Zt'   ''      'final'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_d2_num_<CASE>.h5'  '/d2_num' 'wtmean'   0.0  '/spath_ps_d2_num'   'b_Zt'   ''      'pre_sal' ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_d2_num_<CASE>.h5'  '/d2_num' 'wtmean'   0.0  '/spath_s_d2_num'    'b_Zt'   ''      'sal'     ''       'ge' 10   }
%
%        { 'TsAveragedData/hist_spath_tracer1_<CASE>.h5' '/tracer1' 'wtmean'  0.0  '/spath_tracer1_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer1_<CASE>.h5' '/tracer1' 'wtmean'  0.0  '/spath_i_tracer1'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer1_<CASE>.h5' '/tracer1' 'wtmean'  0.0  '/spath_m_tracer1'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer1_<CASE>.h5' '/tracer1' 'wtmean'  0.0  '/spath_f_tracer1'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer1_<CASE>.h5' '/tracer1' 'wtmean'  0.0  '/spath_ps_tracer1'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer1_<CASE>.h5' '/tracer1' 'wtmean'  0.0  '/spath_s_tracer1'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_tracer2_<CASE>.h5' '/tracer2' 'wtmean'  0.0  '/spath_tracer2_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer2_<CASE>.h5' '/tracer2' 'wtmean'  0.0  '/spath_i_tracer2'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer2_<CASE>.h5' '/tracer2' 'wtmean'  0.0  '/spath_m_tracer2'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer2_<CASE>.h5' '/tracer2' 'wtmean'  0.0  '/spath_f_tracer2'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer2_<CASE>.h5' '/tracer2' 'wtmean'  0.0  '/spath_ps_tracer2'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer2_<CASE>.h5' '/tracer2' 'wtmean'  0.0  '/spath_s_tracer2'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/spath_trdust1_diff_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/spath_i_trdust1_diff'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/spath_m_trdust1_diff'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/spath_f_trdust1_diff'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/spath_ps_trdust1_diff'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_trdust1_diff_<CASE>.h5' '/trdust1_diff' 'wtmean'  0.0  '/spath_s_trdust1_diff'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/spath_trdust2_diff_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/spath_i_trdust2_diff'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/spath_m_trdust2_diff'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/spath_f_trdust2_diff'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/spath_ps_trdust2_diff'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_trdust2_diff_<CASE>.h5' '/trdust2_diff' 'wtmean'  0.0  '/spath_s_trdust2_diff'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/spath_dust_adv_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/spath_i_dust_adv'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/spath_m_dust_adv'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/spath_f_dust_adv'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/spath_ps_dust_adv'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_dust_adv_<CASE>.h5' '/dust_adv' 'wtmean'  0.0  '/spath_s_dust_adv'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/spath_i_dust_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/spath_m_dust_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/spath_f_dust_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/spath_ps_dust_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/spath_s_dust_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_tracer_mass_<CASE>.h5' '/tracer_mass' 'wtmean'  0.0  '/spath_i_tracer_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer_mass_<CASE>.h5' '/tracer_mass' 'wtmean'  0.0  '/spath_m_tracer_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer_mass_<CASE>.h5' '/tracer_mass' 'wtmean'  0.0  '/spath_f_tracer_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer_mass_<CASE>.h5' '/tracer_mass' 'wtmean'  0.0  '/spath_ps_tracer_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_tracer_mass_<CASE>.h5' '/tracer_mass' 'wtmean'  0.0  '/spath_s_tracer_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sp_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_tcond_dust_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_tcond_i_dust_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_tcond_m_dust_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_tcond_f_dust_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_tcond_ps_dust_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_tcond_s_dust_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sp_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_clear_dust_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_clear_i_dust_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_clear_m_dust_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_clear_f_dust_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_clear_ps_dust_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_dust_mass_<CASE>.h5' '/dust_mass' 'wtmean'  0.0  '/sp_clear_s_dust_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_dust_mass_<CASE>.h5'   '/dust_mass'   'wtmean'  0.0  '/spath_dust_mass_ts'    'b_ZT'   ''      ''        ''       'ge' 0.01  }
%        { 'TsAveragedData/hist_spath_tracer_mass_<CASE>.h5' '/tracer_mass' 'wtmean'  0.0  '/spath_tracer_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01  }
%
%        { 'TsAveragedData/hist_sal_dust_mass_<CASE>.h5'   '/dust_mass'   'wtmean'  0.0  '/sal_dust_mass_ts'    'b_ZT'   ''      ''        ''       'ge' 0.01  }
%        { 'TsAveragedData/hist_sal_d1_mass_<CASE>.h5'     '/d1_mass'     'wtmean'  0.0  '/sal_d1_mass_ts'      'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sal_d2_mass_<CASE>.h5'     '/d2_mass'     'wtmean'  0.0  '/sal_d2_mass_ts'      'b_ZT'   ''      ''        ''       'ge' 0.01  }
%        { 'TsAveragedData/hist_sal_tracer_mass_<CASE>.h5' '/tracer_mass' 'wtmean'  0.0  '/sal_tracer_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01  }
%
%        { 'TsAveragedData/hist_sal_ar_dust_mass_<CASE>.h5'   '/dust_mass'   'wtmean'  0.0  '/sal_ar_dust_mass_ts'    'b_ZT'   ''      ''        ''       'ge' 0.01  }
%        { 'TsAveragedData/hist_sal_ar_d1_mass_<CASE>.h5'     '/d1_mass'     'wtmean'  0.0  '/sal_ar_d1_mass_ts'      'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sal_ar_d2_mass_<CASE>.h5'     '/d2_mass'     'wtmean'  0.0  '/sal_ar_d2_mass_ts'      'b_ZT'   ''      ''        ''       'ge' 0.01  }
%        { 'TsAveragedData/hist_sal_ar_tracer_mass_<CASE>.h5' '/tracer_mass' 'wtmean'  0.0  '/sal_ar_tracer_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01  }
%      }
%      'DIAGS/hist_meas_ts_dust_<CASE>.h5'
%    }
%
%    % CCN
%    {
%      'CCN Tsavg'
%      {
%        % sample regions in SAL
%        { 'TsAveragedData/hist_spath_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/spath_i_ccn_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/spath_m_ccn_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/spath_f_ccn_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/spath_ps_ccn_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/spath_s_ccn_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ccn_num_<CASE>.h5'  '/ccn_num'  'wtmean'  0.0  '/spath_ccn_num_ts'   'b_ZT'   ''      ''        ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_ccn_num_<CASE>.h5'  '/ccn_num'  'wtmean'  0.0  '/spath_i_ccn_num'    'b_Zt'   ''      'init'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_ccn_num_<CASE>.h5'  '/ccn_num'  'wtmean'  0.0  '/spath_m_ccn_num'    'b_Zt'   ''      'mid'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_ccn_num_<CASE>.h5'  '/ccn_num'  'wtmean'  0.0  '/spath_f_ccn_num'    'b_Zt'   ''      'final'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_ccn_num_<CASE>.h5'  '/ccn_num'  'wtmean'  0.0  '/spath_ps_ccn_num'   'b_Zt'   ''      'pre_sal' ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_ccn_num_<CASE>.h5'  '/ccn_num'  'wtmean'  0.0  '/spath_s_ccn_num'    'b_Zt'   ''      'sal'     ''       'ge' 10   }
%
%        { 'TsAveragedData/hist_sp_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_tcond_ccn_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_tcond_i_ccn_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_tcond_m_ccn_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_tcond_f_ccn_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_tcond_ps_ccn_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_tcond_s_ccn_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sp_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_clear_ccn_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_clear_i_ccn_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_clear_m_ccn_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_clear_f_ccn_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_clear_ps_ccn_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sp_clear_s_ccn_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/spath_ccn_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sal_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sal_ccn_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sal_ar_ccn_mass_<CASE>.h5' '/ccn_mass' 'wtmean'  0.0  '/sal_ar_ccn_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_ts_ccn_<CASE>.h5'
%    }
%
%    % Regenerated aerosols
%    {
%      'Regen Tsavg'
%      {
%        % sample regions in SAL
%        { 'TsAveragedData/hist_spath_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/spath_i_ra_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/spath_m_ra_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/spath_f_ra_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/spath_ps_ra_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/spath_s_ra_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sp_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_tcond_ra_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_tcond_i_ra_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_tcond_m_ra_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_tcond_f_ra_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_tcond_ps_ra_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_tcond_s_ra_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sp_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_clear_ra_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_clear_i_ra_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_clear_m_ra_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_clear_f_ra_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_clear_ps_ra_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sp_clear_s_ra_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/spath_ra1_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/spath_i_ra1_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/spath_m_ra1_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/spath_f_ra1_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/spath_ps_ra1_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra1_mass_<CASE>.h5' '/ra1_mass' 'wtmean'  0.0  '/spath_s_ra1_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra1_num_<CASE>.h5'  '/ra1_num'  'wtmean'  0.0  '/spath_ra1_num_ts'   'b_ZT'   ''      ''        ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_ra1_num_<CASE>.h5'  '/ra1_num'  'wtmean'  0.0  '/spath_i_ra1_num'    'b_Zt'   ''      'init'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_ra1_num_<CASE>.h5'  '/ra1_num'  'wtmean'  0.0  '/spath_m_ra1_num'    'b_Zt'   ''      'mid'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_ra1_num_<CASE>.h5'  '/ra1_num'  'wtmean'  0.0  '/spath_f_ra1_num'    'b_Zt'   ''      'final'    ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_ra1_num_<CASE>.h5'  '/ra1_num'  'wtmean'  0.0  '/spath_ps_ra1_num'   'b_Zt'   ''      'pre_sal' ''       'ge' 10   }
%        { 'TsAveragedData/hist_spath_ra1_num_<CASE>.h5'  '/ra1_num'  'wtmean'  0.0  '/spath_s_ra1_num'    'b_Zt'   ''      'sal'     ''       'ge' 10   }
%
%        { 'TsAveragedData/hist_spath_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/spath_ra2_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/spath_i_ra2_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/spath_m_ra2_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/spath_f_ra2_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/spath_ps_ra2_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra2_mass_<CASE>.h5' '/ra2_mass' 'wtmean'  0.0  '/spath_s_ra2_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra2_num_<CASE>.h5'  '/ra2_num'  'wtmean'  0.0  '/spath_ra2_num_ts'   'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra2_num_<CASE>.h5'  '/ra2_num'  'wtmean'  0.0  '/spath_i_ra2_num'    'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra2_num_<CASE>.h5'  '/ra2_num'  'wtmean'  0.0  '/spath_m_ra2_num'    'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra2_num_<CASE>.h5'  '/ra2_num'  'wtmean'  0.0  '/spath_f_ra2_num'    'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra2_num_<CASE>.h5'  '/ra2_num'  'wtmean'  0.0  '/spath_ps_ra2_num'   'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_ra2_num_<CASE>.h5'  '/ra2_num'  'wtmean'  0.0  '/spath_s_ra2_num'    'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/spath_ra_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sal_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sal_ra_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sal_ar_ra_mass_<CASE>.h5' '/ra_mass' 'wtmean'  0.0  '/sal_ar_ra_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_ts_ra_<CASE>.h5'
%    }
%
%    % "aero" (dust + regen) measurements
%    {
%      'AERO Tsavg'
%      {
%        % sample regions in SAL
%        { 'TsAveragedData/hist_spath_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/spath_i_aero_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/spath_m_aero_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/spath_f_aero_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/spath_ps_aero_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/spath_s_aero_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sp_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_tcond_aero_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_tcond_i_aero_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_tcond_m_aero_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_tcond_f_aero_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_tcond_ps_aero_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_tcond_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_tcond_s_aero_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sp_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_clear_aero_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_clear_i_aero_mass'   'b_Zt'   ''      'init'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_clear_m_aero_mass'   'b_Zt'   ''      'mid'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_clear_f_aero_mass'   'b_Zt'   ''      'final'    ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_clear_ps_aero_mass'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_sp_clear_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sp_clear_s_aero_mass'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/spath_aero_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sal_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sal_aero_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_sal_ar_aero_mass_<CASE>.h5' '/aero_mass' 'wtmean'  0.0  '/sal_ar_aero_mass_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_ts_aero_<CASE>.h5'
%    }
%
%    % Theta-E
%    {
%      'Theta-E Tsavg'
%      {
%        % sample regions in SAL
%        { 'TsAveragedData/hist_spath_theta_e_<CASE>.h5'  '/theta_e' 'wtmean'  0.0  '/spath_theta_e_ts'   'b_ZT'   ''      ''        ''       'ge' 0 }
%        { 'TsAveragedData/hist_smaxcp_theta_e_<CASE>.h5' '/theta_e' 'wtmean'  0.0  '/smaxcp_theta_e_ts'  'b_ZT'   ''      ''        ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_ts_theta_e_<CASE>.h5'
%    }
%
%    % Theta
%    {
%      'Theta Tsavg'
%      {
%        % sample regions in SAL
%        { 'TsAveragedData/hist_spath_theta_<CASE>.h5'  '/theta' 'wtmean'  0.0  '/spath_theta_ts'  'b_ZT'   ''      ''        ''       'ge' 0 }
%        { 'TsAveragedData/hist_smaxcp_theta_<CASE>.h5' '/theta' 'wtmean'  0.0  '/smaxcp_theta_ts' 'b_ZT'   ''      ''        ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_ts_theta_<CASE>.h5'
%    }
%
%    % Cold pools
%    {
%      'Cold Pools Tsavg'
%      {
%        % storm regions
%        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/lead_cpools_ts'      'b__T'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/lead_ps_cpools_ts'   'b__T'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/lead_s_cpools_ts'    'b__T'   ''      'sal'     ''       'ge' -100 }
%        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_lead_cpools'    'B__t'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_lead_ps_cpools' 'B__t'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_lead_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_lead_s_cpools'  'B__t'   ''      'sal'     ''       'ge' -100 }
%
%        { 'TsAveragedData/hist_core_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/core_cpools_ts'      'b__T'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_core_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/core_ps_cpools_ts'   'b__T'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_core_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/core_s_cpools_ts'    'b__T'   ''      'sal'     ''       'ge' -100 }
%        { 'TsAveragedData/hist_core_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_core_cpools'    'B__t'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_core_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_core_ps_cpools' 'B__t'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_core_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_core_s_cpools'  'B__t'   ''      'sal'     ''       'ge' -100 }
%
%        { 'TsAveragedData/hist_rb_cpools_<CASE>.h5'   '/cpools' 'wtmean'  0.0  '/rb_cpools_ts'        'b__T'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_rb_cpools_<CASE>.h5'   '/cpools' 'wtmean'  0.0  '/rb_ps_cpools_ts'     'b__T'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_rb_cpools_<CASE>.h5'   '/cpools' 'wtmean'  0.0  '/rb_s_cpools_ts'      'b__T'   ''      'sal'     ''       'ge' -100 }
%        { 'TsAveragedData/hist_rb_cpools_<CASE>.h5'   '/cpools' 'wtmean'  0.0  '/hist_rb_cpools'      'B__t'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_rb_cpools_<CASE>.h5'   '/cpools' 'wtmean'  0.0  '/hist_rb_ps_cpools'   'B__t'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_rb_cpools_<CASE>.h5'   '/cpools' 'wtmean'  0.0  '/hist_rb_s_cpools'    'B__t'   ''      'sal'     ''       'ge' -100 }
%
%%        { 'TsAveragedData/hist_env_cpools_<CASE>.h5'  '/cpools' 'wtmean'  0.0  '/env_cpools_ts'       'b__T'   ''      ''        ''       'ge' -100 }
%%        { 'TsAveragedData/hist_env_cpools_<CASE>.h5'  '/cpools' 'wtmean'  0.0  '/env_ps_cpools_ts'    'b__T'   ''      'pre_sal' ''       'ge' -100 }
%%        { 'TsAveragedData/hist_env_cpools_<CASE>.h5'  '/cpools' 'wtmean'  0.0  '/env_s_cpools_ts'     'b__T'   ''      'sal'     ''       'ge' -100 }
%%        { 'TsAveragedData/hist_env_cpools_<CASE>.h5'  '/cpools' 'wtmean'  0.0  '/hist_env_cpools'     'B__t'   ''      ''        ''       'ge' -100 }
%%        { 'TsAveragedData/hist_env_cpools_<CASE>.h5'  '/cpools' 'wtmean'  0.0  '/hist_env_ps_cpools'  'B__t'   ''      'pre_sal' ''       'ge' -100 }
%%        { 'TsAveragedData/hist_env_cpools_<CASE>.h5'  '/cpools' 'wtmean'  0.0  '/hist_env_s_cpools'   'B__t'   ''      'sal'     ''       'ge' -100 }
%
%        % sample regions in SAL
%        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/spath_cpools_ts'      'b__T'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/spath_ps_cpools_ts'   'b__T'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/spath_s_cpools_ts'    'b__T'   ''      'sal'     ''       'ge' -100 }
%        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_spath_cpools'    'B__t'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_spath_ps_cpools' 'B__t'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_spath_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_spath_s_cpools'  'B__t'   ''      'sal'     ''       'ge' -100 }
%
%        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/smaxcp_cpools_ts'      'b__T'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/smaxcp_ps_cpools_ts'   'b__T'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/smaxcp_s_cpools_ts'    'b__T'   ''      'sal'     ''       'ge' -100 }
%        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_smaxcp_cpools'    'B__t'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_smaxcp_ps_cpools' 'B__t'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_smaxcp_cpools_<CASE>.h5' '/cpools' 'wtmean'  0.0  '/hist_smaxcp_s_cpools'  'B__t'   ''      'sal'     ''       'ge' -100 }
%
%      }
%      'DIAGS/hist_meas_ts_cpools_<CASE>.h5'
%    }
%
%    % Precip rate
%    {
%      'Precip Rate Tsavg'
%      {
%        % sample region
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
%        { 'TsAveragedData/hist_spath_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/spath_pcprate_ts'      'b__T'   ''      ''        ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sal_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/sal_pcprate_ts'      'b__T'   ''      ''        ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sal_ar_pcprate_<CASE>.h5' '/pcprate' 'wtmean'  0.0  '/sal_ar_pcprate_ts'      'b__T'   ''      ''        ''       'ge' 0.001 }
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

%    % Vapor
%    {
%      'Vapor Tsavg'
%      {
%        % sample regions in SAL
%        { 'TsAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_vapor_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_ps_vapor'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_lead_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/lead_s_vapor'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_spath_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/spath_vapor_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/spath_ps_vapor'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_spath_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/spath_s_vapor'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%
%        { 'TsAveragedData/hist_smaxcp_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/smaxcp_vapor_ts'  'b_ZT'   ''      ''        ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_smaxcp_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/smaxcp_ps_vapor'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.01 }
%        { 'TsAveragedData/hist_smaxcp_vapor_<CASE>.h5' '/vapor' 'wtmean'  0.0  '/smaxcp_s_vapor'   'b_Zt'   ''      'sal'     ''       'ge' 0.01 }
%      }
%      'DIAGS/hist_meas_ts_vapor_<CASE>.h5'
%    }
%
%    % Temperature
%    {
%      'Temperature Tsavg'
%      {
%        % sample regions in SAL
%        { 'TsAveragedData/hist_spath_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/spath_tempc_ts'  'b_ZT'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_spath_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/spath_ps_tempc'  'b_Zt'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_spath_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/spath_s_tempc'   'b_Zt'   ''      'sal'     ''       'ge' -100 }
%
%        { 'TsAveragedData/hist_smaxcp_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/smaxcp_tempc_ts'  'b_ZT'   ''      ''        ''       'ge' -100 }
%        { 'TsAveragedData/hist_smaxcp_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/smaxcp_ps_tempc'  'b_Zt'   ''      'pre_sal' ''       'ge' -100 }
%        { 'TsAveragedData/hist_smaxcp_tempc_<CASE>.h5' '/tempc' 'wtmean'  0.0  '/smaxcp_s_tempc'   'b_Zt'   ''      'sal'     ''       'ge' -100 }
%      }
%      'DIAGS/hist_meas_ts_tempc_<CASE>.h5'
%    }
%
%    % Relative humidity
%    {
%      'RH Tsavg'
%      {
%        % sample regions in SAL
%        { 'TsAveragedData/hist_lead_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/lead_relhum_ts'  'b_ZT'   ''      ''        ''       'ge' 0 }
%        { 'TsAveragedData/hist_lead_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/lead_ps_relhum'  'b_Zt'   ''      'pre_sal' ''       'ge' 0 }
%        { 'TsAveragedData/hist_lead_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/lead_s_relhum'   'b_Zt'   ''      'sal'     ''       'ge' 0 }
%
%        { 'TsAveragedData/hist_spath_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/spath_relhum_ts'  'b_ZT'   ''      ''        ''       'ge' 0 }
%        { 'TsAveragedData/hist_spath_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/spath_ps_relhum'  'b_Zt'   ''      'pre_sal' ''       'ge' 0 }
%        { 'TsAveragedData/hist_spath_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/spath_s_relhum'   'b_Zt'   ''      'sal'     ''       'ge' 0 }
%
%        { 'TsAveragedData/hist_smaxcp_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/smaxcp_relhum_ts'  'b_ZT'   ''      ''        ''       'ge' 0 }
%        { 'TsAveragedData/hist_smaxcp_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/smaxcp_ps_relhum'  'b_Zt'   ''      'pre_sal' ''       'ge' 0 }
%        { 'TsAveragedData/hist_smaxcp_relhum_<CASE>.h5' '/relhum' 'wtmean'  0.0  '/smaxcp_s_relhum'   'b_Zt'   ''      'sal'     ''       'ge' 0 }
%      }
%      'DIAGS/hist_meas_ts_relhum_<CASE>.h5'
%    }
%
%    % Dust embedded inside specific hydrometeor species
%    {
%      'Dust In Hydrometeors Tsavg'
%      {
%        { 'TsAveragedData/hist_spath_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/spath_dust_cloud_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/spath_i_dust_cloud'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/spath_m_dust_cloud'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/spath_f_dust_cloud'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/spath_ps_dust_cloud'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/spath_s_dust_cloud'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/spath_dust_rain_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/spath_i_dust_rain'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/spath_m_dust_rain'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/spath_f_dust_rain'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/spath_ps_dust_rain'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/spath_s_dust_rain'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/spath_dust_pris_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/spath_i_dust_pris'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/spath_m_dust_pris'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/spath_f_dust_pris'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/spath_ps_dust_pris'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/spath_s_dust_pris'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/spath_dust_snow_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/spath_i_dust_snow'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/spath_m_dust_snow'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/spath_f_dust_snow'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/spath_ps_dust_snow'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/spath_s_dust_snow'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/spath_dust_aggr_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/spath_i_dust_aggr'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/spath_m_dust_aggr'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/spath_f_dust_aggr'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/spath_ps_dust_aggr'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/spath_s_dust_aggr'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/spath_dust_graup_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/spath_i_dust_graup'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/spath_m_dust_graup'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/spath_f_dust_graup'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/spath_ps_dust_graup'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/spath_s_dust_graup'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/spath_dust_hail_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/spath_i_dust_hail'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/spath_m_dust_hail'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/spath_f_dust_hail'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/spath_ps_dust_hail'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/spath_s_dust_hail'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/spath_i_dust_hydro'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/spath_m_dust_hydro'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/spath_f_dust_hydro'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/spath_ps_dust_hydro'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/spath_s_dust_hydro'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sp_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/sp_tcond_dust_cloud_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/sp_tcond_i_dust_cloud'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/sp_tcond_m_dust_cloud'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/sp_tcond_f_dust_cloud'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/sp_tcond_ps_dust_cloud'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_cloud_<CASE>.h5' '/dust_cloud' 'wtmean'  0.0  '/sp_tcond_s_dust_cloud'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sp_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/sp_tcond_dust_rain_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/sp_tcond_i_dust_rain'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/sp_tcond_m_dust_rain'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/sp_tcond_f_dust_rain'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/sp_tcond_ps_dust_rain'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_rain_<CASE>.h5' '/dust_rain' 'wtmean'  0.0  '/sp_tcond_s_dust_rain'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sp_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/sp_tcond_dust_pris_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/sp_tcond_i_dust_pris'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/sp_tcond_m_dust_pris'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/sp_tcond_f_dust_pris'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/sp_tcond_ps_dust_pris'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_pris_<CASE>.h5' '/dust_pris' 'wtmean'  0.0  '/sp_tcond_s_dust_pris'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sp_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/sp_tcond_dust_snow_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/sp_tcond_i_dust_snow'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/sp_tcond_m_dust_snow'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/sp_tcond_f_dust_snow'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/sp_tcond_ps_dust_snow'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_snow_<CASE>.h5' '/dust_snow' 'wtmean'  0.0  '/sp_tcond_s_dust_snow'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sp_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/sp_tcond_dust_aggr_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/sp_tcond_i_dust_aggr'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/sp_tcond_m_dust_aggr'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/sp_tcond_f_dust_aggr'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/sp_tcond_ps_dust_aggr'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_aggr_<CASE>.h5' '/dust_aggr' 'wtmean'  0.0  '/sp_tcond_s_dust_aggr'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sp_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/sp_tcond_dust_graup_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/sp_tcond_i_dust_graup'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/sp_tcond_m_dust_graup'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/sp_tcond_f_dust_graup'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/sp_tcond_ps_dust_graup'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_graup_<CASE>.h5' '/dust_graup' 'wtmean'  0.0  '/sp_tcond_s_dust_graup'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sp_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/sp_tcond_dust_hail_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/sp_tcond_i_dust_hail'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/sp_tcond_m_dust_hail'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/sp_tcond_f_dust_hail'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/sp_tcond_ps_dust_hail'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_hail_<CASE>.h5' '/dust_hail' 'wtmean'  0.0  '/sp_tcond_s_dust_hail'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sp_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/sp_tcond_dust_hydro_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/sp_tcond_i_dust_hydro'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/sp_tcond_m_dust_hydro'   'b_Zt'   ''      'mid'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/sp_tcond_f_dust_hydro'   'b_Zt'   ''      'final'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/sp_tcond_ps_dust_hydro'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_sp_tcond_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/sp_tcond_s_dust_hydro'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/spath_dust_hydro_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sal_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/sal_dust_hydro_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_sal_ar_dust_hydro_<CASE>.h5' '/dust_hydro' 'wtmean'  0.0  '/sal_ar_dust_hydro_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%      }
%      'DIAGS/hist_meas_ts_dust_hydro_<CASE>.h5'
%    }
%
%    % Dust embedded inside specific hydrometeor species that also acted as IFN
%    {
%      'Dust As IFN In Hydrometeors Tsavg'
%      {
%        { 'TsAveragedData/hist_spath_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/spath_dustifn_cloud_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/spath_i_dustifn_cloud'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/spath_ps_dustifn_cloud'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_cloud_<CASE>.h5' '/dustifn_cloud' 'wtmean'  0.0  '/spath_s_dustifn_cloud'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/spath_dustifn_rain_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/spath_i_dustifn_rain'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/spath_ps_dustifn_rain'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_rain_<CASE>.h5' '/dustifn_rain' 'wtmean'  0.0  '/spath_s_dustifn_rain'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/spath_dustifn_pris_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/spath_i_dustifn_pris'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/spath_ps_dustifn_pris'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_pris_<CASE>.h5' '/dustifn_pris' 'wtmean'  0.0  '/spath_s_dustifn_pris'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/spath_dustifn_snow_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/spath_i_dustifn_snow'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/spath_ps_dustifn_snow'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_snow_<CASE>.h5' '/dustifn_snow' 'wtmean'  0.0  '/spath_s_dustifn_snow'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/spath_dustifn_aggr_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/spath_i_dustifn_aggr'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/spath_ps_dustifn_aggr'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_aggr_<CASE>.h5' '/dustifn_aggr' 'wtmean'  0.0  '/spath_s_dustifn_aggr'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/spath_dustifn_graup_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/spath_i_dustifn_graup'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/spath_ps_dustifn_graup'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_graup_<CASE>.h5' '/dustifn_graup' 'wtmean'  0.0  '/spath_s_dustifn_graup'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/spath_dustifn_hail_ts'  'b_ZT'   ''      ''        ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/spath_i_dustifn_hail'   'b_Zt'   ''      'init'    ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/spath_ps_dustifn_hail'  'b_Zt'   ''      'pre_sal' ''       'ge' 0.001 }
%        { 'TsAveragedData/hist_spath_dustifn_hail_<CASE>.h5' '/dustifn_hail' 'wtmean'  0.0  '/spath_s_dustifn_hail'   'b_Zt'   ''      'sal'     ''       'ge' 0.001 }
%
%        { 'TsAveragedData/hist_spath_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/spath_dustifn_hydro_ts'  'b_ZT'   ''      ''        ''       'ge' 1e-6 }
%        { 'TsAveragedData/hist_spath_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/spath_i_dustifn_hydro'   'b_Zt'   ''      'init'    ''       'ge' 1e-6 }
%        { 'TsAveragedData/hist_spath_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/spath_ps_dustifn_hydro'  'b_Zt'   ''      'pre_sal' ''       'ge' 1e-6 }
%        { 'TsAveragedData/hist_spath_dustifn_hydro_<CASE>.h5' '/dustifn_hydro' 'wtmean'  0.0  '/spath_s_dustifn_hydro'   'b_Zt'   ''      'sal'     ''       'ge' 1e-6 }
%      }
%      'DIAGS/hist_meas_ts_dustifn_hydro_<CASE>.h5'
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
  
      icount = 0;
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

        % skip this profile if input file is missing (eg., NODUST cases don't
        % have d1_num, d1_mass, etc. extractions)
        if (exist(InFile, 'file') ~= 2)
          fprintf('      Input files does not exist, skipping: %s (%s)\n', InFile, Vname);
          continue
        else
          icount = icount + 1;
        end

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

        % Simple test for determining if vertical coordinate is height or pressure
        if (length(Z) > 1)
          IsPress = (Z(1) > Z(2));
        else
          IsPress = false
        end
  
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
  
            case 'gt'
              B1 = find(B > SelectVal, 1, 'first');
  
            case 'le'
              B2 = find(B <= SelectVal, 1, 'last');

            case 'lt'
              B2 = find(B < SelectVal, 1, 'last');
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
          if (IsPress)
            % leave in mb
            ZconvFact = 1.0;
          else
            % convert m to km
            ZconvFact = 1.0e-3;
          end
          switch(ZinVar)
            case 'x'
              H = X .* ZconvFact;
            case 'y'
              H = Y .* ZconvFact;
            case 'z'
              H = Z .* ZconvFact;
            case 't'
              H = T .* ZconvFact;
          end

          Z1 = 1;
          Z2 = length(H);
          switch(Zrange)
            case { 'maxlev' 'minlev' 'sfc' }
              if (IsPress)
                Z1 = find(H <= Psfc, 1, 'first');
                Z2 = find(H >= Ptop, 1, 'last');
              else
                Z1 = find(H >= Zsfc, 1, 'first');
                Z2 = find(H <= Ztop, 1, 'last');
              end
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
              ST = (X .* TimeScale) - TimeOffset;
            case 'y'
              ST = (Y .* TimeScale) - TimeOffset;
            case 'z'
              ST = (Z .* TimeScale) - TimeOffset;
            case 't'
              ST = (T .* TimeScale) - TimeOffset;
          end
          T1 = 1;
          T2 = length(ST);
          switch(Trange)
            case 'init'
              T1 = find(ST >= TstartInit, 1, 'first');
              T2 = find(ST <= TendInit,   1, 'last');
  
            case 'mid'
              T1 = find(ST >= TstartMid, 1, 'first');
              T2 = find(ST <= TendMid,   1, 'last');
  
            case 'final'
              T1 = find(ST >= TstartFinal, 1, 'first');
              T2 = find(ST <= TendFinal,   1, 'last');
  
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
        if (icount == 1)
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
