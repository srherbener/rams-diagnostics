function [ ] = GenTsdHdaMeas()

  % make sure output directory exists
  Ddir = 'DIAGS';  % coordinate this with output file names below
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % list of simulation cases
  CaseList = {
    'TSD_SAL_DUST'
%    'TSD_SAL_NODUST'
%    'TSD_NONSAL_DUST'
%    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  % Temporal ranges
  TstartInit = 0;
  TendInit = 2;

  TstartPreSal = 10;
  TendPreSal   = 30;

  TstartSal = 40;
  TendSal   = 60;

  % Height ranges
  Zsfc = 0;
  Ztop = 1.5; % km

  % Input is the output from tsavg hda function, based upon either 2D or 3D field
  %   2D - input will be of the form (h,t)
  %   3D - input will be of the form (h,z,t)
  %
  %   where: h - hda data, size = 2
  %               1 - sum
  %               2 - count
  %               where average = sum/count
  %          z - height
  %          t - time
  %
  % The input form and reduction is specified by a 4 letter string. The string follows
  % the convention of tsavg always writing a 4 dimension variable of the form (x,y,z,t).
  % The reduction spec uses 4 letters that correspond to the xyzt dimension order of the input.
  % If a lower case string is used, then that means to reduce that dimension. If an underscore
  % is used, that means that this input dimension (one of x,y,z,t) is not used.
  %
  % Reduction for the h dimension:
  %       h  -> average (sum/count)
  %       H  -> sum
  %
  % Example, to read in a 3D field from tsavg, and reduce to a vertical profile of
  % average values, use '_hZt'.
  % This means:
  %
  %   input dim     quatity     reduce
  %       x          unused       ---
  %       y          hda          avg
  %       z          height       no
  %       t          time         yes
  %
  % To produce a time series of vertical profiles of sum data from tsavg 3D output,
  % use '_HZT'. This means:
  %
  %   input dim     quatity     reduce
  %       x          unused       ---
  %       y          hda          sum
  %       z          height       no
  %       t          time         no
  %
  % To produce a time series of average data from 2D tsavg output, use '_h_T'. This means:
  %
  %   input dim     quatity     reduce
  %       x          unused       ---
  %       y          hda          avg
  %       z          unused       ---
  %       t          time         no
  %
  %
  % Another notion is to avoid combining counts along the z axis, i.e., preserve
  % levels.
  %
  % Range specs for summing counts:
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
  %     <in_file> <var_name> <out_var_name> <reduce_spec> <t_range> <z_range>
  %       ranges are defined above
  %       <reduce_spec> describes input dimensions (and order) with upper case letters on the dims that are to be reduced.
  %           <reduce_spec> = 4-character string that corresponds to 'xyzt' that describes which dimensions are hda sum and count data, height and time,
  %                           as well as sepcifiying which dimensions are to be reduced (lower case).
  %                For example, '_hZT' means that the input has three dimensions (after squeezing), hda is y, height is z, and time is t
  %                making the input data (h,z,t). The lower case 'h' says to create averages (sum/count) and the result will be a matrix
  %                of the form (z,t)
  %
  %                '_HZt' means hda is y, height is z and time is t, create sums, reduce t, leaving the output (z).

  MeasSets = {

    % Dust
    {
      'Dust Tsavg'
      {
        % region in SAL that is in the storm path
        { 'TsAveragedData/hda_spath_d1_mass_<CASE>.h5' '/d1_mass' '/spath_d1_mass_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_spath_d1_mass_<CASE>.h5' '/d1_mass' '/spath_i_d1_mass'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_spath_d1_mass_<CASE>.h5' '/d1_mass' '/spath_ps_d1_mass'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_spath_d1_mass_<CASE>.h5' '/d1_mass' '/spath_s_d1_mass'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_spath_d2_mass_<CASE>.h5' '/d2_mass' '/spath_d2_mass_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_spath_d2_mass_<CASE>.h5' '/d2_mass' '/spath_i_d2_mass'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_spath_d2_mass_<CASE>.h5' '/d2_mass' '/spath_ps_d2_mass'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_spath_d2_mass_<CASE>.h5' '/d2_mass' '/spath_s_d2_mass'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_spath_dust_mass_<CASE>.h5' '/dust_mass' '/spath_dust_mass_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_spath_dust_mass_<CASE>.h5' '/dust_mass' '/spath_i_dust_mass'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_spath_dust_mass_<CASE>.h5' '/dust_mass' '/spath_ps_dust_mass'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_spath_dust_mass_<CASE>.h5' '/dust_mass' '/spath_s_dust_mass'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_spath_dust_hydro_<CASE>.h5' '/dust_hydro' '/spath_dust_hydro_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_spath_dust_hydro_<CASE>.h5' '/dust_hydro' '/spath_i_dust_hydro'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_spath_dust_hydro_<CASE>.h5' '/dust_hydro' '/spath_ps_dust_hydro'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_spath_dust_hydro_<CASE>.h5' '/dust_hydro' '/spath_s_dust_hydro'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_spath_dust_cloud_<CASE>.h5' '/dust_cloud' '/spath_dust_cloud_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_spath_dust_cloud_<CASE>.h5' '/dust_cloud' '/spath_i_dust_cloud'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_spath_dust_cloud_<CASE>.h5' '/dust_cloud' '/spath_ps_dust_cloud'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_spath_dust_cloud_<CASE>.h5' '/dust_cloud' '/spath_s_dust_cloud'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_spath_dust_rain_<CASE>.h5' '/dust_rain' '/spath_dust_rain_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_spath_dust_rain_<CASE>.h5' '/dust_rain' '/spath_i_dust_rain'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_spath_dust_rain_<CASE>.h5' '/dust_rain' '/spath_ps_dust_rain'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_spath_dust_rain_<CASE>.h5' '/dust_rain' '/spath_s_dust_rain'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_spath_dust_pris_<CASE>.h5' '/dust_pris' '/spath_dust_pris_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_spath_dust_pris_<CASE>.h5' '/dust_pris' '/spath_i_dust_pris'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_spath_dust_pris_<CASE>.h5' '/dust_pris' '/spath_ps_dust_pris'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_spath_dust_pris_<CASE>.h5' '/dust_pris' '/spath_s_dust_pris'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_spath_dust_snow_<CASE>.h5' '/dust_snow' '/spath_dust_snow_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_spath_dust_snow_<CASE>.h5' '/dust_snow' '/spath_i_dust_snow'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_spath_dust_snow_<CASE>.h5' '/dust_snow' '/spath_ps_dust_snow'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_spath_dust_snow_<CASE>.h5' '/dust_snow' '/spath_s_dust_snow'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_spath_dust_aggr_<CASE>.h5' '/dust_aggr' '/spath_dust_aggr_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_spath_dust_aggr_<CASE>.h5' '/dust_aggr' '/spath_i_dust_aggr'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_spath_dust_aggr_<CASE>.h5' '/dust_aggr' '/spath_ps_dust_aggr'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_spath_dust_aggr_<CASE>.h5' '/dust_aggr' '/spath_s_dust_aggr'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_spath_dust_graup_<CASE>.h5' '/dust_graup' '/spath_dust_graup_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_spath_dust_graup_<CASE>.h5' '/dust_graup' '/spath_i_dust_graup'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_spath_dust_graup_<CASE>.h5' '/dust_graup' '/spath_ps_dust_graup'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_spath_dust_graup_<CASE>.h5' '/dust_graup' '/spath_s_dust_graup'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_spath_dust_hail_<CASE>.h5' '/dust_hail' '/spath_dust_hail_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_spath_dust_hail_<CASE>.h5' '/dust_hail' '/spath_i_dust_hail'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_spath_dust_hail_<CASE>.h5' '/dust_hail' '/spath_ps_dust_hail'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_spath_dust_hail_<CASE>.h5' '/dust_hail' '/spath_s_dust_hail'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_spath_dust_mass_<CASE>.h5'   '/dust_mass'   '/spath_sum_dust_mass'    '_HZT'   ''        '' }
        { 'TsAveragedData/hda_spath_dust_hydro_<CASE>.h5'  '/dust_hydro'  '/spath_sum_dust_hydro'   '_HZT'   ''        '' }
        { 'TsAveragedData/hda_spath_tracer_mass_<CASE>.h5' '/tracer_mass' '/spath_sum_tracer_mass'  '_HZT'   ''        '' }
        { 'TsAveragedData/hda_spath_dust_sfc_<CASE>.h5'    '/dust_sfc'    '/spath_sum_dust_sfc'     '_H_T'   ''        '' }

        % SAL sample region (large)
        { 'TsAveragedData/hda_sal_dust_mass_<CASE>.h5'   '/dust_mass'   '/sal_sum_dust_mass'   '_HZT'   ''        '' }
        { 'TsAveragedData/hda_sal_dust_hydro_<CASE>.h5'  '/dust_hydro'  '/sal_sum_dust_hydro'  '_HZT'   ''        '' }
        { 'TsAveragedData/hda_sal_tracer_mass_<CASE>.h5' '/tracer_mass' '/sal_sum_tracer_mass' '_HZT'   ''        '' }
        { 'TsAveragedData/hda_sal_dust_sfc_<CASE>.h5'    '/dust_sfc'    '/sal_sum_dust_sfc'    '_H_T'   ''        '' }

        % region of storm
        { 'TsAveragedData/hda_all_d1_mass_<CASE>.h5' '/d1_mass' '/all_d1_mass_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_all_d1_mass_<CASE>.h5' '/d1_mass' '/all_i_d1_mass'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_all_d1_mass_<CASE>.h5' '/d1_mass' '/all_ps_d1_mass'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_all_d1_mass_<CASE>.h5' '/d1_mass' '/all_s_d1_mass'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_all_d2_mass_<CASE>.h5' '/d2_mass' '/all_d2_mass_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_all_d2_mass_<CASE>.h5' '/d2_mass' '/all_i_d2_mass'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_all_d2_mass_<CASE>.h5' '/d2_mass' '/all_ps_d2_mass'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_all_d2_mass_<CASE>.h5' '/d2_mass' '/all_s_d2_mass'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_all_dust_mass_<CASE>.h5' '/dust_mass' '/all_dust_mass_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_all_dust_mass_<CASE>.h5' '/dust_mass' '/all_i_dust_mass'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_all_dust_mass_<CASE>.h5' '/dust_mass' '/all_ps_dust_mass'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_all_dust_mass_<CASE>.h5' '/dust_mass' '/all_s_dust_mass'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_all_dust_hydro_<CASE>.h5' '/dust_hydro' '/all_dust_hydro_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_all_dust_hydro_<CASE>.h5' '/dust_hydro' '/all_i_dust_hydro'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_all_dust_hydro_<CASE>.h5' '/dust_hydro' '/all_ps_dust_hydro'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_all_dust_hydro_<CASE>.h5' '/dust_hydro' '/all_s_dust_hydro'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_all_dust_cloud_<CASE>.h5' '/dust_cloud' '/all_dust_cloud_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_all_dust_cloud_<CASE>.h5' '/dust_cloud' '/all_i_dust_cloud'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_all_dust_cloud_<CASE>.h5' '/dust_cloud' '/all_ps_dust_cloud'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_all_dust_cloud_<CASE>.h5' '/dust_cloud' '/all_s_dust_cloud'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_all_dust_rain_<CASE>.h5' '/dust_rain' '/all_dust_rain_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_all_dust_rain_<CASE>.h5' '/dust_rain' '/all_i_dust_rain'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_all_dust_rain_<CASE>.h5' '/dust_rain' '/all_ps_dust_rain'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_all_dust_rain_<CASE>.h5' '/dust_rain' '/all_s_dust_rain'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_all_dust_pris_<CASE>.h5' '/dust_pris' '/all_dust_pris_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_all_dust_pris_<CASE>.h5' '/dust_pris' '/all_i_dust_pris'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_all_dust_pris_<CASE>.h5' '/dust_pris' '/all_ps_dust_pris'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_all_dust_pris_<CASE>.h5' '/dust_pris' '/all_s_dust_pris'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_all_dust_snow_<CASE>.h5' '/dust_snow' '/all_dust_snow_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_all_dust_snow_<CASE>.h5' '/dust_snow' '/all_i_dust_snow'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_all_dust_snow_<CASE>.h5' '/dust_snow' '/all_ps_dust_snow'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_all_dust_snow_<CASE>.h5' '/dust_snow' '/all_s_dust_snow'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_all_dust_aggr_<CASE>.h5' '/dust_aggr' '/all_dust_aggr_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_all_dust_aggr_<CASE>.h5' '/dust_aggr' '/all_i_dust_aggr'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_all_dust_aggr_<CASE>.h5' '/dust_aggr' '/all_ps_dust_aggr'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_all_dust_aggr_<CASE>.h5' '/dust_aggr' '/all_s_dust_aggr'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_all_dust_graup_<CASE>.h5' '/dust_graup' '/all_dust_graup_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_all_dust_graup_<CASE>.h5' '/dust_graup' '/all_i_dust_graup'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_all_dust_graup_<CASE>.h5' '/dust_graup' '/all_ps_dust_graup'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_all_dust_graup_<CASE>.h5' '/dust_graup' '/all_s_dust_graup'   '_hZt'   'sal'     '' }

        { 'TsAveragedData/hda_all_dust_hail_<CASE>.h5' '/dust_hail' '/all_dust_hail_ts'  '_hZT'   ''        '' }
        { 'TsAveragedData/hda_all_dust_hail_<CASE>.h5' '/dust_hail' '/all_i_dust_hail'   '_hZt'   'init'    '' }
        { 'TsAveragedData/hda_all_dust_hail_<CASE>.h5' '/dust_hail' '/all_ps_dust_hail'  '_hZt'   'pre_sal' '' }
        { 'TsAveragedData/hda_all_dust_hail_<CASE>.h5' '/dust_hail' '/all_s_dust_hail'   '_hZt'   'sal'     '' }
      }
      'DIAGS/hda_meas_ts_dust_<CASE>.h5'
    }

%    % CCN
%    {
%      'CCN Tsavg'
%      {
%        % region in SAL that is in the storm path
%        { 'TsAveragedData/hda_spath_ccn_mass_<CASE>.h5' '/ccn_mass' '/spath_ccn_mass_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_ccn_mass_<CASE>.h5' '/ccn_mass' '/spath_i_ccn_mass'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_ccn_mass_<CASE>.h5' '/ccn_mass' '/spath_ps_ccn_mass'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_ccn_mass_<CASE>.h5' '/ccn_mass' '/spath_s_ccn_mass'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_ccn_mass_<CASE>.h5' '/ccn_mass' '/spath_sum_ccn_mass'  '_HZT'   ''        '' }
%
%        % SAL sample region (large)
%        { 'TsAveragedData/hda_sal_ccn_mass_<CASE>.h5'    '/ccn_mass'    '/sal_sum_ccn_mass'    '_HZT'   ''        '' }
%
%        % storm region
%        { 'TsAveragedData/hda_all_ccn_mass_<CASE>.h5' '/ccn_mass' '/all_ccn_mass_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_ccn_mass_<CASE>.h5' '/ccn_mass' '/all_i_ccn_mass'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_ccn_mass_<CASE>.h5' '/ccn_mass' '/all_ps_ccn_mass'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_ccn_mass_<CASE>.h5' '/ccn_mass' '/all_s_ccn_mass'   '_hZt'   'sal'     '' }
%      }
%      'DIAGS/hda_meas_ts_ccn_<CASE>.h5'
%    }
%
%    % Regen
%    {
%      'Regen Tsavg'
%      {
%        % region in SAL that is in the storm path
%        { 'TsAveragedData/hda_spath_ra1_mass_<CASE>.h5' '/ra1_mass' '/spath_ra1_mass_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_ra1_mass_<CASE>.h5' '/ra1_mass' '/spath_i_ra1_mass'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_ra1_mass_<CASE>.h5' '/ra1_mass' '/spath_ps_ra1_mass'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_ra1_mass_<CASE>.h5' '/ra1_mass' '/spath_s_ra1_mass'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_ra2_mass_<CASE>.h5' '/ra2_mass' '/spath_ra2_mass_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_ra2_mass_<CASE>.h5' '/ra2_mass' '/spath_i_ra2_mass'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_ra2_mass_<CASE>.h5' '/ra2_mass' '/spath_ps_ra2_mass'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_ra2_mass_<CASE>.h5' '/ra2_mass' '/spath_s_ra2_mass'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_ra_mass_<CASE>.h5' '/ra_mass' '/spath_ra_mass_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_ra_mass_<CASE>.h5' '/ra_mass' '/spath_i_ra_mass'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_ra_mass_<CASE>.h5' '/ra_mass' '/spath_ps_ra_mass'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_ra_mass_<CASE>.h5' '/ra_mass' '/spath_s_ra_mass'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_ra_mass_<CASE>.h5' '/ra_mass' '/spath_sum_ra_mass'  '_HZT'   ''        '' }
%
%        % SAL sample region (large)
%        { 'TsAveragedData/hda_sal_ra_mass_<CASE>.h5'     '/ra_mass'     '/sal_sum_ra_mass'     '_HZT'   ''        '' }
%
%        % storm region
%        { 'TsAveragedData/hda_all_ra1_mass_<CASE>.h5' '/ra1_mass' '/all_ra1_mass_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_ra1_mass_<CASE>.h5' '/ra1_mass' '/all_i_ra1_mass'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_ra1_mass_<CASE>.h5' '/ra1_mass' '/all_ps_ra1_mass'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_ra1_mass_<CASE>.h5' '/ra1_mass' '/all_s_ra1_mass'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_all_ra2_mass_<CASE>.h5' '/ra2_mass' '/all_ra2_mass_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_ra2_mass_<CASE>.h5' '/ra2_mass' '/all_i_ra2_mass'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_ra2_mass_<CASE>.h5' '/ra2_mass' '/all_ps_ra2_mass'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_ra2_mass_<CASE>.h5' '/ra2_mass' '/all_s_ra2_mass'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_all_ra_mass_<CASE>.h5' '/ra_mass' '/all_ra_mass_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_ra_mass_<CASE>.h5' '/ra_mass' '/all_i_ra_mass'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_ra_mass_<CASE>.h5' '/ra_mass' '/all_ps_ra_mass'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_ra_mass_<CASE>.h5' '/ra_mass' '/all_s_ra_mass'   '_hZt'   'sal'     '' }
%      }
%      'DIAGS/hda_meas_ts_ra_<CASE>.h5'
%    }
%
%    % All aerosols
%    {
%      'Aero Tsavg'
%      {
%        % region in SAL that is in the storm path
%        { 'TsAveragedData/hda_spath_aero_mass_<CASE>.h5' '/aero_mass' '/spath_aero_mass_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_aero_mass_<CASE>.h5' '/aero_mass' '/spath_i_aero_mass'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_aero_mass_<CASE>.h5' '/aero_mass' '/spath_ps_aero_mass'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_aero_mass_<CASE>.h5' '/aero_mass' '/spath_s_aero_mass'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_aero_hydro_<CASE>.h5' '/aero_hydro' '/spath_aero_hydro_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_aero_hydro_<CASE>.h5' '/aero_hydro' '/spath_i_aero_hydro'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_aero_hydro_<CASE>.h5' '/aero_hydro' '/spath_ps_aero_hydro'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_aero_hydro_<CASE>.h5' '/aero_hydro' '/spath_s_aero_hydro'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_aero_cloud_<CASE>.h5' '/aero_cloud' '/spath_aero_cloud_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_aero_cloud_<CASE>.h5' '/aero_cloud' '/spath_i_aero_cloud'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_aero_cloud_<CASE>.h5' '/aero_cloud' '/spath_ps_aero_cloud'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_aero_cloud_<CASE>.h5' '/aero_cloud' '/spath_s_aero_cloud'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_aero_rain_<CASE>.h5' '/aero_rain' '/spath_aero_rain_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_aero_rain_<CASE>.h5' '/aero_rain' '/spath_i_aero_rain'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_aero_rain_<CASE>.h5' '/aero_rain' '/spath_ps_aero_rain'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_aero_rain_<CASE>.h5' '/aero_rain' '/spath_s_aero_rain'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_aero_pris_<CASE>.h5' '/aero_pris' '/spath_aero_pris_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_aero_pris_<CASE>.h5' '/aero_pris' '/spath_i_aero_pris'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_aero_pris_<CASE>.h5' '/aero_pris' '/spath_ps_aero_pris'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_aero_pris_<CASE>.h5' '/aero_pris' '/spath_s_aero_pris'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_aero_snow_<CASE>.h5' '/aero_snow' '/spath_aero_snow_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_aero_snow_<CASE>.h5' '/aero_snow' '/spath_i_aero_snow'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_aero_snow_<CASE>.h5' '/aero_snow' '/spath_ps_aero_snow'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_aero_snow_<CASE>.h5' '/aero_snow' '/spath_s_aero_snow'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_aero_aggr_<CASE>.h5' '/aero_aggr' '/spath_aero_aggr_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_aero_aggr_<CASE>.h5' '/aero_aggr' '/spath_i_aero_aggr'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_aero_aggr_<CASE>.h5' '/aero_aggr' '/spath_ps_aero_aggr'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_aero_aggr_<CASE>.h5' '/aero_aggr' '/spath_s_aero_aggr'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_aero_graup_<CASE>.h5' '/aero_graup' '/spath_aero_graup_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_aero_graup_<CASE>.h5' '/aero_graup' '/spath_i_aero_graup'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_aero_graup_<CASE>.h5' '/aero_graup' '/spath_ps_aero_graup'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_aero_graup_<CASE>.h5' '/aero_graup' '/spath_s_aero_graup'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_aero_hail_<CASE>.h5' '/aero_hail' '/spath_aero_hail_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_spath_aero_hail_<CASE>.h5' '/aero_hail' '/spath_i_aero_hail'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_spath_aero_hail_<CASE>.h5' '/aero_hail' '/spath_ps_aero_hail'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_spath_aero_hail_<CASE>.h5' '/aero_hail' '/spath_s_aero_hail'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_spath_aero_mass_<CASE>.h5' '/aero_mass' '/spath_sum_aero_mass'  '_HZT'   ''        '' }
%
%        % SAL sample region (large)
%        { 'TsAveragedData/hda_sal_aero_mass_<CASE>.h5'   '/aero_mass'   '/sal_sum_aero_mass'   '_HZT'   ''        '' }
%
%        % storm region
%        { 'TsAveragedData/hda_all_aero_mass_<CASE>.h5' '/aero_mass' '/all_aero_mass_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_aero_mass_<CASE>.h5' '/aero_mass' '/all_i_aero_mass'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_aero_mass_<CASE>.h5' '/aero_mass' '/all_ps_aero_mass'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_aero_mass_<CASE>.h5' '/aero_mass' '/all_s_aero_mass'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_all_aero_hydro_<CASE>.h5' '/aero_hydro' '/all_aero_hydro_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_aero_hydro_<CASE>.h5' '/aero_hydro' '/all_i_aero_hydro'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_aero_hydro_<CASE>.h5' '/aero_hydro' '/all_ps_aero_hydro'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_aero_hydro_<CASE>.h5' '/aero_hydro' '/all_s_aero_hydro'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_all_aero_cloud_<CASE>.h5' '/aero_cloud' '/all_aero_cloud_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_aero_cloud_<CASE>.h5' '/aero_cloud' '/all_i_aero_cloud'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_aero_cloud_<CASE>.h5' '/aero_cloud' '/all_ps_aero_cloud'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_aero_cloud_<CASE>.h5' '/aero_cloud' '/all_s_aero_cloud'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_all_aero_rain_<CASE>.h5' '/aero_rain' '/all_aero_rain_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_aero_rain_<CASE>.h5' '/aero_rain' '/all_i_aero_rain'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_aero_rain_<CASE>.h5' '/aero_rain' '/all_ps_aero_rain'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_aero_rain_<CASE>.h5' '/aero_rain' '/all_s_aero_rain'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_all_aero_pris_<CASE>.h5' '/aero_pris' '/all_aero_pris_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_aero_pris_<CASE>.h5' '/aero_pris' '/all_i_aero_pris'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_aero_pris_<CASE>.h5' '/aero_pris' '/all_ps_aero_pris'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_aero_pris_<CASE>.h5' '/aero_pris' '/all_s_aero_pris'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_all_aero_snow_<CASE>.h5' '/aero_snow' '/all_aero_snow_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_aero_snow_<CASE>.h5' '/aero_snow' '/all_i_aero_snow'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_aero_snow_<CASE>.h5' '/aero_snow' '/all_ps_aero_snow'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_aero_snow_<CASE>.h5' '/aero_snow' '/all_s_aero_snow'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_all_aero_aggr_<CASE>.h5' '/aero_aggr' '/all_aero_aggr_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_aero_aggr_<CASE>.h5' '/aero_aggr' '/all_i_aero_aggr'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_aero_aggr_<CASE>.h5' '/aero_aggr' '/all_ps_aero_aggr'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_aero_aggr_<CASE>.h5' '/aero_aggr' '/all_s_aero_aggr'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_all_aero_graup_<CASE>.h5' '/aero_graup' '/all_aero_graup_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_aero_graup_<CASE>.h5' '/aero_graup' '/all_i_aero_graup'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_aero_graup_<CASE>.h5' '/aero_graup' '/all_ps_aero_graup'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_aero_graup_<CASE>.h5' '/aero_graup' '/all_s_aero_graup'   '_hZt'   'sal'     '' }
%
%        { 'TsAveragedData/hda_all_aero_hail_<CASE>.h5' '/aero_hail' '/all_aero_hail_ts'  '_hZT'   ''        '' }
%        { 'TsAveragedData/hda_all_aero_hail_<CASE>.h5' '/aero_hail' '/all_i_aero_hail'   '_hZt'   'init'    '' }
%        { 'TsAveragedData/hda_all_aero_hail_<CASE>.h5' '/aero_hail' '/all_ps_aero_hail'  '_hZt'   'pre_sal' '' }
%        { 'TsAveragedData/hda_all_aero_hail_<CASE>.h5' '/aero_hail' '/all_s_aero_hail'   '_hZt'   'sal'     '' }
%      }
%      'DIAGS/hda_meas_ts_aero_<CASE>.h5'
%    }
%
%    % Condensate
%    {
%      'Total Condendsate Tsavg'
%      {
%        % SAL sample region (large)
%        { 'TsAveragedData/hda_sal_tcond_<CASE>.h5'    '/tcond'    '/sal_sum_tcond_mass'    '_HZT'   ''        '' }
%      }
%      'DIAGS/hda_meas_ts_cond_<CASE>.h5'
%    }
%
%    % accum precip
%    {
%      'Accumulated Precip Tsavg'
%      {
%        % SAL sample region (large)
%        { 'TsAveragedData/hda_sal_accpcp_<CASE>.h5'    '/accpcp'    '/sal_sum_accpcp_mass'    '_H_T'   ''        '' }
%      }
%      'DIAGS/hda_meas_ts_precip_<CASE>.h5'
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
        OutVname  = MeasList{imeas}{3};
        Rspec     = num2cell(MeasList{imeas}{4});
        Trange    = MeasList{imeas}{5};
        Zrange    = MeasList{imeas}{6};

        % skip this profile set if doing dust and on a NODUST case
        if ((~isempty(regexp(Case, 'NODUST'))) && ...
            ((~isempty(regexp(Vname, '/d[12]_num'))) || ...
             (~isempty(regexp(Vname, '/d[12]_mass'))) || ...
             (~isempty(regexp(Vname, '/dust_'))) || ...
             (~isempty(regexp(Vname, '/dustifn_'))) || ...
             (~isempty(regexp(Vname, '/tracer[12]'))) || ...
             (~isempty(regexp(Vname, '/trdust[12]_diff')))))
          continue
        else
          icount = icount + 1;
        end

        InFile = regexprep(Ftemplate, '<CASE>', Case);

        % Parse the reduction spec. Always has 4 characters that correspond to x, y, z, t.
        Hindex = 0;
        Zindex = 0;
        Tindex = 0;

        Hreduce = false;
        Zreduce = false;
        Treduce = false;

        InForm = {};
        Ndims = 0;
        InDims = { 'x' 'y' 'z' 't' }';
        for i = 1:length(Rspec)
          switch(Rspec{i})
            case { 'h' 'H' }
              HinVar = InDims{i};
              Ndims = Ndims + 1;
              Hindex = Ndims;
              InForm{Ndims} = 'h';
              if (strcmp(Rspec{i}, 'h'))
                Hreduce = true;
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
        clear H1;
        clear H2;
        clear Z1
        clear Z2
        clear T1;
        clear T2;

        fprintf('        Selection:\n');

        SelectSpec = '';
        % HDA
        if (Hindex > 0)
          switch(HinVar)
            case 'x'
              H = X;
            case 'y'
              H = Y;
            case 'z'
              H = Z;
            case 't'
              H = T;
          end

          if (Hreduce)
            fprintf('          Hda: Average\n');
          else
            fprintf('          Hda: Sum\n');
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
          if (strcmp(SelectSpec, ''))
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
            case 'init'
              T1 = find(ST >= TstartInit, 1, 'first');
              T2 = find(ST <= TendInit,   1, 'last');
  
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
          if (strcmp(SelectSpec, ''))
            SelectSpec = sprintf('%d:%d', T1, T2);
          else
            SelectSpec = sprintf('%s,%d:%d', SelectSpec, T1, T2);
          end
        end

        % Trim down the input data according to the selection indices
        SelectCmd = sprintf('HDATA(1,%s)', SelectSpec);
        SUM = squeeze(eval(SelectCmd));
        SelectCmd = sprintf('HDATA(2,%s)', SelectSpec);
        COUNT = squeeze(eval(SelectCmd));
        fprintf('\n');

        % Adjust indices since we just eliminated the h dimension
        if (Zindex > Hindex)
          Zindex = Zindex - 1;
        end
        if (Tindex > Hindex)
          Tindex = Tindex - 1;
        end
        Hindex = 0;
        Ndims = Ndims - 1;

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
        %    1. Accumulate sum data (dim t)
        %    2. Reduce hda values (dim h)
        %    3. Select z level (dim z)
        %
        % Keep original dimensions intact until finished so that [RHZT]index vars
        % can remain constant.
        fprintf('        Reduction:\n');

        if (Treduce)
          fprintf('          Summing temporal data\n');
          SUM   = squeeze(sum(SUM, Tindex));
          COUNT = squeeze(sum(COUNT, Tindex));

          % Adjust indices, number of dimensions
          if (Zindex > Tindex)
            Zindex = Zindex - 1;
          end
          Tindex = 0;
          Ndims = Ndims - 1;
        end

        % We've already eliminated the h dimension and adjusted the z and t indices
        if (Hreduce)
          fprintf('          Reducing hda to averages\n');
          MEAS = SUM ./ COUNT;
        else
          fprintf('          Reducing hda to sums\n');
          MEAS = SUM;
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

        % if writing out sums (Hreduce is false), then also write out the counts
        if (~Hreduce)
          OutVnameCounts = sprintf('%s_counts', OutVname);
          fprintf('      Writing: %s (%s)\n', OutFile, OutVnameCounts)

          h5create(OutFile, OutVnameCounts, OutSize);
          h5write(OutFile, OutVnameCounts, COUNT);
          AttachDimensionsXyzt(OutFile, OutVnameCounts, DimOrder, Xname, Yname, Zname, Tname);
        end
        fprintf('\n');

      end % measurements
    end % sets
  end % cases
end % function
