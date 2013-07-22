************************************************************************
* plot out time series of single point data
*
* print white background so the colored contours will stand out better
*
* varNum arg:
*   1 --> supercooled cloud droplets
*   2 --> supercooled cloud droplet diameters
*   3 --> precip rate
*   4 --> vertical velocity in supercooled cloud regions
*   5 --> average cloud droplet concentration near the eyewall region

function main(args)

  gExp        = subwrd(args, 1)
  gDir        = subwrd(args, 2)
  gCase       = subwrd(args, 3)
  varNum      = subwrd(args, 4)
  startTime   = subwrd(args, 5)
  timeInc     = subwrd(args, 6)
  plotTstart  = subwrd(args, 7)
  plotTend    = subwrd(args, 8)
  yMin        = subwrd(args, 9)
  yMax        = subwrd(args, 10)
  pFile       = subwrd(args, 11)

* convert plot times (in hrs) to a time step number 
  gTstart = ((plotTstart - startTime) / timeInc) + 1
  gTend = ((plotTend - startTime) / timeInc) + 1

* set up variable names, titles, etc. according to varNum
  if (varNum = 1)
*   supercooled cloud droplet mass
    varName = 'ts_sc_cloud'
    gVar = 'sc_cloud_ts_joi'
    gTitle = gExp': Total supercooled cloud droplet mass, 'gCase
    yTitle = 'SC Cloud Mass (g)'
    avgLen = 4
  endif
  if (varNum = 2)
*   supercooled cloud droplet mean diameter
    varName = 'ts_sc_cloud_diam'
    gVar = 'sc_cloud_diam_t'
    gTitle = gExp': Mean supercooled droplet diameter, 'gCase
    yTitle = 'SC Cloud Mean Diameter (um)'
    avgLen = 4
  endif
  if (varNum = 3)
*   total precip rate
    varName = 'ts_precipr'
    gVar = 'precipr_ts_join'
    gTitle = gExp': Precipitation rate, 'gCase
    yTitle = 'Precip Rate (kg/hr)'
    avgLen = 4
  endif
  if (varNum = 4)
*   vertical velocity
    varName = 'ts_w_up'
    gVar = 'w_up_ts_join'
    gTitle = gExp': Average w in updraft regions, 'gCase
    yTitle = 'Average w (m/s)'
    avgLen = 4
  endif
  if (varNum = 5)
*   eyewall clound droplet concentration
    varName = 'ts_ew_cloud'
    gVar = 'ew_cloud_ts_joi'
    gTitle = gExp': Average cloud droplet conc. near eyewall, 'gCase
    yTitle = 'Average N (#/cc)'
    avgLen = 4
  endif
  if (varNum = 6)
*   supercooled cloud droplet mass, rates
    varName = 'ts_sc_cloud_rates'
    gVar = 'sc_cloud_dt_joi'
    gTitle = gExp': Total supercooled cloud droplet mass dt, 'gCase
    yTitle = 'SC Cloud Mass Rate (g/s)'
    avgLen = 4
  endif
  if (varNum = 7)
*   supercooled cloud droplet mean diameter, rates
    varName = 'ts_sc_cloud_diam_rates'
    gVar = 'sc_cloud_diam_d'
    gTitle = gExp': Mean supercooled droplet diameter dt, 'gCase
    yTitle = 'SC Cloud Mean Diameter Rates (um/s)'
    avgLen = 4
  endif
  if (varNum = 8)
*   total precip rate, rates
    varName = 'ts_precipr_rates'
    gVar = 'precipr_dt_join'
    gTitle = gExp': Precipitation rate dt, 'gCase
    yTitle = 'Precip Rate Rates (kg/hr/s)'
    avgLen = 4
  endif
  if (varNum = 9)
*   vertical velocity, rates
    varName = 'ts_w_up_rates'
    gVar = 'w_up_dt_join'
    gTitle = gExp': Average w in updraft regions dt, 'gCase
    yTitle = 'Average w (m/s/s)'
    avgLen = 4
  endif
  if (varNum = 10)
*   eyewall clound droplet concentration, rates
    varName = 'ts_ew_cloud_rates'
    gVar = 'ew_cloud_dt_joi'
    gTitle = gExp': Average cloud droplet conc. near eyewall dt, 'gCase
    yTitle = 'Average N (#/cc/s)'
    avgLen = 4
  endif
  if (varNum = 11)
*   total kinetic energy
    varName = 'ts_horiz_ke'
    gVar = 'horiz_ke_ts_joi'
    gTitle = gExp': Total kinetic energy, 'gCase
    yTitle = 'Total KE (J)'
    avgLen = 4
  endif
  if (varNum = 12)
*   storm intensity metric
    varName = 'ts_storm_int'
    gVar = 'storm_int_ts_jo'
    gTitle = gExp': Storm intensity metric, 'gCase
    yTitle = 'Intensity'
    avgLen = 4
  endif
  if (varNum = 13)
*   average ccn concentration
    varName = 'ts_ccn_concen'
    gVar = 'ccnconc_ts_join'
    gTitle = gExp': Avg CCN concentration, 'gCase
    yTitle = 'Avg CCN Conc (#/cc)'
    avgLen = 1
  endif
  if (varNum = 14)
*   average ccn concentration rates
    varName = 'ts_ccn_concen_rates'
    gVar = 'ccnconc_dt_join'
    gTitle = gExp': Avg CCN concentration, dt, 'gCase
    yTitle = 'Avg CCN Conc Rates (#/cc/s)'
    avgLen = 1
  endif
  if (varNum = 15)
*   warm cloud droplet mass
    varName = 'ts_wr_cloud'
    gVar = 'wr_cloud_ts_joi'
    gTitle = gExp': Total warm cloud droplet mass, 'gCase
    yTitle = 'Warm Cloud Mass (g)'
    avgLen = 4
  endif
  if (varNum = 16)
*   warm cloud droplet mass, rates
    varName = 'ts_wr_cloud_rates'
    gVar = 'wr_cloud_dt_joi'
    gTitle = gExp': Total warm cloud droplet mass dt, 'gCase
    yTitle = 'Warm Cloud Mass Rate (g/s)'
    avgLen = 4
  endif
  if (varNum = 17)
*   warm cloud droplet mean diameter
    varName = 'ts_wr_cloud_diam'
    gVar = 'wr_cloud_diam_t'
    gTitle = gExp': Mean warm droplet diameter, 'gCase
    yTitle = 'Warm Cloud Mean Diameter (um)'
    avgLen = 4
  endif
  if (varNum = 18)
*   warm cloud droplet mean diameter, rates
    varName = 'ts_wr_cloud_diam_rates'
    gVar = 'wr_cloud_diam_d'
    gTitle = gExp': Mean warm droplet diameter dt, 'gCase
    yTitle = 'Warm Cloud Mean Diameter Rates (um/s)'
    avgLen = 4
  endif
  if (varNum = 19)
*   warm cloud2 droplet mass
    varName = 'ts_wr_cloud2'
    gVar = 'wr_cloud2_ts_jo'
    gTitle = gExp': Total warm cloud2 droplet mass, 'gCase
    yTitle = 'Warm Cloud2 Mass (g)'
    avgLen = 4
  endif
  if (varNum = 20)
*   warm cloud2 droplet mass, rates
    varName = 'ts_wr_cloud2_rates'
    gVar = 'wr_cloud2_dt_jo'
    gTitle = gExp': Total warm cloud2 droplet mass dt, 'gCase
    yTitle = 'Warm Cloud2 Mass Rate (g/s)'
    avgLen = 4
  endif
  if (varNum = 21)
*   warm cloud2 droplet mean diameter
    varName = 'ts_wr_cloud2_diam'
    gVar = 'wr_cloud2_diam_'
    gTitle = gExp': Mean warm droplet diameter, 'gCase
    yTitle = 'Warm Cloud2 Mean Diameter (um)'
    avgLen = 4
  endif
  if (varNum = 22)
*   warm cloud2 droplet mean diameter, rates
    varName = 'ts_wr_cloud2_diam_rates'
    gVar = 'wr_cloud2_diam_'
    gTitle = gExp': Mean warm droplet diameter dt, 'gCase
    yTitle = 'Warm Cloud2 Mean Diameter Rates (um/s)'
    avgLen = 4
  endif
  if (varNum = 23)
*   eyewall cloud2 droplet concentration
    varName = 'ts_ew_cloud2'
    gVar = 'ew_cloud2_ts_jo'
    gTitle = gExp': Average cloud2 droplet conc. near eyewall, 'gCase
    yTitle = 'Average N (#/cc)'
    avgLen = 4
  endif
  if (varNum = 24)
*   eyewall cloud2 droplet concentration, rates
    varName = 'ts_ew_cloud2_rates'
    gVar = 'ew_cloud2_dt_jo'
    gTitle = gExp': Average cloud2 droplet conc. near eyewall dt, 'gCase
    yTitle = 'Average N (#/cc/s)'
    avgLen = 4
  endif

  gcFile = gDir'/'varName'_'gCase'.ctl'

  xTitle = 'Time'

  say 'Plotting time series of given variable:'
  say '  Experiment = 'gExp
  say '  GRADS control file = 'gcFile
  say '  Start time = 'plotTstart'(hr) --> 'gTstart
  say '  End time = 'plotTend'(hr) --> 'gTend
  say '  Running average length = 'avgLen
  say '  Ymin = 'yMin
  say '  Ymax = 'yMax
  say '  Variable number = 'varNum' --> 'varName
  say '  Output plot (GIF) file = 'pFile

* apply smoothing to the time series
* Apply low pass filter to time series (5 point average)
  t1 = gTstart + avgLen
  t2 = gTend - avgLen
  'reinit'
  'clear'
  'open 'gcFile
  'set t 't1' 't2
  'set grads off'
  'set vrange 'yMin' 'yMax
  'd tloop(ave('gVar',t-'avgLen',t+'avgLen'))'
  'draw title 'gTitle
  'draw xlab 'xTitle
  'draw ylab 'yTitle
  'printim 'pFile' white'

return
