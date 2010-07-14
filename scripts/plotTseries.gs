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

function main(args)

  gExp        = subwrd(args, 1)
  gDir        = subwrd(args, 2)
  gCase       = subwrd(args, 3)
  varNum      = subwrd(args, 4)
  startTime   = subwrd(args, 5)
  timeInc     = subwrd(args, 6)
  plotTstart  = subwrd(args, 7)
  plotTend    = subwrd(args, 8)
  pFile       = subwrd(args, 9)

* convert plot times (in hrs) to a time step number 
  gTstart = ((plotTstart - startTime) / timeInc) + 1
  gTend = ((plotTend - startTime) / timeInc) + 1

* set up variable names, titles, etc. according to varNum
  if (varNum = 1)
*   supercooled cloud droplet mass
    varName = 'ts_sc_cloud'
    gVar = 'sc_cloud_test'
    gTitle = gExp': Total supercooled cloud droplet mass (g), 'gCase
    yTitle = 'SC Cloud Mass (g)'
  endif
  if (varNum = 2)
*   supercooled cloud droplet mean diameter
    varName = 'ts_sc_cloud_diam'
    gVar = 'sc_cloud_diam_t'
    gTitle = gExp': Mean supercooled cloud droplet diameter (um), 'gCase
    yTitle = 'SC Cloud Mean Diameter (um)'
  endif
  if (varNum = 3)
*   total precip rate
    varName = 'ts_precipr'
    gVar = 'precipr_test'
    gTitle = gExp': Total precip rate (kg/hr), 'gCase
    yTitle = 'Precip Rate (kg/hr)'
  endif
  if (varNum = 4)
*   vertical velocity
    varName = 'ts_sc_w'
    gVar = 'sc_w_test'
    gTitle = gExp': Average w (m/s), in supercooled cloud region, 'gCase
    yTitle = 'Average w (m/s)'
  endif

  gcFile = gDir'/'varName'_'gCase'.ctl'

  xTitle = 'Time'

  say 'Plotting time series of given variable:'
  say '  Experiment = 'gExp
  say '  GRADS control file = 'gcFile
  say '  Start time = 'plotTstart'(hr) --> 'gTstart
  say '  End time = 'plotTend'(hr) --> 'gTend
  say '  Variable number = 'varNum' --> 'varName
  say '  Output plot (GIF) file = 'pFile

* apply smoothing to the time series
  'reinit'
  'clear'
  'open 'gcFile
  'set t 'gTstart' 'gTend
  'set grads off'
  'd smth9('gVar')'
  'draw title 'gTitle
  'draw xlab 'xTitle
  'draw ylab 'yTitle
  'printim 'pFile' white'

return
