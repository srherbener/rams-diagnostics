************************************************************************
* plot out horzontal slice data
*
* print white background so the colored contours will stand out better
*
* varNum arg:
*   1 --> column integrated liquid
*   2 --> column integrated ice
*   3 --> GCCN
*   4 --> CCN
*   5 --> surface wind
*   6 --> column integrated supercooled liquid
*   7 --> column integrated warm rain liquid
*   8 --> tangential wind
*   9 --> radial wind

function main(args)

  gExp      = subwrd(args, 1)
  gDir      = subwrd(args, 2)
  gCase     = subwrd(args, 3)
  varNum    = subwrd(args, 4)
  startTime = subwrd(args, 5)
  timeInc   = subwrd(args, 6)
  plotTime  = subwrd(args, 7)
  pFile     = subwrd(args, 8)

* convert plotTime (in hrs) to a timeStep number 
  timeStep = ((plotTime - startTime) / timeInc) + 1
  timeStr = plotTime' hrs'

* set up variable names, titles, etc. according to varNum
  if (varNum = 1)
*   column integrated liquid
    varName = 'cint_liq'
    gVar = 'liquid_colint'
    gClevs = '200.0 400.0 600.0 800.0 1000.0 2000.0 4000.0 6000.0 8000.0 10000.0 20000.0 40000.0 60000.0'
    gCcols = '   9    14     4    11     5    13      3     10      7     12      8       2       6'
    gTitle = gExp': LWP (g/m**2), 'gCase', t 'timeStr
  endif
  if (varNum = 2)
*   column integrated ice
    varName = 'cint_ice'
    gVar = 'ice_colint'
    gClevs = '10.0 20.0 40.0 60.0 80.0 100.0 200.0 400.0 600.0 800.0 1000.0 2000.0 4000.0'
    gCcols = '  9  14    4   11    5   13    3    10     7    12     8     2      6'
    gTitle = gExp': Column integrated ice (g/m**2), 'gCase', t 'timeStr
  endif
  if (varNum = 3)
*   GCCN
    varName = 'hslice_2334_gccn'
    gVar = 'gccnconcen_hsli'
    gClevs = '0.00001 0.0001 0.001 0.01 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6'
    gCcols = '      9     14     4   11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': GCCN (#/cc), z=2300m, 'gCase', t 'timeStr
  endif
  if (varNum = 4)
*   CCN
    varName = 'hslice_2334_ccn'
    gVar = 'ccnconcen_hslic'
    gClevs = '20.0 40.0 60.0 80.0 100.0 200.0 400.0 600.0 800.0 1000.0 2000.0 4000.0 6000.0'
    gCcols = '  9   14     4    11     5    13     3    10      7     12      8      2      6'
    gTitle = gExp': CCN (#/cc), z=2300m, 'gCase', t 'timeStr
  endif
  if (varNum = 5)
*   surface wind
    varName = 'sfc_wind'
    gVar = 'sfcwind_mag'
    gClevs = '15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0 70.0 75.0'
    gCcols = '   9   14    4   11    5   13    3   10    7   12    8    2    6'
    gTitle = gExp': Surface wind speed (m/s), 'gCase', t 'timeStr
  endif
  if (varNum = 6)
*   column integrated supercooled liquid
    varName = 'cint_liq_sc'
    gVar = 'liquid_colint'
    gClevs = '10.0 20.0 40.0 60.0 80.0 100.0 200.0 400.0 600.0 800.0 1000.0 2000.0 4000.0'
    gCcols = '  9  14    4   11    5   13    3    10     7    12     8     2      6'
    gTitle = gExp': Supercooled LWP (g/m**2), 'gCase', t 'timeStr
  endif
  if (varNum = 7)
*   column integrated warm liquid
    varName = 'cint_liq_wr'
    gVar = 'liquid_colint'
    gClevs = '200.0 400.0 600.0 800.0 1000.0 2000.0 4000.0 6000.0 8000.0 10000.0 20000.0 40000.0 60000.0'
    gCcols = '  9  14    4   11    5   13    3    10     7    12     8     2      6'
    gTitle = gExp': Warm LWP (g/m**2), 'gCase', t 'timeStr
  endif
  if (varNum = 8)
*   surface tangential wind
    varName = 'storm_winds_tan'
    gVar = 'wind_speed_t'
    gClevs = '5.0 10.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0'
    gCcols = '  9  14    4   11    5   13    3    10     7    12     8     2      6'
    gTitle = gExp': Tan Wind Speed (m/s), 'gCase', t 'timeStr
  endif
  if (varNum = 9)
*   surface radial wind
    varName = 'storm_winds_rad'
    gVar = 'wind_speed_r'
    gClevs = '-18.0 -15.0 -12.0 -9.0 -6.0 -3.0 0.0 3.0 6.0 9.0 12.0 15.0 18.0'
    gCcols = '  9  14    4   11    5   13    3    10     7    12     8     2      6'
    gTitle = gExp': Rad Wind Speed (m/s), 'gCase', t 'timeStr
  endif

  gcFile = gDir'/'varName'_'gCase'.ctl'

  xTitle = 'Longitude'
  yTitle = 'Latitude'

  say 'Plotting horizontal slice of GRADS data:'
  say '  Experiment = 'gExp
  say '  GRADS control file = 'gcFile
  say '  Variable number = 'varNum' --> 'varName
  say '  Timestep = 'timeStr' --> 'timeStep
  say '  Output plot (GIF) file = 'pFile

  'reinit'
  'clear'
  'open 'gcFile
  'set t 'timeStep
  if ((varNum = 8) | (varNum = 9))
    'set z 1'
  endif
  if (gCase = 'INIT')
    'set lon -41 -39'
    'set lat 14 16'
  endif
  'set gxout shaded'
  'set clevs 'gClevs
  'set ccols 'gCcols
  'set grads off'
  'd 'gVar
  'draw title 'gTitle
  'draw xlab 'xTitle
  'draw ylab 'yTitle
  'cbarn 1.0 1'
  'printim 'pFile' white'


return
