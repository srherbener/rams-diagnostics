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

function main(args)

  gExp     = subwrd(args, 1)
  gDir     = subwrd(args, 2)
  gCase    = subwrd(args, 3)
  varNum   = subwrd(args, 4)
  timeStep = subwrd(args, 5)
  pFile    = subwrd(args, 6)

* convert timestep number to hours
  timeStr = 36.0 + ((timeStep - 1.0) * 0.5)
  timeStr = timeStr' hrs'

* set up variable names, titles, etc. according to varNum
  if (varNum = 1)
*   column integrated liquid
    varName = 'cint_liq'
    gVar = 'liquid_colint'
    gClevs = '80.0 100.0 200.0 400.0 600.0 800.0 1000.0 2000.0 4000.0 6000.0 8000.0 10000.0 20000.0'
    gCcols = '   9    14     4    11     5    13      3     10      7     12      8       2       6'
    gTitle = gExp': LWP (g/m**2), 'gCase', t 'timeStr
  endif
  if (varNum = 2)
*   column integrated ice
    varName = 'cint_ice'
    gVar = 'ice_colint'
    gClevs = '6.0 8.0 10.0 20.0 40.0 60.0 80.0 100.0 200.0 400.0 600.0 800.0 1000.0'
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
    gClevs = '0.0 50.0 100.0 150.0 200.0 300.0 400.0 500.0 1000.0 1500.0 2000.0 2500.0 3000.0'
    gCcols = '  9   14     4    11     5    13     3    10      7     12      8      2      6'
    gTitle = gExp': CCN (#/cc), z=2300m, 'gCase', t 'timeStr
  endif
  if (varNum = 5)
*   surface wind
    varName = 'sfc_wind'
    gVar = 'sfcwind_mag'
    gClevs = '10.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0 70.0'
    gCcols = '   9   14    4   11    5   13    3   10    7   12    8    2    6'
    gTitle = gExp': Surface wind speed (m/s), 'gCase', t 'timeStr
  endif
  if (varNum = 6)
*   column integrated supercooled liquid
    varName = 'cint_liq_sc'
    gVar = 'liquid_colint'
    gClevs = '6.0 8.0 10.0 20.0 40.0 60.0 80.0 100.0 200.0 400.0 600.0 800.0 1000.0'
    gCcols = '  9  14    4   11    5   13    3    10     7    12     8     2      6'
    gTitle = gExp': Supercooled LWP (g/m**2), 'gCase', t 'timeStr
  endif

  gcFile = gDir'/'varName'_'gCase'.ctl'

  xTitle = 'Longitude'
  yTitle = 'Latitude'

  say 'Plotting horizontal slice of GRADS data:'
  say '  Experiment = 'gExp
  say '  GRADS control file = 'gcFile
  say '  Variable number = 'varNum' --> 'varName
  say '  Timestep = 'timeStep' --> 'timeStr
  say '  Output plot (GIF) file = 'pFile

  'reinit'
  'clear'
  'open 'gcFile
  'set t 'timeStep
  'set gxout contour'
  'set clevs 'gClevs
  'set ccols 'gCcols
  'd 'gVar
  'draw title 'gTitle
  'draw xlab 'xTitle
  'draw ylab 'yTitle
  'printim 'pFile' white'

*  'cbarn'

return
