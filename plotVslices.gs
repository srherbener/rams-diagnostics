************************************************************************
* plot out vertical slice data of azimutally averaged data
*
* print white background so the colored contours will stand out better
*
* varNum arg:
*   1 --> vertical velocity
*   2 --> tangential horizontal velocity
*   3 --> radial horizontal velocity
*   4 --> liquid water content
*   5 --> ice content
*   6 --> GCCN
*   7 --> CCN
*   8 --> Theta-e

function main(args)

  gExp     = subwrd(args, 1)
  gDir     = subwrd(args, 2)
  gCase    = subwrd(args, 3)
  varNum   = subwrd(args, 4)
  timeStep = subwrd(args, 5)
  pFile    = subwrd(args, 6)

* set up variable names, titles, etc. according to varNum
  if (varNum = 1)
*   vertical velocity
    varName = 'w'
    gVar = 'w_azavg'
    gClevs = '-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0'
    gCcols = '   9   14    4   11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': AZ: w (m/s), 'gCase', t 'timeStep
  endif
  if (varNum = 2)
*   horiz velocity - tangential
    varName = 'speed_t'
    gVar = 'speed_t_azavg'
    gClevs = '-30.0 -20.0 -10.0 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0'
    gCcols = '9    14     4  11    5   13    3   10    7   12    8    2  6'
    gTitle = gExp': AZ: spd_t (m/s), 'gCase', t 'timeStep
  endif
  if (varNum = 3)
*   horiz velocity - radial
    varName = 'speed_r'
    gVar = 'speed_r_azavg'
    gClevs = '-30.0 -25.0 -20.0 -15.0 -10.0 -5.0 0.0 5.0 10.0 15.0 20.0 25.0 30.0'
    gCcols = '9    14     4    11     5   13   3  10    7   12    8    2    6'
    gTitle = gExp': AZ: spd_r (m/s), 'gCase', t 'timeStep
  endif
  if (varNum = 4)
*   liquid water content
    varName = 'liquid'
    gVar = 'liquid_azavg'
    gClevs = '0.0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8'
    gCcols = '9  14   4  11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': AZ: LWC (g/kg), 'gCase', t 'timeStep
  endif
  if (varNum = 5)
*   ice content
    varName = 'ice'
    gVar = 'ice_azavg'
    gClevs = '0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4'
    gCcols = '9  14   4  11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': AZ: ice (g/kg), 'gCase', t 'timeStep
  endif
  if (varNum = 6)
*   GCCN
    varName = 'gccnconcen'
    gVar = 'gccnconcen_azav'
    gClevs = '0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2'
    gCcols = '9  14   4  11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': AZ: GCCN (#/cc), 'gCase', t 'timeStep
  endif
  if (varNum = 7)
*   CCN
    varName = 'ccnconcen'
    gVar = 'ccnconcen_azavg'
    gClevs = '0.0 50.0 100.0 150.0 200.0 300.0 400.0 500.0 1000.0 1500.0 2000.0 2500.0 3000.0'
    gCcols = '9    14     4     11      5     13      3     10      7     12      8      2      6'
    gTitle = gExp': AZ: CCN (#/cc), 'gCase', t 'timeStep
  endif
  if (varNum = 8)
*   Thetae
    varName = 'thetae'
    gVar = 'thetae_azavg'
    gClevs = '300.0 310.0 320.0 330.0 340.0 350.0 360.0 370.0 380.0 390.0 400.0 410.0 420.0'
    gCcols = '9    14     4    11     5    13     3    10     7    12     8     2     6'
    gTitle = gExp': AZ: Theta-e (K), 'gCase', t 'timeStep
  endif

  gcFile = gDir'/'gCase'_'varName'.ctl'

  xTitle = 'Radius (km)'
  yTitle = 'Height (m)'

  say 'Plotting vertical slice of azimuthal averaged data:'
  say '  Experiment = 'gExp
  say '  GRADS control file = 'gcFile
  say '  Variable number = 'varNum' --> 'varName
  say '  Timestep = 'timeStep
  say '  Output plot (GIF) file = 'pFile

  'reinit'
  'clear'
  'open 'gcFile
  'set lev 0 16000'
  'set t 'timeStep
  'set gxout contour'
  'set clevs 'gClevs
  'set ccols 'gCcols
  'set xlab %.1f'
  'd 'gVar
  'draw title 'gTitle
  'draw xlab 'xTitle
  'draw ylab 'yTitle
  'printim 'pFile' white'

*  'cbarn'

return
