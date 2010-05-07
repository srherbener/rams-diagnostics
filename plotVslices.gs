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
*   9 --> T (deg C)
*  10 --> cloud water content

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
*   vertical velocity
    varName = 'w'
    gVar = 'w_azavg'
    gClevs = '-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5'
    gCcols = '   9   14    4   11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': AZ: w (m/s), 'gCase', t 'timeStr
    zTop = 16000
  endif
  if (varNum = 2)
*   horiz velocity - tangential
    varName = 'speed_t'
    gVar = 'speed_t_azavg'
    gClevs = '10.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0 70.0'
    gCcols = '   9   14    4   11    5   13    3   10    7   12    8    2    6'
    gTitle = gExp': AZ: spd_t (m/s), 'gCase', t 'timeStr
    zTop = 16000
  endif
  if (varNum = 3)
*   horiz velocity - radial
    varName = 'speed_r'
    gVar = 'speed_r_azavg'
    gClevs = '-18.0 -15.0 -12.0 -9.0 -6.0 -3.0 0.0 3.0 6.0 9.0 12.0 15.0 18.0'
    gCcols = '    9    14     4   11    5   13   3  10   7  12    8    2    6'
    gTitle = gExp': AZ: spd_r (m/s), 'gCase', t 'timeStr
    zTop = 16000
  endif
  if (varNum = 4)
*   liquid water content
    varName = 'liquid'
    gVar = 'liquid_azavg'
    gClevs = '0.0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8'
    gCcols = '  9  14   4  11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': AZ: LWC (g/kg), 'gCase', t 'timeStr
    zTop = 10000
  endif
  if (varNum = 5)
*   ice content
    varName = 'ice'
    gVar = 'ice_azavg'
    gClevs = '0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4'
    gCcols = '  9  14   4  11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': AZ: ice (g/kg), 'gCase', t 'timeStr
    zTop = 16000
  endif
  if (varNum = 6)
*   GCCN
    varName = 'gccnconcen'
    gVar = 'gccnconcen_azav'
    gClevs = '0.00001 0.0001 0.001 0.01 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6'
    gCcols = '      9     14     4   11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': AZ: GCCN (#/cc), 'gCase', t 'timeStr
    zTop = 16000
  endif
  if (varNum = 7)
*   CCN
    varName = 'ccnconcen'
    gVar = 'ccnconcen_azavg'
    gClevs = '0.0 50.0 100.0 150.0 200.0 300.0 400.0 500.0 1000.0 1500.0 2000.0 2500.0 3000.0'
    gCcols = '  9   14     4    11     5    13     3    10      7     12      8      2      6'
    gTitle = gExp': AZ: CCN (#/cc), 'gCase', t 'timeStr
    zTop = 16000
  endif
  if (varNum = 8)
*   Thetae
    varName = 'thetae'
    gVar = 'thetae_azavg'
    gClevs = '348.0 350.0 352.0 354.0 356.0 358.0 360.0 362.0 364.0 366.0 368.0 370.0 372.0'
    gCcols = '    9    14     4    11     5    13     3    10     7    12     8     2     6'
    gTitle = gExp': AZ: Theta-e (K), 'gCase', t 'timeStr
    zTop = 2000
  endif
  if (varNum = 9)
*   Temp (deg C)
    varName = 'tempc'
    gVar = 'tempc_azavg'
    gClevs = '18.0 20.0 21.0 22.0 22.5 23.0 23.5 24.0 24.5 25.0 25.5 26.0 26.5'
    gCcols = '   9   14    4   11    5   13    3   10    7   12    8    2    6'
    gTitle = gExp': AZ: T (C), 'gCase', t 'timeStr
    zTop = 2000
  endif
  if (varNum = 10)
*   relative humidity
    varName = 'relhum'
    gVar = 'relhum_azavg'
    gClevs = '80.0 82.0 84.0 86.0 88.0 90.0 92.0 94.0 96.0 97.0 98.0 99.0 100.0'
    gCcols = '   9   14    4   11    5   13    3   10    7   12    8    2    6'
    gTitle = gExp': AZ: RH (%), 'gCase', t 'timeStr
    zTop = 2000
  endif
  if (varNum = 11)
*   pressure
    varName = 'press'
    gVar = 'press_azavg'
    gClevs = '740.0 760.0 780.0 800.0 820.0 840.0 860.0 880.0 900.0 920.0 940.0 960.0 980.0'
    gCcols = '    9    14     4    11     5    13     3    10     7    12     8     2     6'
    gTitle = gExp': AZ: Pressure (mb), 'gCase', t 'timeStr
    zTop = 2000
  endif
  if (varNum = 12)
*   cloud water content
    varName = 'cloud'
    gVar = 'cloud_azavg'
    gClevs = '0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2'
    gCcols = '  9  14   4  11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': AZ: CWC (g/kg), 'gCase', t 'timeStr
    zTop = 10000
  endif
  if (varNum = 13)
*   cloud droplet number concentration
    varName = 'cloudconcen_cm3'
    gVar = 'cloudconcen_cm3'
*     gClevs = '0.0 200.0 400.0 600.0 800.0 1000.0 1200.0 1400.0 1600.0 1800.0 2000.0 2200.0 2400.0'
    gClevs = '0.0 20.00 40.00 60.00 80.00 100.00 120.00 140.00 160.00 180.00 200.00 220.00 240.00'
    gCcols = '  9  14   4  11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': AZ: Cloud Droplet Number (#/cc), 'gCase', t 'timeStr
*    zTop = 10000
    zTop = 4000
  endif
  if (varNum = 14)
*   cloud2 droplet number concentration
    varName = 'cloud2concen_cm'
    gVar = 'cloud2concen_cm'
    gClevs = '0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2'
    gCcols = '  9  14   4  11   5  13   3  10   7  12   8   2   6'
    gTitle = gExp': AZ: Cloud2 Droplet Number (#/cc), 'gCase', t 'timeStr
*     zTop = 10000
    zTop = 4000
  endif

  gcFile = gDir'/'gCase'_'varName'.ctl'

  xTitle = 'Radius (km)'
  yTitle = 'Height (m)'

  say 'Plotting vertical slice of azimuthal averaged data:'
  say '  Experiment = 'gExp
  say '  GRADS control file = 'gcFile
  say '  Variable number = 'varNum' --> 'varName
  say '  Timestep = 'timeStep' --> 'timeStr
  say '  Output plot (GIF) file = 'pFile

* The normal plot area in landscape is:
*     X: 1.5 to 10.5
*     Y: 1 to 8
* Need to push up the bottom of the plot to make room
* for the colorbar, so use parea to reset Y: 2 to 8
  'reinit'
  'clear'
  'open 'gcFile
  'set lev 0 'zTop
  'set t 'timeStep
  'set gxout shaded'
  'set clevs 'gClevs
  'set ccols 'gCcols
  'set xlab %.1f'
  'set grads off'
  'set parea 1.5 10.5 2 8'
  'd 'gVar
  'draw title 'gTitle
  'draw xlab 'xTitle
  'draw ylab 'yTitle
  'cbarn 1.0 0'
  'printim 'pFile' white'


return
