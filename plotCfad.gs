************************************************************************
* plot out cfad data
*
* print white background so the colored contours will stand out better
*
* varNum arg:
*   1 --> vertical velocity

function main(args)

  gExp     = subwrd(args, 1)
  gDir     = subwrd(args, 2)
  gCase    = subwrd(args, 3)
  varNum   = subwrd(args, 4)
  rBand    = subwrd(args, 5)
  timeStep = subwrd(args, 6)
  pFile    = subwrd(args, 7)

* convert timestep number to hours
  timeStr = 36.0 + ((timeStep - 1.0) * 0.5)
  timeStr = timeStr' hrs'

* convert rband to radius value
  numRbands = 6
  totalRadius = 210.0
  rbandInc = totalRadius / numRbands
  rbandStart = (rBand - 1) * rbandInc
  rbandEnd = rBand * rbandInc

* set up variable names, titles, etc. according to varNum
  if (varNum = 1)
*   vertical velocity
    varName = 'w'
    gVar = 'w_cfad'
*    gClevs = '0.0 2.0 4.0 6.0 8.0 10.0 15.0 20.0 25.0 30.0 35.0 40.0 50.0'
    gClevs = '0.0 5.0 10.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0'
    gCcols = '  9  14    4   11    5   13    3   10    7   12    8    2    6'
    gTitle = gExp': CFAD: w (%), 'gCase', r: 'rbandStart'-'rbandEnd' km, t: 'timeStr
    zTop = 10000
  endif

  gcFile = gDir'/cfad_'varName'_'gCase'.ctl'

  xTitle = 'W (m/s)'
  yTitle = 'Height (m)'

  say 'Plotting CFAD:'
  say '  Experiment = 'gExp
  say '  GRADS control file = 'gcFile
  say '  Variable number = 'varNum' --> 'varName
  say '  Radial band number = 'rBand' --> 'rbandStart'-'rbandEnd
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
  'set x 'rBand
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
