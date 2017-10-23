* script to form radial wind speed values
* assumes that u, v are in the grads1* files

function main(args)

  sCase    = subwrd(args, 1)
  sDir     = subwrd(args, 2)
  timeStep = subwrd(args, 3)
  sLat     = subwrd(args, 4)
  sLon     = subwrd(args, 5)
  pDir     = subwrd(args, 6)

  gcFile = 'TC_SEED_'sCase'/GRADS/grads1TC_SEED_'sCase'-AS-1998-08-22-120000-g3.ctl'

* convert timestep number to hours
  timeStr = 36.0 + ((timeStep - 1.0) * 0.5)
  timeStr = timeStr' hrs'

  say 'Plotting radial wind speed'
  say '  Case: 'sCase
  say '  Dir: 'sDir
  say '  Storm lat: 'sLat
  say '  Storm lon: 'sLon
  say '  Time step: 'timeStep
  say '  Time string: 'timeStr
  say '  GRADS file: 'gcFile
  say '  Plot dir: 'pDir
  
  'reinit'
  'clear'
  'open 'gcFile
  'set gxout contour'

* set for rval calculation
*  x,y,z covering the entire data space
  'set t 'timeStep
  'set z 1 39'
  'set x 1 207'
  'set y 1 201'
  
* radial speed relative to storm center given by sLat, sLon
  'define rval = cos(atan2(v,u)-atan2((lat-'sLat'),(lon-'sLon')))*mag(u,v)'
  
* reset to 2D slice around vortex on surface for drawing plot
  'set z 1'
  
  'd rval'
  'draw title Radial Wind Speed: z 1, ts 'timeStr
* get the x,y (page coord) from the given lat,lon
  'query w2xy 'sLon' 'sLat
  sX = subwrd(result,3)
  sY = subwrd(result,6)
  'draw mark 1 'sX' 'sY' 0.5'
* dump out gif file with white background
  pFile = pDir'/rspd_'sCase'_z1_t'timeStep'.gif'
  'printim 'pFile' white'
  
  'clear'
  'set z 2'
  'd rval'
  'query w2xy 'sLon' 'sLat
  sX = subwrd(result,3)
  sY = subwrd(result,6)
  'draw mark 1 'sX' 'sY' 0.5'
  'draw title Radial Wind Speed:z 2, ts 'timeStr
  pFile = pDir'/rspd_'sCase'_z2_t'timeStep'.gif'
  'printim 'pFile' white'
  
return
