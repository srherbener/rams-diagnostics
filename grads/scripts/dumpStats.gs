************************************************************************
* dump out stats for a given GRADS variable in the "set gxout stat" format
*
* varNum arg:
*   1 --> GRADS control file
*   2 --> GRADS var name

function main(args)

  gFile = subwrd(args, 1)
  gVar  = subwrd(args, 2)

  say 'Dumping stats for GRADS variable:'
  say '  GRADS control file = 'gFile
  say '  GRADS variable = 'gVar
  say ''

  'reinit'
  'clear'
  say 'Reading GRADS control file: 'gFile
  'open 'gFile
* get the dimensions
  'q file'
  gDims = sublin(result,5)
  xMax  = subwrd(gDims,3)
  yMax  = subwrd(gDims,6)
  zMax  = subwrd(gDims,9)
  tMax  = subwrd(gDims,12)

  say '  xMax = 'xMax
  say '  yMax = 'yMax
  say '  zMax = 'zMax
  say '  tMax = 'tMax
  say ''

* loop through z and t values, keep x and y fixed to cover entire horizontal
* plane at each z and t
  'set x 1 'xMax
  'set y 1 'yMax
  it = 1
  while (it <= tMax)
    iz = 1
    while (iz <= zMax)
* print out the min, max, mean, std deviation values
      'set z 'iz
      'set t 'it
      'set gxout stat'
      'd 'gVar

      minMaxLine = sublin(result,8)
      meanLine   = sublin(result,12)
      stdDevLine = sublin(result,14)
      gMin       = subwrd(minMaxLine,4)
      gMax       = subwrd(minMaxLine,5)
      gMean      = subwrd(meanLine,2)
      gSdev      = subwrd(stdDevLine,2)

      say 'Stats:'
      say '  t: 'it
      say '  z: 'iz
      say '  Min: 'gMin
      say '  Max: 'gMax
      say '  Mean: 'gMean
      say '  StdDev: 'gSdev

      iz = iz + 1
    endwhile
    it = it + 1
  endwhile

return
