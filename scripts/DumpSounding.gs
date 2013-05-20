************************************************************************
* write a sounding into an ascii file
*   each line is a pressure
*   the entries on a line contain:
*       P T Td RH WindSpeed WindDirection
*
* Input arguments:
*    1 --> Pressure file
*    2 --> Temperature file
*    3 --> Dew point file
*    4 --> Relative Humidity file
*    5 --> u (zonal winds) file
*    6 --> v (meridional winds) file
*    7 --> Longitude
*    8 --> Latitude
*    9 --> Zmin
*   10 --> Zmax
*   11 --> Timestep

function main(args)

  Pfile    = subwrd(args, 1)
  Tfile    = subwrd(args, 2)
  TdFile   = subwrd(args, 3)
  RhFile   = subwrd(args, 4)
  Ufile    = subwrd(args, 5)
  Vfile    = subwrd(args, 6)
  Xlon     = subwrd(args, 7)
  Ylat     = subwrd(args, 8)
  Zmin     = subwrd(args, 9)
  Zmax     = subwrd(args, 10)
  Tstep    = subwrd(args, 11)

  'reinit'
  'clear'

  say 'Extracting Sounding:'
  say '  Pressure file: 'Pfile
  say '  Temperature file: 'Tfile
  say '  Dew point file: 'TdFile
  say '  Relative humidity file: 'RhFile
  say '  U file: 'Ufile
  say '  V file: 'Vfile
  say ''
  say '  Lon: 'Xlon
  say '  Lat: 'Ylat
  say '  Zmin: 'Zmin
  say '  Zmax: 'Zmax
  say '  Time step: 'Tstep
  say ''
  
* Pressure    -> dfile 1
* Temperature -> dfile 2
* Dew Point   -> dfile 3
* RH          -> dfile 4
* U           -> dfile 5
* V           -> dfile 6

  'sdfopen 'Pfile
  'sdfopen 'Tfile
  'sdfopen 'TdFile
  'sdfopen 'RhFile
  'sdfopen 'Ufile
  'sdfopen 'Vfile

* Sounding will be taken from one lon, lat location at one time step.
* Only z varies. Use q dims to find the indices for Zmin and Zmax. Then
* use these indices to loop through the levels and dump out the sounding
* data.
  'set dfile 1'
  'set lon 'Xlon
  'set lat 'Ylat
  'set t 'Tstep

  'set lev 'Zmin
  'q dims'
  zspec = sublin(result, 4)
  Zstart = subwrd(zspec, 9)
  'set lev 'Zmax
  'q dims'
  zspec = sublin(result, 4)
  Zend = subwrd(zspec, 9)
  
  say '  Index for Zmin: 'Zstart
  say '  Index for Zmax: 'Zend
  say ''

* Walk through each level and dump out the sounding infomration
*   Set z
*   Grab variables from their respective files
*   Dump out a line with each value for that level
*      Mark the line so the caller can easily identify it
  i = Zstart
  while (i <= Zend)

    'set dfile 1'
    'set z 'i
    'q dims'
    zspec = sublin(result, 4)
    Zlev = subwrd(zspec, 6)

    if (Zlev >= 0)
      say '  Level = 'Zlev

      'set dfile 1'
      'd press'
      P = subwrd(result, 4)

      'set dfile 2'
      'd 'tempc
      T = subwrd(result, 4)

      'set dfile 3'
      'd 'dewptc
      Td = subwrd(result, 4)

      'set dfile 4'
      'd 'relhum
      RH = subwrd(result, 4)

      'set dfile 5'
      'd 'u
      Uwind = subwrd(result, 4)

      'set dfile 6'
      'd 'v
      Vwind = subwrd(result, 4)

      say 'SND: 'Zlev' 'P' 'T' 'Td' 'RH' 'Uwind' 'Vwind
    endif

    i = i + 1
  endwhile

return
