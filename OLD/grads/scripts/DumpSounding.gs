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
*    5 --> Cloud file
*    6 --> Ice file
*    7 --> Rain file
*    8 --> u (zonal winds) file
*    9 --> v (meridional winds) file
*   10 --> Longitude
*   11 --> Latitude
*   12 --> Zmin
*   13 --> Zmax
*   14 --> Timestep

function main(args)

  Pfile    = subwrd(args, 1)
  Tfile    = subwrd(args, 2)
  TdFile   = subwrd(args, 3)
  RhFile   = subwrd(args, 4)
  Cfile    = subwrd(args, 5)
  Ifile    = subwrd(args, 6)
  Rfile    = subwrd(args, 7)
  Ufile    = subwrd(args, 8)
  Vfile    = subwrd(args, 9)
  Xlon     = subwrd(args, 10)
  Ylat     = subwrd(args, 11)
  Zmin     = subwrd(args, 12)
  Zmax     = subwrd(args, 13)
  Tstep    = subwrd(args, 14)

  'reinit'
  'clear'

  say 'Extracting Sounding:'
  say '  Pressure file: 'Pfile
  say '  Temperature file: 'Tfile
  say '  Dew point file: 'TdFile
  say '  Relative humidity file: 'RhFile
  say '  Cloud file: 'Cfile
  say '  Ice file: 'Ifile
  say '  Rain file: 'Rfile
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
  'sdfopen 'Cfile
  'sdfopen 'Ifile
  'sdfopen 'Rfile
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
*
* Temp and dew point come of of REVU in deg C, convert to K --> add 273.15
* Cloud, Rain and Ice come out of REVU in g/kg, convert to kg/kg --> divide by 1000
*
* When displaying variables, shut off warnings and display the variable twice. It
* seems that GRADS insists upon writing its warning about chunk sizes at least once
* even when warnings are shut off.
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
      'set warn off'
      'd 'press
      'd 'press
      P = subwrd(result, 4)

      'set dfile 2'
      'set warn off'
      'd 'tempc
      'd 'tempc
      T = subwrd(result, 4) + 273.15

      'set dfile 3'
      'set warn off'
      'd 'dewptc
      'd 'dewptc
      Td = subwrd(result, 4) + 273.15

      'set dfile 4'
      'set warn off'
      'd 'relhum
      'd 'relhum
      RH = subwrd(result, 4)

      'set dfile 5'
      'set warn off'
      'd 'cloud
      'd 'cloud
      Cloud = subwrd(result, 4) / 1000

      'set dfile 6'
      'set warn off'
      'd 'ice
      'd 'ice
      Ice = subwrd(result, 4) / 1000

      'set dfile 7'
      'set warn off'
      'd 'rain
      'd 'rain
      Rain = subwrd(result, 4) / 1000

      'set dfile 8'
      'set warn off'
      'd 'u
      'd 'u
      Uwind = subwrd(result, 4)

      'set dfile 9'
      'set warn off'
      'd 'v
      'd 'v
      Vwind = subwrd(result, 4)

      say 'SND: 'Zlev' 'P' 'T' 'Td' 'RH' 'Cloud' 'Ice' 'Rain' 'Uwind' 'Vwind
    endif

    i = i + 1
  endwhile

return
