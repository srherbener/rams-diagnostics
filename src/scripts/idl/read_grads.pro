; IDL procedure to read in GRADS data
;
; GRADS files come in pairs: the control file and the data file
;   control file - description of the format in the data file
;   data file - binary (unformatted) data
;
; Inputs:
;   ctl_file - GRADS control file (description of data in data_file)
;   data_file - GRADS data file
;
; Outputs:
;   num[xyztv] - number of x, y, z, t (time), v (variables) points
;   [xyztv]vals - 1D arrays with values of x,y,z,t,v data (along dimension)
;   gdata - 5D data array (x,y,z,t,v)
;      holds value of variable v at time t at spatial location (x,y,z)
;      x is fastest changing index, v the slowest
;

; put this in so we can index arrays beyond 32,768 entries
Compile_Opt DEFINT32

pro read_grads, ctl_file, data_file, numx, numy, numz, numt, numv, $
  xvals, yvals, zvals, tvals, vvals, gdata

  ; first read in the data description, then use the sizes
  ; to allocate the array for the data
  ; GRADS always uses 4D data
  ;    (x,y,z,t) --> (lon, lat, height, time)

  ; string
  line = ''

  print, 'Reading grads control file: ', ctl_file
  print, ''

  get_lun, lun           ; grab a free lun
  openr, lun, ctl_file

  inVarsSection = 0
  iv = 0

  while not eof(lun) do begin
    readf, lun, line
    fields = strsplit(line, ' ', /extract)

    if (strmatch(fields[0], "XDEF") or strmatch(fields[0], "xdef")) then begin
      numx = fields[1]
      xvals = fltarr(numx)
      genArrVals, numx, fields, xvals
    endif else $
    if (strmatch(fields[0], "YDEF") or strmatch(fields[0], "ydef")) then begin
      numy = fields[1]
      yvals = fltarr(numy)
      genArrVals, numy, fields, yvals
    endif else $
    if (strmatch(fields[0], "ZDEF") or strmatch(fields[0], "zdef")) then begin
      numz = fields[1]
      zvals = fltarr(numz)
      genArrVals, numz, fields, zvals
    endif else $
    if (strmatch(fields[0], "TDEF") or strmatch(fields[0], "tdef")) then begin
      numt = fields[1]
      tvals = intarr(numt)
      for i = 0, numt-1 do begin
        tvals[i] = i
      endfor
    endif else if (strmatch(fields[0], "VARS") or strmatch(fields[0], "vars")) then begin
      numv = fields[1]
      vvals = strarr(numv)
      inVarsSection = 1     
    endif else if (strmatch(fields[0], "ENDVARS") or strmatch(fields[0], "endvars")) then begin
      inVarsSection = 0
    endif else if (inVarsSection eq 1) then begin
      vvals[iv] = fields[0]
      iv = iv + 1
    endif
  endwhile
  free_lun, lun

  ; allocate the data array and read in the values from the data array

  gdata = fltarr(numx,numy,numz,numt,numv)

  print, 'Reading grads data file: ', data_file
  print, '  numx: ', numx
  print, '  numy: ', numy
  print, '  numz: ', numz
  print, '  numt: ', numt
  print, '  numv: ', numv
  print, ''

  get_lun, lun
  openr, lun, data_file
  readu, lun, gdata
  free_lun, lun

end

;*****************************************************************
; genArrVals()
;
; This procedure will generate the values of an array given
; the number of elements, and the GRADS description (fields) of
; that array.
;
; Args
;  1. num - number of values in the array
;  2. fields - the GRADS control file line split into fields
;              describing the array
;  3. array - output array
;
; Output
;  fill up the output array with values according to the instructions
;  on the GRADS control line
;

pro genArrVals, num, fields, array

  ; The GRADS control line has two forms:
  ;  XDEF <numpts> LINEAR <start> <incr>
  ;  XDEF <numpts> LEVELS <val1> <val2> ... <valn>

  if (strmatch(fields[2], "LINEAR") or strmatch(fields[2], "linear")) then begin
    start = fields[3]
    incr = fields[4]
    for i = 0, num-1 do begin
      array[i] = start + (float(i) * incr)
    endfor
  endif else if (strmatch(fields[2], "LEVELS") or strmatch(fields[2], "levels")) then begin
    ; just copy the values out of the fields array
    for i = 0, num-1 do begin
      array[i] = fields[i+3]
    endfor
  endif else begin
    print, "ERROR: unknown GRADS control: ", fields[2]
  endelse
end
