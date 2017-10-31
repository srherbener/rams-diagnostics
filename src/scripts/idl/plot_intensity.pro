; IDL procedure to generate plots
;
; This routine depends on compiling the read_grads and plot_utils routines before running.
;

;**********************************************************
; plot_intensity
;
; plot the time series of area coverage of hurricane force winds
;
; Args
;  1. plot_file - file to place output plot
;

pro plot_intensity, gExp, cases, startTime, timeInc, plot_file, f_yrange=force_yrange

  ; put this in so we can index arrays beyond 32,768 entries
  Compile_Opt DEFINT32

  ; Assume that the GRADS files are always called sfc_wind_<case>.[ctl|gra]
;  ncases = 5
;  cases = [ 'R0100', 'R0500', 'R1000', 'R1500', 'R2000' ]
;  cases = [ 'CLNM', 'R0100_R0500', 'R0100_R1000', 'R0100_R1500', 'R0100_R2000' ]
  ;c_colors = [ 48, 96, 160, 224, 240 ]

  use_yrange = 0
  if (n_elements(force_yrange) eq 2) then begin
    use_yrange = 1
  endif
  ncases = n_elements(cases)
  c_colors = [ 1, 1, 1, 1, 1 ]
  c_lstyles = [ 0, 2, 3, 4, 5 ]

  ; read in all cases - assume that dimensions of GRADS data will be the same
  ; for each case so keep overwriting the num*, *vals vars
  ; read the first file to get the dimensions, allocate the arrays and then
  ; read in the rest of the data
  
  gctl_file = string( 'GRADS/sfc_wind_', cases[0], '.ctl')
  gdat_file = string( 'GRADS/sfc_wind_', cases[0], '.gra')

  read_grads, gctl_file, gdat_file, numx, numy, numz, numt, numv, $
    xvals, yvals, zvals, tvals, vvals, gdata

  sfc_wind_data = fltarr(ncases,numx,numy,numz,numt,numv)
  hfw_coverage = fltarr(ncases,numt) ;hfw -> hurricane force wind

  sfc_wind_data[0,*,*,*,*,*] = gdata

  for ic = 1, (ncases-1) do begin
    gctl_file = string( 'GRADS/sfc_wind_', cases[ic], '.ctl')
    gdat_file = string( 'GRADS/sfc_wind_', cases[ic], '.gra')

    read_grads, gctl_file, gdat_file, numx, numy, numz, numt, numv, $
      xvals, yvals, zvals, tvals, vvals, gdata

    sfc_wind_data[ic,*,*,*,*,*] = gdata
  endfor

help, sfc_wind_data
help, gdata

  ; convert the surface wind data into weighted coverage data
  ; run through each time point and each x,y location of the sfc wind data
  ; weight each count by the Saffir-Simpson category number
  ;
  ;   Wind speed   Weight
  ;    <  33 m/s -> 0  (not hurricane force)
  ;    33-42 m/s -> 1  (Category 1 hurricane)
  ;    43-49 m/s -> 2  (etc.)
  ;    50-58 m/s -> 3
  ;    59-69 m/s -> 4
  ;    >= 70 m/s -> 5
  ;

  iv = 0 ; should be just one variable
  iz = 0 ; should be just one z level

  num_hpts = 0.0
  num_hpts = float(numy) * float(numx)
help, num_hpts

  for ic = 0, (ncases-1) do begin
    for it = 0, (numt -1) do begin
      hfw_coverage[ic, it] = 0.0
      ncat0 = 0
      ncat1 = 0
      ncat2 = 0
      ncat3 = 0
      ncat4 = 0
      ncat5 = 0
      for ix = 0, (numx -1) do begin
        for iy = 0, (numy -1) do begin
          swind = sfc_wind_data[ic, ix, iy, iz, it, iv]
          if (swind ge 70.0) then begin
            ncat5 = ncat5 + 1
          endif else if (swind ge 59.0) then begin
            ncat4 = ncat4 + 1
          endif else if (swind ge 50.0) then begin
            ncat3 = ncat3 + 1
          endif else if (swind ge 43.0) then begin
            ncat2 = ncat2 + 1
          endif else if (swind ge 33.0) then begin
            ncat1 = ncat1 + 1
          endif else begin
            ncat0 = ncat0 + 1
          endelse
        endfor
      endfor
; linear weighting
      hfw_coverage[ic, it] = (ncat1 * 1.0) + (ncat2 * 2.0) + (ncat3 * 3.0) + $
         (ncat4 * 4.0) + (ncat5 * 5.0)
; power weighting (square)
;      hfw_coverage[ic, it] = (ncat1 * 1.0) + (ncat2 * 4.0) + (ncat3 * 9.0) + $
;         (ncat4 * 16.0) + (ncat5 * 25.0)
; power weighting (cube)
;      hfw_coverage[ic, it] = (ncat1 * 1.0) + (ncat2 * 8.0) + (ncat3 * 27.0) + $
;         (ncat4 * 64.0) + (ncat5 * 125.0)

      ; weight by the number of horizontal grid cell points
      hfw_coverage[ic, it] = hfw_coverage[ic, it] / num_hpts

      print, 'Case, timestep: ', cases[ic], (it+1)
      print, '  counts(0,1,2,3,4,5,tot): ', ncat0, ncat1, ncat2, ncat3, ncat4, ncat5, $
          (ncat0 + ncat1 + ncat2 + ncat3 + ncat4 + ncat5)
      print, ''
    endfor
  endfor

  ; plot the curves
  plabel = string(gExp, ': Sfc Wind Coverage Metric')
  xlabel = 'Time (hr)'
  ylabel = 'Coverage metric'

  ; convert the time step numbers to time in hours
  timevals = fltarr(numt)
  for it = 0, (numt - 1) do begin
    timevals[it] = startTime + (float(tvals[it]-1) * timeInc)
  endfor

  if (use_yrange eq 0) then begin
    plot_multi_curves, ncases, npoints, timevals, hfw_coverage, $
    cases, plabel, xlabel, ylabel, c_colors, c_lstyles, plot_file
  endif else begin
    plot_multi_curves, ncases, npoints, timevals, hfw_coverage, $
    cases, plabel, xlabel, ylabel, c_colors, c_lstyles, plot_file, f_yrange=force_yrange
  endelse

end
