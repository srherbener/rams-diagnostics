; IDL procedure to generate plots
;
; This routine depends on compiling the read_grads and plot_utils routines before running.
;

;**********************************************************
; plot_h_integral
;
; plot the time series of horizontal integral
;
; Args
;  1. plot_file - file to place output plot
;

pro plot_h_integral, gExp, gVar, gUnits, cases, startTime, timeInc, delta_x, delta_y, gctl_file, gdat_file, plot_file, f_yrange=force_yrange, z_level=user_z_level

  ; put this in so we can index arrays beyond 32,768 entries
  Compile_Opt DEFINT32

  ncases = 1

  use_yrange = 0
  if (n_elements(force_yrange) eq 2) then begin
    use_yrange = 1
  endif

  c_colors = [ 1, 1, 1, 1, 1 ]
  c_lstyles = [ 0, 2, 3, 4, 5 ]

  ; read in the GRADS data
  
  read_grads, gctl_file, gdat_file, numx, numy, numz, numt, numv, $
    xvals, yvals, zvals, tvals, vvals, gdata

  var_data = fltarr(numx,numy,numz,numt,numv)
  var_h_int = fltarr(1,numt)

  var_data[*,*,*,*] = gdata

help, var_data
help, gdata

  iz = 0 ; level nearest surface
  if (n_elements(user_z_level) eq 1) then begin
    iz = user_z_level
  endif

print, 'Z level = ', iz

  for it = 0, (numt -1) do begin
    ; sum up over the horizontal slice
    var_h_int[0,it] = 0.0
    for iy = 0, (numy -1) do begin
      for ix = 0, (numx -1) do begin
          var_h_int[0,it] = var_h_int[0,it] + var_data[ix,iy,iz,it]
      endfor
    endfor

    ; multiply by area of slice, this can be done once outside the
    ; x,y loops only if the size of all the grids are the same (which is
    ; the case here)
    ; 
    ;   doing precipitation rate which comes in mm/hr
    ;   multiplying by area of slice (m**2) will yield kg of water per hr
    ;     -> density of water ~1000 kg/m**3
    ;     -> mass of water = area * x mm/hr * density of water
    ;           = m**2 * 0.001m * 1000 kg/m**3 /hr = kg/hr
    ;           = area of grid cell * precip_rate = mass of water per hr

    var_h_int[0,it] = delta_x * delta_y * var_h_int[0,it]
  endfor

  ; plot the curves
  plabel = string(gExp, ': Horiz Integral of: ', gVar, ', Z level = ', iz)
  xlabel = 'Time (hr)'
  ylabel = string('Horiz Integral ', gVar, ' (', gUnits, ')')

  ; convert the time step numbers to time in hours
  timevals = fltarr(numt)
  for it = 0, (numt - 1) do begin
    timevals[it] = startTime + (float(tvals[it]-1) * timeInc)
  endfor

  if (use_yrange eq 0) then begin
    plot_multi_curves, ncases, numt, timevals, var_h_int, $
    cases, plabel, xlabel, ylabel, c_colors, c_lstyles, plot_file
  endif else begin
    plot_multi_curves, ncases, numt, timevals, var_h_int, $
    cases, plabel, xlabel, ylabel, c_colors, c_lstyles, plot_file, f_yrange=force_yrange
  endelse

end
