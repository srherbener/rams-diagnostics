; IDL procedure to generate plots
;
; This routine depends on compiling the read_grads and plot_utils routines before running.
;

;**********************************************************
; plot_max
;
; plot the time series of maximum values
;
; Args
;  1. plot_file - file to place output plot
;

pro plot_max, gExp, gVar, gUnits, cases, startTime, timeInc, gctl_file, gdat_file, plot_file, f_yrange=force_yrange, plot_min=used_plot_min, levs=set_levels

  ; put this in so we can index arrays beyond 32,768 entries
  Compile_Opt DEFINT32

  ncases = 1

  plot_min_vals = 0
  if (n_elements(used_plot_min) eq 1) then begin
    plot_min_vals = used_plot_min
  endif

  use_yrange = 0
  if (n_elements(force_yrange) eq 2) then begin
    use_yrange = 1
  endif

  user_levels = 0
  user_start_lev = 0
  user_end_lev = 0
  if (n_elements(set_levels) eq 2) then begin
    user_levels = 1
    user_start_lev = set_levels[0]
    user_end_lev = set_levels[1]
  endif

  c_colors = [ 1, 1, 1, 1, 1 ]
  c_lstyles = [ 0, 2, 3, 4, 5 ]

  ; read in the GRADS data
  
  read_grads, gctl_file, gdat_file, numx, numy, numz, numt, numv, $
    xvals, yvals, zvals, tvals, vvals, gdata

  var_data = fltarr(numx,numy,numz,numt,numv)
  var_extreme = fltarr(1,numt)

  var_data[*,*,*,*] = gdata

help, var_data
help, gdata

  if (user_levels eq 1) then begin
    k_start = user_start_lev
    k_end = user_end_lev
  endif else begin
    k_start = 0
    k_end = numz - 1
  endelse

print, 'k_start = ', k_start
print, 'k_end = ', k_end

  for it = 1, (numt -1) do begin
    if (plot_min_vals eq 0) then begin
      var_extreme[0,it] = 0.0
    endif else begin
      var_extreme[0,it] = 100000000.0
    endelse
  endfor

  for it = 0, (numt -1) do begin
    for iz = k_start, k_end do begin
      for iy = 0, (numy -1) do begin
        for ix = 0, (numx -1) do begin
          if (plot_min_vals eq 0) then begin
            if (var_data[ix,iy,iz,it] gt var_extreme[0,it]) then begin
               var_extreme[0,it] = var_data[ix,iy,iz,it]
            endif
          endif else begin
            if (var_data[ix,iy,iz,it] lt var_extreme[0,it]) then begin
               var_extreme[0,it] = var_data[ix,iy,iz,it]
            endif
          endelse
        endfor
      endfor
    endfor
  endfor

  ; plot the curves
  xlabel = 'Time (hr)'
  if (plot_min_vals eq 0) then begin
    plabel = string(gExp, ': Max value of: ', gVar)
    ylabel = string('Max ', gVar, ' (', gUnits, ')')
  endif else begin
    plabel = string(gExp, ': Min value of: ', gVar)
    ylabel = string('Min ', gVar, ' (', gUnits, ')')
  endelse

  ; convert the time step numbers to time in hours
  timevals = fltarr(numt)
  for it = 0, (numt - 1) do begin
    timevals[it] = startTime + (float(tvals[it]-1) * timeInc)
  endfor

  if (use_yrange eq 0) then begin
    plot_multi_curves, ncases, numt, timevals, var_extreme, $
    cases, plabel, xlabel, ylabel, c_colors, c_lstyles, plot_file
  endif else begin
    plot_multi_curves, ncases, numt, timevals, var_extreme, $
    cases, plabel, xlabel, ylabel, c_colors, c_lstyles, plot_file, f_yrange=force_yrange
  endelse

end
