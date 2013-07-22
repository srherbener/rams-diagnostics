; IDL procedures to plot GRADS data
;
; Routines:
;   plot_radial - this routine will plot a given variable as a function of
;                 radius for a given z level and time point
;
;   plot_tseries - this routine will plot a given variable as a function of
;                  time for a given z level and radius point
;
;   plot_profile - this routine will plot a vertical profile of a given variable
;                  for a given radius and time point
;
; These routines all depend on compiling the read_grads routine before running.
;

; put this in so we can index arrays beyond 32,768 entries
Compile_Opt DEFINT32

;**********************************************************
; plot_radial
;
; Plot a variable as a function of radius, at a given
; z level and time
;
; Args
;  1. varname - name of variable to grab from azavg data
;  2. zlev - z level, integer index for the height value to use
;            (the hieghts are defined in the GRADS control file)
;  3. time - time point, integer index for the time point
;            (again see GRADS control file)
;  4. plot_file - name of output PS file for resulting plot
;

pro plot_radial, varname, zlev, time, plot_file

  ; Assume that the GRADS files are always called <case>_<varname>.[ctl|gra]
;  cases = [ 'CLNM', 'R0100_R0500', 'R0100_R1000', 'R0100_R1500', 'R0100_R2000' ]
  cases = [ 'CCNB_R0100', 'CCNB_R0500', 'CCNB_R1000', 'CCNB_R1500', 'CCNB_R2000' ]
  ncases = 5
  npoints = 50
  ;c_colors = [ 48, 96, 160, 224, 240 ]
  c_colors = [ 1, 1, 1, 1, 1 ]
  c_lstyles = [ 0, 2, 3, 4, 5 ]
  azdata = fltarr(ncases,npoints)

  ; read in all cases - assume that dimensions of GRADS data will be the same
  ; for each case so keep overwriting the num*, *vals vars
  for i = 0, (ncases-1) do begin
    gctl_file = string( cases[i], '_', varname, '.ctl')
    gdat_file = string( cases[i], '_', varname, '.gra')

    read_grads, gctl_file, gdat_file, numx, numy, numz, numt, numv, $
      xvals, yvals, zvals, tvals, vvals, az_data

    azdata[i,*] = az_data[*, 0, zlev, time]
  endfor

  ; plot the curves
  case varname of
    'press'    : begin 
                   title_name = 'Pressure'
                   ylabel = 'Pressure (mb)'
                 end
    'speed_t'  : begin 
                   title_name = 'Tangential Wind Speed'
                   ylabel = 'Wind Speed (m/s)'
                 end
    else       : begin
                   title_name = varname
                   ylabel = varname
                 end
  endcase

  plabel = string('Azavg ', title_name, ' - Level: ', zvals[zlev], 'm, Timestep: ', tvals[time])
  xlabel = 'Radius From Storm Center (km)'

  plot_multi_curves, ncases, npoints, xvals, azdata, $
  cases, plabel, xlabel, ylabel, c_colors, c_lstyles, plot_file
end

;**********************************************************
; plot_tseries
;
; plot a time series of a varaiable given a z level
; and radius
;
; Args
;  1. varname - name of variable to grab from azavg data
;  1. zlev - index of height level
;  2. radius - index of radius from storm center
;  3. plot_file - output post script plot file
;

pro plot_tseries, varname, zlev, radius, plot_file

  ; Assume that the GRADS files are always called <case>_speed_t.[ctl|gra]
;  cases = [ 'CLNM', 'R0100_R0500', 'R0100_R1000', 'R0100_R1500', 'R0100_R2000' ]
  cases = [ 'CCNB_R0100', 'CCNB_R0500', 'CCNB_R1000', 'CCNB_R1500', 'CCNB_R2000' ]
  ncases = 5
  npoints = 73
  ;c_colors = [ 48, 96, 160, 224, 240 ]
  c_colors = [ 1, 1, 1, 1, 1 ]
  c_lstyles = [ 0, 2, 3, 4, 5 ]
  azdata = fltarr(ncases,npoints)

  ; read in all cases - assume that dimensions of GRADS data will be the same
  ; for each case so keep overwriting the num*, *vals vars
  for i = 0, (ncases-1) do begin
    gctl_file = string( cases[i], '_', varname, '.ctl')
    gdat_file = string( cases[i], '_', varname, '.gra')

    read_grads, gctl_file, gdat_file, numx, numy, numz, numt, numv, $
      xvals, yvals, zvals, tvals, vvals, az_data

    azdata[i,*] = az_data[radius, 0, zlev, 0:72]
  endfor

  ; plot the curves
  case varname of
    'press'    : begin 
                   title_name = 'Pressure'
                   ylabel = 'Pressure (mb)'
                 end
    'speed_t'  : begin 
                   title_name = 'Tangential Wind Speed'
                   ylabel = 'Wind Speed (m/s)'
                 end
    else       : begin
                   title_name = varname
                   ylabel = varname
                 end
  endcase

  plabel = string('Azavg ', title_name, ' - Level: ', zvals[zlev], 'm, Radius: ', xvals[radius], 'km')
  xlabel = 'Timestep number'

  plot_multi_curves, ncases, npoints, tvals, azdata, $
  cases, plabel, xlabel, ylabel, c_colors, c_lstyles, plot_file
  
end

;**********************************************************
; plot_profile
;
; Plot the vertical ice profile
;
; Args
;  1. varname - name of variable from azavg data
;  2. radius - distance from storm center, integer index
;            (the radius values are defined in the GRADS control file)
;  3. time - time point, integer index for the time point
;            (again see GRADS control file)
;  4. plot_file - name of output PS file for resulting plot
;

pro plot_profile, varname, radius, time, plot_file

;  cases = [ 'CLNM', 'R0100_R0500', 'R0100_R1000', 'R0100_R1500', 'R0100_R2000' ]
  cases = [ 'CCNB_R0100', 'CCNB_R0500', 'CCNB_R1000', 'CCNB_R1500', 'CCNB_R2000' ]
  ncases = 5
  npoints = 39
  ;c_colors = [ 48, 96, 160, 224, 240 ]
  c_colors = [ 1, 1, 1, 1, 1 ]
  c_lstyles = [ 0, 2, 3, 4, 5 ]
  prof_data = fltarr(ncases,npoints)
  temp_data = fltarr(ncases,npoints)
  avg_temps = fltarr(npoints)

  ; read in all cases - assume that dimensions of GRADS data will be the same
  ; for each case so keep overwriting the num*, *vals vars
  for i = 0, (ncases-1) do begin
    ; Assume that the GRADS files are always called <case>_<varname>.[ctl|gra]
    ; Also, need the 0 degC level as a reference: <case>_tempc.[ctl|gra]
    gctl_file = string( cases[i], '_', varname, '.ctl')
    gdat_file = string( cases[i], '_', varname, '.gra')

    tempc_ctl_file = string( cases[i], '_tempc.ctl')
    tempc_dat_file = string( cases[i], '_tempc.gra')

    ; the ice data and tempc data should have the exact same grid
    ; so just overwrite num* and *vals arrays
    read_grads, gctl_file, gdat_file, numx, numy, numz, numt, numv, $
      xvals, yvals, zvals, tvals, vvals, az_data
    read_grads, tempc_ctl_file, tempc_dat_file, numx, numy, numz, numt, numv, $
      xvals, yvals, zvals, tvals, vvals, az_tempc

    prof_data[i,*] = az_data[radius, 0, *, time]
    temp_data[i,*] = az_tempc[radius, 0, *, time]
  endfor

  case varname of
    'ice'        : begin 
                     title_name = 'Ice'
                     xlabel = 'Ice Mixing Ratio (g/kg)'
                   end
    'liquid'     : begin 
                     title_name = 'Liquid Water'
                     xlabel = 'Liquid Water Mixing Ratio (g/kg)'
                   end
    'ccnconcen'  : begin 
                     title_name = 'CCN'
                     xlabel = 'CCN Concentration (#/cm3)'
                   end
    'gccnconcen' : begin 
                     title_name = 'GCCN'
                     xlabel = 'GCCN Concentration (#/cm3)'
                   end
    else         : begin
                     title_name = varname
                     xlabel = varname
                   end
  endcase
  plabel = string('Azavg ', title_name, ' Profile - Radius: ', xvals[radius], 'km, Timestep: ', tvals[time])
  ylabel = 'Height (m)'

  ; create average temperatures for each z level for plotting a zero degree isotherm
  ; reference on the profile plot
  for i = 0, (npoints-1) do begin
    ; ith element is mean of the ith points for all cases
    avg_temps[i] = mean(temp_data[*,i])
  endfor

  plot_multi_curves, ncases, npoints, zvals, prof_data, $
  cases, plabel, xlabel, ylabel, c_colors, c_lstyles, plot_file, $
  vertical=1, add_zero_isotherm=avg_temps

end

;**********************************************************
; plot_multi_curves
;
; This routine will take a 2D array of data [case, curve]
; and plot the curves (one for each case) on the same plot.
;
; Args
;  1. ncases - number of cases
;  2. npoints - number of points in each curve
;  3. xvals - 1D array, [points], x values
;  4. pdata - 2D array, [cases, points], y values
;  5, cases - 1D array, [cases], names of cases
;  6. plabel - main title for plot
;  7. xlabel - x axis label
;  8. ylabel - y axis label
;  9. colors - 1D array, [cases], colors for plot data for each case
; 10. lstyles - 1D array, [cases], line styles for plot data for each case
; 11. plot_file - post script file to dump plot into
;
; Keyword options
;  vertical - if 1, then use xvals as y values on the plot,
;                   and the different curves in pdata as x values
;                   (to make vertical profile plots)
;             not 1, plot horizontally
;

pro plot_multi_curves, ncases, npoints, xvals, pdata, $
  cases, plabel, xlabel, ylabel, colors, lstyles, plot_file, $
  vertical=plot_vertically, add_zero_isotherm=temps, f_yrange=force_yrange

  ; set the keyword opts
  if (n_elements(plot_vertically) eq 0) then begin
    ; keyword 'vertical' was not used
    plot_vertically = 0
  endif

  if (plot_vertically eq 1) then begin
    ; vertical profile plot
    xrange = [ min(pdata), max(pdata) ]
    yrange = [ min(xvals), max(xvals) ]
  endif else begin
    ; horizontal plot
    xrange = [ min(xvals), max(xvals) ]
    yrange = [ min(pdata), max(pdata) ]
  endelse

  ; if keyword yrange was used properly it will contain 2 elements
  if (n_elements(force_yrange) eq 2) then begin
    yrange = force_yrange
  endif

;  set_plot, 'ps' ; do this before loadct so that x server is not required
;  loadct, 13  ; rainbow, each 16 indecies are one color with different shades
;  device, color=1, filename=plot_file, xsize=8.0, ysize=10.0, /inches, $
;      bits_per_pixel=8, /encapsulated

  ;set up the z buffer, all the drawing commands will create image in
  ; z buffer, then the z buffer image can be converted to a GIF file
   
  set_plot, 'Z' ; do this before loadct so that x server is not required
  device, set_pixel_depth=24; use colors
  device, decomposed=0; use a color table
  device, z_buffering=0 ; 2D plotting
  device, set_resolution=[1024, 768]
   
  loadct, 39  ; rainbow+white
  ;loadct, 0 ; gray scale
  
  ; clear the buffer
  erase

  plot, xrange, yrange, color=1, title=plabel, xtitle=xlabel, ytitle=ylabel, $
        yrange=yrange, position=[ 0.20, 0.30, 0.95, 0.85 ], /nodata, $
        background=255, charsize=2.0, thick=2.0

  for i = 0, (ncases-1) do begin
    if (plot_vertically eq 1) then begin
      ; vertical profile plot
      oplot, pdata[i,*], xvals, color=colors[i], linestyle=lstyles[i], thick=2.0
    endif else begin
      ; horizontal plot
      oplot, xvals, pdata[i,*], color=colors[i], linestyle=lstyles[i], thick=2.0
    endelse

    y_legend = 0.20 - (float(i) * 0.03)
    plots, [0.10, 0.15], [y_legend, y_legend], color=colors[i], linestyle=lstyles[i], $
       /normal, thick=2.0
    xyouts, 0.17, (y_legend-0.005), cases[i], color=colors[i], $
       /normal, charsize=2.0
  endfor

  if (n_elements(temps) gt 0) then begin
    ; assume used only for vertical profile plots
    ; xvals actually holds the z levels
    plot_zero_isotherm, xrange, npoints, xvals, temps, 1
  endif
 
  ;read the image into an array, and then dump the array to a GIF file
  image = color_quan(TVRD(TRUE=1), 1, r, g, b)
  write_gif, plot_file, image, r, g, b

;  device, /close_file

end

;*****************************************************************
; plot_zero_isotherm
;
; This routine will take an array (z levels) of temperatures (degC)
; and find the level that represents the zero degree isotherm. Then
; it will plot it on the current plot.
;
; Args
;  1. xrange - 2 element array containing the x min and max values
;  2. numz - number of levels in the array
;  3. zvals - array of z levels
;  4. tempc - array of temperatures (degC)
;               entries correspond to entries in zvals
;  5. pcolor - color for zero degC isotherm line
;

pro plot_zero_isotherm, xrange, numz, zvals, tempc, pcolor
  ; find the height where the temperature (degC) becomes negative, record
  ; the 0 degC level as the last level before temp becomes negative
  level_0 = 0.0
  for i = 0, (numz-1) do begin
    if (tempc[i] ge 0.0) then begin
      level_0 = zvals[i]
    endif
  endfor

  text_x = (mean(xrange) + max(xrange)) / 2.0

  plots, xrange, [level_0, level_0], linestyle=2, color=pcolor
  xyouts, text_x, level_0, 'zero deg. C level'
end

;******************************************************************
; plot_single_curve
;
; This routine will plot a single curve.
; The image is placed in the Z buffer and then written to a GIF
; file at the end.
;
; Args
;   1. xvals - array of independent var values
;   2. yvals - array of dependent var values
;   3. pfile - name of file to put the GIF image into
;
;   Notes number of points in xvays and yvals must match.
;

pro plot_single_curve, xvals, yvals, plabel, xlabel, ylabel, color, lstyle, pfile

  ;set up the z buffer, all the drawing commands will create image in
  ; z buffer, then the z buffer image can be converted to a GIF file

  set_plot, 'Z' ; do this before loadct so that x server is not required
  device, set_pixel_depth=24, decomposed=0 ; use colors, with a color table

  loadct, 39  ; rainbow+white
  ;loadct, 0 ; gray scale

  ; clear the buffer
  erase

  ; plot the frame first in black, then plot the data
  xrange = [ min(xvals), max(xvals) ]
  yrange = [ min(yvals), max(yvals) ]
  plot, xrange, yrange, color=1, title=plabel, xtitle=xlabel, ytitle=ylabel, $
        yrange=yrange, background=255, /nodata

  oplot, xvals, yvals, color=color, linestyle=lstyle

  ;read the image into an array, and then dump the array to a GIF file
  image = color_quan(TVRD(TRUE=1), 1, r, g, b)
  write_gif, pfile, image, r, g, b

end

;******************************************************************
; plot_azavg_vslices
;
; This routine will generate the vertical slice plots for the
; CCNB experiments
;

pro plot_azavg_vslices, expname
  cases = [ 'CCNB_R0100', 'CCNB_R0500', 'CCNB_R1000', 'CCNB_R1500', 'CCNB_R2000' ]
  ncases = 5

  vars = [ 'w', 'speed_t', 'speed_r', 'thetae', 'liquid', 'ice', 'gccnconcen', 'ccnconcen' ]
  nvars = 8

  gdir = 'AzAveragedData'

  case expname of 
    'GCCN_ON' : begin 
    times = [ 1, 10, 20, 30, 40, 50, 60, 70 ]
    ntimes = 8
    end
    'GCCN_OFF' : begin 
    times = [ 1, 10, 20, 30, 40, 50, 60 ]
    ntimes = 7
    end
  endcase

  for ic = 0, (ncases-1) do begin
    for iv = 0, (nvars-1) do begin
      for it = 0, (ntimes-1) do begin
        time_string = string(times[it], format='(i2.2)')
        pfile = string(gdir, '/vslice_plots/', cases[ic], '_', vars[iv], '_', time_string, '.gif')
        plot_vslice, gdir, cases[ic], vars[iv], times[it], pfile
      endfor
    endfor
  endfor
end

;******************************************************************
; plot_vslice
;
; This routine will create a contour plot of a vertical slice data
; set in grads.
;

pro plot_vslice, gdir, gcase, varname, timestep, pfile

  nclevels = 10
  nz_clip = 25 ; approx 16km (top of storm)

  ; piece together the GRADS file names
  gc_file = string(gdir, '/', gcase, '_', varname, '.ctl')
  gd_file = string(gdir, '/', gcase, '_', varname, '.gra')

  ; read in the grads data
  read_grads, gc_file, gd_file, numx, numy, numz, numt, numv, $
     xvals, yvals, zvals, tvals, vvals, gdata

  ; extract the particular vertical slice out of the data
  vs_data = fltarr(numx, nz_clip)
  vs_heights = fltarr(nz_clip)
  nx = numx - 1
  nz = nz_clip - 1
  vs_data[*,*] = gdata[0:nx,0,0:nz,timestep]
  vs_heights = zvals[0:nz] / 1000.0 ; convert to km

print, 'Data min: ',min(vs_data)
print, 'Data max: ',max(vs_data)

  ; set up the labels, etc.
  case varname of
    'w' : begin 
      var_descrip = 'vertical velocity (m/s)'
      c_start = -0.5
      c_end = 4.0
      end
    'speed_t' : begin 
      var_descrip = 'horizontal tangential velocity (m/s)'
      c_start = -10.0
      c_end = 80.0
      end
    'speed_r' : begin 
      var_descrip = 'horizontal radial velocity (m/s)'
      c_start = -20.0
      c_end = 25.0
      end
    'thetae' : begin 
      var_descrip = 'equivalent potential temperature (K)'
      c_start = 320.0
      c_end = 410.0
      end
    'liquid' : begin 
      var_descrip = 'liquid water content (g/kg)'
      c_start = 0.0
      c_end = 4.5
      end
    'ice' : begin 
      var_descrip = 'ice content (g/kg)'
      c_start = 0.0
      c_end = 1.8
      end
    'gccconcen' : begin 
      var_descrip = 'GCCN concentration (#/cc)'
      c_start = 0.0
      c_end = 1.8
      end
    'ccconcen' : begin 
      var_descrip = 'CCN concentration (#/cc)'
      c_start = 0.0
      c_end = 4500.0
      end
    else : begin
      var_descrip = varname
      c_start = min(vs_data)
      c_end = max(vs_data)
      end
  endcase

  time_string = string(timestep, format='(i2.2)')
  ptitle = string('Azimuthally averaged ', var_descrip, ', ', gcase, ', ts ', time_string)
  xtitle = 'Radius (km)'
  ytitle = 'Height (km)'
  c_colors=[16,  40,  64,  88, 112, 136, 160, 184, 208, 232]

  c_levels = fltarr(nclevels)
  c_inc = (c_end - c_start)/float(nclevels-1)
  c_lev = c_start
  for i = 0, (nclevels-1) do begin
    c_levels[i] = c_lev
    c_lev = c_lev + c_inc
  endfor

print, 'ptitle: ', ptitle
  ; create the contour plot
  plot_contour, vs_data, xvals, vs_heights, ptitle, xtitle, ytitle, $
    c_colors, c_levels, pfile
  
end

;******************************************************************
; plot_contour
;
; This routine will create a contour plot
;
; Args
;   1. 2D array with data to plot
;   2. output GIF file
;

pro plot_contour, pdata, xvals, yvals, ptitle, xtitle, ytitle, $
    c_colors, c_levels, pfile

  ;set up the z buffer, all the drawing commands will create image in
  ; z buffer, then the z buffer image can be converted to a GIF file

  set_plot, 'Z' ; do this before loadct so that x server is not required
  device, set_pixel_depth=24; use colors
  device, decomposed=0; use a color table
  device, z_buffering=0 ; 2D plotting
  device, set_resolution=[1024, 768]

  loadct, 39  ; rainbow+white
  ;loadct, 0 ; gray scale
  
  ; clear the buffer
  erase

  contour, pdata, xvals, yvals, color=1, background=255, /follow, $
      title=ptitle, xtitle=xtitle, ytitle=ytitle, $
      c_colors=c_colors, levels=c_levels

  ;read the image into an array, and then dump the array to a GIF file
  image = color_quan(TVRD(TRUE=1), 1, r, g, b)
  write_gif, pfile, image, r, g, b

end
