; IDL procedure to generate plots
;
; This routine depend on compiling the read_grads and plot_utils routines before running.
;

; put this in so we can index arrays beyond 32,768 entries
Compile_Opt DEFINT32

;**********************************************************
; do_plots
;
; Args
;  1. plot_dir - directory to place output plots
;

pro do_plots, plot_dir

  times = [ 0, 10, 20, 30, 40, 50, 60, 70 ]
  radii = [ 5, 15, 25, 35, 45 ]

  for it = 0, 7 do begin
    ; plot surface wind and pressure
    plot_file = string( format='(a,"/SfcWind.",i2.2,".ps")', plot_dir, times[it] )
    print, 'Plot wind --> ', plot_file
    plot_radial, 'speed_t', 0, times[it], plot_file

    plot_file = string( format='(a,"/SfcPress.",i2.2,".ps")', plot_dir, times[it] )
    print, 'Plot pressure --> ', plot_file
    plot_radial, 'press', 0, times[it], plot_file

    ; plot ice, liquid water, temp, and ccn profiles at various radii
    print, 'Plot profiles:'
    for ir = 0, 4 do begin
      plot_file = string( format='(a,"/ProfIce.", i2.2,".",i2.2,".ps")', $
         plot_dir, radii[ir], times[it] )
      print, '  ', ir, ' --> ', plot_file
      plot_profile, 'ice', radii[ir], times[it], plot_file

      plot_file = string( format='(a,"/ProfLiq.", i2.2,".",i2.2,".ps")', $
         plot_dir, radii[ir], times[it] )
      print, '  ', ir, ' --> ', plot_file
      plot_profile, 'liquid', radii[ir] , times[it], plot_file

      plot_file = string( format='(a,"/ProfTemp.", i2.2,".",i2.2,".ps")', $
         plot_dir, radii[ir], times[it] )
      print, '  ', ir, ' --> ', plot_file
      plot_profile, 'tempc', radii[ir] , times[it], plot_file

      plot_file = string( format='(a,"/ProfCcn.", i2.2,".",i2.2,".ps")', $
         plot_dir, radii[ir], times[it] )
      print, '  ', ir, ' --> ', plot_file
      plot_profile, 'ccnconcen', radii[ir] , times[it], plot_file

      plot_file = string( format='(a,"/ProfW.", i2.2,".",i2.2,".ps")', $
         plot_dir, radii[ir], times[it] )
      print, '  ', ir, ' --> ', plot_file
      plot_profile, 'w', radii[ir] , times[it], plot_file
    endfor
  endfor

  ; plot time series of surface wind and pressure
  print, 'Plot time series:'
  for ir = 0, 4 do begin
    plot_file = string( format='(a,"/TsSfcWind.", i2.2,".ps")', $
       plot_dir, radii[ir] )
    print, '  ', ir, ' --> ', plot_file
    plot_tseries, 'speed_t', 0, radii[ir], plot_file

    plot_file = string( format='(a,"/TsSfcPress.", i2.2,".ps")', $
       plot_dir, radii[ir] )
    print, '  ', ir, ' --> ', plot_file
    plot_tseries, 'press', 0, radii[ir], plot_file
  endfor

end
