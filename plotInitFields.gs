************************************************************************
* plot out sample fields from the init runs (first 36 hours)
*
* naming scheme for files and variables
* SET 1
*  grads1TC_CCNB_INIT-AS-1998-08-22-120000-g1.ctl
*    vars: u, v, w, press
*
* SET 3
*  grads3TC_CCNB_INIT-AS-1998-08-22-120000-g1.ctl
*    vars: ccnconcen
*

* vars
last_x = 79
last_y = 78
last_z = 39
last_t = 33

cent_lat = 15.0
cent_lon = -40.0

title_sp = "Surface pressure [mb], 36 hrs"
title_w = "w [m/s], 1000m, 36 hrs"
title_ccn = "CCN concentration [#/cc], 2500m, 36 hrs"

'reinit'
'open grads1TC_CCNB_INIT-AS-1998-08-22-120000-g1.ctl '
'open grads3TC_CCNB_INIT-AS-1998-08-22-120000-g1.ctl '

* Make surface plots of pressure, horizontal wind
* Make plot of w at 1000m
* Make plot of ccnconcen at 2500m

* surface pressure plot at last time point
'set dfile 1'
'set grads off'
'set lon -41.0 -39.0'
'set lat 14.0 16.0'
'set lev 0'
'set t 'last_t
'set gxout contour'
'clear'
'd press'
'draw title 'title_sp
'printim press_sfc.gif'

* w at 1000m
'set dfile 1'
'set grads off'
'set lon -41.0 -39.0'
'set lat 14.0 16.0'
'set lev 1000'
'set t 'last_t
'set gxout contour'
'clear'
'd w'
'draw title 'title_w
'printim w_1000.gif'

* ccnconcen at 2500m
'set dfile 2'
'set grads off'
'set lon -41.0 -39.0'
'set lat 14.0 16.0'
'set lev 2500'
'set t 'last_t
'set gxout shaded'
'clear'
'd ccnconcen'
'draw title 'title_ccn
'printim ccn_2500.gif'

'quit'
