************************************************************************
* plot out CCN fields from test sensitivity run (R2000.TEST)
*
* naming scheme for files and variables
* files containing results of scheme with grid3 northernmost row being set
*  grads4TC_CCNB_R2000.TEST-AS-1998-08-22-120000-g[123].ctl
*    vars: ccnconcen, ifnconc, dn0, gccnconcen  
*
* files containgin results of scheme without grid3 points being set
*  GRID3_OFF/grads4TC_CCNB_R2000.TEST-AS-1998-08-22-120000-g[123].ctl
*    vars: ccnconcen, ifnconc, dn0, gccnconcen  
*
* print white background so the colored contours will stand out better
* setting background to white doesn't seem to work with gif format,
* but does work with jpeg

* vars
title_ccn_g1_e = "CCN [#/cc], grid1, 2500m, 38 hrs"
title_ccn_g1_l = "CCN [#/cc], grid1, 2500m, 47.5 hrs"
title_ccn_g2_e = "CCN [#/cc], grid2, 2500m, 38 hrs"
title_ccn_g2_l = "CCN [#/cc], grid2, 2500m, 47.5 hrs"
title_ccn_g3_e = "CCN [#/cc], grid3, 2500m, 38 hrs"
title_ccn_g3_l = "CCN [#/cc], grid3, 2500m, 47.5 hrs"

early_t = 5
later_t = 24

g1_x = 79
g1_y = 78
g2_x = 104
g2_y = 101
g3_x = 207
g3_y = 201

'reinit'
* don't change the order of these without changing the 'set dfile #' commands below
'open grads4TC_CCNB_R2000.TEST-AS-1998-08-22-120000-g1.ctl '
'open grads4TC_CCNB_R2000.TEST-AS-1998-08-22-120000-g2.ctl '
'open grads4TC_CCNB_R2000.TEST-AS-1998-08-22-120000-g3.ctl '
'open GRID3_OFF/grads4TC_CCNB_R2000.TEST-AS-1998-08-22-120000-g1.ctl '
'open GRID3_OFF/grads4TC_CCNB_R2000.TEST-AS-1998-08-22-120000-g2.ctl '
'open GRID3_OFF/grads4TC_CCNB_R2000.TEST-AS-1998-08-22-120000-g3.ctl '

* Make plots of ccnconcen at 2500m

****************** setting northermost row of grid3 *******************************
'set dfile 1'
'set grads off'
'set x 1 'g1_x
'set y 1 'g1_y
'set lev 2500'
'set t 'early_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g1_e
'printim ccn_2500.t5.g3on.g1.jpg white'

'set dfile 1'
'set grads off'
'set x 1 'g1_x
'set y 1 'g1_y
'set lev 2500'
'set t 'later_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g1_l
'printim ccn_2500.t24.g3on.g1.jpg white'

'set dfile 2'
'set grads off'
'set x 1 'g2_x
'set y 1 'g2_y
'set lev 2500'
'set t 'early_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g2_e
'printim ccn_2500.t5.g3on.g2.jpg white'

'set dfile 2'
'set grads off'
'set x 1 'g2_x
'set y 1 'g2_y
'set lev 2500'
'set t 'later_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g2_l
'printim ccn_2500.t24.g3on.g2.jpg white'

'set dfile 3'
'set grads off'
'set x 1 'g3_x
'set y 1 'g3_y
'set lev 2500'
'set t 'early_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g3_e
'printim ccn_2500.t5.g3on.g3.jpg white'

'set dfile 3'
'set grads off'
'set x 1 'g3_x
'set y 1 'g3_y
'set lev 2500'
'set t 'later_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g3_l
'printim ccn_2500.t24.g3on.g3.jpg white'

****************** not setting grid3 *******************************
'set dfile 4'
'set grads off'
'set x 1 'g1_x
'set y 1 'g1_y
'set lev 2500'
'set t 'early_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g1_e
'printim ccn_2500.t5.g3off.g1.jpg white'

'set dfile 4'
'set grads off'
'set x 1 'g1_x
'set y 1 'g1_y
'set lev 2500'
'set t 'later_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g1_l
'printim ccn_2500.t24.g3off.g1.jpg white'

'set dfile 5'
'set grads off'
'set x 1 'g2_x
'set y 1 'g2_y
'set lev 2500'
'set t 'early_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g2_e
'printim ccn_2500.t5.g3off.g2.jpg white'

'set dfile 5'
'set grads off'
'set x 1 'g2_x
'set y 1 'g2_y
'set lev 2500'
'set t 'later_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g2_l
'printim ccn_2500.t24.g3off.g2.jpg white'

'set dfile 6'
'set grads off'
'set x 1 'g3_x
'set y 1 'g3_y
'set lev 2500'
'set t 'early_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g3_e
'printim ccn_2500.t5.g3off.g3.jpg white'

'set dfile 6'
'set grads off'
'set x 1 'g3_x
'set y 1 'g3_y
'set lev 2500'
'set t 'later_t
'set gxout contour'
'clear'
'd ccnconcen'
'draw title 'title_ccn_g3_l
'printim ccn_2500.t24.g3off.g3.jpg white'

'quit'
