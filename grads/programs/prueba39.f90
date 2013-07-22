!###############################################################################
!   
!###############################################################################
implicit none
!
!
! ==============================================================================
! =============================> parameters  <===++=============================
! ==============================================================================
!
integer, parameter :: nexp = 1    !!!6
real,  parameter :: delta = 10.  !<--- used by cinta
!!
!---> begining time 
! 1=24/00:00Z, 36h
! 3=24/00:03Z, 39h
! 6=24/00:00Z, 42h
! 9=24/00:00Z, 45h
!
integer, parameter :: i_start= 1
integer,  parameter :: i_end = 37
!
! ==============================================================================
! ==============================================================================
integer,  parameter :: nt = 37
integer,  parameter :: nx = 207 
integer,  parameter :: ny = 201  
integer,  parameter :: nz = 39 
real,  parameter :: dx = 1500.,  dy = 1500.
!==============================
!
!-------------------------> model levels
real,dimension(nz) :: hdata      
data hdata/147.639,457.235,786.955,1138.108,1512.085,1910.370,2334.544,2786.290,3267.399,3779.781,4325.467,   &
4906.623,5525.554,6184.715,6886.722,7634.360,8430.595,9278.585,10181.694,11146.951,12147.144,13147.6,14147.6, &
15147.6,16147.6,17147.6,18147.6,19147.6,20147.6,21147.6,22147.6,23147.6,24147.6,25147.6,26147.6, 27147.6, &
28147.6, 29147.6, 30147.609/
! ==============================================================================
integer :: i,j,k,it,ivar, i_exp, i_aer, j_sta,i_aux, k_ser, iit, ii, jj, kk, i_cond
!
!---> inputs
real , dimension(nx,ny,nz) :: var
real , dimension(nx,ny,nt) :: acum, rate
real , dimension(nx,ny,2) :: var2
real , dimension(nx,ny,nz,3) :: var3a, var3b
real , dimension(nx,ny,nz,7) :: var7a, var7b
!
real , dimension(nx,ny,nz) ::   w, w_f
real , dimension(nx,ny,nz) ::   tc, qv, densi
real , dimension(nx,ny,nz) ::   qc,qr,qi,qs,qa,qg,qh,qsc,qia 
real , dimension(nx,ny,nz) ::   ccn,nc,dc,nsc,dsc,dr,ni,di,da
!
!---> outputs
real, dimension(nt,nexp,3) ::v_scmax,v_iamax,v_rmax,v_hmax,wmax,wmax_h,maxrat,averat, sclwp,umodmax
real, dimension(nt,nexp,3) ::qiamax,qrmax,qhmax,qscmax,nscmax,dscmax,qimax,qamax,dimax,damax,areaw
!
!---> auxiliares
!
real , dimension(nx,ny) :: dist
real :: ncount1,ncount2,ncount3,aux1,aux2,aux3, aux4
real , dimension(nx,ny) :: cinta, afuera,todo, condition
!
character*70 :: file_in5,file_in1,file_in2,file_in3,file_in4,file_out1,file_out 
real :: ccn_low, fimax, fvmax, dist0, umod
!############
!

do it=1,nt
do i_exp= 1, nexp
do i_cond=1, 3
  v_scmax(it,i_exp,i_cond) =0.
  v_iamax(it,i_exp,i_cond) =0.
  v_rmax(it,i_exp,i_cond) =0.
  v_hmax(it,i_exp,i_cond) =0.
  !
  wmax(it, i_exp,i_cond)= 0. 
  wmax_h(it, i_exp,i_cond)= 0. 
  areaw(it,i_exp,i_cond) = 0.
  !
  qiamax(it,i_exp,i_cond)= 0.
  qrmax(it,i_exp,i_cond)=  0.
  qhmax(it,i_exp,i_cond)=  0.
  qscmax(it,i_exp,i_cond) =   0.
  nscmax(it,i_exp,i_cond) =   0.
  dscmax(it,i_exp,i_cond) =   0.
  qiamax(it,i_exp,i_cond) = 0.
  qimax(it,i_exp,i_cond)  = 0.
  qamax(it,i_exp,i_cond)  = 0.
  dimax(it,i_exp,i_cond)  = 0.
  damax(it,i_exp,i_cond)  = 0.
  maxrat(it,i_exp,i_cond)=0.
  averat(it,i_exp,i_cond)=0.
  !
  qscmax(it,i_exp,i_cond) =  0.
  nscmax(it,i_exp,i_cond) =  0.
  dscmax(it,i_exp,i_cond) =  0.
  qiamax(it,i_exp,i_cond) =  0.
  qimax(it,i_exp,i_cond)  =  0.
  qamax(it,i_exp,i_cond)  =  0.
  dimax(it,i_exp,i_cond)  =  0.
  damax(it,i_exp,i_cond)  =  0.
  !
  sclwp(it,i_exp,i_cond)  =  0.
enddo
enddo
enddo
!
!
!######################################################################################################
                    do i_exp= 1 , nexp  !!!! <------- experiment
!######################################################################################################
!
! ---> circle
!
if (i_exp==1) dist0= 60.
if (i_exp==2) dist0= 60.
if (i_exp==3) dist0= 60.
!
!****************************************************************************
!
!
if (i_exp==1)  file_in1 = '/media/disk/DHS/temp/00.gra'
if (i_exp==2)  file_in1 = '/media/disk/DHS/temp/39b.4000.gra'
if (i_exp==3)  file_in1 = '/media/disk/DHS/temp/39b.6000.gra'
if (i_exp==4)  file_in1 = '/media/disk/DHS/temp/39b.8000.gra'
if (i_exp==5)  file_in1 = '/media/disk/DHS/temp/39b.d6k.gra'
if (i_exp==6)  file_in1 = '/media/disk/DHS/temp/39b.d8k.gra'
!
!------------
!
!
if (i_exp==1)  file_in2 = '/media/disk/DHS/micro/00_q.gra'
if (i_exp==2)  file_in2 = '/media/disk/DHS/micro/39b_q.4000.gra'
if (i_exp==3)  file_in2 = '/media/disk/DHS/micro/39b_q.6000.gra'
if (i_exp==4)  file_in2 = '/media/disk/DHS/micro/39b_q.8000.gra'
if (i_exp==5)  file_in2 = '/media/disk/DHS/micro/39b_q.d6k.gra'
if (i_exp==6)  file_in2 = '/media/disk/DHS/micro/39b_q.d8k.gra'
!
!
!------------
!
!
if (i_exp==1)  file_in3 = '/media/disk/DHS/uvw/00.gra'
if (i_exp==2)  file_in3 = '/media/disk/DHS/uvw/39b.4000.gra'
if (i_exp==3)  file_in3 = '/media/disk/DHS/uvw/39b.6000.gra'
if (i_exp==4)  file_in3 = '/media/disk/DHS/uvw/39b.8000.gra'
if (i_exp==5)  file_in3 = '/media/disk/DHS/uvw/39b.d6k.gra'
if (i_exp==6)  file_in3 = '/media/disk/DHS/uvw/39b.d8k.gra'
!
!
!------------
!
!
if (i_exp==1)  file_in4 = '/media/disk/DHS/precip/00.gra'
if (i_exp==2)  file_in4 = '/media/disk/DHS/precip/39b.4000.gra'
if (i_exp==3)  file_in4 = '/media/disk/DHS/precip/39b.6000.gra'
if (i_exp==4)  file_in4 = '/media/disk/DHS/precip/39b.8000.gra'
if (i_exp==5)  file_in4 = '/media/disk/DHS/precip/39b.d6k.gra'
if (i_exp==6)  file_in4 = '/media/disk/DHS/precip/39b.d8k.gra'
!'
!
!------------
!
if (i_exp==1)  file_in5 = '/media/disk/DHS/micro/00_c.gra'
if (i_exp==2)  file_in5 = '/media/disk/DHS/micro/39b_c.4000.gra'
if (i_exp==3)  file_in5 = '/media/disk/DHS/micro/39b_c.6000.gra'
if (i_exp==4)  file_in5 = '/media/disk/DHS/micro/39b_c.8000.gra'
if (i_exp==5)  file_in5 = '/media/disk/DHS/micro/39b_c.d6k.gra'
if (i_exp==6)  file_in5 = '/media/disk/DHS/micro/39b_c.d8k.gra'
!
!****************************************************************************
!
print*, file_in1
print*, file_in2
print*, file_in3
print*, file_in4
print*, file_in5
!
!*********************************************************************************************
!1. ----> define conditions
!*********************************************************************************************
!
do i=1, nx 
do j=1, ny
  !
  dist(i,j)= (real(i)-real(nx)/2.)**2 + (real(j)-real(ny)/2.)**2 
  !
  dist(i,j)= sqrt(dist(i,j))
  !
  if ((dist(i,j)> (dist0-delta)).and.(dist(i,j)<(dist0+delta))) then
         cinta(i,j)= 1.
  else
         cinta(i,j)= 0.
  endif
  !
  if (dist(i,j)> (dist0-delta)) then
        afuera(i,j) = 1. 
  else
        afuera(i,j) = 0. 
  endif

  todo(i,j)=1
  !
enddo
enddo
!
!
! ############################################################################################
! ############################################################################################
                                       do it = i_start, i_end
! ############################################################################################
! ############################################################################################
!
!*********************************************************************************************
!3. ----> read temperature, vapor, and density.
!*********************************************************************************************
!
! ----------------------------------------------------
if (it==i_start) open (unit=1,file=file_in1 ,status='old',form='unformatted',access='direct',recl=3*nz*nx*ny*4)
!! ----------------------------------------------------
! ----------------------------------------------------
!
read(1, rec=it) ((((var3a(i,j,k,ivar),i=1,nx),j=1,ny),k=1,nz), ivar=1,3)
do i=1, nx 
do j=1, ny
do k=1, nz
  tc(i,j,k) = var3a(i,j,k,1)
  qv(i,j,k) = var3a(i,j,k,2)
  densi(i,j,k) =var3a(i,j,k,3)
enddo
enddo
enddo
!
! ----------------------------------------------------
if (it==i_end)  close (1)
! ----------------------------------------------------
!
!*********************************************************************************************
!3. ----> mixing ratios
!*********************************************************************************************
!
! ----------------------------------------------------
if (it==i_start) open (unit=2,file=file_in2 ,status='old',form='unformatted',access='direct',recl=7*nz*nx*ny*4)
!! ----------------------------------------------------
! ----------------------------------------------------
!

read(2, rec=it) ((((var7a(i,j,k,ivar),i=1,nx),j=1,ny),k=1,nz), ivar=1,7)
!
do i=1, nx 
do j=1, ny
do k=1, nz
   !
   qc(i,j,k)= var7a(i,j,k,1)
   qr(i,j,k)= var7a(i,j,k,2)
   qi(i,j,k)= var7a(i,j,k,3) 
   qs(i,j,k)= var7a(i,j,k,4) 
   qa(i,j,k)= var7a(i,j,k,5) 
   qg(i,j,k)= var7a(i,j,k,6) 
   qh(i,j,k)= var7a(i,j,k,7) 
   qia(i,j,k)=  qi(i,j,k)+qa(i,j,k)
   !
   if (tc(j,j,k)<-0.5) then
      qsc(i,j,k) =  qc(i,j,k)
   else
      qsc(i,j,k) = 0.
   endif
   !
enddo
enddo
enddo
!
do i=1, nx 
do j=1, ny
      if (cinta(i,j)==1.) then
          do k=1, nz-1
              v_scmax(it,i_exp,1) = v_scmax(it,i_exp,1) + 0.001 * qsc(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
              v_iamax(it,i_exp,1) = v_iamax(it,i_exp,1) + 0.001 * qia(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
              v_rmax(it,i_exp,1) = v_rmax(it,i_exp,1) + 0.001 * qr(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
              v_hmax(it,i_exp,1) = v_hmax(it,i_exp,1) + 0.001 * qh(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
          enddo
      endif
      !
      if (afuera(i,j)==1.) then
          do k=1, nz-1
              v_scmax(it,i_exp,2) = v_scmax(it,i_exp,2) + 0.001 * qsc(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
              v_iamax(it,i_exp,2) = v_iamax(it,i_exp,2) + 0.001 * qia(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
              v_rmax(it,i_exp,2) = v_rmax(it,i_exp,2) + 0.001 * qr(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
              v_hmax(it,i_exp,2) = v_hmax(it,i_exp,2) + 0.001 * qh(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
          enddo
      endif
      !

      if (todo(i,j)==1.) then
          !
          aux4=0.
          !
          do k=1, nz-1
              !
              aux4 = aux4  + 0.001 * qsc(i,j,k)*densi(i,j,k)*(hdata(k+1)-hdata(k))
              !
              v_scmax(it,i_exp,3) = v_scmax(it,i_exp,3) + 0.001 * qsc(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
              v_iamax(it,i_exp,3) = v_iamax(it,i_exp,3) + 0.001 * qia(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
              v_rmax(it,i_exp,3) = v_rmax(it,i_exp,3) + 0.001 * qr(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
              v_hmax(it,i_exp,3) = v_hmax(it,i_exp,3) + 0.001 * qh(i,j,k)*densi(i,j,k)*dx*dy* (hdata(k+1)-hdata(k))
          enddo
          sclwp(it,i_exp,3)= max( sclwp(it,i_exp,3),aux4)  
      endif
      !
enddo
enddo 
!
! ----------------------------------------------------
if (it==i_end)  close (2)
! ----------------------------------------------------
!
!*********************************************************************************************
!5. ----> read updrafts
!*********************************************************************************************
!
! ----------------------------------------------------
if (it==i_start) open (unit=3,file=file_in3 ,status='old',form='unformatted',access='direct',recl=3*nz*nx*ny*4)
!! ----------------------------------------------------
! ----------------------------------------------------
!
read(3, rec=it) ((((var3b(i,j,k,ivar),i=1,nx),j=1,ny),k=1,nz), ivar=1,3)
!
do i=1, nx
do j=1, ny
    aux1=0.
    aux2=0.
    aux3=0.
    !
    do k=1, nz
      !
      umod = sqrt(var3b(i,j,k,1)**2 + var3b(i,j,k,2)**2)
      !
      if (umod>umodmax(it,i_exp,1)) then
           umodmax(it,i_exp,1) = umod
           umodmax(it,i_exp,2) = umod
           umodmax(it,i_exp,3) = umod
      endif
      !
      if ((cinta(i,j)==1.).and.(var3b(i,j,k,3)>wmax(it,i_exp,1))) then
           wmax(it,i_exp,1) = var3b(i,j,k,3)
           wmax_h(it,i_exp,1) = hdata(k) 
      endif
      !
      if ((afuera(i,j)==1.).and.(var3b(i,j,k,3)>wmax(it,i_exp,2))) then
           wmax(it,i_exp,2) = var3b(i,j,k,3)
           wmax_h(it,i_exp,2) = hdata(k) 
      endif
      !
      if ((todo(i,j)==1.).and.(var3b(i,j,k,3)>wmax(it,i_exp,3))) then
           wmax(it,i_exp,3) = var3b(i,j,k,3)
           wmax_h(it,i_exp,3) = hdata(k) 
      endif
      !
    enddo
    !
    do k=1, nz
       if ((cinta(i,j)==1.).and.(var3b(i,j,k,3)>5.))  aux1=1.
       if ((afuera(i,j)==1.).and.(var3b(i,j,k,3)>5.)) aux1=2.
       if ((todo(i,j)==1.).and.(var3b(i,j,k,3)>5.))   aux1=3.
    enddo
    !
    if (aux1==1.) areaw(it,i_exp,1) = areaw(it,i_exp,1) + dx/1000. * dy/1000.
    if (aux2==1.) areaw(it,i_exp,2) = areaw(it,i_exp,2) + dx/1000. * dy/1000.
    if (aux3==1.) areaw(it,i_exp,3) = areaw(it,i_exp,3) + dx/1000. * dy/1000.
    !
enddo
enddo
! 
! ----------------------------------------------------
if (it==i_end)  close (3)
! ----------------------------------------------------
!
!
!*********************************************************************************************
!6. ----> read ppp
!*********************************************************************************************
!
! ----------------------------------------------------
if (it==i_start) open (unit=4,file=file_in4,status='old',form='UNFORMATTED',access='direct',recl=2*ny*nx*4)
!! ----------------------------------------------------
!
ncount1= 0.
ncount2= 0.
ncount3= 0.
!
read(4, rec=it) (((var2(i,j,ivar),i=1,nx),j=1,ny), ivar=1,2)
!
do i=1, nx
do j=1, ny
    acum(i,j,it)= var2(i,j,1)
    rate(i,j,it)= var2(i,j,2)
    !
    if (cinta(i,j) ==1.) then
       maxrat(it,i_exp,1)=max(maxrat(it,i_exp,1),rate(i,j,it))
       ncount1=ncount1+1.
       averat(it,i_exp,1)=averat(it,i_exp,1)+rate(i,j,it)
    endif
    !
    if (afuera(i,j) ==1.) then
       maxrat(it,i_exp,2)=max(maxrat(it,i_exp,2),rate(i,j,it))
       ncount2=ncount2+1.
       averat(it,i_exp,2)=averat(it,i_exp,2)+rate(i,j,it)
    endif
    !
    if (todo(i,j) ==1.) then
       maxrat(it,i_exp,3)=max(maxrat(it,i_exp,3),rate(i,j,it))
       ncount3=ncount3+1.
       averat(it,i_exp,3)=averat(it,i_exp,3)+rate(i,j,it)
    endif
    !
enddo
enddo
!
averat(it,i_exp,1)=averat(it,i_exp,1)/ncount1
averat(it,i_exp,2)=averat(it,i_exp,2)/ncount2
averat(it,i_exp,3)=averat(it,i_exp,3)/ncount3
! 
! ----------------------------------------------------
if (it==i_end)  close (4)
! ----------------------------------------------------
!
!*********************************************************************************************
!7. ----> concen's
!*********************************************************************************************
!!
! ----------------------------------------------------
if (it==i_start) open (unit=5,file=file_in5 ,status='old',form='unformatted',access='direct',recl=7*nz*nx*ny*4)
!! ----------------------------------------------------
! ----------------------------------------------------
!
!
read(5, rec=it) ((((var7b(i,j,k,ivar),i=1,nx),j=1,ny),k=1,nz), ivar=1,7)
!
do i=1, nx 
do j=1, ny
do k=1, nz
 !
 ccn(i,j,k)= var7b(i,j,k,1)
 nc(i,j,k)= var7b(i,j,k,2)
 dc(i,j,k)= var7b(i,j,k,3) 
 dr(i,j,k)= var7b(i,j,k,4)
 ni(i,j,k)= var7b(i,j,k,5) 
 di(i,j,k)= var7b(i,j,k,6) 
 da(i,j,k)= var7b(i,j,k,7) 
 !
 if (cinta(i,j)==1.) then
 !
      if (qsc(i,j,k)>qscmax(it,i_exp,1)) then
         qscmax(it,i_exp,1) =  qsc(i,j,k)
         nscmax(it,i_exp,1) =  nc(i,j,k) 
         dscmax(it,i_exp,1) =  dc(i,j,k) 
      endif
      !
      if ((qia(i,j,k)>qiamax(it,i_exp,1)).and.(k>=27) ) then
         qiamax(it,i_exp,1) = qia(i,j,k) 
         qimax(it,i_exp,1)  =  qi(i,j,k)
         qamax(it,i_exp,1)  =  qa(i,j,k)
         dimax(it,i_exp,1)  =  di(i,j,k)
         damax(it,i_exp,1)  =  da(i,j,k)
      endif 
      !
      if (qr(i,j,k)>qrmax(it,i_exp,1)) then
         qrmax(it,i_exp,1) = qr(i,j,k) 
      endif 
      if (qh(i,j,k)>qhmax(it,i_exp,1)) then
         qhmax(it,i_exp,1) = qh(i,j,k) 
      endif 
    !
  endif !(cinta)
  !
 if (afuera(i,j)==1.) then
 !
      if (qsc(i,j,k)>qscmax(it,i_exp,2)) then
         qscmax(it,i_exp,2) =  qsc(i,j,k)
         nscmax(it,i_exp,2) =  nc(i,j,k) 
         dscmax(it,i_exp,2) =  dc(i,j,k) 
      endif
      !
      if ((qia(i,j,k)>qiamax(it,i_exp,2))) then
         qiamax(it,i_exp,2) = qia(i,j,k) 
         qimax(it,i_exp,2)  =  qi(i,j,k)
         qamax(it,i_exp,2)  =  qa(i,j,k)
         dimax(it,i_exp,2)  =  di(i,j,k)
         damax(it,i_exp,2)  =  da(i,j,k)
      endif 
      !
      if (qr(i,j,k)>qrmax(it,i_exp,2)) then
         qrmax(it,i_exp,2) = qr(i,j,k) 
      endif 
      if (qh(i,j,k)>qhmax(it,i_exp,2)) then
         qhmax(it,i_exp,2) = qh(i,j,k) 
      endif 
    !
  endif !(afuera)
  !
 if (todo(i,j)==1.) then
 !
      if (qsc(i,j,k)>qscmax(it,i_exp,3)) then
         qscmax(it,i_exp,3) =  qsc(i,j,k)
         nscmax(it,i_exp,3) =  nc(i,j,k) 
         dscmax(it,i_exp,3) =  dc(i,j,k) 
      endif
      !
      if ((qia(i,j,k)>qiamax(it,i_exp,3)).and.(k>=27) ) then
         qiamax(it,i_exp,3) = qia(i,j,k) 
         qimax(it,i_exp,3)  =  qi(i,j,k)
         qamax(it,i_exp,3)  =  qa(i,j,k)
         dimax(it,i_exp,3)  =  di(i,j,k)
         damax(it,i_exp,3)  =  da(i,j,k)
      endif 
      !
      if (qr(i,j,k)>qrmax(it,i_exp,3)) then
         qrmax(it,i_exp,3) = qr(i,j,k) 
      endif 
      if (qh(i,j,k)>qhmax(it,i_exp,3)) then
         qhmax(it,i_exp,3) = qh(i,j,k) 
      endif 
    !
  endif !(todo)
  !
enddo
enddo
enddo
! 
! ----------------------------------------------------
if (it==i_end)  close (5)
! ----------------------------------------------------
!
! ############################################################################################
! ############################################################################################
                                           enddo
! ############################################################################################
! ############################################################################################

!##########################################
enddo !!!! <------ experiment
!##########################################
!
!
!
do i_exp=1, nexp
!
if (i_exp==1)  file_out = '00_39.gra'
if (i_exp==2)  file_out = '39b_4.gra'
if (i_exp==3)  file_out = '39b_6.gra'
if (i_exp==4)  file_out = '39b_8.gra'
if (i_exp==5)  file_out = '39b_12.gra'
if (i_exp==6)  file_out = '39b_16.gra'
!
 open(unit=i_exp,file=file_out,status='unknown',form='UNFORMATTED', access='sequential',position='append')
  do it= i_start, i_end
     !
     i_cond = 3 
       !
       print*, it, umodmax(it,i_exp,i_cond)

       write(i_exp) wmax(it,i_exp,i_cond) 
       write(i_exp) wmax_h(it,i_exp,i_cond) 
       write(i_exp) qrmax(it,i_exp,i_cond) 
       write(i_exp) v_scmax(it,i_exp,i_cond) 
       write(i_exp) qscmax(it,i_exp,i_cond) 
       write(i_exp) sclwp(it,i_exp,3) 
       write(i_exp) averat(it,i_exp,i_cond)
       write(i_exp) umodmax(it,i_exp,i_cond)
       !
       !
     !
  enddo ! time
  !  
enddo !<--- output file
!
close (i_exp)





end 
! ############################################################################################
! ############################################################################################
! ############################################################################################





