implicit none
!********************************************************************
integer,  parameter :: n_exp = 9
!
!
real, parameter :: delta_temp= 4.5
!
!********************************************************************
integer,  parameter :: nt = 37
integer,  parameter :: nx = 207 
integer,  parameter :: ny = 201  
integer,  parameter :: nz = 39 
real,  parameter :: dx = 1500.,  dy = 1500.
!
integer :: i, j, k, i_exp, it, ivar
real , dimension(nx,ny,nz,3) :: var3
!
character*70 :: file_in5,file_in1,file_in2,file_in3,file_in4,file_out1,file_out !
!
!====> PROPIAS
real , dimension(nx,ny,nt) :: var
real :: a_med, a_max, g_med, diff_max, aux1
real, dimension(nt,n_exp) :: tcmed
real, dimension(n_exp) :: aux
real , dimension(nx,ny) ::  mask_big
integer :: t_big
!
!********************************************************************
print*
print*
print*
print*
print*
print*
!*************
do i_exp=1, n_exp
!*************
!
if (i_exp==1)  file_in1 = '/media/disk/DHS/temp/00.gra'
if (i_exp==2)  file_in1 = '/media/disk/DHS/temp/36a.8000.gra'
if (i_exp==3)  file_in1 = '/media/disk/DHS/temp/36b.8000.gra'
if (i_exp==4)  file_in1 = '/media/disk/DHS/temp/42b.8000.gra'
if (i_exp==5)  file_in1 = '/media/disk/DHS/temp/39b.8000.gra'
if (i_exp==6)  file_in1 = '/media/disk/DHS/temp/36ab.8000.gra'
if (i_exp==7)  file_in1 = '/media/disk/DHS/temp/42ab.8000.gra'
if (i_exp==8) file_in1 = '/media/disk/DHS/temp/39ab.8000.gra'
if (i_exp==9)  file_in1 = '/media/disk/DHS/temp/39b.6000.gra'


!
!
open (unit=1,file=file_in1 ,status='old',form='unformatted',access='direct',recl=3*nz*nx*ny*4)
diff_max =0.
do it =1 , nt
       open (unit=1,file=file_in1 ,status='old',form='unformatted',access='direct',recl=3*nz*nx*ny*4)
       read(1, rec=it) ((((var3(i,j,k,ivar),i=1,nx),j=1,ny),k=1,nz), ivar=1,3)
       !
       ! tcmed = ?
       !
       tcmed(it,i_exp) = 0.
       do i=1, nx 
       do j=1, ny
           tcmed(it,i_exp) = tcmed(it,i_exp) + var3(i,j,1,1)
       enddo
       enddo
       tcmed(it,i_exp) = tcmed(it,i_exp)/real(nx*ny) 
       !
       ! anomaly = ?
       !
       do i=1, nx 
       do j=1, ny
           !
           if ((var3(i,j,1,1) - tcmed(it,i_exp))< -delta_temp) then  !<--- relative to each experiment's ave 
!!!!!!!!!!!if ((var3(i,j,1,1) - tcmed(it,1))< -delta_temp) then      !<--- relative to control's ave 

                var(i,j,it) = 1.
                diff_max=min (diff_max,(var3(i,j,1,1) - tcmed(it,i_exp))    )
           else
                var(i,j,it) = 0.
           endif
           !
       enddo
       enddo
!
enddo ! time
       !
       !>>> PDF_a(nx, ny, ntimes, delta_xy, var, ntinit, ntend, icellmin, a_med, a_max, mask_big,t_big) <<<    
       !
       if (i_exp ==1) then
       print*,'******************************************************************' 
       print*, '              results for DT=', delta_temp
       endif
       !
       print*,'******************************************************************' 
       print*, file_in1
       !
       call PDF_a(nx, ny,   nt,    1.5,      var,   1,     nt,      4,    a_med, a_max, mask_big,t_big)  
       !
       aux1=0
       g_med = 0.
       do i=1, nx 
       do j=1, ny
            if (mask_big(i,j)==1.) then
                  g_med = g_med + var3(i,j,1,1)
                  aux1 = aux1 + 1.
            endif
       enddo
       enddo
       g_med = g_med/aux1
       !
       print*,i_exp, a_med, a_max, t_big
!*************
enddo !exp
!*************
!
print*,'******************************************************************' 
Do it=1, nt
   !print*, it, tcmed(it,2)-tcmed(it,1), tcmed(it,3)-tcmed(it,1), tcmed(it,4)-tcmed(it,1)
enddo
!
!
do i=1, n_exp
    aux(i)=0.
    do it=13, 37
       aux(i) = aux(i) + tcmed(it,i)
    enddo
    aux(i)=aux(i)/24.
    print*, i,  aux(i)-aux(1)
enddo
!!
end
!########################################################################################
!########################################################################################
subroutine PDF_a(nx, ny, ntimes, delta_xy, var, ntinit, ntend, icellmin, a_med, a_max, mask_big,t_big)  
implicit none
integer :: nx, ny, ntimes,ntinit, ntend, icellmin,t_big
real , dimension(nx,ny,ntimes) :: var
real :: delta_xy
!
real, dimension(nx,ny) :: group, mask_big           ! cell state and mask_big for biggest group
integer  ncount                                 ! count of cells of any group
integer  ncell                                  ! number of center cells checked 
integer  ncount_max
!
real, dimension(nx*ny,ntinit:ntend) :: ngroup ! Group counters 
!
real, dimension(nx*ny) :: bin, dist, pdf
integer    ilist (nx*ny)                         ! list of x's within one cloud
integer    jlist (nx*ny)                         ! list of y's within one cloud
!
!--> aux
integer  ::  aux2
!
integer ::    ntot,is,js,l,m,ivar,it, dummy, todo,i, j, k, ncat 
real :: ntot_groups, a_med, s_dist, a_max
!
!     Parameters for diplay 
integer    coarser
integer,   parameter :: deltacoarser=1
real,parameter :: lowpercent=0.01
real :: percentage, accum
!
!********************************************************************
!
ncat = nx*ny
!
!---->Initialize group counters for nx*ny categories
do k=1,ncat
   !for number of category
   dist(k)=0.
   do m=ntinit, ntend 
       !for time
       ngroup(k,m)=0.
   enddo
enddo
!
!---->Initialize counter for biggest group
ncount_max=0



!********************************************************************
!********************************************************************
                      do it= ntinit, ntend
!********************************************************************
!********************************************************************
!
!---->define  group for time it 
do i=1,nx
do j=1,ny
   !
   group(i,j) = var(i,j, it)
   !
enddo
enddo
!
!--------------------> FINDS GROUPS ---------------------------------
!
!----> Finds a cell to start a group 
!
1000 dummy=dummy  
todo = 1
do i=1, nx
do j=1, ny
  if (group(i,j).eq.1) then
       !----> initialize lists for the group and counters
       do k=1, ncat
                ilist(k)=0
                jlist(k)=0
       enddo
       ncell=0
       ncount=1
       !-----> First cell of the group   
       ilist(1)=i
       jlist(1)=j
       !----->Takes out first cell of this group 
       group(i,j)=2
       todo = 0
       goto 1001
   endif
enddo
enddo
!
if (todo == 1) goto 3333
!
!----> Begins new group -
! Initialize counter for new group
1001 dummy=dummy   
ncell=ncell + 1
is=ilist(ncell)
js=jlist(ncell)
!
!-----> Check neighbors count, change state and cells to the list
! 
!----->East
if (is<nx) then
if (group(is+1,js).eq.1) then
             ncount=ncount+1
             group(is+1,js)=2
             ilist(ncount)=is+1
             jlist(ncount)=js
endif
endif
!----->West
if (is>1) then
if (group(is-1,js).eq.1) then
             ncount=ncount+1
             group(is-1,js)=2
             ilist(ncount)=is-1
             jlist(ncount)=js
endif
endif
!----->North
if (js<ny) then
if (group(is,js+1).eq.1) then
             ncount=ncount+1
             group(is,js+1)=2
             ilist(ncount)=is
             jlist(ncount)=js+1
endif
endif
!----->South
if (js>1) then
if (group(is,js-1).eq.1) then
             ncount=ncount+1
             group(is,js-1)=2
             ilist(ncount)=is
             jlist(ncount)=js-1
endif
endif
! -------------------------------------------------------------------
! -------------------------------------------------------------------
if (ncount>ncell) then
!------->Go to check neighbors of the next cell in the list
  goto 1001
else
!-----> Counts one group of "ncount" cells in this cloud
! and go finds new group
  ngroup(ncount,it)=ngroup(ncount,it)+1
  !---> 
  if (ncount > ncount_max) then
     !
     do i=1, nx
     do j=1, nx
         mask_big(i,j)=0.
     enddo
     enddo
     !
     do k=1, ncount
         mask_big(ilist(k),jlist(k))= 1. 
     enddo
     !
     t_big = it
  endif
  !
  goto 1000
endif
! 
!---------------------  Begins new cloud -----------------------
! Initialize counter for new cloud
2001  ncell=ncell+1
is=ilist(ncell)
js=jlist(ncell)

if (ncount>ncell) then
!------->Go to check neighbors of the next cell in the list
         goto 2001
endif
!------------------------------------------------------------
!
! -------------------------------------------------------------------
! -------------------------------------------------------------------
!

!
!--------------------------------------->temporal sumation <----------------------------------
3333  dummy =  dummy
do k=icellmin, ncat
  dist(k)=dist(k)+ngroup(k, it)
enddo
do k=1, icellmin
  dist(k)=0.
enddo        
!********************************************************************
                         enddo ! <----- temporal loop
!********************************************************************
!
!
!---> count groups
ntot_groups=0
do k=icellmin, ncat
   ntot_groups= ntot_groups + dist(k)
enddo
!
!---> Bins
do k=icellmin, ncat
   bin(k)= (k-1)*delta_xy*delta_xy
enddo
!
!---> Average value
a_med = 0.
!
do k=icellmin, ncat
   a_med= a_med + dist(k) * bin(k)
enddo
a_med= a_med/ntot_groups
!
!---> convert dist to a density fuction
do k=icellmin, ncat
   pdf(k)= dist(k)/(delta_xy*delta_xy)
enddo
!
!---> PDF
s_dist=0.
do k=icellmin, ncat
   pdf(k)= dist(k)/ntot_groups
enddo
!
!---> Maximum area
a_max = 0.
k = ncat
1234 if (pdf(k)> 0.) then
          a_max = bin(k)
     else
          k = k-1
          goto 1234
     endif
!
goto 5678
!
! --------------------------------------------------------------------------------------------------
OPEN(1,FILE= 'cld_a.txt',STATUS='UNKNOWN')
123      Format( F9.4,F9.4,F9.4)
OPEN(24,FILE='cld_a.gra',STATUS='UNKNOWN',FORM='UNFORMATTED')
!----------------------------------------------------------------------------------------------------
!                                 displays  every deltacourser cells (0.1 km2)
       print*,'******************************************************************'
       print*,'            INTERVAL [Km2]        PERCENTAGE     ACC % '
 
        aux2=0
!       counts all clouds


        do k=icellmin, ncat
           aux2= aux2+dist(k)
        enddo
!
accum=0 
do k=icellmin,400, deltacoarser
      coarser=0
!     counts per category 
      !
      do m=k,k+deltacoarser-1
                coarser= coarser+dist(m)
      enddo
      percentage=100.*coarser/(1.*aux2)
      !
      if (percentage>lowpercent) then
                accum=accum+percentage
                print*,(k-1)*delta_xy*delta_xy,'-',(k+deltacoarser-1)*delta_xy*delta_xy,percentage, accum



                write(1,123)(k+deltacoarser-1)*delta_xy*delta_xy,percentage, accum

                WRITE(24) percentage
                WRITE(24) accum

      endif
enddo
CLOSE(24)  

5678 return 
end