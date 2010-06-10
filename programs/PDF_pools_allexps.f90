include 'gDataTypes.h'
use GfileTypes

implicit none
!********************************************************************
!
real, parameter :: delta_temp= 4.5
!
!********************************************************************
!
integer, parameter :: MAX_STR_LEN = 256
integer, parameter :: MAX_NUM_EXP = 25
integer, parameter :: MAX_TIME_STEPS = 500
integer, parameter :: MAX_X_DIM = 500
integer, parameter :: MAX_Y_DIM = 500

integer :: iargc, Ierror
integer :: Nexp, Nt, Nx, Ny, Nz, Nvars
integer :: i, j, k, i_exp, it
integer :: t_big

character (len=MAX_STR_LEN) :: ConfigFile
character (len=MAX_STR_LEN) :: TempcVar
character (len=MAX_STR_LEN), dimension(1:MAX_NUM_EXP) :: GradsControlFiles

real, dimension(:,:,:,:), allocatable :: TempcData
real :: a_med, a_max, g_med, diff_max, aux1
real, dimension(1:MAX_X_DIM,1:MAX_Y_DIM) ::  mask_big
real, dimension(1:MAX_X_DIM,1:MAX_Y_DIM,1:MAX_TIME_STEPS) :: var
real, dimension(1:MAX_TIME_STEPS,1:MAX_NUM_EXP) :: tcmed
real, dimension(1:MAX_NUM_EXP) :: aux

type (GradsDataDescription), dimension(1:MAX_NUM_EXP) :: GdataDescrip
type (GradsVarLocation) :: TempcLoc

!********************************************************************
! Read in arugumets and configuration.
!********************************************************************

if (iargc() .ne. 1) then
  write (0,*) 'ERROR: must supply exactly one argument'
  write (0,*) ''
  write (0,*) 'Usage PDF_pools_allexps <config_file>'
  write (0,*) '   <config_file>: file containing description of input data, format:'
  write (0,*) '                  <n_exp>'
  write (0,*) '                  <var_name>'
  write (0,*) '                  <GRADS control file> <-- first one is control'
  write (0,*) '                  <GRADS control file>'
  write (0,*) '                  ...'
  stop
endif

call getarg(1, ConfigFile)

write (0,*) 'Reading configuration file: ', trim(ConfigFile)
open (unit=10, file=ConfigFile, status='old', form='formatted')
read (10,'(i)') Nexp
read (10,'(a)') TempcVar
if (Nexp .gt. MAX_NUM_EXP) then
  write (0,*) 'ERROR: Number of experiments, ', Nexp, ', exceeds maximum allowed, ', MAX_NUM_EXP
  stop
else
  do i_exp=1, Nexp
    read(10,'(a)') GradsControlFiles(i_exp)
  enddo
endif

write (0,*) 'PDF pools all experiments:'
write (0,*) '  Number of experiments: ', Nexp
write (0,*) '  Variable name: ', trim(TempcVar)
write (0,*) '  GRADS control files:'
do i_exp=1, Nexp
  write (0,*) '    ', trim(GradsControlFiles(i_exp))
enddo
write (0,*) ''

!********************************************************
! For each experiment:
!    Read in the surface temperature data
!    Calculate the mean temperature at each time step
!    Mark the cold pools at each time step
!    Calculate the PDF of cold pools at each time step
!********************************************************

do i_exp=1, Nexp

  write (0,*) 'Experiment: ', i_exp
  write (0,*) ''
  write (0,*) 'Reading GRADS control file: ', trim(GradsControlFiles(i_exp))
  call ReadGradsCtlFile(GradsControlFiles(i_exp), GdataDescrip(1))
  write (0,*) ''
  call CheckDataDescrip_VS(GdataDescrip, 1, Nx, Ny, Nz, Nt, Nvars, TempcLoc, TempcVar)
  write (0,*) ''

  write (0,*) 'Gridded data information:'
  write (0,*) '  Number of x (longitude) points:      ', Nx
  write (0,*) '  Number of y (latitude) points:       ', Ny
  write (0,*) '  Number of z (vertical level) points: ', Nz
  write (0,*) '  Number of t (time) points:           ', Nt
  write (0,*) '  Total number of grid variables       ', Nvars
  write (0,*) ''
  write (0,*) 'Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (0,*) ''
  write (0,*) 'Location of tempc variable in GRADS data file (file number, var number):'
  write (0,'(a20,i3,a2,i3,a1)') 'tempc: (', TempcLoc%Fnum, ', ', TempcLoc%Vnum, ')'
  write (0,*) ''

  allocate (TempcData(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (0,*) 'ERROR: Data array allocation failed'
    stop
  endif

  call ReadGradsData(GdataDescrip, TempcVar, TempcLoc, TempcData, Nx, Ny, Nz, Nt)
 
  ! At this point TempcData contains the temperature for one whole
  ! experirmental run - all x,y,z,t values. We want just the surface temp
  ! at each time step -> set k=1 for this loop.

  ! Mean sfc temp
  k = 1
  diff_max =0.
  do it = 1, Nt
    tcmed(it,i_exp) = 0.
    do i=1, Nx 
      do j=1, Ny
        tcmed(it,i_exp) = tcmed(it,i_exp) + TempcData(i,j,k,it)
      enddo
    enddo
    tcmed(it,i_exp) = tcmed(it,i_exp)/real(Nx*Ny) 

    ! Mark cold pools. var() holds an image of the cold pools, like pixels, you get a 1 where the
    ! temp is colder than the mean temp by delta_temp degrees C, otherwise you get a 0.
    do i=1, Nx 
      do j=1, Ny
        if ((TempcData(i,j,k,it) - tcmed(it,i_exp))< -delta_temp) then
          var(i,j,it) = 1.
          diff_max=min (diff_max,(TempcData(i,j,k,it) - tcmed(it,i_exp))    )
        else
          var(i,j,it) = 0.
        endif
      enddo
    enddo
  enddo ! it = 1, Nt

  ! Calculate the PDF. i_exp =1 is the control experiment
  ! The call to PDF_a is:
  !>>> PDF_a(nx, ny, ntimes, delta_xy, var, ntinit, ntend, icellmin, a_med, a_max, mask_big,t_big) <<<    
  !
  if (i_exp ==1) then
    print*,'******************************************************************' 
    print*, '              results for DT=', delta_temp
  endif
  !
  print*,'******************************************************************' 
  print*, trim(GradsControlFiles(i_exp))
  !
  call PDF_a(MAX_X_DIM, MAX_Y_DIM, MAX_TIME_STEPS, Nx, Ny, Nt, 1.5, var, 1, Nt, 4, a_med, a_max, mask_big, t_big)  
  !
  ! This code looks like it produces data that is not used
  ! Development of this code must be in progress at this point
  aux1=0
  g_med = 0.
  do i=1, Nx 
    do j=1, Ny
      if (mask_big(i,j)==1.) then
        g_med = g_med + TempcData(i,j,k,Nt)  !<--- this is on last time step from loop above
                                             !     is this really what we want?
        aux1 = aux1 + 1.
      endif
    enddo
  enddo
  g_med = g_med/aux1
  !
  print*,i_exp, a_med, a_max, t_big

  deallocate (TempcData, stat=Ierror)
  if (Ierror .ne. 0) then
    write (0,*) 'ERROR: Data array memory de-allocation failed'
    stop
  endif
enddo  ! i_exp = 1, Nexp

!*************
!
print*,'******************************************************************' 
do it=1, Nt
   !print*, it, tcmed(it,2)-tcmed(it,1), tcmed(it,3)-tcmed(it,1), tcmed(it,4)-tcmed(it,1)
enddo
!
!
do i=1, Nexp
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
subroutine PDF_a(max_x_dim, max_y_dim, max_time_steps, nx, ny, ntimes, delta_xy, var, ntinit, ntend, icellmin, a_med, a_max, mask_big,t_big)  
implicit none
integer :: nx, ny, ntimes,ntinit, ntend, icellmin,t_big
integer :: max_x_dim, max_y_dim, max_time_steps
real , dimension(max_x_dim,max_y_dim,max_time_steps) :: var
real :: delta_xy
!
real, dimension(max_x_dim,max_y_dim) :: mask_big
real, dimension(nx,ny) :: group                 ! cell state and mask_big for biggest group
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
open(unit=11,file='cld_a.txt',status='REPLACE')
123      format(f9.4,f9.4,f9.4)
open(unit=12,file='cld_a.gra',status='REPLACE',form='UNFORMATTED')
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



                write(11,123)(k+deltacoarser-1)*delta_xy*delta_xy,percentage, accum

                write(12) percentage
                write(12) accum

      endif
enddo
close(24)  

5678 return 
end
