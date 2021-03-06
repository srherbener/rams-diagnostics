!***************************************************************
! Program to take a vertical slice out of the GRADS data
!
! This program will read in GRADS data from a RAMS simulation
! extract a vertical slice along the given line segment
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
!   3. name of variable to extract data for
!   4. coordinates of line segment
!
! Output
!   The output will be two files which make a GRADS data set. One
!   file is the control file and the other is the binary data.
!

program main
  use gdata_utils
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarToSlice

  type (GradsControlFiles) :: GctlFiles

  ! Data arrays: need one the var we are doing the slice on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of the var in the GRADS input files
  ! data files: the first index is the file number, the second index is the
  ! var number
  type (GradsVar) :: Gvar, VarSlice
  real, dimension(:), allocatable :: Xcoords, Ycoords
  real :: Xstart, Xinc, Zstart, Zinc

  integer :: i

  integer :: ix, iy, iz, it
  integer :: ix1, ix2, iy1, iy2, NumPoints

  ! Get the command line argument
  call GetMyArgs(Infiles, OfileBase, VarToSlice, ix1, iy1, ix2, iy2, NumPoints)
  call String2List(Infiles, ':', GctlFiles%Fnames, MaxFiles, GctlFiles%Nfiles, 'input files')

  write (*,*) 'Extracting vertical slice of data out of GRADS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, GctlFiles%Nfiles
    write (*,*) '  ', i, ': ', trim(GctlFiles%Fnames(i))
  end do
  write (*,*) '  Output file base name: ', trim(OfileBase)
  write (*,*) '  Variable that is being extracted: ', trim(VarToSlice)
  write (*,*) '  Line segment defining slice location: '
  write (*,*) '    (', ix1, ', ', iy1, ') (', ix2, ', ', iy2, ')'
  write (*,*) '  Number of points to divide slice line segment into: ', NumPoints
  write (*,*) ''

  ! Read the GRADS data description files and collect the information about the data
  call ReadGradsCtlFiles(GctlFiles)

  call InitGvarFromGdescrip(GctlFiles, Gvar, VarToSlice)

  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', Gvar%Nx
  write (*,*) '  Number of y (latitude) points:           ', Gvar%Ny
  write (*,*) '  Number of z (vertical level) points:     ', Gvar%Nz
  write (*,*) '  Number of t (time) points:               ', Gvar%Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Gvar%Nx*Gvar%Ny*Gvar%Nz*Gvar%Nt
  write (*,*) ''

  write (*,*) 'Locations of variables in GRADS data (file, var number):'
  write (*,'(a17,a3,a,a2,i3,a1)') trim(VarToSlice), ': (', trim(Gvar%DataFile), ', ', Gvar%Vnum, ')'
  write (*,*) ''


  ! Read in the data for the vars using the description and location information
  call ReadGradsData(Gvar)

  ! Divide up the vslice line segment
  allocate(Xcoords(1:NumPoints), Ycoords(1:NumPoints))
  call DivideLineSegment(ix1, iy1, ix2, iy2, NumPoints, Xcoords, Ycoords)
!do i = 1, NumPoints
!  write (*,*) 'DEBUG: (i,x,y): ', i, Xcoords(i), Ycoords(i)
!end do

  ! Cut out the vertical slice
  call InitGradsVar(VarSlice, VarToSlice, NumPoints, 1, Gvar%Nz, Gvar%Nt, &
                    1.0, 1.0, 0.0, 1.0, Gvar%Zcoords, Gvar%Tstart, Gvar%Tinc, &
                    Gvar%UndefVal, '<NONE>', 0, 0)
  call CreateVslice(NumPoints, Gvar, VarSlice, Xcoords, Ycoords)

  ! Write out the slice in GRADS format
  call WriteGrads(VarSlice, OfileBase, 'vslice')

  stop
end

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   OfileBase - output GRADS file, base name for two files
!   VarToSlice - RAMS variable to do the column integration on
!   CoordSpec - end points of line segment defining location of vertical slice
!   NumPoints - number of points to extract along vslice line segment
!

subroutine GetMyArgs(Infiles, OfileBase, VarToSlice, ix1, iy1, ix2, iy2, NumPoints)
  use gdata_utils
  implicit none

  integer, parameter :: NumReqCoords=4

  character (len=*) :: Infiles, OfileBase, VarToSlice
  integer :: ix1, iy1, ix2, iy2, NumPoints

  character (len=128), dimension(1:NumReqCoords) :: CoordList
  integer :: Ncoords

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 5) then
    write (*,*) 'ERROR: must supply exactly 5 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: vslice <in_data_files> <out_data_file> <variable> <coords>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <variable>: name of variable that will be sliced out of data'
    write (*,*) '        <coords>: Coordinates of line segment defining the slice'
    write (*,*) '                  Colon separated list: <x1>:<y1>:<x2>:<y2>'
    write (*,*) '                  Use integer grid indecies'
    write (*,*) '        <num_points>: number of points to divide slice line segment into'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)
  call getarg(3, VarToSlice)

  ! coordinate spec
  call getarg(4, arg)
  call String2List(arg, ':', CoordList, NumReqCoords, Ncoords, 'coordinates')
  if (Ncoords .ne. NumReqCoords) then
    write (*,*) 'ERROR: Must specify exactly 4 coordinates'
    stop
  endif

  read (CoordList(1), '(i)') ix1
  read (CoordList(2), '(i)') iy1
  read (CoordList(3), '(i)') ix2
  read (CoordList(4), '(i)') iy2

  call getarg(5, arg)
  read (arg, '(i)') NumPoints
  if (NumPoints .le. 1) then
    write (*,*) 'ERROR: <num_points> must be > 1'
    stop
  endif

  return
end subroutine

!*********************************************************************
! DivideLineSegment()
!
! This routine will take the line segment end points and the number
! of points to divide it up into and create the lists of x,y coordinates
! of all the points.
!

subroutine DivideLineSegment(ix1, iy1, ix2, iy2, NumPoints, Xcoords, Ycoords)

  implicit none

  integer :: ix1, iy1, ix2, iy2
  integer :: NumPoints
  real, dimension(1:NumPoints) :: Xcoords, Ycoords

  real :: Xinc, Yinc
  integer :: i

  ! Figure out the x and y distances between adjacent points (Xinc, Yinc),
  ! Then assign the given end points to the ends of the output arrays
  ! (Xcoords, Ycoords) and fill in the middles using the Xinc, Yinc values.

  Xinc = float(ix2 - ix1) / float(NumPoints - 1)
  Yinc = float(iy2 - iy1) / float(NumPoints - 1)

  Xcoords(1) = float(ix1)
  Ycoords(1) = float(iy1)
  Xcoords(NumPoints) = float(ix2)
  Ycoords(NumPoints) = float(iy2)

  do i = 2, (NumPoints - 1)
    Xcoords(i) = Xcoords(i - 1) + Xinc
    Ycoords(i) = Ycoords(i - 1) + Yinc
  end do

  return
end subroutine

!**************************************************************************
! CreateSlice()
!
! This routine will cut out the slice from the GRADS data.
!
subroutine CreateVslice(NumPoints, Gvar, VarSlice, Xcoords, Ycoords)
  use gdata_utils
  implicit none

  real ::  BilinInterp

  integer :: NumPoints
  type (GradsVar) :: Gvar, VarSlice
  real, dimension(1:NumPoints) :: Xcoords, Ycoords

  integer ix, iy, iz, it, ip

  do it = 1, Gvar%Nt
    do iz = 1, Gvar%Nz
      ! Go through the x,y pairs (Xcoords,Ycoords) and interpolate
      ! from the original data (Gvar)
      do ip = 1, NumPoints
        VarSlice%Vdata(it,iz,ip,1) = BilinInterp(NumPoints,Gvar,Xcoords(ip),Ycoords(ip),iz,it)
!write(*,*) 'DEBUG(CV): (ip,iz,it,bi) ', ip,iz,it,VarSlice(it,iz,ip,1)
      end do
    end do
  end do

  return
end subroutine

!******************************************************************************
! BilinInterp()
!
! This function performs a bilinear interpolation of the horizontal data (specified
! by iz,it setting) given the grid coordinates in x,y
!

real function BilinInterp(NumPoints,Gvar,x,y,Izloc,Itloc)
  use gdata_utils
  implicit none

  integer :: NumPoints, Izloc, Itloc
  type (GradsVar) :: Gvar
  real :: x, y

  integer :: ix1, iy1, ix2, iy2
  real :: f1, f2, f3, f4, Xprop, Yprop

  ! Read x and y out of Xcoords, Ycoords - these will be in grid coordinates, so
  ! find the four surrounding (integer) coordinates that will index into Var.

  ix1 = int(x)
  iy1 = int(y)

  ! If ix1 is at the end of the Var array (equal to Nx) then set ix2 equal to ix1
  ! which effectively extends the boundary value of Var to the right. Ditto for iy2.
  ! Code before the call to this routine forces the endpoints in Xcoords, Ycoords
  ! to match exactly on grid points in Var. This will prevent ix1, iy1 from going too
  ! low in value (to the left or below the grid in Var), but not prevent ix2, iy2 from
  ! going too high (to the right or above the grid in Var).
  ix2 = ix1 + 1
  if (ix2 .gt. Gvar%Nx) then
    ix2 = ix1
  endif
  iy2 = iy1 + 1
  if (iy2 .gt. Gvar%Ny) then
    iy2 = iy1
  endif

  ! Grab f1 .. f4 starting at the lower left of the box (defined by ix1,iy2 ix2,iy2)
  ! and moving around counter clockwise.
  f1 = Gvar%Vdata(Itloc,Izloc,ix1,iy1)
  f2 = Gvar%Vdata(Itloc,Izloc,ix2,iy1)
  f3 = Gvar%Vdata(Itloc,Izloc,ix2,iy2)
  f4 = Gvar%Vdata(Itloc,Izloc,ix1,iy2)

  ! Find the proportion of x between ix1 and ix2, ditto for y.
  ! If ix1 equals ix2, then set Xprop to zero.

  if (ix1 .ne. ix2) then
    Xprop = (x - float(ix1)) / float(ix2 - ix1)
  else
    Xprop = 0.0
  endif
  if (iy1 .ne. iy2) then
    Yprop = (y - float(iy1)) / float(iy2 - iy1)
  else
    Yprop = 0.0
  endif

!write (*,*) 'DEBUG(BI): (x,y,ix1,iy1,ix2,iy2,f1,f2,f3,f4): ', x, y, ix1, iy1, ix2, iy2, f1, f2, f3, f4

  ! Do the interpolation of the point between the four points from Var (f1 .. f4)
  BilinInterp = (((1.0-Xprop) * (1.0-Yprop) * f1) + &
                 (Xprop       * (1.0-Yprop) * f2) + &
                 (Xprop       * Yprop       * f3) + &
                 ((1.0-Xprop) * Yprop       * f4))

  return
end function
