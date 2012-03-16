module azavg_utils
!*********************************************************
! DATA
!*********************************************************

contains
!*********************************************************
! SUBROUTINES
!*********************************************************

!***************************************************************
! Utilities for the azimuthial averaging and CFAD code.
!
! This program will read in GRADS data from a RAMS simulation, find
! the storm center and perform azimuthial averaging on the given
! quantity.
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
!   3. number of radial bands to split data into
!   4. RAMS quantity to perform the averaging on
!
! Output
!   The output will be two files which make a GRADS data set. One
!   file is the control file and the other is the binary data.
!

!**********************************************************************
! RecordStormCenter()
!
! This routine will record the storm center for each time step
!

subroutine RecordStormCenter(Press, StmIx, StmIy, MinP)
  use gdata_utils
  implicit none

  type (GradsVar) :: Press
  integer, dimension(:), allocatable :: StmIx, StmIy
  real, dimension(:), allocatable :: MinP

  integer :: it
  integer :: Ierror

  allocate (StmIx(1:Press%Nt), StmIy(1:Press%Nt), MinP(1:Press%Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: StmIx, StmIy, MinP data array memory allocation failed'
    stop
  end if

  do it = 1, Press%Nt
    call FindStormCenter(Press, it, StmIx(it), StmIy(it), MinP(it))
  end do

  return
end subroutine

!**********************************************************************
! FindStormCenter
!
! This routine will locate the storm center using the simple hueristic
! of the center being where the minimum surface pressure exists.
!
! Argument it holds the time step that you want to analyze. (iStmCtr, jStmCtr)
! hold the grid position of the minumum pressure value on the first vertical
! level (iz = 1).
!

subroutine FindStormCenter(Press, iTime, ixStmCtr, iyStmCtr, MinP)
  use gdata_utils
  implicit none

  type (GradsVar) :: Press
  integer :: iTime, ixStmCtr, iyStmCtr
  real :: MinP

  integer :: ix, iy, iz, it

  iz = 1
  it = iTime
  MinP = 1e10 ! ridiculously large pressure
  ixStmCtr = 0
  iyStmCtr = 0
  
  do ix = 1, Press%Nx
    do iy = 1, Press%Ny
      if (Press%Vdata(it,iz,ix,iy) .lt. MinP) then
        MinP = Press%Vdata(it,iz,ix,iy)
        ixStmCtr = ix
        iyStmCtr = iy
      end if
    end do
  end do
  
  return
end subroutine

!*************************************************************************
! AzimuthalAverage
!
! This routine will do the azimuthal averaging. It will create a 4-D array
! for its output (to keep GRADS write routine consistent/general) even though it
! really is 3-D data (radial band, z level, time). It will put the radial bands 
! in the x dimension and make the y dimension 1 long.
!
! This routine will step though each z and t value and compute the azimuthal
! average of the x-y plane at each of these points. Then it will store that
! in the output array so that you get azimuthal averaging for every original
! z and t point.
!

subroutine AzimuthalAverage(NumRbands, W, StmIx, StmIy, MinP, Avar, AzAvg, &
          Xcoords, Ycoords, RadialDist, RbandInc, WfilterMin, WfilterMax, UndefVal)

  use gdata_utils
  implicit none

  integer :: NumRbands
  type (GradsVar) :: W, Avar, AzAvg
  integer, dimension(1:W%Nt) :: StmIx, StmIy
  real, dimension(1:W%Nt) :: MinP
  real, dimension(1:W%Nx) :: Xcoords
  real, dimension(1:W%Ny) :: Ycoords
  real :: RadialDist, RbandInc, WfilterMin, WfilterMax, UndefVal

  integer :: ix, iy, iz, it
  real, dimension(1:NumRbands) :: Rcounts
  integer :: ir, iRband
  real :: DeltaX, DeltaY, Radius
  real, dimension(1:W%Nt,1:W%Nx,1:W%Ny) :: WmaxVals, WminVals
  real :: Wmax, Wmin

  ! Mask the input data (Avar) with the W data, ie if the corresponding
  ! W data is outside the interval defined by WfilterMin, WfilterMax,
  ! use the Avar data in the averaging; otherwise skip the Avar data
  !
  ! The storm center is taken to be the min pressure location of the
  ! x-y plane on the surface (iz .eq. 1)

  ! We want to filter where convection is likely taking place. Approximating
  ! this by filtering on larger w values. However, the larger w values won't
  ! persist through the entire vertical column. Handle this by finding the 
  ! min and max w at each time step in each column - and use these values
  ! to do the filtering.

  do it = 1, W%Nt
    do ix = 1, W%Nx
      do iy = 1, W%Ny
        Wmax = W%Vdata(it,1,ix,iy)
        Wmin = W%Vdata(it,1,ix,iy)
        do iz = 2, W%Nz
          if (W%Vdata(it,iz,ix,iy) .gt. Wmax) then
            Wmax = W%Vdata(it,iz,ix,iy)
          end if
          if (W%Vdata(it,iz,ix,iy) .lt. Wmin) then
            Wmin = W%Vdata(it,iz,ix,iy)
          end if
        end do
        WmaxVals(it,ix,iy) = Wmax
        WminVals(it,ix,iy) = Wmin
      end do
    end do
  end do

  write (*,*) 'Averaging Data:'

  do it = 1, W%Nt
    if (modulo(it,10) .eq. 0) then
    write (*,*) '  Timestep: ', it
    write (*,'(a,i3,a,i3,a,g,a,g,a)') '    Storm Center: (', StmIx(it), ', ', StmIy(it), ') --> (', &
          Xcoords(StmIx(it)), ', ', Ycoords(StmIy(it)), ')'
    write (*,*) '    Minimum Pressue: ', MinP(it)
    end if
    do iz = 1, Avar%Nz

      ! For the averaging
      do ir = 1, NumRbands
        Rcounts(ir) = 0.0
        AzAvg%Vdata(it, iz, ir, 1) = 0.0
      end do

      ! Get the averages for this x-y plane
      do iy = 1, W%Ny
        do ix = 1, W%Nx
          ! Only keep the points where W max or W min are outside the filter interval
          if (((WmaxVals(it,ix,iy) .gt. WfilterMax) .or. &
               (WminVals(it,ix,iy) .lt. WfilterMin)) .and. &
              (Avar%Vdata(it,iz,ix,iy) .ne. UndefVal)) then
             DeltaX = Xcoords(ix)-Xcoords(StmIx(it))
             DeltaY = Ycoords(iy)-Ycoords(StmIy(it))
             Radius = sqrt(DeltaX*DeltaX + DeltaY*DeltaY)
             iRband = int(Radius / RbandInc) + 1
             ! iRband will go from 1 to n, but n might extend beyond the last radial
             ! band due to approximations made in deriving RbandInc
             if (iRband .gt. NumRbands) then
               iRband = NumRbands
             end if

             AzAvg%Vdata(it,iz,iRband,1) = AzAvg%Vdata(it,iz,iRband,1) + Avar%Vdata(it,iz,ix,iy)
             Rcounts(iRband) = Rcounts(iRband) + 1.0
          end if
        end do
      end do

      do ir = 1, NumRbands
        ! If we didn't put anything into an AzAvg slot then set it
        ! to the undefined value so the data isn't biased by trying to chose
        ! a default value.
        if (Rcounts(ir) .ne. 0.0) then
          AzAvg%Vdata(it,iz,ir,1) = AzAvg%Vdata(it,iz,ir,1) / Rcounts(ir)
        else
          AzAvg%Vdata(it,iz,ir,1) = UndefVal
        end if
      end do

    end do
  end do
  write (*,*) ''

  return
end subroutine

!*******************************************************************************
! ConvertHorizVelocity()
!
! This routine will convert the horizontal velocity vectors (described in U and V)
! into tangential or radial components given the storm center location.
!

subroutine ConvertHorizVelocity(U, V, StmIx, StmIy, Xcoords, Ycoords, Avar, DoTangential)
  use gdata_utils
  implicit none

  type (GradsVar) :: U, V, Avar
  integer, dimension(1:U%Nt) :: StmIx, StmIy
  real, dimension(1:U%Nx) :: Xcoords
  real, dimension(1:U%Ny) :: Ycoords
  logical :: DoTangential

  integer :: ix, iy, iz, it
  real :: StmX, StmY         ! x,y location relative to the storm center
  real :: Theta, Phi, Alpha  ! angle values, in radians
  real :: WindMag, WindX, WindY
  
  do it = 1, U%Nt
    do iz = 1, U%Nz
      do iy = 1, U%Ny
        StmY = Ycoords(iy) - Ycoords(StmIy(it))
        do ix = 1, U%Nx
          StmX = Xcoords(ix) - Xcoords(StmIx(it))

          WindX = U%Vdata(it,iz,ix,iy)
          WindY = V%Vdata(it,iz,ix,iy)

          Theta = atan2(StmY, StmX) ! Angle of radius line from horizontal
          Phi = atan2(WindY, WindX) ! Angle of wind vector from horizontal
          Alpha = Phi - Theta       ! Angle of wind vector from raduis line
                                    !    radius line is the line from storm center through
                                    !    the point (StmX, StmY)

          WindMag = sqrt(WindX**2 + WindY**2)
          
          if (DoTangential) then
            ! Tangential wind
            Avar%Vdata(it,iz,ix,iy) = WindMag * sin(Alpha)
          else
            ! Radial wind
            Avar%Vdata(it,iz,ix,iy) = WindMag * cos(Alpha)
          end if
        end do
      end do
    end do
  end do
  
  return
end subroutine

!******************************************************************
! ConvertGridCoords()
!
! This routine will convert the longitude, latitude angle values
! in the input GRADS var to a flat plane (x and y length values)
!
subroutine ConvertGridCoords(Gvar, Xcoords, Ycoords)
  use gdata_utils
  implicit none

  real, parameter :: RadiusEarth = 6378.1  ! km
  real, parameter :: PI = 3.14159

  type (GradsVar) :: Gvar
  real, dimension(:), allocatable :: Xcoords, Ycoords

  integer :: ix, iy
  integer :: Ierror
  real :: ConvDeg2Rad;
  real :: DeltaX, DeltaY
  real :: DeltaT, DeltaP  ! Theta is longitude angle, Phi is latitude angle, radians
  real :: RadiusX
  real :: PhiN, Phi1, Phi2, Theta1, Theta2

  ! The Xcoords in the Gvar structure are longitude angles in degrees
  ! The Ycoords in the Gvar structure are latitude angles in degrees
  !
  ! To convert, figure out what DeltaX and DeltaY are according to the longitude, latitude
  ! angles. Call the lower left point of the grid (0,0) and just add in the delta values to
  ! get the remaining coordinate values. Put the x,y values in km.
  !
  ! Technically, the DeltaX values change for each unique latitude angle since DeltaX depends
  ! on the arm perpendicular to the axis of rotation of the Earth (not the center of Earth).
  ! However, since the storms are near the equator this arm doesn't change much in length.
  ! Approximate the spherical surface with a rectangle using the delta x that is an average
  ! of the deltax you find at the southernmost latidute and northernmost latitude. This will
  ! greatly simplify the code by allowing the use of just one list of x coordinates for the 
  ! entire grid.
  !
  !   DeltaX = (RadiusEarth * cos(Phi)) * DeltaTheta
  !   DeltaY = RadiusEarth * DeltaPhi
  !   angle values are in radians
  allocate (Xcoords(1:Gvar%Nx), Ycoords(1:Gvar%Ny), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Xcoords, Ycoords data array memory allocation failed'
    stop
  end if

  ! write a warning if the grid is located far away from the equator
  if (abs(Gvar%Ycoords(Gvar%Ny)) .gt. 23.0) then
    write (*,*) 'Warning: extent of grid goes outside tropics, this code assumes grid is near equator'
    write (*,*) ''
  end if

  ! convert to radians
  ConvDeg2Rad = (2.0 * PI) / 360.0
  Theta1 = Gvar%Xcoords(1) * ConvDeg2Rad
  Theta2 = Gvar%Xcoords(2) * ConvDeg2Rad

  Phi1 = Gvar%Ycoords(1) * ConvDeg2Rad
  Phi2 = Gvar%Ycoords(2) * ConvDeg2Rad
  PhiN = Gvar%Ycoords(Gvar%Ny) * ConvDeg2Rad

  DeltaT = Theta2 - Theta1
  DeltaP = Phi2 - Phi1

  ! average of the arms at south lat and north lat
  RadiusX = (RadiusEarth * 0.5) * (cos(Phi1) + cos(PhiN))
  DeltaX = RadiusX * DeltaT
  DeltaY = RadiusEarth * DeltaP

  ! DeltaX and DeltaY are in units of RadiusEarth (km for now)
  do ix = 1, Gvar%Nx
    Xcoords(ix) = real(ix - 1) * DeltaX
  end do

  do iy = 1, Gvar%Ny
    Ycoords(iy) = real(iy - 1) * DeltaY
  end do

  return
end subroutine

!*****************************************************************************
! InsideCylVol()
!
! This function will see if a given grid location falls within a cylindrical
! volume relative to the storm center.
!

logical function InsideCylVol(Nx, Ny, Nz, Nt, Ix, Iy, Iz, It, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ,  StmIx, StmIy, Xcoords, Ycoords, Zcoords)
  implicit none

  real, parameter :: PI = 3.141592654

  integer :: Nx, Ny, Nz, Nt, Ix, Iy, Iz, It
  real :: MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real, dimension(1:Nz) :: Zcoords

  real :: dX, dY, Radius, Phi, Z

  dX = Xcoords(Ix) - Xcoords(StmIx(It))
  dY = Ycoords(Iy) - Ycoords(StmIy(It))
  Radius = sqrt(dX*dX + dY*dY)
  ! atan2 returns a value in radians between -PI and +PI
  ! convert to value between 0 and 2*PI
  Phi = atan2(dY, dX)
  if (Phi .lt. 0.0) then
    Phi = Phi + (2.0 * PI)
  end if
  Z = Zcoords(Iz)

  InsideCylVol = .false.
  if (((Radius .ge. MinR) .and. (Radius .le. MaxR)) .and. &
      ((Phi .ge. MinPhi) .and. (Phi .le. MaxPhi)) .and. &
      ((Z .ge. MinZ) .and. (Z .le. MaxZ))) then
    InsideCylVol = .true.
  end if
  return
end function

end module azavg_utils
