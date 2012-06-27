module azavg_utils
!*********************************************************
! DATA
!*********************************************************

contains
!*********************************************************
! SUBROUTINES
!*********************************************************

!***************************************************************
! Utilities for the azimuthial averaging code.
!***************************************************************

!**********************************************************************
! RecordStormCenter()
!
! This routine will record the storm center for each time step
!

subroutine RecordStormCenter(Nx, Ny, Nz, Nt, Press, StmIx, StmIy, MinP)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(Nx,Ny,Nt) :: Press
  integer, dimension(Nt) :: StmIx, StmIy
  real, dimension(Nt) :: MinP

  integer :: it

  do it = 1, Nt
    call FindStormCenter(Nx, Ny, Nz, Nt, Press, it, StmIx(it), StmIy(it), MinP(it))
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

subroutine FindStormCenter(Nx, Ny, Nz, Nt, Press, iTime, ixStmCtr, iyStmCtr, MinP)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(Nx,Ny,Nt) :: Press
  integer :: iTime, ixStmCtr, iyStmCtr
  real :: MinP

  integer :: ix, iy, it

  it = iTime
  MinP = 1e10 ! ridiculously large pressure
  ixStmCtr = 0
  iyStmCtr = 0

  do ix = 1, Nx
    do iy = 1, Ny
      if (Press(ix,iy,it) .lt. MinP) then
        MinP = Press(ix,iy,it)
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

subroutine AzimuthalAverage(Nx, Ny, Nz, Nt, AvarNz, NumRbands, W, StmIx, StmIy, MinP, Avar, AzAvg, &
          Xcoords, Ycoords, RadialDist, RbandInc, WfilterMin, WfilterMax, UndefVal)

  use rhdf5_utils
  implicit none

  integer :: Nx, Ny, Nz, Nt
  integer :: AvarNz
  integer :: NumRbands
  real, dimension(Nx,Ny,Nz,Nt) :: W
  real, dimension(Nx,Ny,AvarNz,Nt) :: Avar
  real, dimension(NumRbands,AvarNz,Nt) :: AzAvg
  integer, dimension(Nt) :: StmIx, StmIy
  real, dimension(Nt) :: MinP
  real, dimension(Nx) :: Xcoords
  real, dimension(Ny) :: Ycoords
  real :: RadialDist, RbandInc, WfilterMin, WfilterMax, UndefVal

  integer :: ix, iy, iz, it
  real, dimension(NumRbands) :: Rcounts
  integer :: ir, iRband
  real :: DeltaX, DeltaY, Radius
  real, dimension(Nx,Ny,Nt) :: WmaxVals, WminVals
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

  do it = 1, Nt
    do ix = 1, Nx
      do iy = 1, Ny
        Wmax = W(ix,iy,1,it)
        Wmin = W(ix,iy,1,it)
        do iz = 2, Nz
          if (W(ix,iy,iz,it) .gt. Wmax) then
            Wmax = W(ix,iy,iz,it)
          end if
          if (W(ix,iy,iz,it) .lt. Wmin) then
            Wmin = W(ix,iy,iz,it)
          end if
        end do
        WmaxVals(ix,iy,it) = Wmax
        WminVals(ix,iy,it) = Wmin
      end do
    end do
  end do

  write (*,*) 'Averaging Data:'

  do it = 1, Nt
    if (modulo(it,10) .eq. 0) then
    write (*,*) '  Timestep: ', it
    write (*,'(a,i3,a,i3,a,g,a,g,a)') '    Storm Center: (', StmIx(it), ', ', StmIy(it), ') --> (', &
          Xcoords(StmIx(it)), ', ', Ycoords(StmIy(it)), ')'
    write (*,*) '    Minimum Pressue: ', MinP(it)
    end if
    do iz = 1, AvarNz
      ! For the averaging
      do ir = 1, NumRbands
        Rcounts(ir) = 0.0
        AzAvg(ir,iz,it) = 0.0
      end do

      ! Get the averages for this x-y plane
      do iy = 1, Ny
        do ix = 1, Nx
          ! Only keep the points where W max or W min are outside the filter interval
          if (((WmaxVals(ix,iy,it) .gt. WfilterMax) .or. &
               (WminVals(ix,iy,it) .lt. WfilterMin)) .and. &
              (Avar(ix,iy,iz,it) .ne. UndefVal)) then

             DeltaX = Xcoords(ix)-Xcoords(StmIx(it))
             DeltaY = Ycoords(iy)-Ycoords(StmIy(it))
             Radius = sqrt(DeltaX*DeltaX + DeltaY*DeltaY)
             iRband = int(Radius / RbandInc) + 1
             ! iRband will go from 1 to n, but n might extend beyond the last radial
             ! band due to approximations made in deriving RbandInc
             if (iRband .gt. NumRbands) then
               iRband = NumRbands
             end if

             AzAvg(iRband,iz,it) = AzAvg(iRband,iz,it) + Avar(ix,iy,iz,it)
             Rcounts(iRband) = Rcounts(iRband) + 1.0
          end if
        end do
      end do

      do ir = 1, NumRbands
        ! If we didn't put anything into an AzAvg slot then set it
        ! to the undefined value so the data isn't biased by trying to chose
        ! a default value.
        if (Rcounts(ir) .ne. 0.0) then
          AzAvg(ir,iz,it) = AzAvg(ir,iz,it) / Rcounts(ir)
        else
          AzAvg(ir,iz,it) = UndefVal
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

subroutine ConvertHorizVelocity(Nx, Ny, Nz, Nt, U, V, StmIx, StmIy, Xcoords, Ycoords, Avar, DoTangential)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(Nx,Ny,Nz,Nt) :: U, V, Avar
  integer, dimension(Nt) :: StmIx, StmIy
  real, dimension(Nx) :: Xcoords
  real, dimension(Ny) :: Ycoords
  logical :: DoTangential

  integer :: ix, iy, iz, it
  real :: StmX, StmY         ! x,y location relative to the storm center
  real :: Theta, Phi, Alpha  ! angle values, in radians
  real :: WindMag, WindX, WindY
  
  do it = 1, Nt
    do iz = 1, Nz
      do iy = 1, Ny
        StmY = Ycoords(iy) - Ycoords(StmIy(it))
        do ix = 1, Nx
          StmX = Xcoords(ix) - Xcoords(StmIx(it))

          WindX = U(ix,iy,iz,it)
          WindY = V(ix,iy,iz,it)

          Theta = atan2(StmY, StmX) ! Angle of radius line from horizontal
          Phi = atan2(WindY, WindX) ! Angle of wind vector from horizontal
          Alpha = Phi - Theta       ! Angle of wind vector from raduis line
                                    !    radius line is the line from storm center through
                                    !    the point (StmX, StmY)

          WindMag = sqrt(WindX**2 + WindY**2)
          
          if (DoTangential) then
            ! Tangential wind
            Avar(ix,iy,iz,it) = WindMag * sin(Alpha)
          else
            ! Radial wind
            Avar(ix,iy,iz,it) = WindMag * cos(Alpha)
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
! in the input var to a flat plane (x and y length values)
!
subroutine ConvertGridCoords(Nx, Ny, Nz, Nt, Lon, Lat, Xcoords, Ycoords)
  implicit none

  real, parameter :: RadiusEarth = 6378.1  ! km
  real, parameter :: PI = 3.14159

  integer :: Nx, Ny, Nz, Nt
  real, dimension(Nx) :: Lon
  real, dimension(Ny) :: Lat
  real, dimension(Nx) :: Xcoords
  real, dimension(Ny) :: Ycoords

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

  ! write a warning if the grid is located far away from the equator
  if ((abs(Lat(Ny)) .gt. 23.0) .or. (abs(Lat(1)) .gt. 23.0)) then
    write (*,*) 'Warning: extent of grid goes outside tropics, this code assumes grid is near equator'
    write (*,*) ''
  end if

  ! convert to radians
  ConvDeg2Rad = (2.0 * PI) / 360.0
  Theta1 = Lon(1) * ConvDeg2Rad
  Theta2 = Lon(2) * ConvDeg2Rad

  Phi1 = Lat(1) * ConvDeg2Rad
  Phi2 = Lat(2) * ConvDeg2Rad
  PhiN = Lat(Ny) * ConvDeg2Rad

  DeltaT = Theta2 - Theta1
  DeltaP = Phi2 - Phi1

  ! average of the arms at south lat and north lat
  RadiusX = (RadiusEarth * 0.5) * (cos(Phi1) + cos(PhiN))
  DeltaX = RadiusX * DeltaT
  DeltaY = RadiusEarth * DeltaP

  ! DeltaX and DeltaY are in units of RadiusEarth (km for now)
  do ix = 1, Nx
    Xcoords(ix) = real(ix - 1) * DeltaX
  end do

  do iy = 1, Ny
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
