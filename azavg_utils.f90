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

subroutine RecordStormCenter(Nx, Ny, Nz, Nt, Press, StmIx, StmIy, MinP)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: Press
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nt) :: MinP

  integer it

  do it = 1, Nt
    call FindStormCenter(Press, Nx, Ny, Nz, Nt, it, StmIx(it), StmIy(it), MinP(it))
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

subroutine FindStormCenter(Press, Nx, Ny, Nz, Nt, iTime, ixStmCtr, iyStmCtr, MinP)
  implicit none

  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: Press
  integer :: Nx, Ny, Nz, Nt
  integer :: iTime, ixStmCtr, iyStmCtr
  real :: MinP

  integer :: ix, iy, iz, it

  iz = 1
  it = iTime
  MinP = 1e10 ! ridiculously large pressure
  ixStmCtr = 0
  iyStmCtr = 0
  
  do ix = 1, Nx
    do iy = 1, Ny
      if (Press(ix,iy,iz,it) .lt. MinP) then
        MinP = Press(ix,iy,iz,it)
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

subroutine AzimuthalAverage(Nx, Ny, Nz, Nt, NumRbands, W, StmIx, StmIy, MinP, Avar, AzAvg, &
          Xcoords, Ycoords, RadialDist, RbandInc, WfilterMin, WfilterMax, UndefVal)

  implicit none

  integer :: Nx, Ny, Nz, Nt, NumRbands
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: W, Avar
  real, dimension(1:NumRbands,1:1,1:Nz,1:Nt) :: AzAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nt) :: MinP
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real :: RadialDist, RbandInc, WfilterMin, WfilterMax, UndefVal

  integer :: ix, iy, iz, it
  real, dimension(1:NumRbands) :: Rcounts
  integer :: ir, iRband
  real :: DeltaX, DeltaY, Radius
  real, dimension(1:Nx,1:Ny,1:Nt) :: WmaxVals, WminVals
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
    do iz = 1, Nz

      ! For the averaging
      do ir = 1, NumRbands
        Rcounts(ir) = 0.0
        AzAvg(ir, 1, iz, it) = 0.0
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

             AzAvg(iRband, 1, iz, it) = AzAvg(iRband, 1, iz, it) + Avar(ix, iy, iz, it)
             Rcounts(iRband) = Rcounts(iRband) + 1.0
          end if
        end do
      end do

      do ir = 1, NumRbands
        ! If we didn't put anything into an AzAvg slot then set it
        ! to the undefined value so the data isn't biased by trying to chose
        ! a default value.
        if (Rcounts(ir) .ne. 0.0) then
          AzAvg(ir, 1, iz, it) = AzAvg(ir, 1, iz, it) / Rcounts(ir)
        else
          AzAvg(ir, 1, iz, it) = UndefVal
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
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: U, V, Avar
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
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
! in the input GRADS grid to a flat plane (x and y length values)
!
subroutine ConvertGridCoords(Nx, Ny, GdataDescrip , Xcoords, Ycoords)
  use GfileTypes
  implicit none

  real, parameter :: RadiusEarth = 6378.1  ! km
  real, parameter :: PI = 3.14159

  integer :: Nx, Ny
  type (GradsDataDescription) :: GdataDescrip
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords

  integer :: ix, iy
  real :: ConvDeg2Rad;
  real :: DeltaX, DeltaY
  real :: DeltaT, DeltaP  ! Theta is longitude angle, Phi is latitude angle, radians
  real :: RadiusX
  real :: PhiN, Phi1, Phi2, Theta1, Theta2

  ! The Xcoords in the GdataDescrip structure are longitude angles in degrees
  ! The Ycoords in the GdataDescrip structure are latitude angles in degrees
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

  if (abs(GdataDescrip%Ycoords(Ny)) .gt. 23.0) then
    write (*,*) 'Warning: extent of grid goes outside tropics, this code assumes grid is near equator'
    write (*,*) ''
  end if

  ! convert to radians
  ConvDeg2Rad = (2.0 * PI) / 360.0
  Theta1 = GdataDescrip%Xcoords(1) * ConvDeg2Rad
  Theta2 = GdataDescrip%Xcoords(2) * ConvDeg2Rad

  Phi1 = GdataDescrip%Ycoords(1) * ConvDeg2Rad
  Phi2 = GdataDescrip%Ycoords(2) * ConvDeg2Rad
  PhiN = GdataDescrip%Ycoords(Ny) * ConvDeg2Rad

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
