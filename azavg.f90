!***************************************************************
! Program to do azimuthial averaging
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

include 'gDataTypes.h'

program main
  use GfileTypes
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  integer NumRbands
  real :: Wthreshold
  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarToAvg
  logical :: DoHorizVel, DoTangential

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles
  integer :: Nfiles

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars
  type (GradsOutDescription) :: GoutDescrip

  ! Data arrays: need one for w (vertical velocity), press (pressure)
  ! and the var we are doing the averaging on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: U, V, W, Press, Avar, AzAvg
  integer, dimension(:), allocatable :: StmIx, StmIy
  real, dimension(:), allocatable :: MinP, Xcoords, Ycoords
  type (GradsVarLocation) :: Wloc, Uloc, Vloc, PressLoc, VarLoc

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it

  real :: DeltaX, DeltaY, RadialDist, RbandInc
  real :: TestData, TestGridX, TestGridY

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBase, NumRbands, Wthreshold, VarToAvg, DoHorizVel, DoTangential)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Calculating azimuthal average for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  Number of radial bands: ', NumRbands
  write (*,*) '  W threshold :           ', Wthreshold
  if (DoHorizVel) then
    if (DoTangential) then
      write (*,*) '  RAMS variable that is being averaged: Tangential Horizontal Velocity'
    else
      write (*,*) '  RAMS variable that is being averaged: Radial Horizontal Velocity'
    end if
  else
    write (*,*) '  RAMS variable that is being averaged: ', trim(VarToAvg)
  end if
  write (*,*) ''

  if (trim(VarToAvg) /= 'test') then
    ! Not running a test so grab the data from GRADS input files and perform the averaging.

    ! Read the GRADS data description files and collect the information about the data
    do i = 1, Nfiles
      write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
      call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
    end do
    write (*,*) ''
  
    call CheckDataDescrip(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, Uloc, Vloc, Wloc, PressLoc, VarLoc, VarToAvg, DoHorizVel)
  
    write (*,*) 'Gridded data information:'
    write (*,*) '  Number of x (longitude) points:          ', Nx
    write (*,*) '  Number of y (latitude) points:           ', Ny
    write (*,*) '  Number of z (vertical level) points:     ', Nz
    write (*,*) '  Number of t (time) points:               ', Nt
    write (*,*) '  Total number of grid variables:          ', Nvars
    write (*,*) ''
    write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
    write (*,*) ''
  
    write (*,*) 'Locations of variables in GRADS data (file number, var number):'
    write (*,'(a20,i3,a2,i3,a1)') 'w: (', Wloc%Fnum, ', ', Wloc%Vnum, ')'
    write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'
    if (DoHorizVel) then
      write (*,'(a20,i3,a2,i3,a1)') 'speed - u: (', Uloc%Fnum, ', ', Uloc%Vnum, ')'
      write (*,'(a20,i3,a2,i3,a1)') 'speed - v: (', Vloc%Fnum, ', ', Vloc%Vnum, ')'
    else
      write (*,'(a17,a3,i3,a2,i3,a1)') trim(VarToAvg), ': (', VarLoc%Fnum, ', ', VarLoc%Vnum, ')'
    end if
    write (*,*) ''
  
    ! Allocate the data arrays and read in the data from the GRADS data files
    allocate (U(1:Nx,1:Ny,1:Nz,1:Nt), V(1:Nx,1:Ny,1:Nz,1:Nt), W(1:Nx,1:Ny,1:Nz,1:Nt), &
              Press(1:Nx,1:Ny,1:Nz,1:Nt), Avar(1:Nx,1:Ny,1:Nz,1:Nt), &
              StmIx(1:Nt), StmIy(1:Nt), MinP(1:Nt), Xcoords(1:Nx), Ycoords(1:Ny), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory allocation failed'
      stop
    end if
  
    ! Convert the GRADS grid coordinates from longitude, latitude to flat plane (x and y).
    call ConvertGridCoords(Nx, Ny, GdataDescrip(1), Xcoords, Ycoords)
  
    write (*,*) 'Horzontal Grid Coordinate Info:'
    write (*,*) '  X Range (min lon, max lon) --> (min x, max x): '
    write (*,*) '    ', GdataDescrip(1)%Xcoords(1), GdataDescrip(1)%Xcoords(Nx), Xcoords(1), Xcoords(Nx)
    write (*,*) '  Y Range: '
    write (*,*) '    ', GdataDescrip(1)%Ycoords(1), GdataDescrip(1)%Ycoords(Ny), Ycoords(1), Ycoords(Ny)
    write (*,*) ''
  
    ! Figure out how big to make each radial band. Assume that the storm center stays near the center of
    ! the grid domain. Take the diagonal distance of the domain, cut it in half and then break this up
    ! into NumRbands sections of equal length.
  
    DeltaX = Xcoords(Nx) - Xcoords(1)
    DeltaY = Ycoords(Ny) - Ycoords(1)
    RadialDist = sqrt(DeltaX*DeltaX + DeltaY*DeltaY) / 2.0
    RbandInc = RadialDist / real(NumRbands)
  
    write (*,*) 'Radial band information:'
    write (*,*) '  Delta x of domain:     ', DeltaX
    write (*,*) '  Delta y of domain:     ', DeltaY
    write (*,*) '  Radial distance:       ', RadialDist
    write (*,*) '  Radial band increment: ', RbandInc
    write (*,*) ''
  
    ! Read in the data for the vars using the description and location information
    call ReadGradsData(GdataDescrip, 'w', Wloc, W, Nx, Ny, Nz, Nt)
    call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
    call RecordStormCenter(Nx, Ny, Nz, Nt, Press, StmIx, StmIy, MinP)
    if (DoHorizVel) then
      call ReadGradsData(GdataDescrip, 'u', Uloc, U, Nx, Ny, Nz, Nt)
      call ReadGradsData(GdataDescrip, 'v', Vloc, V, Nx, Ny, Nz, Nt)
      call ConvertHorizVelocity(Nx, Ny, Nz, Nt, U, V, StmIx, StmIy, &
                      Xcoords, Ycoords, Avar, DoTangential)
    else
      call ReadGradsData(GdataDescrip, VarToAvg, VarLoc, Avar, Nx, Ny, Nz, Nt)
    end if
  
    ! Allocate the output array and compute the azimuthal averaging
    allocate (AzAvg(1:NumRbands,1:1,1:Nz,1:Nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Ouput data array memory allocation failed'
      stop
    end if
  else
    ! Performing a test, load up the variables with data that will produce a known result and
    ! finish out the program (averaging and output).

    ! Create one horizontal slice, 100 x 100 tiles, one z level, 5 time points.
    !   Make the x,y increments 1km
    ! Move the storm center around a little bit but keep it near the center.
    ! Make 10 radial bands, 

    Nx = 101
    Ny = 101
    Nz = 1
    Nt = 5

    NumRbands = 10
    TestGridX = 100.0
    TestGridY = 100.0

    allocate (U(1:Nx,1:Ny,1:Nz,1:Nt), V(1:Nx,1:Ny,1:Nz,1:Nt), W(1:Nx,1:Ny,1:Nz,1:Nt), &
              Press(1:Nx,1:Ny,1:Nz,1:Nt), Avar(1:Nx,1:Ny,1:Nz,1:Nt), &
              StmIx(1:Nt), StmIy(1:Nt), MinP(1:Nt), Xcoords(1:Nx), Ycoords(1:Ny), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: TEST: Data array memory allocation failed'
      stop
    end if
    allocate (AzAvg(1:NumRbands,1:1,1:Nz,1:Nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Ouput data array memory allocation failed'
      stop
    end if

    RadialDist = sqrt(TestGridX*TestGridX + TestGridY*TestGridY) / 2.0
    RbandInc = RadialDist / real(NumRbands)

    ! load up the coordinates
    DeltaX = TestGridX / real(Nx - 1)
    DeltaY = TestGridY / real(Ny - 1)
    do ix = 1, Nx
      Xcoords(ix) = real(ix-1) * DeltaX
    end do
    do iy = 1, Ny
      Ycoords(iy) = real(iy-1) * DeltaY
    end do
    do iz = 1, Nz
      GdataDescrip(1)%Zcoords(iz) = iz
    end do
    GdataDescrip(1)%Tstart = '00:00Z24aug1998'
    GdataDescrip(1)%Tinc =  '01mn'

    ! Undefined value
    GdataDescrip(1)%UndefVal = 10.0e30

    do it = 1, Nt
      do iz = 1, Nz
        do iy = 1, Ny
          do ix = 1, Nx

            ! Storm center drifts from roughly the center of the grid
            StmIx(it) = int(Nx/2) + it
            StmIy(it) = int(Ny/2) - it
            MinP(it) = 980.0 - real(it)  ! Just make it up, it's only used in diagnostic msg
                                         ! in AzimuthalAverage().

            ! Fill up W with the it value so you can see if screening works with different
            ! Wthreshold values
            W(ix,iy,iz,it) = real(it)

            ! Make the data match the radial band number after the averaging takes place.
            !   int(sqrt(DeltaX*DeltaX + DeltaY*DeltaY) / RbandInc) + 1) gives you the radial band number
            !   just multiply it by it so you see increasing slope for successive time steps
            DeltaX = Xcoords(StmIx(it)) - Xcoords(ix)
            DeltaY = Ycoords(StmIy(it)) - Ycoords(iy)
            Avar(ix,iy,iz,it) = real((int(sqrt(DeltaX*DeltaX + DeltaY*DeltaY) / RbandInc) + 1) * it)

          end do 
        end do 
      end do 
    end do 


  end if

  ! Do the averageing. Note that you can pick any index in GdataDescrip below since this data
  ! has been (superficially) checked out to be consistent.
  call AzimuthalAverage(Nx, Ny, Nz, Nt, NumRbands, W, StmIx, StmIy, MinP, Avar, AzAvg, &
          Xcoords, Ycoords, RadialDist, RbandInc, Wthreshold, GdataDescrip(1)%UndefVal)

  call BuildGoutDescrip(NumRbands, 1, Nz, Nt, AzAvg, OfileBase, GdataDescrip(1)%UndefVal, VarToAvg, &
          0.0, RbandInc, 0.0, 1.0, GdataDescrip(1)%Zcoords, GdataDescrip(1)%Tstart, &
          GdataDescrip(1)%Tinc, GoutDescrip, 'azavg')

  call WriteGrads(GoutDescrip, AzAvg)

  ! Clean up
  deallocate (U, V, W, Press, Avar, StmIx, StmIy, MinP, AzAvg, stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory de-allocation failed'
    stop
  end if

  stop
end

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   OfileBase - output GRADS file, base name for two files
!   NumRbands - number of radial bands to split data into
!   VarToAvg - RAMS variable to do the averaging on
!

subroutine GetMyArgs(Infiles, OfileBase, NumRbands, Wthreshold, VarToAvg, DoHorizVel, DoTangential)
  implicit none

  integer :: NumRbands
  real :: Wthreshold
  character (len=*) :: Infiles, OfileBase, VarToAvg
  logical :: DoHorizVel, DoTangential

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 5) then
    write (*,*) 'ERROR: must supply exactly 5 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_file> <num_radial_bands> <w_threshold> <var_to_average>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <num_radial_bands>: number of bands to split data into'
    write (*,*) '        <w_threshold>: only take values of the variable where w is greater than or equal to this number'
    write (*,*) '        <var_to_average>: name of RAMS variable to do the averaging on'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)

  call getarg(3, arg)       !this is a string
  read (arg, '(i)') NumRbands !convert to integer

  call getarg(4, arg)
  read (arg, '(f)') Wthreshold

  call getarg(5, VarToAvg)
  if (VarToAvg == 'speed_t') then
    DoHorizVel = .true.
    DoTangential = .true.
  else
    if (VarToAvg == 'speed_r') then
      DoHorizVel = .true.
      DoTangential = .false.
    else
      DoHorizVel = .false.
      DoTangential = .false.
    end if
  end if

  return
end subroutine

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
          Xcoords, Ycoords, RadialDist, RbandInc, Wthreshold, UndefVal)

  implicit none

  integer :: Nx, Ny, Nz, Nt, NumRbands
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: W, Avar
  real, dimension(1:NumRbands,1:1,1:Nz,1:Nt) :: AzAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nt) :: MinP
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real :: RadialDist, RbandInc, Wthreshold, UndefVal

  integer :: ix, iy, iz, it
  real, dimension(1:NumRbands) :: Rcounts
  integer :: ir, iRband
  real :: DeltaX, DeltaY, Radius

  ! Mask the input data (Avar) with the W data, ie if the corresponding
  ! W data is >= to the Wthreshold use the Avar data in the averaging;
  ! otherwise skip the Avar data
  !
  ! The storm center is taken to be the min pressure location of the
  ! x-y plane on the surface (iz .eq. 1)

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
        AzAvg(ir, 1, Nz, Nt) = 0.0
      end do

      ! Get the averages for this x-y plane
      do iy = 1, Ny
        do ix = 1, Nx
          ! Only keep the points where W meets the threshold
          if ((W(ix,iy,iz,it) .ge. Wthreshold) .and. (Avar(ix,iy,iz,it) .ne. UndefVal)) then
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
        ! if we didn't put anything into an AzAvg slot then it should be zero from
        ! from the intialization above, and we can keep it that way for the result.
        if (Rcounts(ir) .ne. 0.0) then
          AzAvg(ir, 1, iz, it) = AzAvg(ir, 1, iz, it) / Rcounts(ir)
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
