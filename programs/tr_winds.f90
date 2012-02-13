!***************************************************************
! Program to convert u,v wind values to tangential,radial wind
! values relative to the storm center.
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

program main
  use gdata_utils
  use azavg_utils
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBaseTan, OfileBaseRad
  character (len=LittleString) :: VarName

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
  real, dimension(:,:,:,:), allocatable :: U, V, Press, Rad, Tang
  integer, dimension(:), allocatable :: StmIx, StmIy
  real, dimension(:), allocatable :: MinP, Xcoords, Ycoords
  type (GradsVarLocation) :: Ploc, Uloc, Vloc

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it

  real :: Xstart, Xinc, Ystart, Yinc

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBaseTan, OfileBaseRad)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Calculating storm tangential and radial winds for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name, (tangential winds):  ', trim(OfileBaseTan)
  write (*,*) '  Output file base name, (radial winds):  ', trim(OfileBaseRad)
  write (*,*) ''

  ! Read the GRADS data description files and collect the information about the data
  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''

  call CheckDataDescrip_TR(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, Ploc, Uloc, Vloc)

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
  write (*,'(a20,i3,a2,i3,a1)') 'pressure: (', Ploc%Fnum, ', ', Ploc%Vnum, ')'
  write (*,'(a20,i3,a2,i3,a1)') 'speed - u: (', Uloc%Fnum, ', ', Uloc%Vnum, ')'
  write (*,'(a20,i3,a2,i3,a1)') 'speed - v: (', Vloc%Fnum, ', ', Vloc%Vnum, ')'
  write (*,*) ''

  ! Allocate the data arrays and read in the data from the GRADS data files
  allocate (U(1:Nx,1:Ny,1:Nz,1:Nt), V(1:Nx,1:Ny,1:Nz,1:Nt), Press(1:Nx,1:Ny,1:Nz,1:Nt), &
            Rad(1:Nx,1:Ny,1:Nz,1:Nt), Tang(1:Nx,1:Ny,1:Nz,1:Nt), &
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

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(GdataDescrip, 'press', Ploc, Press, Nx, Ny, Nz, Nt)
  call RecordStormCenter(Nx, Ny, Nz, Nt, Press, StmIx, StmIy, MinP)
  call ReadGradsData(GdataDescrip, 'u', Uloc, U, Nx, Ny, Nz, Nt)
  call ReadGradsData(GdataDescrip, 'v', Vloc, V, Nx, Ny, Nz, Nt)
  call ConvertHorizVelocity(Nx, Ny, Nz, Nt, U, V, StmIx, StmIy, &
                  Xcoords, Ycoords, Tang, .true.)
  call ConvertHorizVelocity(Nx, Ny, Nz, Nt, U, V, StmIx, StmIy, &
                  Xcoords, Ycoords, Rad, .false.)

  ! Want lat, lon values so grab them from the description arrays
  Xstart = GdataDescrip(1)%Xcoords(1)
  Xinc = GdataDescrip(1)%Xcoords(2) - GdataDescrip(1)%Xcoords(1)
  Ystart = GdataDescrip(1)%Ycoords(1)
  Yinc = GdataDescrip(1)%Ycoords(2) - GdataDescrip(1)%Ycoords(1)
  VarName = 'wind_speed'

  ! Write the tangential wind values
  call BuildGoutDescrip(Nx, Ny, Nz, Nt, Tang, OfileBaseTan, GdataDescrip(1)%UndefVal, VarName, &
          Xstart, Xinc, Ystart, Yinc, GdataDescrip(1)%Zcoords, GdataDescrip(1)%Tstart, &
          GdataDescrip(1)%Tinc, GoutDescrip, 'tan_winds')
  call WriteGrads(GoutDescrip, Tang)

  ! Write the radial wind values
  call BuildGoutDescrip(Nx, Ny, Nz, Nt, Rad, OfileBaseRad, GdataDescrip(1)%UndefVal, VarName, &
          Xstart, Xinc, Ystart, Yinc, GdataDescrip(1)%Zcoords, GdataDescrip(1)%Tstart, &
          GdataDescrip(1)%Tinc, GoutDescrip, 'rad_winds')
  call WriteGrads(GoutDescrip, Rad)

  ! Clean up
  deallocate (U, V, Press, Tang, Rad, StmIx, StmIy, MinP, stat=Ierror)
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
!

subroutine GetMyArgs(Infiles, OfileBaseTan, OfileBaseRad)
  implicit none

  character (len=*) :: Infiles, OfileBaseTan, OfileBaseRad

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 3) then
    write (*,*) 'ERROR: must supply exactly 3 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: tr_winds <in_data_files> <out_file_tan> <out_file_rad>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_file_tan>: tangential wind output file, GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <out_file_rad>: radial wind output file, GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBaseTan)
  call getarg(3, OfileBaseRad)

  return
end subroutine
