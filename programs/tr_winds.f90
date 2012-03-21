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

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBaseTan, OfileBaseRad
  character (len=LittleString) :: VarName

  type (GradsControlFiles) :: GctlFiles

  ! Data arrays: need one for w (vertical velocity), press (pressure)
  ! and the var we are doing the averaging on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  type (GradsVar) :: U, V, Press, Rad, Tang
  integer, dimension(:), allocatable :: StmIx, StmIy
  real, dimension(:), allocatable :: MinP, Xcoords, Ycoords

  integer :: i

  integer :: ix, iy, iz, it

  real :: Xstart, Xinc, Ystart, Yinc

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBaseTan, OfileBaseRad)
  call String2List(Infiles, ':', GctlFiles%Fnames, MaxFiles, GctlFiles%Nfiles, 'input files')

  write (*,*) 'Calculating storm tangential and radial winds for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, GctlFiles%Nfiles
    write (*,*) '  ', i, ': ', trim(GctlFiles%Fnames(i))
  end do
  write (*,*) '  Output file base name, (tangential winds):  ', trim(OfileBaseTan)
  write (*,*) '  Output file base name, (radial winds):  ', trim(OfileBaseRad)
  write (*,*) ''

  ! Read the GRADS data description files and collect the information about the data
  call ReadGradsCtlFiles(GctlFiles)

  call InitGvarFromGdescrip(GctlFiles, U, 'u')
  call InitGvarFromGdescrip(GctlFiles, V, 'v')
  call InitGvarFromGdescrip(GctlFiles, Press, 'press')

  if (.not. (GvarDimsMatch(U, V, .false.) .and. GvarDimsMatch(U, Press, .false.))) then
    write (*,*) 'ERROR: dimensions of u, v and press do not match'
    stop
  endif

  ! At this point the dimensions of U, V, Press match so use U as a representative
  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', U%Nx
  write (*,*) '  Number of y (latitude) points:           ', U%Ny
  write (*,*) '  Number of z (vertical level) points:     ', U%Nz
  write (*,*) '  Number of t (time) points:               ', U%Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', U%Nx*U%Ny*U%Nz*U%Nt
  write (*,*) ''

  write (*,*) 'Locations of variables in GRADS data (file, var number):'
  write (*,'(a20,a,a2,i3,a1)') 'speed - u: (', trim(U%DataFile), ', ', U%Vnum, ')'
  write (*,'(a20,a,a2,i3,a1)') 'speed - v: (', trim(V%DataFile), ', ', V%Vnum, ')'
  write (*,'(a20,a,a2,i3,a1)') 'press: (', trim(Press%DataFile), ', ', Press%Vnum, ')'
  write (*,*) ''

  ! Convert the GRADS grid coordinates from longitude, latitude to flat plane (x and y).
  call ConvertGridCoords(U, Xcoords, Ycoords)

  write (*,*) 'Horzontal Grid Coordinate Info:'
  write (*,*) '  X Range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', U%Xcoords(1), U%Xcoords(U%Nx), Xcoords(1), Xcoords(U%Nx)
  write (*,*) '  Y Range: '
  write (*,*) '    ', U%Ycoords(1), U%Ycoords(U%Ny), Ycoords(1), Ycoords(U%Ny)
  write (*,*) ''

  ! initialize the output variables
  call InitGradsVar(Tang, 'wind_speed' , U%Nx, U%Ny, U%Nz, U%Nt, &
                    U%Xstart, U%Xinc, U%Ystart, U%Yinc, U%Zcoords, U%Tstart, U%Tinc, &
                    U%UndefVal, '<NONE>', 0, 0)
  call InitGradsVar(Rad, 'wind_speed' , U%Nx, U%Ny, U%Nz, U%Nt, &
                    U%Xstart, U%Xinc, U%Ystart, U%Yinc, U%Zcoords, U%Tstart, U%Tinc, &
                    U%UndefVal, '<NONE>', 0, 0)


  ! Read in the data for the vars using the description and location information
  call ReadGradsData(U)
  call ReadGradsData(V)
  call ReadGradsData(Press)

  call RecordStormCenter(Press, StmIx, StmIy, MinP)
  call ConvertHorizVelocity(U, V, StmIx, StmIy, Xcoords, Ycoords, Tang, .true.)
  call ConvertHorizVelocity(U, V, StmIx, StmIy, Xcoords, Ycoords, Rad, .false.)

  call WriteGrads(Tang, OfileBaseTan, 'tan_winds')
  call WriteGrads(Rad,  OfileBaseRad, 'rad_winds')

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
