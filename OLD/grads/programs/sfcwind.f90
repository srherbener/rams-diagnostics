!***************************************************************
! Program to calculate magnitude of the surface wind
!
! This program will read in GRADS data from a RAMS simulation, and
! generate the column integrated value of a given variable.
! Friction at the real surface in the simulation will tend to slow
! down the wind speed which tend to make all the storms look
! similar (dampens out the structure of the wind speed vs radius). To
! help mitigate that effect, take the wind speeds from the second z level
! instead of the first.
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
!   3. variable used for column integration
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
  character (len=LittleString) :: VarToInt

  type (GradsControlFiles):: GctlFiles

  integer :: Nx, Ny, Nz, Nt

  ! Data arrays: need one for dn0 (density) and the var we are doing the
  ! column integration on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  type (GradsVar) :: U, V, SfcWind
  real, dimension(1:1) :: DummyZcoords
  real :: Xstart, Xinc, Ystart, Yinc

  integer :: i

  integer :: ix, iy, iz, it

  ! Get the command line argument
  call GetMyArgs(Infiles, OfileBase)
  call String2List(Infiles, ':', GctlFiles%Fnames, MaxFiles, GctlFiles%Nfiles, 'input files')

  write (*,*) 'Calculating magnitude of surface wind for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, GctlFiles%Nfiles
    write (*,*) '  ', i, ': ', trim(GctlFiles%Fnames(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) ''


  ! Read the GRADS data description files and collect the information about the data
  call ReadGradsCtlFiles(GctlFiles)

  call InitGvarFromGdescrip(GctlFiles, U, 'u')
  call InitGvarFromGdescrip(GctlFiles, V, 'v')

  if (.not. (GvarDimsMatch(U, V, .false.))) then
    write (*,*) 'ERROR: dimensions of u, and v do not match'
    stop
  endif

  ! U and V have matching dimensions if we got to here, so use U as a representative
  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', U%Nx
  write (*,*) '  Number of y (latitude) points:           ', U%Ny
  write (*,*) '  Number of z (vertical level) points:     ', U%Nz
  write (*,*) '  Number of t (time) points:               ', U%Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', U%Nx*U%Ny*U%Nz*U%Nt
  write (*,*) ''

  write (*,*) 'Locations of variables in GRADS data (file, var number):'
  write (*,'(a20,a,a2,i3,a1)') 'u: (', trim(U%DataFile), ', ', U%Vnum, ')'
  write (*,'(a20,a,a2,i3,a1)') 'v: (', trim(V%DataFile), ', ', V%Vnum, ')'
  write (*,*) ''

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(U)
  call ReadGradsData(V)

  DummyZcoords(1) = 0.0
  call InitGradsVar(SfcWind, 'sfcwind', U%Nx, U%Ny, 1, U%Nt, &
                    U%Xstart, U%Xinc, U%Ystart, U%Yinc, DummyZcoords, U%Tstart, U%Tinc, &
                    U%UndefVal, '<NONE>', 0, 0)


  ! Generate a surface wind magnitude each time step
  !   mag = sqrt (u**2 + v**2)
  iz = 2 ! use second level to help avoid damping from surface friction
  do it = 1, U%Nt
    do ix = 1, U%Nx
      do iy = 1, U%Ny
        SfcWind%Vdata(it,1,ix,iy) = sqrt(U%Vdata(it,iz,ix,iy)**2 + V%Vdata(it,iz,ix,iy)**2)
      enddo
    enddo
  enddo

  call WriteGrads(SfcWind, OfileBase, 'sfc_wind')

  stop
end

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   OfileBase - output GRADS file, base name for two files
!   VarToInt - RAMS variable to do the column integration on
!

subroutine GetMyArgs(Infiles, OfileBase)
  implicit none

  character (len=*) :: Infiles, OfileBase

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 2) then
    write (*,*) 'ERROR: must supply exactly 2 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: sfcwind <in_data_files> <out_data_file>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)

  return
end subroutine

