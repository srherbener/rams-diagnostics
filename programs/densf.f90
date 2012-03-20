!***************************************************************
! Program to do "density functions"
!
! This program will read in GRADS data from a RAMS simulation, and
! generate a density function (histogram) 
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
!   3. variable used for density function
!   4-9. specs for a cylindrical volume from which to select data
!          min/max radius, min/max phi, min/max z
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
  character (len=LittleString) :: VarName

  type (GradsControlFiles) :: GctlFiles

  integer :: Nbins
  real :: MinR, MaxR
  real :: MinPhi, MaxPhi
  real :: MinZ, MaxZ

  ! Data arrays: need one for dn0 (density) and the var we are doing the
  ! column integration on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  type (GradsVar) :: Gvar, DensFunc
  real, dimension(1:1) :: DummyZcoords

  integer :: i

  integer :: ix, iy, iz, it

  ! Get the command line argument
  call GetMyArgs(Infiles, OfileBase, VarName, Nbins, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ)
  call String2List(Infiles, ':', GctlFiles%Fnames, MaxFiles, GctlFiles%Nfiles, 'input files')

  write (*,*) 'Calculating density function for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, GctlFiles%Nfiles
    write (*,*) '  ', i, ': ', trim(GctlFiles%Fnames(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  RAMS variable: ', trim(VarName)
  write (*,*) ''


  ! Read the GRADS data description files and collect the information about the data
  call ReadGradsCtlFiles(GctlFiles)

  call InitGvarFromGdescrip(GctlFiles, Gvar, VarName)

  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', Gvar%Nx
  write (*,*) '  Number of y (latitude) points:           ', Gvar%Ny
  write (*,*) '  Number of z (vertical level) points:     ', Gvar%Nz
  write (*,*) '  Number of t (time) points:               ', Gvar%Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Gvar%Nx*Gvar%Ny*Gvar%Nz*Gvar%Nt
  write (*,*) ''

  write (*,*) 'Locations of variables in GRADS data (file, var number):'
  write (*,'(a17,a3,a,a2,i3,a1)') trim(VarName), ': (', trim(Gvar%DataFile), ', ', Gvar%Vnum, ')'
  write (*,*) ''

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(Gvar)

  ! WIP: Initialize DensFunc, fill it in, write it out

  stop
end

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   OfileBase - output GRADS file, base name for two files
!   VarName - RAMS variable to do the column integration on
!

subroutine GetMyArgs(Infiles, OfileBase, VarName, Nbins, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ)
  implicit none

  character (len=*) :: Infiles, OfileBase, VarName
  integer :: Nbins, MinZ, MaxZ
  real :: MinR, MaxR, MinPhi, MaxPhi

  integer :: iargc
  character (len=128) :: arg

  logical :: BadArgs

  if (iargc() .ne. 10) then
    write (*,*) 'ERROR: must supply exactly 10 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: densf <in_data_files> <out_data_file> <var> <nbins> <min_r> <max_r> <min_phi> <max_phi> <min_z> <max_z>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <var>: name of GRADS variable to do the density function calculation on'
    write (*,*) '        <nbins>: number of bins to split up the data into'
    write (*,*) ''
    write (*,*) '        The following args are used to select data'
    write (*,*) '          Select data inside cylindrical (r,phi,z) volume:'
    write (*,*) '            <min_r> <max_r> (in km)'
    write (*,*) '            <max_phi> <max_phi> (in radians)'
    write (*,*) '            <min_z> <max_z> (in m)'

    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)

  call getarg(3, VarName)

  call getarg(4, arg)
  read(arg, '(i)') Nbins

  call getarg(5, arg)
  read(arg, '(f)') MinR

  call getarg(6, arg)
  read(arg, '(f)') MaxR

  call getarg(7, arg)
  read(arg, '(f)') MinPhi

  call getarg(8, arg)
  read(arg, '(f)') MaxPhi

  call getarg(9, arg)
  read(arg, '(i)') MinZ

  call getarg(10, arg)
  read(arg, '(i)') MaxZ

  BadArgs = .false.

  if (MinR .ge. MaxR) then
    write (*,*) 'ERROR: <min_r> must be < <max_r>'
    BadArgs = .true.
  endif

  if (MinPhi .ge. MaxPhi) then
    write (*,*) 'ERROR: <min_phi> must be < <max_phi>'
    BadArgs = .true.
  endif

  if (MinZ .ge. MaxZ) then
    write (*,*) 'ERROR: <min_z> must be < <max_z>'
    BadArgs = .true.
  endif

  if (BadArgs) then
    stop
  endif

  return
end subroutine

