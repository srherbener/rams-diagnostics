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
  use GfileTypes
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarName

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles
  integer :: Nfiles

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars
  type (GradsOutDescription) :: GoutDescrip

  ! Data arrays: need one for dn0 (density) and the var we are doing the
  ! column integration on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: Gvar, DensFunc
  real, dimension(1:1) :: DummyZcoords
  real :: Xstart, Xinc, Ystart, Yinc
  type (GradsVarLocation) :: VarLoc

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it

  ! Get the command line argument
  call GetMyArgs(Infiles, OfileBase, VarName)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Calculating density function for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  RAMS variable: ', trim(VarName)
  write (*,*) ''


  ! Read the GRADS data description files and collect the information about the data
  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''

  call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, VarLoc, VarName)

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
  write (*,'(a17,a3,i3,a2,i3,a1)') trim(VarName), ': (', VarLoc%Fnum, ', ', VarLoc%Vnum, ')'
  write (*,*) ''

  ! Allocate the data arrays and read in the data from the GRADS data files
  allocate (Gvar(1:Nx,1:Ny,1:Nz,1:Nt), DensFunc(1:Nx, 1:Ny, 1:1, 1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(GdataDescrip, VarName, VarLoc, Gvar, Nx, Ny, Nz, Nt)


!  DummyZcoords(1) = 0.0
!  Xstart = GdataDescrip(1)%Xcoords(1)
!  Xinc = GdataDescrip(1)%Xcoords(2) - Xstart
!  Ystart = GdataDescrip(1)%Ycoords(1)
!  Yinc = GdataDescrip(1)%Ycoords(2) - Ystart
!  call BuildGoutDescrip(Nx, Ny, 1, Nt, DensFunc, OfileBase, GdataDescrip(1)%UndefVal, &
!          VarName, Xstart, Xinc, Ystart, Yinc, DummyZcoords, &
!          GdataDescrip(1)%Tstart, GdataDescrip(1)%Tinc, GoutDescrip, 'colint')
!
!  call WriteGrads(GoutDescrip, DensFunc)

  ! Clean up
  deallocate (Gvar, DensFunc, stat=Ierror)
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
!   VarName - RAMS variable to do the column integration on
!

subroutine GetMyArgs(Infiles, OfileBase, VarName, Nbins, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ)
  implicit none

  character (len=*) :: Infiles, OfileBase, VarName
  integer :: Nbins, MinZ, MaxZ
  real :: MinR, MaxR, MinPhi, MaxPhi

  integer :: iargc
  character (len=128) :: arg

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

  return
end subroutine

