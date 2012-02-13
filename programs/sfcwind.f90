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
  integer, parameter :: MaxFiles=10

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarToInt

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
  real, dimension(:,:,:,:), allocatable :: U, V, SfcWind
  real, dimension(1:1) :: DummyZcoords
  real :: Xstart, Xinc, Ystart, Yinc
  type (GradsVarLocation) :: Uloc, Vloc

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it

  ! Get the command line argument
  call GetMyArgs(Infiles, OfileBase)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Calculating magnitude of surface wind for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) ''


  ! Read the GRADS data description files and collect the information about the data
  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''

  call CheckDataDescrip_SW(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, Uloc, Vloc)

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
  write (*,'(a20,i3,a2,i3,a1)') 'u: (', Uloc%Fnum, ', ', Uloc%Vnum, ')'
  write (*,'(a20,i3,a2,i3,a1)') 'v: (', Vloc%Fnum, ', ', Vloc%Vnum, ')'
  write (*,*) ''

  ! Allocate the data arrays and read in the data from the GRADS data files
  allocate (U(1:Nx,1:Ny,1:Nz,1:Nt), V(1:Nx,1:Ny,1:Nz,1:Nt), &
            SfcWind(1:Nx, 1:Ny, 1:1, 1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(GdataDescrip, 'u', Uloc, U, Nx, Ny, Nz, Nt)
  call ReadGradsData(GdataDescrip, 'v', Vloc, V, Nx, Ny, Nz, Nt)

  ! Generate a surface wind magnitude each time step
  !   mag = sqrt (u**2 + v**2)
  iz = 2 ! use second level to help avoid damping from surface friction
  do it = 1, Nt
    do ix = 1, Nx
      do iy = 1, Ny
        SfcWind(ix,iy,1,it) = sqrt(U(ix,iy,iz,it)**2 + V(ix,iy,iz,it)**2)
      enddo
    enddo
  enddo

  DummyZcoords(1) = 0.0
  Xstart = GdataDescrip(1)%Xcoords(1)
  Xinc = GdataDescrip(1)%Xcoords(2) - Xstart
  Ystart = GdataDescrip(1)%Ycoords(1)
  Yinc = GdataDescrip(1)%Ycoords(2) - Ystart
  call BuildGoutDescrip(Nx, Ny, 1, Nt, SfcWind, OfileBase, GdataDescrip(1)%UndefVal, &
          'sfcwind', Xstart, Xinc, Ystart, Yinc, DummyZcoords, &
          GdataDescrip(1)%Tstart, GdataDescrip(1)%Tinc, GoutDescrip, 'sfc_wind')

  call WriteGrads(GoutDescrip, SfcWind)

  ! Clean up
  deallocate (U, V, SfcWind, stat=Ierror)
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

