!***************************************************************
! Program to extract a horizontal slice out of GRADS data for
! a variable
!
! This program will read in GRADS data from a RAMS simulation, and
! generate the column integrated value of a given variable.
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
!   3. variable
!
! Output
!   The output will be two files which make a GRADS data set. One
!   file is the control file and the other is the binary data.
!

program main
  use gdata_mod
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarToSlice
  integer :: Zlevel

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles
  integer :: Nfiles

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars
  type (GradsOutDescription) :: GoutDescrip

  ! Data arrays: need one for the var we are extracting
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: Var, VarSlice
  real, dimension(1:1) :: Zcoords
  real :: Xstart, Xinc, Ystart, Yinc
  type (GradsVarLocation) :: VarLoc

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it

  ! Get the command line argument
  call GetMyArgs(Infiles, OfileBase, VarToSlice, Zlevel)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Extracting horizontal slice out of GRADS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  VarToSlice: ', trim(VarToSlice)
  write (*,*) '  Zlevel: ', Zlevel
  write (*,*) ''


  ! Read the GRADS data description files and collect the information about the data
  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''

  call CheckDataDescrip_VS(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, VarLoc, VarToSlice)

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
  write (*,'(a17,a3,i3,a2,i3,a1)') trim(VarToSlice), ': (', VarLoc%Fnum, ', ', VarLoc%Vnum, ')'
  write (*,*) ''

  ! Allocate the data arrays and read in the data from the GRADS data files
  allocate (Var(1:Nx,1:Ny,1:Nz,1:Nt), VarSlice(1:Nx, 1:Ny, 1:1, 1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(GdataDescrip, VarToSlice, VarLoc, Var, Nx, Ny, Nz, Nt)

  ! Cut out the horiz slice
  do it = 1, Nt
    do ix = 1, Nx
      do iy = 1, Ny
        VarSlice(ix,iy,1,it) = Var(ix,iy,Zlevel,it)
      enddo
    enddo
  enddo

  Zcoords(1) = GdataDescrip(1)%Zcoords(Zlevel)
  Xstart = GdataDescrip(1)%Xcoords(1)
  Xinc = GdataDescrip(1)%Xcoords(2) - Xstart
  Ystart = GdataDescrip(1)%Ycoords(1)
  Yinc = GdataDescrip(1)%Ycoords(2) - Ystart
  call BuildGoutDescrip(Nx, Ny, 1, Nt, VarSlice, OfileBase, GdataDescrip(1)%UndefVal, &
          VarToSlice, Xstart, Xinc, Ystart, Yinc, Zcoords, &
          GdataDescrip(1)%Tstart, GdataDescrip(1)%Tinc, GoutDescrip, 'hslice')

  call WriteGrads(GoutDescrip, VarSlice)

  ! Clean up
  deallocate (Var, VarSlice, stat=Ierror)
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

subroutine GetMyArgs(Infiles, OfileBase, VarToSlice, Zlevel)
  implicit none

  character (len=*) :: Infiles, OfileBase, VarToSlice
  integer :: Zlevel

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 4) then
    write (*,*) 'ERROR: must supply exactly 4 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: sfcwind <in_data_files> <out_data_file>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <variable>: name of variable'
    write (*,*) '        <z_level>: z level for horizontal slice in grid coordinate'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)
  call getarg(3, VarToSlice)

  call getarg(4, arg)
  read (arg, '(i)') Zlevel

  return
end subroutine

