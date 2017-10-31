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
  use gdata_utils
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarToSlice
  integer :: Zlevel

  type (GradsControlFiles) :: GctlFiles

  integer :: Nx, Ny, Nz, Nt

  ! Data arrays: need one for the var we are extracting
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  type (GradsVar) :: Gvar, VarSlice
  real, dimension(1:1) :: Zcoords
  real :: Xstart, Xinc, Ystart, Yinc

  integer :: i

  integer :: ix, iy, iz, it

  ! Get the command line argument
  call GetMyArgs(Infiles, OfileBase, VarToSlice, Zlevel)
  call String2List(Infiles, ':', GctlFiles%Fnames, MaxFiles, GctlFiles%Nfiles, 'input files')

  write (*,*) 'Extracting horizontal slice out of GRADS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, GctlFiles%Nfiles
    write (*,*) '  ', i, ': ', trim(GctlFiles%Fnames(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  VarToSlice: ', trim(VarToSlice)
  write (*,*) '  Zlevel: ', Zlevel
  write (*,*) ''


  ! Read the GRADS data description files and collect the information about the data
  call ReadGradsCtlFiles(GctlFiles)

  call InitGvarFromGdescrip(GctlFiles, Gvar, VarToSlice)

  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', Gvar%Nx
  write (*,*) '  Number of y (latitude) points:           ', Gvar%Ny
  write (*,*) '  Number of z (vertical level) points:     ', Gvar%Nz
  write (*,*) '  Number of t (time) points:               ', Gvar%Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Gvar%Nx*Gvar%Ny*Gvar%Nz*Gvar%Nt
  write (*,*) ''

  write (*,*) 'Locations of variables in GRADS data (file, var number):'
  write (*,'(a17,a3,a,a2,i3,a1)') trim(VarToSlice), ': (', trim(Gvar%DataFile), ', ', Gvar%Vnum, ')'
  write (*,*) ''

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(Gvar)

  ! Initialize the output var
  Zcoords(1) = Gvar%Zcoords(Zlevel)
  call InitGradsVar(VarSlice, VarToSlice, Gvar%Nx, Gvar%Ny, 1, Gvar%Nt, &
                    Gvar%Xstart, Gvar%Xinc, Gvar%Ystart, Gvar%Yinc, Zcoords, Gvar%Tstart, Gvar%Tinc, &
                    Gvar%UndefVal, '<NONE>', 0, 0)

  ! Cut out the horiz slice
  do it = 1, Gvar%Nt
    do ix = 1, Gvar%Nx
      do iy = 1, Gvar%Ny
        VarSlice%Vdata(it,1,ix,iy) = Gvar%Vdata(it,Zlevel,ix,iy)
      enddo
    enddo
  enddo

  call WriteGrads(VarSlice, OfileBase, 'hslice')

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
    write (*,*) 'USAGE: hslice <in_data_files> <out_data_file>'
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

