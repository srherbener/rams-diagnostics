!***************************************************************
! Program to take a list of files and join the GRADS data within
! them together.
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
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

  type (GradsControlFiles), dimension(1:MaxFiles) :: GctlFiles
  integer :: MinLevel, MaxLevel
  character (len=MediumString), dimension(1:MaxFiles) :: InFileList
  integer :: Nfiles

  integer :: Nx, Ny, Nz, Nt

  ! Data arrays
  ! OutVar dims: t, z, x, y
  ! Invars, array where each element has a data array with dims: t, z, x, y
  type (GradsVar), dimension(1:MaxFiles) :: InVars
  type (GradsVar) :: OutVar

  integer :: i

  integer :: ix, iy, iz, it
  integer :: OutTime
  
  real :: Xstart, Xinc, Ystart, Yinc
  logical :: BadData

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBase, VarName)
  call String2List(Infiles, ':', InFileList, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Joining GRADS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(InFileList(i))
  end do
  write (*,*) '  GRADS variable:', trim(VarName)
  write (*,*) ''

  ! For this function we need to locate the variable in every input file since
  ! we want to join these together. Use arrayed versions of the GradsControlFiles type
  ! and GradVar type. Walk though each individual input file given and collect the variable
  ! information from just that particular file. Store each result in the elements of the
  ! arrayed versions of GradsControlFiles and GradsVar types.
  do i = 1, Nfiles
    GctlFiles(i)%Fnames(1) = InFileList(i)
    GctlFiles(i)%Nfiles = 1
    call ReadGradsCtlFiles(GctlFiles(i))
    call InitGvarFromGdescrip(GctlFiles(i), InVars(i), VarName)

    ! Do some dimension checking
    if (i .eq. 1) then
      Nx = InVars(i)%Nx
      Ny = InVars(i)%Ny
      Nz = InVars(i)%Nz
      Nt = InVars(i)%Nt
    else
      BadData = .false.
      if (Nx .ne. InVars(i)%Nx) then
        write (*,*) 'ERROR: number of x points in GRADS control files do not match'
        BadData = .true.
      endif
      if (Ny .ne. InVars(i)%Ny) then
        write (*,*) 'ERROR: number of y points in GRADS control files do not match'
        BadData = .true.
      endif
      if (Nz .ne. InVars(i)%Nz) then
        write (*,*) 'ERROR: number of z points in GRADS control files do not match'
        BadData = .true.
      endif
      if (BadData) then
        write (*,*) 'ERROR:  Control file: ', trim(InFileList(i))
        stop
      endif

      ! Accumulate the number of time steps so we know how many time steps will exist
      ! in the output data
      Nt = Nt + InVars(i)%Nt
    endif
  enddo

  write (*,*) 'Locations of variables in GRADS data (file, var number):'
  do i = 1, Nfiles
    write (*,'(a17,a3,a,a2,i3,a1)') trim(VarName), ': (', trim(InVars(i)%DataFile), ', ', InVars(i)%Vnum, ')'
  end do
  write (*,*) 'Variable dimensions after joining files:'
  write (*,*) '  Nx: ', Nx
  write (*,*) '  Ny: ', Ny
  write (*,*) '  Nz: ', Nz
  write (*,*) '  Nt: ', Nt
  write (*,*) ''
  flush(6)

  ! Initialize the output variable, go through the list of input variables and read
  ! in their data and append it to the output variable
  ! We have verified that the x, y, z dimensions of the input data in all GRADS input
  ! files match so it's okay use InVars(1) as a representative.
  call InitGradsVar(OutVar, VarName, Nx, Ny, Nz, Nt, &
                    InVars(1)%Xstart, InVars(1)%Xinc, InVars(1)%Ystart, InVars(1)%Yinc, &
                    InVars(1)%Zcoords, InVars(1)%Tstart, InVars(1)%Tinc, &
                    InVars(1)%UndefVal, '<NONE>', 0, 0)

  OutTime = 0
  do i = 1, Nfiles
    call ReadGradsData(InVars(i))
    
    do it = 1, InVars(i)%Nt
      OutTime = OutTime + 1
      do iz = 1, InVars(i)%Nz
        do ix = 1, InVars(i)%Nx
          do iy = 1, InVars(i)%Ny
            OutVar%Vdata(OutTime,iz,ix,iy) = InVars(i)%Vdata(it,iz,ix,iy)
          end do
        end do
      end do
    end do
  end do

  !Write out the joined data
  call WriteGrads(OutVar, OfileBase, 'join')

  stop
end

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   OfileBase - output GRADS file, base name for two files
!

subroutine GetMyArgs(Infiles, OfileBase, VarName)
  implicit none

  integer, parameter :: MAX_ITEMS = 5

  character (len=*) :: Infiles, OfileBase, VarName

  integer :: iargc
  character (len=128) :: arg
  character (len=128), dimension(1:MAX_ITEMS) :: ArgList
  integer :: Nitems

  logical :: BadArgs

  if (iargc() .ne. 3) then
    write (*,*) 'ERROR: must supply exactly 3 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_file> <variable_name>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this program will tag on .ctl, .gra suffixes'
    write (*,*) '        <variable_name>: GRADS variable name, this program will look for this variable in all of the <in_data_files>'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)
  call getarg(3, VarName)

  return
end subroutine
