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

  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarName

  type (GradsControlFiles) :: GctlFiles
  integer :: MinLevel, MaxLevel

  ! Data arrays
  ! OutVar dims: t, z, x, y
  ! InVar dims: t, z, x, y
  type (GradsVar) :: OutVar, InVar

  integer :: i

  integer :: ix, iy, iz, it
  integer :: OutTime
  integer :: OutNt
  
  ! Get the command line arguments
  call GetMyArgs(GctlFiles%Fnames(1), OfileBase, VarName)
  GctlFiles%Nfiles = 1

  write (*,*) 'Doubling points in the GRADS data:'
  write (*,*) '  GRADS input control file: ', trim(GctlFiles%Fnames(1))
  write (*,*) '  GRADS variable: ', trim(VarName)
  write (*,*) ''
  flush(6)

  ! Read the GRADS data description file and collect the information about the data
  call ReadGradsCtlFiles(GctlFiles)

  call InitGvarFromGdescrip(GctlFiles, InVar, VarName)
  OutNt = (2*InVar%Nt) - 1  ! will be insterting Nt-1 new points inbetween the existing Nt points

  write (*,*) 'Locations of variables in GRADS data (file, var number):'
  write (*,'(a17,a3,a,a2,i3,a1)') trim(VarName), ': (', trim(InVar%DataFile), ', ', InVar%Vnum, ')'
  write (*,*) ''
  flush(6)

  write (*,*) 'Variable dimensions (after doubling time points):'
  write (*,*) '  Nx: ', InVar%Nx
  write (*,*) '  Ny: ', InVar%Ny
  write (*,*) '  Nz: ', InVar%Nz
  write (*,*) '  Nt: ', InVar%Nt
  write (*,*) '  Output Nt: ', OutNt
  write (*,*) ''
  flush(6)

  call InitGradsVar(OutVar, VarName, InVar%Nx, InVar%Ny, InVar%Nz, OutNt, &
                    InVar%Xstart, InVar%Xinc, InVar%Ystart, InVar%Yinc, InVar%Zcoords, &
                    InVar%Tstart, InVar%Tinc, InVar%UndefVal, '<NONE>', 0, 0)

  call ReadGradsData(InVar)

  ! Not quite doubling points, instead we are filling in new time points in
  ! between all the existing points in the input data (so if there are n
  ! input points, then we end up with 2n-1 output points).
  !
  ! The pattern we want for creating the output data is:
  !
  !     input time step             output time step
  !            1  copy of 1   ------->    1
  !               avg of 1,2  ------->    2
  !            2  copy of 2   ------->    3
  !               avg of 2,3  ------->    4
  !            3  copy of 3   ------->    5
  !            .     .           .        .
  !            .     .           .        .
  !            .     .           .        .
  !
  ! Algorithm:
  !   1. if on first input time step:
  !        - set output time step to 1
  !        - copy data at first intput time step to first output time step
  !   2. if not on first time step
  !        - increment output time step by 2
  !        - place average of data at the current and previous input time steps
  !          into the previous output time step
  !        - place data at current intput time step to current output time step

  do it = 1, InVar%Nt
    if (it .eq. 1) then
      OutTime = 1
    else
      OutTime = OutTime + 2
    end if
    do iz = 1, InVar%Nz
      do iy = 1, InVar%Ny
        do ix = 1, InVar%Nx
          if (it .eq. 1) then
            OutVar%Vdata(OutTime,iz,ix,iy) = InVar%Vdata(it,iz,ix,iy)
          else
            OutVar%Vdata((OutTime-1),iz,ix,iy) = (InVar%Vdata(it,iz,ix,iy) + InVar%Vdata((it-1),iz,ix,iy)) / 2.0
            OutVar%Vdata(OutTime,iz,ix,iy) = InVar%Vdata(it,iz,ix,iy)
          end if
        end do
      end do
    end do
  end do

  call WriteGrads(OutVar, VarName, 'dpts')

  stop
end

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infile - input GRADS file
!   OfileBase - output GRADS file, base name for two files
!

subroutine GetMyArgs(Infile, OfileBase, VarName)
  implicit none

  integer, parameter :: MAX_ITEMS = 5

  character (len=*) :: Infile, OfileBase, VarName

  integer :: iargc
  character (len=128) :: arg
  character (len=128), dimension(1:MAX_ITEMS) :: ArgList
  integer :: Nitems

  logical :: BadArgs

  if (iargc() .ne. 3) then
    write (*,*) 'ERROR: must supply exactly 3 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_data_file> <out_data_file> <averaging_function>'
    write (*,*) '        <in_data_file>: GRADS format, control file'
    write (*,*) '        <out_data_file>: GRADS format, this program will tag on .ctl, .gra suffixes'
    write (*,*) '        <variable_name>: GRADS variable name, this program will look for this variable in <in_data_file>'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infile)
  call getarg(2, OfileBase)
  call getarg(3, VarName)

  return
end subroutine

