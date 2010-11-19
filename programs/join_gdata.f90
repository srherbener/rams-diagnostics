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

include 'gDataTypes.h'

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
  integer :: Nfiles, MinLevel, MaxLevel

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars
  type (GradsOutDescription) :: GoutDescrip

  ! Data arrays
  ! OutData dims: x, y, z, t
  ! VarData dims: x, y, z, t
  real, dimension(:,:,:,:), allocatable :: OutData, VarData
  type (GradsVarLocation), dimension(1:MaxFiles) :: VarLocs

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it
  integer :: outTime
  
  real :: Xstart, Xinc, Ystart, Yinc

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBase, VarName)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Joining GRADS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  GRADS variable:', trim(VarName)
  write (*,*) ''

  ! Read the GRADS data description files and collect the information about the data
  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''
  flush(6)

  ! Check the data description for consistency and locate the variables in the GRADS control files
  call CheckDataDescripJoin(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, VarLocs, VarName)
  do i = 1, Nfiles
    write (*,*) 'Variable location for file: ', trim(GradsCtlFiles(i))
    write (*,*) '  Fnum: ', VarLocs(i)%Fnum
    write (*,*) '  Vnum: ', VarLocs(i)%Vnum
  end do
  write (*,*) 'Variable dimensions after joining files:'
  write (*,*) '  Nx: ', Nx
  write (*,*) '  Ny: ', Ny
  write (*,*) '  Nz: ', Nz
  write (*,*) '  Nt: ', Nt
  write (*,*) ''
  flush(6)

  !Read in the data and copy into the output array as you go
  allocate (OutData(1:Nx,1:Ny,1:Nz,1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  outTime = 0
  do i = 1, Nfiles
    allocate(VarData(1:GdataDescrip(i)%nx,1:GdataDescrip(i)%ny,1:GdataDescrip(i)%nz, &
             1:GdataDescrip(i)%nt), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: VarData array memory allocation failed'
      stop
    end if

    call ReadGradsData(GdataDescrip, VarName, VarLocs(i), VarData, GdataDescrip(i)%nx, &
         GdataDescrip(i)%ny, GdataDescrip(i)%nz, GdataDescrip(i)%nt)
    
    do it = 1, GdataDescrip(i)%nt
      outTime = outTime + 1
      do iz = 1, GdataDescrip(i)%nz
        do iy = 1, GdataDescrip(i)%ny
          do ix = 1, GdataDescrip(i)%nx
            OutData(ix,iy,iz,outTime) = VarData(ix,iy,iz,it)
          end do
        end do
      end do
    end do

    deallocate (VarData, stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: VarData array memory de-allocation failed'
      stop
    end if
  end do

  !Write out the joined data

  Xstart = GdataDescrip(1)%Xcoords(1)
  if (Nx .gt. 1) then
    Xinc = GdataDescrip(1)%Xcoords(2)-GdataDescrip(1)%Xcoords(1)
  else
    Xinc = 1.0
  end if
  Ystart = GdataDescrip(1)%Ycoords(1)
  if (Ny .gt. 1) then
    Yinc = GdataDescrip(1)%Ycoords(2)-GdataDescrip(1)%Ycoords(1)
  else
    Yinc = 1.0
  end if

  call BuildGoutDescrip(Nx, Ny, Nz, Nt, OutData, OfileBase, GdataDescrip(1)%UndefVal, VarName, &
          Xstart, Xinc, Ystart, Yinc, GdataDescrip(1)%Zcoords, GdataDescrip(1)%Tstart, &
          GdataDescrip(1)%Tinc, GoutDescrip, 'join')

  call WriteGrads(GoutDescrip, OutData)

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
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_file> <averaging_function>'
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

!*************************************************************************
! CheckDataDescripJoin()
!
! This routine will check the consistency of the data description obtained
! from the GRADS control files. The number of x, y, z, t points should match
! in all files. Also, the total number of variables are calculated and returned
! in Nvars.
!
! This routein will check one variable at a time, so it needs to be called
! multiple times for multiple variables.
!

subroutine CheckDataDescripJoin(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, VarLocs, VarName)
  use GfileTypes
  implicit none

  type (GradsDataDescription), dimension(*) :: GdataDescrip
  integer Nx, Ny, Nz, Nt, Nfiles
  integer Nvars
  type (GradsVarLocation), dimension(*) :: VarLocs
  character (len=*) :: VarName
  logical :: DoHorizVel

  ! i -> file number, j -> var number
  integer i, j
  logical BadData, GotNz

  BadData = .false.
  GotNz = .false.
  ! Make Nt be the total time steps (after joining the data)
  do i = 1, Nfiles
    VarLocs(i)%Fnum = 0
    VarLocs(i)%Vnum = 0

    ! Make sure nx, ny and for all vars (2D and 3D) match
    if (i .eq. 1) then
      Nx = GdataDescrip(i)%nx
      Ny = GdataDescrip(i)%ny
      Nt = GdataDescrip(i)%nt
   
      Nvars = GdataDescrip(i)%nvars
    else
      if (GdataDescrip(i)%nx .ne. Nx) then
        write (0,*) 'ERROR: number of x points in GRADS control files do not match'
        BadData = .true.
      end if
      if (GdataDescrip(i)%ny .ne. Ny) then
        write (0,*) 'ERROR: number of y points in GRADS control files do not match'
        BadData = .true.
      end if

      Nt = Nt + GdataDescrip(i)%nt
      Nvars = Nvars + GdataDescrip(i)%nvars
    end if

    ! allow for a 2D var (nz eq 1), but make sure all 3D vars have matching nz
    if (GotNz) then
      if ((GdataDescrip(i)%nz .ne. 1) .and. (GdataDescrip(i)%nz .ne. Nz)) then
        write (0,*) 'ERROR: number of z points in GRADS control files do not match'
        BadData = .true.
      end if
    else
      ! grab the first nz that is not equal to 1 (first 3D var)
      if ((GdataDescrip(i)%nz .ne. 1) .and. (.not. GotNz)) then
        Nz = GdataDescrip(i)%nz
        GotNz = .true.
      end if
    end if

    ! Record the location of the given variable in each of the files
    do j =1, GdataDescrip(i)%nvars
      if (GdataDescrip(i)%VarNames(j) .eq. VarName) then
        VarLocs(i)%Fnum = i
        VarLocs(i)%Vnum = j
      end if
    end do
  end do

  ! If the entire list was 2D variables, then the code above will have not recorded an
  ! Nz value (GotNz is still .false.). Set Nz to '1' in this case
  if (.not. GotNz) then
    Nz = 1
  end if

  ! Check to see if you got the variable
  do i = 1, Nfiles
    if (VarLocs(i)%Fnum .eq. 0) then
      write (0,*) 'ERROR: cannot find grid var "', trim(VarName), '" in the GRADS data file number:', i
      BadData = .true.
    end if
  end do

  if (BadData) then
    stop
  end if
end subroutine
