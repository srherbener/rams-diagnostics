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
  use GfileTypes
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarName

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles
  integer :: MinLevel, MaxLevel

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
  call GetMyArgs(GradsCtlFiles(1), OfileBase, VarName)

  write (*,*) 'Joining GRADS data:'
  write (*,*) '  GRADS input control file:', trim(GradsCtlFiles(1))
  write (*,*) '  GRADS variable:', trim(VarName)
  write (*,*) ''

  ! Read the GRADS data description file and collect the information about the data
  write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(1))
  call ReadGradsCtlFile(GradsCtlFiles(1), GdataDescrip(1))
  write (*,*) ''
  flush(6)

  ! Check the data description for consistency and locate the variables in the GRADS control files
  call CheckDataDescripDouble(GdataDescrip, 1, Nx, Ny, Nz, Nt, Nvars, VarLocs, VarName)
  write (*,*) 'Variable location for file: ', trim(GradsCtlFiles(1))
  write (*,*) '  Fnum: ', VarLocs(1)%Fnum
  write (*,*) '  Vnum: ', VarLocs(1)%Vnum
  write (*,*) 'Variable dimensions (after doubling time points):'
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

  allocate(VarData(1:GdataDescrip(1)%nx,1:GdataDescrip(1)%ny,1:GdataDescrip(1)%nz, &
           1:GdataDescrip(1)%nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: VarData array memory allocation failed'
    stop
  end if

  call ReadGradsData(GdataDescrip, VarName, VarLocs(1), VarData, GdataDescrip(1)%nx, &
       GdataDescrip(1)%ny, GdataDescrip(1)%nz, GdataDescrip(1)%nt)

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

  do it = 1, GdataDescrip(1)%nt
    if (it .eq. 1) then
      outTime = 1
    else
      outTime = outTime + 2
    end if
    do iz = 1, GdataDescrip(1)%nz
      do iy = 1, GdataDescrip(1)%ny
        do ix = 1, GdataDescrip(1)%nx
          if (it .eq. 1) then
            OutData(ix,iy,iz,outTime) = VarData(ix,iy,iz,it)
          else
            OutData(ix,iy,iz,(outTime-1)) = (VarData(ix,iy,iz,it) + VarData(ix,iy,iz,(it-1))) / 2.0
            OutData(ix,iy,iz,outTime) = VarData(ix,iy,iz,it)
          end if
        end do
      end do
    end do
  end do

  deallocate (VarData, stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: VarData array memory de-allocation failed'
    stop
  end if

  !Write out the new data

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
          GdataDescrip(1)%Tinc, GoutDescrip, 'dpts')

  call WriteGrads(GoutDescrip, OutData)

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

!*************************************************************************
! CheckDataDescripDouble()
!
! This routine will check the consistency of the data description obtained
! from the GRADS control files. The number of x, y, z, t points should match
! in all files. Also, the total number of variables are calculated and returned
! in Nvars.
!
! This routein will check one variable at a time, so it needs to be called
! multiple times for multiple variables.
!

subroutine CheckDataDescripDouble(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, VarLocs, VarName)
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

  ! Make Nt be the total time steps (after doubling the time steps)
  ! Actually adding Nt-1 time steps which fit in between the existing
  ! intput time steps.
  Nt = (2*Nt) - 1

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
