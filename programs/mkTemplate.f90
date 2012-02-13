!***************************************************************
! Program to find "interesting" regions in the storm
!
! This program will read in GRADS data from a RAMS simulation, and
! form groupings based on how a given quantity relates to a 
! given threshold.
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
  use GfileTypes
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarName, ThresholdDir

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles
  integer :: Nfiles

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars
  type (GradsOutDescription) :: GoutDescrip

  ! Data arrays: need one for var we are creating templates from
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of the var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: Var, Template
  real, dimension(1:1) :: DummyZcoords
  real :: Xstart, Xinc, Ystart, Yinc
  type (GradsVarLocation) :: VarLoc

  real :: VarThreshold
  integer :: SampleZ

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it

  ! Get the command line argument
  call GetMyArgs(Infiles, OfileBase, VarName, VarThreshold, ThresholdDir, SampleZ)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Calculating column integration for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  RAMS variable that is being templated: ', trim(VarName)
  write (*,*) '  Threshold value for determining which grid cells get included in the template: ', VarThreshold
  write (*,*) '  Direction variable value must be from threshold value for grid cell to be selected: ', ThresholdDir
  write (*,*) '  Sample Z: ', SampleZ
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

  if ((SampleZ .gt. Nz) .or. (SampleZ .lt. 1)) then
    write (*,*) 'ERROR: <sample_z> is out of range, must pick a value from 1 to ', Nz
    stop
  end if

  ! Allocate the data arrays and read in the data from the GRADS data files
  allocate (Var(1:Nx,1:Ny,1:Nz,1:Nt), Template(1:Nx,1:Ny,1:1,1:Nt), &
            stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(GdataDescrip, VarName , VarLoc, Var, Nx, Ny, Nz, Nt)

  ! Create the template
  !   mark regions where the variable value passes the threshold criteria
  !   group regions (give each region a unique id) according to locality (regions
  !     near each other in successive time steps are the same region)
  !
  ! In the Template array, the cell values have the following meaning:
  !   0 --> outside regions of interest
  !  -1 --> region of interest, but hasn't been assigned a unique group id yet
  !  >0 --> unique group ids  
  call MarkTemplateRegions(Nx, Ny, Nz, Nt, Var, Template, VarThreshold, ThresholdDir, SampleZ)
write (*,*) 'DEBUG: after MarkTemplateRegions' 
call DebugDumpTemplate(Nx, Ny, Nt, Template, 24)
  call GroupTemplateRegions(Nx, Ny, Nt, Template)
write (*,*) 'DEBUG: after GroupTemplateRegions' 
call DebugDumpTemplate(Nx, Ny, Nt, Template, 24)

  ! Write out result in GRADS format
  DummyZcoords(1) = 0.0
  Xstart = GdataDescrip(1)%Xcoords(1)
  Xinc = GdataDescrip(1)%Xcoords(2) - Xstart
  Ystart = GdataDescrip(1)%Ycoords(1)
  Yinc = GdataDescrip(1)%Ycoords(2) - Ystart
  call BuildGoutDescrip(Nx, Ny, 1, Nt, Template, OfileBase, GdataDescrip(1)%UndefVal, &
          VarName, Xstart, Xinc, Ystart, Yinc, DummyZcoords, &
          GdataDescrip(1)%Tstart, GdataDescrip(1)%Tinc, GoutDescrip, 'template')

  call WriteGrads(GoutDescrip, Template)

  ! Clean up
  deallocate (Var, Template, stat=Ierror)
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
!   VarName - RAMS variable which is used to build the template
!   VarThreshold - Threshold value for labeling the var location (cell)
!                  as part of the template.
!

subroutine GetMyArgs(Infiles, OfileBase, VarName, VarThreshold, ThresholdDir, SampleZ)
  implicit none

  character (len=*) :: Infiles, OfileBase, VarName, ThresholdDir
  real :: VarThreshold
  integer :: SampleZ

  integer :: iargc
  character (len=128) :: arg
  logical :: BadArgs

  if (iargc() .ne. 6) then
    write (*,*) 'ERROR: must supply exactly 6 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: mkTemplate <in_data_files> <out_data_file> <var_name> <var_threshold> <threshold_direction> <sample_z>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <var_name>: name of RAMS variable which is used for building the template'
    write (*,*) '        <var_threshold>: Threshold value for determining which cells become part of the template'
    write (*,*) '        <threshold_direction>: Direction past threshold that specifies how to select a grid cell'
    write (*,*) '           "up" --> select when variable value is greater than threshold'
    write (*,*) '           "down" --> select when variable value is less than threshold'
    write (*,*) '        <sample_z>: level to sample variable values (GRADS z number)'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)

  call getarg(2, OfileBase)

  call getarg(3, VarName)

  call getarg(4, arg)
  read (arg, '(f)') VarThreshold

  call getarg(5, ThresholdDir)

  call getarg(6, arg)
  read (arg, '(i)') SampleZ

  BadArgs = .false.

  if ((ThresholdDir .ne. 'up') .and. (ThresholdDir .ne. 'down')) then
    write (*,*) 'ERROR: <threshold_direction> must be one of "up" or "down"'
    BadArgs = .true.
  end if

  if (BadArgs) then
    stop
  end if

  return
end subroutine

!**********************************************************************
! MarkTemplateRegions()
!
! This routine will inspect the variable values and select grid cells based
! on the specified threshold
!

subroutine MarkTemplateRegions(Nx, Ny, Nz, Nt, Var, Template, VarThreshold, ThresholdDir, SampleZ)
  implicit none

  integer :: Nx, Ny, Nz, Nt, SampleZ
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: Var
  real, dimension(1:Nx,1:Ny,1:1,1:Nt) :: Template
  real :: VarThreshold
  character (len=*) :: ThresholdDir

  integer :: ix, iy, iz, it

  iz = SampleZ

  do it = 1, Nt
    do ix = 1, Nx
      do iy = 1, Ny
        if (ThresholdDir .eq. "up") then
           if (Var(ix,iy,iz,it) .gt. VarThreshold) then
             Template(ix,iy,1,it) = -1.0
           else
             Template(ix,iy,1,it) = 0.0
           end if
        else
           ! ThresholdDir is "down"
           if (Var(ix,iy,iz,it) .lt. VarThreshold) then
             Template(ix,iy,1,it) = -1.0
           else
             Template(ix,iy,1,it) = 0.0
           end if
        end if
      end do
    end do
  end do

  return
end subroutine

!**************************************************************************
! GroupTemplateRegions()
!
! This routine will group regions selected from MarkTemplateRegions() and
! assign unique ids to each group.
!
! A group is a contiguous collection of cells with a nonzero marker (initially,
! cells selected by MarkTemplateRegions() ).

subroutine GroupTemplateRegions(Nx, Ny, Nt, Template)
  implicit none

  integer :: Nx, Ny, Nt
  real, dimension(1:Nx,1:Ny,1:1,1:Nt) :: Template

  integer :: ix,iy,it
  real :: nextId, neighborId

  do it = 1, Nt
    nextId = 1.0
    do ix = 1, Nx
      do iy = 1, Ny
        if (Template(ix,iy,1,it) < 0.0) then
          call GetNeighborId(Nx, Ny, Nt, Template, ix, iy, it, neighborId)
          if (neighborId > 0.0) then
            ! found an assigned id in a neighbor cell
            Template(ix,iy,1,it) = neighborId
          else
            ! no assigned ids in neighbor cells, grab the next available id
            Template(ix,iy,1,it) = nextId
            nextId = nextId + 1.0
          end if
        end if
      end do
    end do
  end do

  return
end subroutine

!********************************************************************
! GetNeighborId()
!
! This routine will walk through the neighbor cells to a given cell
! and look at their ids. It will return -1 if there are no >0 values
! in the neighbor cells to signify that a new id is needed. If there
! are >0 values in neighbor cells, then the first one of these is
! returned.

subroutine GetNeighborId(Nx, Ny, Nt, Template, CurIx, CurIy, CurIt, neighborId)
  implicit none

  integer :: Nx, Ny, Nt, CurIx, CurIy, CurIt
  real, dimension(1:Nx,1:Ny,1:1,1:Nt) :: Template
  real :: neighborId

  integer :: ix, iy, it
  integer :: StartIx, StartIy, EndIx, EndIy

  ! make sure we don't run off the edges of the array
  StartIx = CurIx - 1
  EndIx = CurIx + 1
  if (StartIx .lt. 1) then
    StartIx = 1
  end if
  if (EndIx .gt. Nx) then
    EndIx = Nx
  end if

  StartIy = CurIy - 1
  EndIy = CurIy + 1
  if (StartIy .lt. 1) then
    StartIy = 1
  end if
  if (EndIy .gt. Ny) then
    EndIy = Ny
  end if

  it = CurIt
  neighborId = -1.0
  NLOOP: do ix = StartIx, EndIx
           do iy = StartIy, EndIy
             if (Template(ix,iy,1,it) > 0.0) then
               neighborId = Template(ix,iy,1,it)
               exit NLOOP
             end if
           end do
         end do NLOOP

  return
end subroutine

!**********************************************************************
! DebugDumpTemplate()
!
! Debug routine to dump out contents of template array at a given time
! step.
!

subroutine DebugDumpTemplate(Nx, Ny, Nt, Template, CurIt)
  implicit none

  integer :: Nx, Ny, Nt, CurIt
  real, dimension(1:Nx,1:Ny,1:1,1:Nt) :: Template

  integer :: ix, iy, it

  ! make the rows come out top to bottom in the output
  write (*,*) 'DEBUG: Contents of Template array at time: ', CurIt
  write (*,*) ''
  it = CurIt
  do iy = Ny, 1, -1
    do ix = 1, Nx
      write (*,*) ix, iy, Template(ix,iy,1,it)
    end do
  end do
  write (*,*) ''

  return
end subroutine
