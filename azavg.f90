!***************************************************************
! Program to do azimuthial averaging
!
! This program will read in GRADS data from a RAMS simulation, find
! the storm center and perform azimuthial averaging on the given
! quantity.
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
!   3. number of radial bands to split data into
!   4. RAMS quantity to perform the averaging on
!
! Output
!   The output will be two files which make a GRADS data set. One
!   file is the control file and the other is the binary data.
!

module GfileTypes
  integer, parameter :: MaxString=256
  integer, parameter :: MaxCoords=1000
  integer, parameter :: MaxVars=50
  integer, parameter :: InUnit=8
  integer, parameter :: OutUnit=10
  integer, parameter :: BinRecFactor=4
  
  type GradsDataDescription
    character (len=MaxString) :: DataFile
    real :: UndefVal
    real, dimension(1:MaxCoords) :: Xcoords
    real, dimension(1:MaxCoords) :: Ycoords
    real, dimension(1:MaxCoords) :: Zcoords
    real, dimension(1:MaxCoords) :: Tcoords
    character (len=MaxString), dimension(1:MaxVars) :: VarNames
    character (len=MaxString) :: Tstart, Tinc
    integer :: nx
    integer :: ny
    integer :: nz
    integer :: nt
    integer :: nvars
  end type GradsDataDescription

  type GradsVarLocation
    integer :: Fnum
    integer :: Vnum
  end type GradsVarLocation

  type GradsOutDescription
    character (len=MaxString) :: CtlFile
    character (len=MaxString) :: DataFile
    character (len=MaxString) :: Title
    character (len=MaxString) :: VarName
    real :: UndefVal
    real :: Xstart, Xinc
    real :: Ystart, Yinc
    real, dimension(1:MaxCoords) :: Zcoords
    character (len=MaxString) :: Tstart, Tinc
    integer :: nx
    integer :: ny
    integer :: nz
    integer :: nt
  end type GradsOutDescription
end module GfileTypes

program main
  use GfileTypes
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  integer NumRbands
  real :: Wthreshold
  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarToAvg
  logical :: DoHorizVel, DoTangential

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles
  integer :: Nfiles

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars
  type (GradsOutDescription) :: GoutDescrip

  ! Data arrays: need one for w (vertical velocity), press (pressure)
  ! and the var we are doing the averaging on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: U, V, W, Press, Avar, AzAvg
  integer, dimension(:), allocatable :: StmIx, StmIy
  real, dimension(:), allocatable :: MinP
  type (GradsVarLocation) :: Wloc, Uloc, Vloc, PressLoc, VarLoc

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it

  real :: DeltaX, DeltaY, RadialDist, RbandInc

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBase, NumRbands, Wthreshold, VarToAvg, DoHorizVel, DoTangential)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Calculating azimuthal average for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  Number of radial bands: ', NumRbands
  write (*,*) '  W threshold :           ', Wthreshold
  if (DoHorizVel) then
    if (DoTangential) then
      write (*,*) '  RAMS variable that is being averaged: Tangential Horizontal Velocity'
    else
      write (*,*) '  RAMS variable that is being averaged: Radial Horizontal Velocity'
    end if
  else
    write (*,*) '  RAMS variable that is being averaged: ', trim(VarToAvg)
  end if
  write (*,*) ''

  ! Read the GRADS data description files and collect the information about the data

  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''

  call CheckDataDescrip(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, Uloc, Vloc, Wloc, PressLoc, VarLoc, VarToAvg, DoHorizVel)

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
  write (*,'(a20,i3,a2,i3,a1)') 'w: (', Wloc%Fnum, ', ', Wloc%Vnum, ')'
  write (*,'(a20,i3,a2,i3,a1)') 'press: (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'
  if (DoHorizVel) then
    write (*,'(a20,i3,a2,i3,a1)') 'speed - u: (', Uloc%Fnum, ', ', Uloc%Vnum, ')'
    write (*,'(a20,i3,a2,i3,a1)') 'speed - v: (', Vloc%Fnum, ', ', Vloc%Vnum, ')'
  else
    write (*,'(a17,a3,i3,a2,i3,a1)') trim(VarToAvg), ': (', VarLoc%Fnum, ', ', VarLoc%Vnum, ')'
  end if
  write (*,*) ''

  ! Figure out how big to make each radial band. Assume that the storm center stays near the center of
  ! the grid domain. Take the diagonal distance of the domain, cut it in half and then break this up
  ! into NumRbands sections of equal length.

  DeltaX = GdataDescrip(1)%Xcoords(Nx) - GdataDescrip(1)%Xcoords(1)
  DeltaY = GdataDescrip(1)%Ycoords(Ny) - GdataDescrip(1)%Ycoords(1)
  RadialDist = sqrt(DeltaX*DeltaX + DeltaY*DeltaY) / 2.0
  RbandInc = RadialDist / real(NumRbands)

  write (*,*) 'Radial band information:'
  write (*,*) '  Delta x of domain:     ', DeltaX
  write (*,*) '  Delta y of domain:     ', DeltaY
  write (*,*) '  Radial distance:       ', RadialDist
  write (*,*) '  Radial band increment: ', RbandInc
  write (*,*) ''

  ! Allocate the data arrays and read in the data from the GRADS data files
  allocate (U(1:Nx,1:Ny,1:Nz,1:Nt), V(1:Nx,1:Ny,1:Nz,1:Nt), W(1:Nx,1:Ny,1:Nz,1:Nt), &
            Press(1:Nx,1:Ny,1:Nz,1:Nt), Avar(1:Nx,1:Ny,1:Nz,1:Nt), &
            StmIx(1:Nt), StmIy(1:Nt), MinP(1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(GdataDescrip, 'w', Wloc, W, Nx, Ny, Nz, Nt)
  call ReadGradsData(GdataDescrip, 'press', PressLoc, Press, Nx, Ny, Nz, Nt)
  call RecordStormCenter(Nx, Ny, Nz, Nt, Press, StmIx, StmIy, MinP)
  if (DoHorizVel) then
    call ReadGradsData(GdataDescrip, 'u', Uloc, U, Nx, Ny, Nz, Nt)
    call ReadGradsData(GdataDescrip, 'v', Vloc, V, Nx, Ny, Nz, Nt)
    call ConvertHorizVelocity(Nx, Ny, Nz, Nt, U, V, StmIx, StmIy, &
                    GdataDescrip(1)%Xcoords, GdataDescrip(1)%Ycoords, Avar, DoTangential)
  else
    call ReadGradsData(GdataDescrip, VarToAvg, VarLoc, Avar, Nx, Ny, Nz, Nt)
  end if

  ! Allocate the output array and compute the azimuthal averaging
  allocate (AzAvg(1:NumRbands,1:1,1:Nz,1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Ouput data array memory allocation failed'
    stop
  end if

  ! Do the averageing. Note that you can pick any index in GdataDescrip below since this data
  ! has been (superficially) checked out to be consistent.
  call AzimuthalAverage(Nx, Ny, Nz, Nt, NumRbands, W, StmIx, StmIy, MinP, Avar, AzAvg, &
          GdataDescrip(1)%Xcoords, GdataDescrip(1)%Ycoords, RadialDist, RbandInc, Wthreshold)

  call BuildGoutDescrip(NumRbands, 1, Nz, Nt, AzAvg, OfileBase, GdataDescrip(1)%UndefVal, VarToAvg, &
          0.0, RbandInc, 0.0, 1.0, GdataDescrip(1)%Zcoords, GdataDescrip(1)%Tstart, &
          GdataDescrip(1)%Tinc, GoutDescrip)

  call WriteGrads(GoutDescrip, AzAvg)

  ! Clean up
  deallocate (U, V, W, Press, Avar, StmIx, StmIy, MinP, AzAvg, stat=Ierror)
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
!   NumRbands - number of radial bands to split data into
!   VarToAvg - RAMS variable to do the averaging on
!

subroutine GetMyArgs(Infiles, OfileBase, NumRbands, Wthreshold, VarToAvg, DoHorizVel, DoTangential)
  implicit none

  integer :: NumRbands
  real :: Wthreshold
  character (len=*) :: Infiles, OfileBase, VarToAvg
  logical :: DoHorizVel, DoTangential

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 5) then
    write (*,*) 'ERROR: must supply exactly 5 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_file> <num_radial_bands> <w_threshold> <var_to_average>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <num_radial_bands>: number of bands to split data into'
    write (*,*) '        <w_threshold>: only take values of the variable where w is greater than or equal to this number'
    write (*,*) '        <var_to_average>: name of RAMS variable to do the averaging on'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)

  call getarg(3, arg)       !this is a string
  read (arg, '(i)') NumRbands !convert to integer

  call getarg(4, arg)
  read (arg, '(f)') Wthreshold

  call getarg(5, VarToAvg)
  if (VarToAvg == 'speed_t') then
    DoHorizVel = .true.
    DoTangential = .true.
  else
    if (VarToAvg == 'speed_r') then
      DoHorizVel = .true.
      DoTangential = .false.
    else
      DoHorizVel = .false.
      DoTangential = .false.
    end if
  end if

  return
end subroutine

!***********************************************************
! String2List()
!
! This routine will split the input string (colon separated list)
! into the individual strings (between the colons) and load these
! into the output array.
!
! Args
!  1. input list of files, colon separated items (strings)
!  2. character used for separator
!  3. output array
!  4. maximum number of input items
!  5. actual number of input items
!  6. description of what the items are

subroutine String2List(InString, Separator, OutList, MaxItems, NumItems, ItemType)
  implicit none

  character (len=*) InString
  character Separator
  character (len=*), dimension(*) :: OutList
  character (len=*) ItemType
  integer MaxItems, NumItems

  integer i, StrStart, StrEnd
  integer BetweenSep, PrevBetweenSep

  ! walk through the string from start to end
  ! when a colon is encountered, store the string
  ! from the "start" up to the colon in the output array

  NumItems = 0
  PrevBetweenSep = 0  ! treat the beginning of the line as separator
  i = 1
  do
    if (i .gt. len_trim(InString)) exit

    ! Record if looking at the string between the separators
    ! of if looking at the separators.
    ! Remember the previous state of BetweenSep:
    !   If you just changed from separators to string, then this
    !   is the beginning of the string
    !   If you just changed from string to separators, then this
    !   is the end of the string --> then record it in the list

    if (InString(i:i) .eq. Separator) then
      BetweenSep = 0
    else
      BetweenSep = 1
    end if

    if ((PrevBetweenSep .eq. 0) .and. (BetweenSep .eq. 1)) then
      ! Start of next string
      StrStart = i
    else
      if ((PrevBetweenSep .eq. 1) .and. (BetweenSep .eq. 0)) then
        ! End of string, record it in the list
        NumItems = NumItems + 1
        StrEnd = i - 1
        if (NumItems .le. MaxItems) then
          OutList(NumItems) = InString(StrStart:StrEnd)
        else
          write (*,*) 'ERROR: Maximum number of ', ItemType, ' (', MaxItems, ') exceeded'
          stop
        end if
      end if
    end if

    i = i + 1
    PrevBetweenSep = BetweenSep
  end do

  ! may have one more string at this point
  ! treat the end of the line as a separator

  if (BetweenSep .eq. 1) then
    NumItems = NumItems + 1
    StrEnd = i - 1
    if (NumItems .le. MaxItems) then
      OutList(NumItems) = InString(StrStart:StrEnd)
    else
      write (*,*) 'ERROR: Maximum number of ', ItemType, ' (', MaxItems, ') exceeded'
      stop
    end if
  end if

  return
end subroutine

!**********************************************************
! ReadGradsCtlFile()
!
! This routine will read in the GRADS control files
!
! Args
!   1. Name of grads control file
!   2. Grads data description
!

subroutine  ReadGradsCtlFile(GradsCtlFile, GdataDescrip)
  use GfileTypes
  implicit none

  integer, parameter :: MaxLine = 1000
  integer, parameter :: MaxFields = 100

  character (len=*) :: GradsCtlFile
  type (GradsDataDescription) GdataDescrip

  integer InError
  character (len=MaxLine) InLine
  character (len=MaxLine), dimension(1:MaxFields) :: InFields
  integer Nfields

  character (len=MaxString) :: GfileDir

  open(unit=InUnit, file=GradsCtlFile, form='formatted', &
       status='old', action='read', iostat=InError)
  if (InError .ne. 0) then
    write (*,*) 'ERROR: cannot open GRADS control file: ', trim(GradsCtlFile)
    stop
  end if

  InError = 0
  do
    read (unit=InUnit, fmt='(a)', iostat=InError) InLine
    if (InError .ne. 0) exit

    ! split the line into a list of fields separated by white space
    call String2List(InLine, ' ', InFields, MaxFields, Nfields, 'control file fields')

    if ((InFields(1) .eq. 'dset') .or. (InFields(1) .eq. 'DSET')) then
      if (InFields(2)(1:1) .eq. '^') then
        !replace the '^' with the directory portion of input control file
        if (index(GradsCtlFile, '/', .true.) .eq. 0) then
          ! no leading path
          GfileDir = ''
        else
          GfileDir = GradsCtlFile(1:index(GradsCtlFile, '/', .true.))
        end if
        ! Make sure to use trim(), without it the first string GfileDir fills up the
        ! GdataDescrip%DataFile string since these are both declared to be the same length
        GdataDescrip%DataFile = trim(GfileDir) // trim(InFields(2)(2:len(InFields(2))))
      else
        GdataDescrip%DataFile = InFields(2)
      end if
    else
      if ((InFields(1) .eq. 'undef') .or. (InFields(1) .eq. 'UNDEF')) then
        read (InFields(2), '(f)') GdataDescrip%UndefVal
      else
        if ((InFields(1) .eq. 'xdef') .or. (InFields(1) .eq. 'XDEF')) then
          call GenCoords(InFields, Nfields, GdataDescrip%Xcoords, GdataDescrip%nx, MaxCoords, 0)
        else
          if ((InFields(1) .eq. 'ydef') .or. (InFields(1) .eq. 'YDEF')) then
            call GenCoords(InFields, Nfields, GdataDescrip%Ycoords, GdataDescrip%ny, MaxCoords, 0)
          else
            if ((InFields(1) .eq. 'zdef') .or. (InFields(1) .eq. 'ZDEF')) then
              call GenCoords(InFields, Nfields, GdataDescrip%Zcoords, GdataDescrip%nz, MaxCoords, 0)
            else
              if ((InFields(1) .eq. 'tdef') .or. (InFields(1) .eq. 'TDEF')) then
                GdataDescrip%Tstart = InFields(4)
                GdataDescrip%Tinc   = InFields(5)
                call GenCoords(InFields, Nfields, GdataDescrip%Tcoords, GdataDescrip%nt, MaxCoords, 1)
              else
                if ((InFields(1) .eq. 'vars') .or. (InFields(1) .eq. 'VARS')) then
                  read (InFields(2), '(i)') GdataDescrip%nvars
                  if (GdataDescrip%nvars .gt. MaxVars) then
                    write (*,*) 'ERROR: Maximum number of variables (', MaxVars, &
                                ') exceeded in GRADS control file specs'
                    stop
                  end if
                  call GetVarNames(InUnit, GdataDescrip%VarNames, GdataDescrip%nvars, MaxLine, MaxFields)
                end if
              end if
            end if
          end if
        end if
      end if
    end if
  end do

  close(unit=InUnit, status='keep')

  return
end subroutine

!*****************************************************************
! GenCoords()
!
! This routine will read in the dimension data from the given
! string (which is an input line from the GRADS control file); and
! generate the coordinate values for that dimension.
!
! This should be one of the xdef, ydef, zdef, tdef specs.
!

subroutine GenCoords(InFields, Nfields, Coords, Ncoords, MaxCoords, GenDmyCoords)
  implicit none

  character (len=*), dimension(*) :: InFields
  integer Nfields
  real Coords(*)
  integer Ncoords, MaxCoords, GenDmyCoords

  real Xstart, Xincr
  integer i

  read (InFields(2), '(i)') Ncoords 
  if (Ncoords .gt. MaxCoords) then
    write (*,*) 'ERROR: Maximum number of coordinates (', MaxCoords, ') exceeded in GRADS control file specs'
    stop
  end if

  if (GenDmyCoords .ne. 1) then
    ! read in the real coordinate values
    if ((InFields(3) .eq. 'linear') .or. (InFields(3) .eq. 'LINEAR')) then
      read (InFields(4), '(f)') Xstart
      read (InFields(5), '(f)') Xincr

      do i = 1, Ncoords
        Coords(i) = Xstart + (real(i - 1) * Xincr)
      end do

    else
      if ((InFields(3) .eq. 'levels') .or. (InFields(3) .eq. 'LEVELS')) then
        do i = 1, Ncoords
          read(InFields(i+3), '(f)') Coords(i)
        end do
      end if
    end if
  else
    ! Generate dummy coordinate values
    do i = 1, Ncoords
      Coords(i) = real(i)
    end do
  end if

  return
end subroutine

!*****************************************************************************
! GetVarNames()
!
! This routine will read Nvars lines from the file given by InUnit, and record
! the names of the variables found in the VarNames array.
!

subroutine  GetVarNames(InUnit, VarNames, Nvars, MaxLine, MaxFields)
  implicit none

  integer InUnit
  character (len=*), dimension(*) :: VarNames
  integer Nvars
  integer MaxLine
  integer MaxFields

  integer InError
  character (len=MaxLine) InLine
  character (len=MaxLine), dimension(1:MaxFields) :: InFields
  integer Nfields

  integer i

  do i = 1, Nvars
    read (unit=InUnit, fmt='(a)', iostat=InError) InLine
    if (InError .ne. 0) then
      write (*,*) 'ERROR: Reached end of GRADS control file before reading in all (', &
                  Nvars, ') variable names'
      stop
    end if

    call String2List(InLine, ' ', InFields, MaxFields, Nfields, 'control file fields')

    VarNames(i) = InFields(1)
  end do

  return
end subroutine

!*************************************************************************
! CheckDataDescrip()
!
! This routine will check the consistency of the data description obtained
! from the GRADS control files. The number of x, y, z, t points should match
! in all files. Also, the total number of variables are calculated and returned
! in Nvars.
!

subroutine CheckDataDescrip(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, Uloc, Vloc, Wloc, &
           PressLoc, VarLoc, VarName, DoHorizVel)
  use GfileTypes
  implicit none

  type (GradsDataDescription), dimension(*) :: GdataDescrip
  integer Nx, Ny, Nz, Nt, Nfiles
  integer Nvars
  type (GradsVarLocation) :: Uloc, Vloc, Wloc, PressLoc, VarLoc
  character (len=*) :: VarName
  logical :: DoHorizVel

  ! i -> file number, j -> var number
  integer i, j
  logical BadData

  Uloc%Fnum = 0
  Uloc%Vnum = 0
  Vloc%Fnum = 0
  Vloc%Vnum = 0
  Wloc%Fnum = 0
  Wloc%Vnum = 0
  PressLoc%Fnum = 0
  PressLoc%Vnum = 0
  VarLoc%Fnum = 0
  VarLoc%Vnum = 0

  BadData = .false.
  do i = 1, Nfiles
    if (i .eq. 1) then
      Nx = GdataDescrip(i)%nx
      Ny = GdataDescrip(i)%ny
      Nz = GdataDescrip(i)%nz
      Nt = GdataDescrip(i)%nt
   
      Nvars = GdataDescrip(i)%nvars
    else
      if (GdataDescrip(i)%nx .ne. Nx) then
        write (*,*) 'ERROR: number of x points in GRADS control files do not match'
        BadData = .true.
      end if
      if (GdataDescrip(i)%ny .ne. Ny) then
        write (*,*) 'ERROR: number of y points in GRADS control files do not match'
        BadData = .true.
      end if
      if (GdataDescrip(i)%nz .ne. Nz) then
        write (*,*) 'ERROR: number of z points in GRADS control files do not match'
        BadData = .true.
      end if
      if (GdataDescrip(i)%nt .ne. Nt) then
        write (*,*) 'ERROR: number of t points in GRADS control files do not match'
        BadData = .true.
      end if

      Nvars = Nvars + GdataDescrip(i)%nvars
    end if

    do j =1, GdataDescrip(i)%nvars
      if (GdataDescrip(i)%VarNames(j) .eq. 'u') then
        Uloc%Fnum = i
        Uloc%Vnum = j
      end if
      if (GdataDescrip(i)%VarNames(j) .eq. 'v') then
        Vloc%Fnum = i
        Vloc%Vnum = j
      end if
      if (GdataDescrip(i)%VarNames(j) .eq. 'w') then
        Wloc%Fnum = i
        Wloc%Vnum = j
      end if
      if (GdataDescrip(i)%VarNames(j) .eq. 'press') then
        PressLoc%Fnum = i
        PressLoc%Vnum = j
      end if
      if (GdataDescrip(i)%VarNames(j) .eq. VarName) then
        VarLoc%Fnum = i
        VarLoc%Vnum = j
      end if
    end do
  end do

  ! Check to see if you got all three vars (w, press, var)
  if (Wloc%Fnum .eq. 0) then
    write (*,*) 'ERROR: cannot find grid var "w" in the GRADS data files'
    BadData = .true.
  end if
  if (PressLoc%Fnum .eq. 0) then
    write (*,*) 'ERROR: cannot find grid var "press" in the GRADS data files'
    BadData = .true.
  end if
  if (DoHorizVel) then
    if (Uloc%Fnum .eq. 0) then
      write (*,*) 'ERROR: cannot find grid var "u" in the GRADS data files'
      BadData = .true.
    end if
    if (Vloc%Fnum .eq. 0) then
      write (*,*) 'ERROR: cannot find grid var "v" in the GRADS data files'
      BadData = .true.
    end if
  else
    if (VarLoc%Fnum .eq. 0) then
      write (*,*) 'ERROR: cannot find grid var "', trim(VarName), '" in the GRADS data files'
      BadData = .true.
    end if
  end if

  if (BadData) then
    stop
  end if
end subroutine

!********************************************************************
! ReadGradsData()
!
! This routine will use the information in GdataDescrip and VarLoc to
! locate the data for the var in the GRADS data files. This data will
! then be loaded into the VarData array.
!

subroutine ReadGradsData(GdataDescrip, VarName, VarLoc, VarData, Nx, Ny, Nz, Nt)
  use GfileTypes
  implicit none

  type (GradsDataDescription), dimension(*) :: GdataDescrip
  character (len=*) :: VarName
  type (GradsVarLocation) :: VarLoc
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: VarData
  integer Nx, Ny, Nz, Nt

  integer Nvars, VarNum
  character (len=MaxString) :: DataFile
  integer RecLen, Ierror
  integer ix,iy,iz,it
  integer RecNum
  integer NumRecs

  Nvars = GdataDescrip(VarLoc%Fnum)%nvars
  VarNum = VarLoc%Vnum

  ! Each record is a horizontal slice so it has Nx * Ny elements in it
  ! of BinRecFactor size
  RecLen = Nx * Ny * BinRecFactor

  DataFile = GdataDescrip(VarLoc%Fnum)%DataFile
  open (unit=InUnit, file=DataFile, form='unformatted', access='direct', &
        recl=RecLen, status='old', action='read', iostat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Cannot open GRADS data file: ', trim(DataFile)
    stop
  end if
 
  ! The data file contains Nx * Ny * Nz * Nt * Nvars total elements
  ! One record is Nx * Ny elements so the number of records is
  ! Nz * Nt * Nvars
  ! The data records are organized with Z changing fastest, then variable,
  ! then time the slowest. Eg. Nz = 2, Nvars =3, Nt = 4,
  !
  !   Record  Z_level  Variable Time
  !     1       1         1      1
  !     2       2         1      1
  !     3       1         2      1
  !     4       2         2      1
  !     5       1         3      1
  !     6       2         3      1
  !
  !     7       1         1      2
  !     8       2         1      2
  !     9       1         2      2
  !    10       2         2      2
  !    11       1         3      2
  !    12       2         3      2
  !
  !    13       1         1      3
  !    14       2         1      3
  !    15       1         2      3
  !    16       2         2      3
  !    17       1         3      3
  !    18       2         3      3
  !
  !    19       1         1      4
  !    20       2         1      4
  !    21       1         2      4
  !    22       2         2      4
  !    23       1         3      4
  !    24       2         3      4

  NumRecs = Nz * Nt * Nvars

  write (*,*) 'Reading GRADS data file: ', trim(DataFile)
  write (*,*) '  Var: ', trim(VarName)
  write (*,*) '  X points:          ', Nx
  write (*,*) '  Y points:          ', Ny
  write (*,*) '  Z points:          ', Nz
  write (*,*) '  T points:          ', Nt
  write (*,*) '  Number of Vars:    ', Nvars
  write (*,*) '  Record length:     ', RecLen
  write (*,*) '  Number of records: ', NumRecs
  write (*,*) ''

  ! Need to derive the correct record number in the file
  ! from the it and iz values. it cycles by Nz*Nvars records.
  ! Within each it cycle, Nz*(VarNum-1) is the offset from the
  ! start of the cycle to the VarNum position.
  ! The formula becomes:
  !
  !  RecNum = (it-1)*(Nz*Nvars) + (Nz*(VarNum-1)) + iz
  
  NumRecs = 0
  do it = 1, Nt
    do iz = 1, Nz
      RecNum = (it-1)*(Nz*Nvars) + (Nz*(VarNum-1)) + iz
      read(unit=InUnit, rec=RecNum) ((VarData(ix,iy,iz,it), ix = 1, Nx), iy = 1, Ny)

      NumRecs = NumRecs + 1
      if (modulo(NumRecs,100) .eq. 0) then
        write (*,*) '  Reading: ', trim(VarName), ': ', NumRecs
      end if
    end do   
  end do 
  write (*,*) '  Total records read: ', NumRecs
  write (*,*) ''

  close (unit=InUnit, status='keep')

  return
end subroutine

!**************************************************************************
! WriteGrads()
!
! This routine will write out the given data array into the GRADS (input)
! format. This creates two files: 
!   control file - contains description of the data in the data file
!   data file - raw binary values in the data array
!

subroutine WriteGrads(GoutDescrip, AzAvg)
  use GfileTypes
  implicit none

  type (GradsOutDescription) :: GoutDescrip
  real, dimension(1:GoutDescrip%nx,1:GoutDescrip%ny,1:GoutDescrip%nz,1:GoutDescrip%nt) :: AzAvg
  integer :: Nx, Ny, Nz, Nt

  integer :: OutRecLen
  integer :: RecNum
  integer :: Ierror
  integer :: ix
  integer :: iy
  integer :: iz
  integer :: it

  write (*,*) 'Writing out result in GRADS format:'
  write (*,*) '  Control file: ', trim(GoutDescrip%CtlFile)
  write (*,*) '  Data file:    ', trim(GoutDescrip%DataFile)
  write (*,*) '  Total number of data points: ', &
              GoutDescrip%nx * GoutDescrip%ny * GoutDescrip%nz * GoutDescrip%nt
  write (*,*) ''

  ! Control (data description) file
  open (OutUnit, file=GoutDescrip%CtlFile, form='formatted', action='write', status='replace', iostat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: cannot open output GRADS control file for writing: ', trim(GoutDescrip%CtlFile)
    stop
  end if

  write (OutUnit, '(a,a)')          'DSET ', trim(GoutDescrip%DataFile)
  write (OutUnit, '(a,a)')          'TITLE ', trim(GoutDescrip%Title)
  write (OutUnit, '(a,g)')          'UNDEF ', GoutDescrip%UndefVal
  write (OutUnit, '(a)')            'OPTIONS  little_endian'
  write (OutUnit, '(a,i,a,g,g)')    'XDEF ', GoutDescrip%nx, ' LINEAR ', GoutDescrip%Xstart, GoutDescrip%Xinc
  write (OutUnit, '(a,i,a,g,g)')    'YDEF ', GoutDescrip%ny, ' LINEAR ', GoutDescrip%Ystart, GoutDescrip%Yinc
  write (OutUnit, '(a,i,a,100g')    'ZDEF ', GoutDescrip%nz, ' LEVELS ', &
                                    (GoutDescrip%Zcoords(iz), iz = 1, GoutDescrip%nz)
  write (OutUnit, '(a,i,a,a,5x,a)') 'TDEF ', GoutDescrip%nt, ' LINEAR ', &
                                    trim(GoutDescrip%Tstart), trim(GoutDescrip%Tinc)
  write (OutUnit, '(a)')            'VARS 1 '
  write (OutUnit, '(a,i,a)')        trim(GoutDescrip%VarName) , GoutDescrip%nz, ' 99 Azimuthal Averaged Data'
  write (OutUnit, '(a)')            'ENDVARS'

  close (OutUnit, status='keep')

  ! Data file
  ! dimensions from fastest changing to slowest changing are: x, y, z, t
  ! One record is a single horizontal slice -> (nx * ny) data points
  OutRecLen = GoutDescrip%nx * GoutDescrip%ny * BinRecFactor
  open (OutUnit, file=GoutDescrip%DataFile, form='unformatted', access='direct', &
        recl=OutRecLen, action='write', status='replace', iostat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: cannot open output GRADS data file for writing: ', trim(GoutDescrip%DataFile)
    stop
  end if

  RecNum = 1
  do it = 1, GoutDescrip%nt
    do iz = 1, GoutDescrip%nz
      if (modulo(RecNum,100) .eq. 0) then
        write (*,*) '  Writing: ', trim(GoutDescrip%VarName), ': ', RecNum
      end if
      write (OutUnit, rec=RecNum) ((AzAvg(ix,iy,iz,it), ix = 1, GoutDescrip%nx), iy = 1, GoutDescrip%ny)
      RecNum = RecNum + 1
    end do
  end do

  close (OutUnit, status='keep')

  return
end subroutine

!**********************************************************************
! RecordStormCenter()
!
! This routine will record the storm center for each time step
!

subroutine RecordStormCenter(Nx, Ny, Nz, Nt, Press, StmIx, StmIy, MinP)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: Press
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nt) :: MinP

  integer it

  do it = 1, Nt
    call FindStormCenter(Press, Nx, Ny, Nz, Nt, it, StmIx(it), StmIy(it), MinP(it))
  end do

  return
end subroutine

!**********************************************************************
! FindStormCenter
!
! This routine will locate the storm center using the simple hueristic
! of the center being where the minimum surface pressure exists.
!
! Argument it holds the time step that you want to analyze. (iStmCtr, jStmCtr)
! hold the grid position of the minumum pressure value on the first vertical
! level (iz = 1).
!

subroutine FindStormCenter(Press, Nx, Ny, Nz, Nt, iTime, ixStmCtr, iyStmCtr, MinP)
  implicit none

  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: Press
  integer :: Nx, Ny, Nz, Nt
  integer :: iTime, ixStmCtr, iyStmCtr
  real :: MinP

  integer :: ix, iy, iz, it

  iz = 1
  it = iTime
  MinP = 1e10 ! ridiculously large pressure
  ixStmCtr = 0
  iyStmCtr = 0
  
  do ix = 1, Nx
    do iy = 1, Ny
      if (Press(ix,iy,iz,it) .lt. MinP) then
        MinP = Press(ix,iy,iz,it)
        ixStmCtr = ix
        iyStmCtr = iy
      end if
    end do
  end do
  
  return
end subroutine

!*************************************************************************
! AzimuthalAverage
!
! This routine will do the azimuthal averaging. It will create a 4-D array
! for its output (to keep GRADS write routine consistent/general) even though it
! really is 3-D data (radial band, z level, time). It will put the radial bands 
! in the x dimension and make the y dimension 1 long.
!
! This routine will step though each z and t value and compute the azimuthal
! average of the x-y plane at each of these points. Then it will store that
! in the output array so that you get azimuthal averaging for every original
! z and t point.
!

subroutine AzimuthalAverage(Nx, Ny, Nz, Nt, NumRbands, W, StmIx, StmIy, MinP, Avar, AzAvg, &
          Xcoords, Ycoords, RadialDist, RbandInc, Wthreshold)

  implicit none

  integer :: Nx, Ny, Nz, Nt, NumRbands
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: W, Avar
  real, dimension(1:NumRbands,1:1,1:Nz,1:Nt) :: AzAvg
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nt) :: MinP
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  real :: RadialDist, RbandInc, Wthreshold

  integer :: ix, iy, iz, it
  real, dimension(1:NumRbands) :: Rcounts
  integer :: ir, iRband
  real :: DeltaX, DeltaY, Radius

  ! Mask the input data (Avar) with the W data, ie if the corresponding
  ! W data is >= to the Wthreshold use the Avar data in the averaging;
  ! otherwise skip the Avar data
  !
  ! The storm center is taken to be the min pressure location of the
  ! x-y plane on the surface (iz .eq. 1)

  write (*,*) 'Averaging Data:'

  do it = 1, Nt
    if (modulo(it,10) .eq. 0) then
    write (*,*) '  Timestep: ', it
    write (*,'(a,i3,a,i3,a,g,a,g,a)') '    Storm Center: (', StmIx(it), ', ', StmIy(it), ') --> (', &
          Xcoords(StmIx(it)), ', ', Ycoords(StmIy(it)), ')'
    write (*,*) '    Minimum Pressue: ', MinP(it)
    end if
    do iz = 1, Nz

      ! For the averaging
      do ir = 1, NumRbands
        Rcounts(ir) = 0.0
        AzAvg(ir, 1, Nz, Nt) = 0.0
      end do

      ! Get the averages for this x-y plane
      do iy = 1, Ny
        do ix = 1, Nx
          ! Only keep the points where W meets the threshold
          if (W(ix,iy,iz,it) .ge. Wthreshold) then
             DeltaX = Xcoords(ix)-Xcoords(StmIx(it))
             DeltaY = Ycoords(iy)-Ycoords(StmIy(it))
             Radius = sqrt(DeltaX*DeltaX + DeltaY*DeltaY)
             iRband = int(Radius / RbandInc) + 1
             ! iRband will go from 0 to n, but n might extend beyond the last radial
             ! band due to approximations made in deriving RbandInc
             if (iRband .gt. NumRbands) then
               iRband = NumRbands
             end if

             AzAvg(iRband, 1, iz, it) = AzAvg(iRband, 1, iz, it) + Avar(ix, iy, iz, it)
             Rcounts(iRband) = Rcounts(iRband) + 1.0
          end if
        end do
      end do

      do ir = 1, NumRbands
        ! if we didn't put anything into an AzAvg slot then it should be zero from
        ! from the intialization above, and we can keep it that way for the result.
        if (Rcounts(ir) .ne. 0.0) then
          AzAvg(ir, 1, iz, it) = AzAvg(ir, 1, iz, it) / Rcounts(ir)
        end if
      end do

    end do
  end do
  write (*,*) ''

  return
end subroutine

!******************************************************************************
! BuildGoutDescrip()
!
! This routine will fill in the GoutDescrip record for the AzAvg data. This
! record is used for creating the GRADS control file.
!

subroutine BuildGoutDescrip(Nx, Ny, Nz, Nt, AzAvg, OfileBase, UndefVal, VarName, &
          Xstart, Xinc, Ystart, Yinc, Zcoords, Tstart, Tinc, GoutDescrip)

  use GfileTypes
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: AzAvg
  character (len=*) :: OfileBase, VarName
  real :: UndefVal, Xstart, Xinc, Ystart, Yinc
  real, dimension(1:Nz) :: Zcoords
  character (len=*) :: Tstart, Tinc
  type (GradsOutDescription) :: GoutDescrip

  integer :: iz

  GoutDescrip%CtlFile = trim(OfileBase) // '.ctl'
  GoutDescrip%DataFile = trim(OfileBase) // '.gra'
  GoutDescrip%Title = 'Azimuthal Averaged: ' // trim(VarName)
  GoutDescrip%UndefVal = UndefVal
  GoutDescrip%nx = Nx
  GoutDescrip%ny = Ny
  GoutDescrip%nz = Nz
  GoutDescrip%nt = Nt
  GoutDescrip%Xstart = Xstart
  GoutDescrip%Xinc = Xinc
  GoutDescrip%Ystart = Ystart
  GoutDescrip%Yinc = Yinc
  GoutDescrip%VarName = trim(VarName) // '_azavg'
  do iz = 1, Nz
    GoutDescrip%Zcoords(iz) = Zcoords(iz)
  end do
  GoutDescrip%Tstart = Tstart
  GoutDescrip%Tinc = Tinc

  return
end subroutine

!*******************************************************************************
! ConvertHorizVelocity()
!
! This routine will convert the horizontal velocity vectors (described in U and V)
! into tangential or radial components given the storm center location.
!

subroutine ConvertHorizVelocity(Nx, Ny, Nz, Nt, U, V, StmIx, StmIy, Xcoords, Ycoords, Avar, DoTangential)
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: U, V, Avar
  integer, dimension(1:Nt) :: StmIx, StmIy
  real, dimension(1:Nx) :: Xcoords
  real, dimension(1:Ny) :: Ycoords
  logical :: DoTangential

  integer :: ix, iy, iz, it

  do it = 1, Nt
    do iz = 1, Nz
      do iy = 1, Ny
        do ix = 1, Nx
          Avar(ix,iy,iz,it) = U(ix,iy,iz,it)
        end do
      end do
    end do
  end do
  
  return
end subroutine
