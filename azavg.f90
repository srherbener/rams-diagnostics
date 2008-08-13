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

  type GradsDataDescription
    character (len=MaxString) :: DataFile
    real UndefVal
    real, dimension(1:MaxCoords) :: Xcoords
    real, dimension(1:MaxCoords) :: Ycoords
    real, dimension(1:MaxCoords) :: Zcoords
    real, dimension(1:MaxCoords) :: Tcoords
    character (len=MaxString), dimension(1:MaxVars) :: VarNames
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
end module GfileTypes

program main
  use GfileTypes

  integer, parameter :: BinRecFactor=4
  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  integer NumRbands
  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarToAvg

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles
  integer Nfiles

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer Nx, Ny, Nz, Nt, Nvars

  ! Data arrays: need one for w (vertical velocity), press (pressure)
  ! and the var we are doing the averaging on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: W, Press, Avar
  type (GradsVarLocation) :: Wloc, PressLoc, VarLoc

  ! Index vars for the grid data
  !   i -> x -> longitude
  !   j -> y -> latitude
  !   k -> z -> vertical levels
  !   l -> var -> gridded variables
  !   m -> t -> time
  integer i, j, k, l, m
  integer Ierror

!  integer DoutUnit, CoutUnit
!  character *128 DoutFile, CoutFile
!  integer OutRecLen
!  integer i,j,k,nt
!  real Xstart, Xinc, Ystart, Yinc, Zstart, Zinc
!  character *20 Tstart, Tinc
!  real OutData(1:MaxI,1:MaxJ,1:MaxK,1:MaxT)
!  real UndefVal;

  ! Get the command line arguments
  call GetMyArgs(Infiles, OfileBase, NumRbands, VarToAvg)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Calculating azimuthal average for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name: ', trim(OfileBase)
  write (*,*) '  Number of radial bands: ', NumRbands
  write (*,*) '  RAMS variable that is being averaged: ', trim(VarToAvg)
  write (*,*) ''

  ! Read the GRADS data description files and collect the information about the data

  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''

  call CheckDataDescrip(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, Wloc, PressLoc, VarLoc, VarToAvg)

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
  write (*,*) '  w (vertical velocity): (', Wloc%Fnum, ', ', Wloc%Vnum, ')'
  write (*,*) '  press (pressure): (', PressLoc%Fnum, ', ', PressLoc%Vnum, ')'
  write (*,*) '  ', trim(VarToAvg), ': (', VarLoc%Fnum, ', ', VarLoc%Vnum, ')'
  write (*,*) ''

  ! Allocate the data array and read in the data from the GRADS data files
  allocate (W(1:Nx,1:Ny,1:Nz,1:Nt),Press(1:Nx,1:Ny,1:Nz,1:Nt), Avar(1:Nx,1:Ny,1:Nz,1:Nt),STAT=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  ! Clean up
  deallocate (W, Press, Avar, STAT=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory de-allocation failed'
    stop
  end if


!   ! Set up the output files, GRADS format (two files)
!   CoutUnit = 8
!   CoutFile = '/home/sherbener/projects/grads/udf/data/testGfile.ctl' 
!   DoutUnit = 10
!   DoutFile = '/home/sherbener/projects/grads/udf/data/testGfile.dat'
! 
!   OutRecLen = MaxI * MaxJ * RecFactor
!   open (CoutUnit, file=CoutFile, form='formatted')
!   open (DoutUnit, file=DoutFile, form='unformatted', access='direct', recl=OutRecLen)
! 
!   !******************************************************************
!   ! Try building a known function and see if get the right picture
!   ! in GRADS. In the binary data file, the order of changing indexes
!   ! go lon, lat, lev, time which are i, j, k, t respectively. These
!   ! correspond to x, y, z, t in the plots.
! 
!   ! Try a 2D plot x-y planes
! 
!   ! Control
!   UndefVal = 1.0e30
!   Xstart = 1.0
!   Xinc = 1.0
!   Ystart = 1.0
!   Yinc = 1.0
!   Zstart = 1.0
!   Zinc = 1.0
!   Tstart = '00:00Z01jan2000'
!   Tinc = '1hr'
! 
!   write (CoutUnit, '(a,a)')       'DSET ', DoutFile
!   write (CoutUnit, '(a)')         'TITLE  1D plot example'
!   write (CoutUnit, '(a,g)')       'UNDEF ', UndefVal
!   write (CoutUnit, '(a)')         'OPTIONS  little_endian'
!   write (CoutUnit, '(a,i,a,g,g)') 'XDEF ', MaxI, ' LINEAR ', Xstart, Xinc
!   write (CoutUnit, '(a,i,a,g,g)') 'YDEF ', MaxJ, ' LINEAR ', Ystart, Yinc
!   write (CoutUnit, '(a,i,a,g,g)') 'ZDEF ', MaxK, ' LINEAR ', Zstart, Zinc
!   write (CoutUnit, '(a,i,a,g,g)') 'TDEF ', MaxT, ' LINEAR ', Tstart, Tinc
!   write (CoutUnit, '(a)')         'VARS 1 '
!   write (CoutUnit, '(a,i,a)')     'f ', MaxK, ' 99 Test 1D data'
!   write (CoutUnit, '(a)')         'ENDVARS'
! 
!   ! Data
!   do j = 1, MaxJ
!     do i = 1, MaxI
!       do nt = 1, MaxT
!         OutData(i,j,1,nt) = real(i * (j+nt))
!       end do
!     end do
!   end do
!   do nt = 1, MaxT
!     write (DoutUnit, rec=nt) ((OutData(i,j,1,nt), i = 1, MaxI), j = 1, MaxJ)
!   end do

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

subroutine GetMyArgs(Infiles, OfileBase, NumRbands, VarToAvg)

  integer NumRbands
  character (len=*) :: Infiles, OfileBase, VarToAvg

  integer iarg
  character (len=128) :: arg

  if (iargc() .ne. 4) then
    write (*,*) 'ERROR: must supply exactly 4 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: azavg <in_data_files> <out_data_file> <num_radial_bands> <var_to_average>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <num_radial_bands>: number of bands to split data into'
    write (*,*) '        <var_to_average>: name of RAMS variable to do the averaging on'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)

  call getarg(3, arg)       !this is a string
  read (arg, '(i)') NumRbands !convert to integer

  call getarg(4, VarToAvg)

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
        if (index(GradsCtlFile, '/', .TRUE.) .eq. 0) then
          ! no leading path
          GfileDir = ''
        else
          GfileDir = GradsCtlFile(1:index(GradsCtlFile, '/', .TRUE.))
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

subroutine CheckDataDescrip(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, Wloc, PressLoc, VarLoc, VarName)

  use GfileTypes

  type (GradsDataDescription), dimension(*) :: GdataDescrip
  integer Nx, Ny, Nz, Nt, Nfiles
  integer Nvars
  type (GradsVarLocation) :: Wloc, PressLoc, VarLoc
  character (len=*) :: VarName

  ! i -> file number, j -> var number
  integer i, j
  logical BadData

  Wloc%Fnum = 0
  Wloc%Vnum = 0
  PressLoc%Fnum = 0
  PressLoc%Vnum = 0
  VarLoc%Fnum = 0
  VarLoc%Vnum = 0

  BadData = .FALSE.
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
        BadData = .TRUE.
      end if
      if (GdataDescrip(i)%ny .ne. Ny) then
        write (*,*) 'ERROR: number of y points in GRADS control files do not match'
        BadData = .TRUE.
      end if
      if (GdataDescrip(i)%nz .ne. Nz) then
        write (*,*) 'ERROR: number of z points in GRADS control files do not match'
        BadData = .TRUE.
      end if
      if (GdataDescrip(i)%nt .ne. Nt) then
        write (*,*) 'ERROR: number of t points in GRADS control files do not match'
        BadData = .TRUE.
      end if

      Nvars = Nvars + GdataDescrip(i)%nvars
    end if

    do j =1, GdataDescrip(i)%nvars
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
    BadData = .TRUE.
  end if
  if (PressLoc%Fnum .eq. 0) then
    write (*,*) 'ERROR: cannot find grid var "press" in the GRADS data files'
    BadData = .TRUE.
  end if
  if (VarLoc%Fnum .eq. 0) then
    write (*,*) 'ERROR: cannot find grid var "', trim(VarName), '" in the GRADS data files'
    BadData = .TRUE.
  end if

  if (BadData) then
    stop
  end if
end subroutine
