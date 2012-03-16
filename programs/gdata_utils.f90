!***************************************************************
! Data types for GRADS data routines
!

module gdata_utils
!**********************************************************
! DATA
!**********************************************************

  integer, parameter :: MaxString=256
  integer, parameter :: MaxCoords=1000
  integer, parameter :: MaxVars=50
  integer, parameter :: MaxFiles=10
  integer, parameter :: InUnit=8
  integer, parameter :: OutUnit=10
  integer, parameter :: BinRecFactor=4
  
  type GradsDataDescription
    character (len=MaxString) :: DataFile
    real :: UndefVal
    real :: Xstart
    real :: Xinc
    real :: Ystart
    real :: Yinc
    real, dimension(1:MaxCoords) :: Zlevels
    character (len=MaxString) :: Tstart
    character (len=MaxString) :: Tinc
    character (len=MaxString), dimension(1:MaxVars) :: VarNames
    integer, dimension(1:MaxVars) :: VarNz
    integer :: Nx
    integer :: Ny
    integer :: Nz
    integer :: Nt
    integer :: Nvars
  end type GradsDataDescription

  type GradsControlFiles
    integer :: Nfiles
    character (len=MaxString), dimension(1:MaxFiles) :: Fnames
    type (GradsDataDescription), dimension(1:MaxFiles) :: DataDscrs
  end type GradsControlFiles

  type GradsVar
    character (len=MaxString) :: Vname
    integer :: Nx
    integer :: Ny
    integer :: Nz
    integer :: Nt
    real, dimension(:), allocatable :: Xcoords
    real, dimension(:), allocatable :: Ycoords
    real, dimension(:), allocatable :: Zcoords
    real, dimension(:), allocatable :: Tcoords
    real, dimension(:,:,:,:), allocatable :: Vdata
    ! Followin vars are for use by grads input
    ! routines. Vnum is the variable number (index)
    ! within the grads data file (DataFile)
    character (len=MaxString) :: DataFile
    integer :: Vnum
    integer :: Nvars
    ! Following vars are carried along for use by
    ! grads output routines. This is done since not
    ! all vars of this type are associated with an
    ! input control file.
    real :: Xstart
    real :: Xinc
    real :: Ystart
    real :: Yinc
    character (len=MaxString) :: Tstart
    character (len=MaxString) :: Tinc
    real :: UndefVal
  end type GradsVar

contains
!**********************************************************
! SUBROUTINES
!**********************************************************

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
          write (0,*) 'ERROR: Maximum number of ', ItemType, ' (', MaxItems, ') exceeded'
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
      write (0,*) 'ERROR: Maximum number of ', ItemType, ' (', MaxItems, ') exceeded'
      stop
    end if
  end if

  return
end subroutine

!**********************************************************
! ReadGradsCtlFiles()
!
! This routine will read in the GRADS control files

subroutine  ReadGradsCtlFiles(GctlFiles)
  implicit none

  integer, parameter :: MaxLine = 1000
  integer, parameter :: MaxFields = 100

  type (GradsControlFiles) :: GctlFiles

  integer :: InError
  character (len=MaxLine) :: InLine
  character (len=MaxLine), dimension(1:MaxFields) :: InFields
  integer :: i, iz
  integer :: Nfields

  character (len=MaxString) :: GfileDir

  do i = 1, GctlFiles%Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GctlFiles%Fnames(i))

    open(unit=InUnit, file=GctlFiles%Fnames(i), form='formatted', &
         status='old', action='read', iostat=InError)
    if (InError .ne. 0) then
      write (0,*) 'ERROR: cannot open GRADS control file: ', trim(GctlFiles%Fnames(i))
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
          if (index(GctlFiles%Fnames(i), '/', .true.) .eq. 0) then
            ! no leading path
            GfileDir = ''
          else
            GfileDir = GctlFiles%Fnames(i)(1:index(GctlFiles%Fnames(i), '/', .true.))
          end if
          ! Make sure to use trim(), without it the first string GfileDir fills up the
          ! GctlFiles%DataDscrs(i)%DataFile string since these are both declared to be the same length
          GctlFiles%DataDscrs(i)%DataFile = trim(GfileDir) // trim(InFields(2)(2:len(InFields(2))))
        else
          GctlFiles%DataDscrs(i)%DataFile = InFields(2)
        end if
      else
        if ((InFields(1) .eq. 'undef') .or. (InFields(1) .eq. 'UNDEF')) then
          read (InFields(2), '(f)') GctlFiles%DataDscrs(i)%UndefVal
        else if ((InFields(1) .eq. 'xdef') .or. (InFields(1) .eq. 'XDEF')) then
          read (InFields(2), '(i)') GctlFiles%DataDscrs(i)%Nx 
          read (InFields(4), '(f)') GctlFiles%DataDscrs(i)%Xstart
          read (InFields(5), '(f)') GctlFiles%DataDscrs(i)%Xinc
        else if ((InFields(1) .eq. 'ydef') .or. (InFields(1) .eq. 'YDEF')) then
          read (InFields(2), '(i)') GctlFiles%DataDscrs(i)%Ny 
          read (InFields(4), '(f)') GctlFiles%DataDscrs(i)%Ystart
          read (InFields(5), '(f)') GctlFiles%DataDscrs(i)%Yinc
        else if ((InFields(1) .eq. 'zdef') .or. (InFields(1) .eq. 'ZDEF')) then
          read (InFields(2), '(i)') GctlFiles%DataDscrs(i)%Nz
          do iz = 1, GctlFiles%DataDscrs(i)%Nz
            read(InFields(iz+3), '(f)') GctlFiles%DataDscrs(i)%Zlevels(iz)
          enddo
        else if ((InFields(1) .eq. 'tdef') .or. (InFields(1) .eq. 'TDEF')) then
          read (InFields(2), '(i)') GctlFiles%DataDscrs(i)%Nt
          GctlFiles%DataDscrs(i)%Tstart = InFields(4)
          GctlFiles%DataDscrs(i)%Tinc   = InFields(5)
        else if ((InFields(1) .eq. 'vars') .or. (InFields(1) .eq. 'VARS')) then
          read (InFields(2), '(i)') GctlFiles%DataDscrs(i)%Nvars
          if (GctlFiles%DataDscrs(i)%Nvars .gt. MaxVars) then
            write (0,*) 'ERROR: Maximum number of variables (', MaxVars, &
                        ') exceeded in GRADS control file specs'
            stop
          endif
          call GetVarInfo(InUnit, GctlFiles%DataDscrs(i)%VarNames, GctlFiles%DataDscrs(i)%VarNz, GctlFiles%DataDscrs(i)%Nvars, MaxLine, MaxFields)
        endif
      endif
    enddo
  
    close(unit=InUnit, status='keep')
  
    if ((GctlFiles%DataDscrs(i)%Nx .gt. MaxCoords) .or. (GctlFiles%DataDscrs(i)%Ny .gt. MaxCoords) .or. &
        (GctlFiles%DataDscrs(i)%Nz .gt. MaxCoords) .or. (GctlFiles%DataDscrs(i)%Nt .gt. MaxCoords)) then
      write (0,*) 'ERROR: Maximum number of coordinates (', MaxCoords, ') exceeded in GRADS control file: ', trim(GctlFiles%Fnames(i))
      stop
    endif
  enddo
  write (*,*) ''

  return
end subroutine

!*****************************************************************************
! GetVarInfo()
!
! This routine will read Nvars lines from the file (GRADS control file) given
! by InUnit, and record the names and number of Z levels of the variables
! listed in the file.

subroutine  GetVarInfo(InUnit, VarNames, VarNz, Nvars, MaxLine, MaxFields)
  implicit none

  integer InUnit
  character (len=*), dimension(1:Nvars) :: VarNames
  integer, dimension(1:Nvars) :: VarNz
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
      write (0,*) 'ERROR: Reached end of GRADS control file before reading in all (', &
                  Nvars, ') variable names'
      stop
    end if

    call String2List(InLine, ' ', InFields, MaxFields, Nfields, 'control file fields')

    VarNames(i) = InFields(1)
    read(InFields(2), '(i)') VarNz(i)
  end do

  return
end subroutine

!**********************************************************************
! InitGvarFromGdescrip()
!
! This routine will find the GRADS variable description in the set of
! GRADS control file information. Once found the variable description
! is used to fill in all of the data into the GRADS variable except
! for the actual 3D or 2D data field which is deferred so that all
! of the GRADS variables can be checked for consistency before trying
! to read in the data field (can be time consuming).

subroutine InitGvarFromGdescrip(GctlFiles, Gvar, Vname)
  implicit none

  type (GradsControlFiles) :: GctlFiles
  type (GradsVar) :: Gvar
  character (len=*) :: Vname

  ! i -> file number, j -> var number
  integer :: i, j
  integer :: Fnum, Vnum

  Fnum = 0
  Vnum = 0

  do i = 1, GctlFiles%Nfiles
    do j =1, GctlFiles%DataDscrs(i)%Nvars
      if (GctlFiles%DataDscrs(i)%VarNames(j) .eq. Vname) then
        Fnum = i
        Vnum = j
      endif
    enddo
  enddo

  ! If we found the variable, then fill in the dimensions, coordinate data, and allocate
  ! the space for the data.
  if (Fnum .gt. 0) then
    call InitGradsVar(Gvar, Vname, GctlFiles%DataDscrs(Fnum)%Nx, GctlFiles%DataDscrs(Fnum)%Ny, &
                      GctlFiles%DataDscrs(Fnum)%Nz, GctlFiles%DataDscrs(Fnum)%Nt, &
                      GctlFiles%DataDscrs(Fnum)%Xstart, GctlFiles%DataDscrs(Fnum)%Xinc, &
                      GctlFiles%DataDscrs(Fnum)%Ystart, GctlFiles%DataDscrs(Fnum)%Yinc, &
                      GctlFiles%DataDscrs(Fnum)%Zlevels, &
                      GctlFiles%DataDscrs(Fnum)%Tstart, GctlFiles%DataDscrs(Fnum)%Tinc, &
                      GctlFiles%DataDscrs(Fnum)%UndefVal, &
                      GctlFiles%DataDscrs(Fnum)%DataFile, Vnum, GctlFiles%DataDscrs(Fnum)%Nvars)
  else
    write (0,*) 'ERROR: cannot find grid var "', trim(Vname), '" in the GRADS data files'
    stop
  endif

  return
end subroutine

!**********************************************************************
! InitGvarFromGvar()
!
! This routine will copy the set up from GvarSrc into GvarDes.
! Pretty light weight but makes for convenient use from the caller.

subroutine InitGvarFromGvar(GvarSrc, GvarDest, Vname)
  implicit none

  type (GradsVar) :: GvarSrc, GvarDest
  character (len=*) :: Vname

  ! Since there is no input GRADS data file associated with this new variable, set the
  ! DataFile, Vnum and Nvars specs to '<NONE>', 0 and 0 respectively
  call InitGradsVar(GvarDest, Vname, GvarSrc%Nx, GvarSrc%Ny, GvarSrc%Nz, GvarSrc%Nt, &
                    GvarSrc%Xstart, GvarSrc%Xinc, GvarSrc%Ystart, GvarSrc%Yinc, &
                    GvarSrc%Zcoords, GvarSrc%Tstart, GvarSrc%Tinc, &
                    GvarSrc%UndefVal, '<NONE>', 0, 0)

  return
end subroutine

!*******************************************************************
! InitGradsVar()
!
! Generic routine to initialize a GradsVar type.
! All members are filled in except for Vdata which will get filled
! in by various means (read from GRADS data file, convert u,v data
! to tangential wind, etc.) so for now just allocate the space for Vdata.
subroutine InitGradsVar(Gvar, Vname, Nx, Ny, Nz, Nt, Xstart, Xinc, Ystart, Yinc, Zcoords, &
                        Tstart, Tinc, UndefVal, DataFile, Vnum, Nvars)
  implicit none

  type (GradsVar) :: Gvar
  character (len=*) :: Vname
  integer :: Nx
  integer :: Ny
  integer :: Nz
  integer :: Nt
  real :: Xstart
  real :: Xinc
  real :: Ystart
  real :: Yinc
  real, dimension(1:Nz) :: Zcoords
  character (len=*) :: Tstart
  character (len=*) :: Tinc
  real :: UndefVal
  character (len=*) :: DataFile
  integer :: Vnum
  integer :: Nvars

  integer :: i
  integer :: Ierror

  ! Copy in specs, allocate the coord and data arrays, fill in the coord array using
  ! the specs
  Gvar%Vname = Vname
  
  Gvar%Nx = Nx
  Gvar%Ny = Ny
  Gvar%Nz = Nz
  Gvar%Nt = Nt

  Gvar%Xstart = Xstart
  Gvar%Xinc = Xinc
  Gvar%Ystart = Ystart
  Gvar%Yinc = Yinc
  Gvar%Tstart = Tstart
  Gvar%Tinc = Tinc
  Gvar%UndefVal = UndefVal

  Gvar%DataFile = DataFile
  Gvar%Vnum = Vnum
  Gvar%Nvars = Nvars

  allocate(Gvar%Xcoords(1:Nx), Gvar%Ycoords(1:Ny), Gvar%Zcoords(1:Nz), Gvar%Tcoords(1:Nt), &
           Gvar%Vdata(1:Nt,1:Nz,1:Nx,1:Ny), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Coordinate, data array memory allocation failed'
    stop
  endif

  ! fill in the coordinate arrays
  do i = 1, Nx
    Gvar%Xcoords(i) = Xstart + (real(i - 1) * Xinc)
  enddo
  do i = 1, Ny
    Gvar%Ycoords(i) = Ystart + (real(i - 1) * Yinc)
  enddo
  do i = 1, Nz
    Gvar%Zcoords(i) = Zcoords(i)
  enddo
  do i = 1, Nt
    Gvar%Tcoords(i) = real(i)
  enddo

  return
end subroutine

!***********************************************************************
! GvarDimsMatch()
!
! This function will return true/false according to whether or not the
! dimensions of the two given GRADS variables match. When the variables
! are 3D then x, y, z, and t dimensions are checked. If one of the variables
! is 2D then only x, y, and t dimensions are checked.

logical function GvarDimsMatch(Gvar1, Gvar2, Is2D)
  implicit none

  type (GradsVar) :: Gvar1, Gvar2
  logical :: Is2D

  logical :: GcoordsMatch

  if (Is2D) then
    ! 2D var: check that the first entry in the Z coordinates match
    GvarDimsMatch = GcoordsMatch(Gvar1%Xcoords, Gvar1%Nx, Gvar2%Xcoords, Gvar2%Nx) .and. &
                    GcoordsMatch(Gvar1%Ycoords, Gvar1%Ny, Gvar2%Ycoords, Gvar2%Ny) .and. &
                    GcoordsMatch(Gvar1%Zcoords, 1,        Gvar2%Zcoords, 1       ) .and. &
                    GcoordsMatch(Gvar1%Tcoords, Gvar1%Nt, Gvar2%Tcoords, Gvar2%Nt)
  else
    ! 3D var: check that all entries in the X coordinates match
    GvarDimsMatch = GcoordsMatch(Gvar1%Xcoords, Gvar1%Nx, Gvar2%Xcoords, Gvar2%Nx) .and. &
                    GcoordsMatch(Gvar1%Ycoords, Gvar1%Ny, Gvar2%Ycoords, Gvar2%Ny) .and. &
                    GcoordsMatch(Gvar1%Zcoords, Gvar1%Nz, Gvar2%Zcoords, Gvar2%Nz) .and. &
                    GcoordsMatch(Gvar1%Tcoords, Gvar1%Nt, Gvar2%Tcoords, Gvar2%Nt)
  endif

  return
end function

!***********************************************************************
! GcoodsMatch()
!
! This function will return true/false according to whether or not the
! coordinate values in the two given arrays match.

logical function GcoordsMatch(Coords1, Nc1, Coords2, Nc2)
  implicit none

  integer :: Nc1, Nc2
  real, dimension(1:Nc1) :: Coords1
  real, dimension(1:Nc2) :: Coords2

  integer :: i

  if (Nc1 .eq. Nc2) then
    ! lengths of arrays match
    ! check if entries in the arrays match
    GcoordsMatch = .true.
    do i = 1, Nc1
      if (abs(Coords1(i)-Coords2(i)) .gt. 1.0e-6) then
        GcoordsMatch = .false.
        exit
      endif
    enddo
  else
    ! lengths of arrays do not match
    GcoordsMatch = .false.
  endif

  return
end function

!********************************************************************
! ReadGradsData()
!
! This routine will use the information in GradsVar to locate the
! data for the var in the GRADS data files. This data will then be
! loaded into the Vdata array.
!

subroutine ReadGradsData(Gvar)
  implicit none

  type (GradsVar) :: Gvar

  integer RecLen, Ierror
  integer ix,iy,iz,it
  integer RecNum
  integer NumRecs

  ! Each record is a horizontal slice so it has Nx * Ny elements in it
  ! of BinRecFactor size
  RecLen = Gvar%Nx * Gvar%Ny * BinRecFactor

  open (unit=InUnit, file=Gvar%DataFile, form='unformatted', access='direct', &
        recl=RecLen, status='old', action='read', iostat=Ierror)
  if (Ierror .ne. 0) then
    write (0,*) 'ERROR: Cannot open GRADS data file: ', trim(Gvar%DataFile)
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

  NumRecs = Gvar%Nz * Gvar%Nt * Gvar%Nvars

  write (0,*) 'Reading GRADS data file: ', trim(Gvar%DataFile)
  write (0,*) '  Var: ', trim(Gvar%Vname)
  write (0,*) '  X points:          ', Gvar%Nx
  write (0,*) '  Y points:          ', Gvar%Ny
  write (0,*) '  Z points:          ', Gvar%Nz
  write (0,*) '  T points:          ', Gvar%Nt
  write (0,*) '  Number of Vars:    ', Gvar%Nvars
  write (0,*) '  Record length:     ', RecLen
  write (0,*) '  Number of records: ', NumRecs
  write (0,*) ''

  ! Need to derive the correct record number in the file
  ! from the it and iz values. it cycles by Nz*Nvars records.
  ! Within each it cycle, Nz*(Vnum-1) is the offset from the
  ! start of the cycle to the Vnum position.
  ! The formula becomes:
  !
  !  RecNum = (it-1)*(Nz*Nvars) + (Nz*(Vnum-1)) + iz
  
  NumRecs = 0
  do it = 1, Gvar%Nt
    do iz = 1, Gvar%Nz
      RecNum = (it-1)*(Gvar%Nz*Gvar%Nvars) + (Gvar%Nz*(Gvar%Vnum-1)) + iz
      read(unit=InUnit, rec=RecNum) ((Gvar%Vdata(it,iz,ix,iy), ix = 1, Gvar%Nx), iy = 1, Gvar%Ny)

      NumRecs = NumRecs + 1
      if (modulo(NumRecs,100) .eq. 0) then
        write (0,*) '  Reading: ', trim(Gvar%Vname), ': ', NumRecs
      end if
    end do   
  end do 
  write (0,*) '  Total records read: ', NumRecs
  write (0,*) ''

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

subroutine WriteGrads(Gvar, OfileBase, Diagnostic)
  implicit none

  type (GradsVar) :: Gvar
  character (len=*) :: OfileBase
  character (len=*) :: Diagnostic

  integer :: OutRecLen
  integer :: RecNum
  integer :: Ierror
  integer :: ix
  integer :: iy
  integer :: iz
  integer :: it

  character (len=MaxString) :: CtlFile
  character (len=MaxString) :: DataFile
  character (len=MaxString) :: Title

  CtlFile = trim(OfileBase) // '.ctl'
  DataFile = trim(OfileBase) // '.gra'
  Title = 'Diagnostic: ' // trim(Diagnostic)

  write (0,*) 'Writing out result in GRADS format:'
  write (0,*) '  Control file: ', trim(CtlFile)
  write (0,*) '  Data file:    ', trim(DataFile)
  write (0,*) '  Total number of data points: ', &
              Gvar%Nx * Gvar%Ny * Gvar%Nz * Gvar%Nt
  write (0,*) ''

  ! Control (data description) file
  open (OutUnit, file=CtlFile, form='formatted', action='write', status='replace', iostat=Ierror)
  if (Ierror .ne. 0) then
    write (0,*) 'ERROR: cannot open output GRADS control file for writing: ', trim(CtlFile)
    stop
  end if

  write (OutUnit, '(a,a)')          'DSET ^./', trim(DataFile)
  write (OutUnit, '(a,a)')          'TITLE ', trim(Title)
  write (OutUnit, '(a,g)')          'UNDEF ', Gvar%UndefVal
  write (OutUnit, '(a)')            'OPTIONS  little_endian'
  write (OutUnit, '(a,i,a,g,g)')    'XDEF ', Gvar%Nx, ' LINEAR ', Gvar%Xstart, Gvar%Xinc
  write (OutUnit, '(a,i,a,g,g)')    'YDEF ', Gvar%Ny, ' LINEAR ', Gvar%Ystart, Gvar%Yinc
  write (OutUnit, '(a,i,a,100g)')   'ZDEF ', Gvar%Nz, ' LEVELS ', &
                                    (Gvar%Zcoords(iz), iz = 1, Gvar%Nz)
  write (OutUnit, '(a,i,a,a,5x,a)') 'TDEF ', Gvar%Nt, ' LINEAR ', &
                                    trim(Gvar%Tstart), trim(Gvar%Tinc)
  write (OutUnit, '(a)')            'VARS 1 '
  write (OutUnit, '(a,i,a)')        trim(Gvar%Vname) , Gvar%Nz, ' 99 Diagnostic'
  write (OutUnit, '(a)')            'ENDVARS'

  close (OutUnit, status='keep')

  ! Data file
  ! dimensions from fastest changing to slowest changing are: x, y, z, t
  ! One record is a single horizontal slice -> (nx * ny) data points
  OutRecLen = Gvar%Nx * Gvar%Ny * BinRecFactor
  open (OutUnit, file=DataFile, form='unformatted', access='direct', &
        recl=OutRecLen, action='write', status='replace', iostat=Ierror)
  if (Ierror .ne. 0) then
    write (0,*) 'ERROR: cannot open output GRADS data file for writing: ', trim(DataFile)
    stop
  end if

  RecNum = 1
  do it = 1, Gvar%Nt
    do iz = 1, Gvar%Nz
      if (modulo(RecNum,100) .eq. 0) then
        write (0,*) '  Writing: ', trim(Gvar%Vname), ': ', RecNum
      end if
      write (OutUnit, rec=RecNum) ((Gvar%Vdata(it,iz,ix,iy), ix = 1, Gvar%Nx), iy = 1, Gvar%Ny)
      RecNum = RecNum + 1
    end do
  end do

  close (OutUnit, status='keep')

  return
end subroutine

!***********************************************************************
! SetZmHeights()
!
! This routine will set the Zm heights which are used for defining
! the layer thicknesses in the column integration calculations.
! Zm grid points are the heights at the vector data points (whereas Zt
! are the heights of the scalar data points) - staggered grid scheme in RAMS.
! Zm heights are defined by the RAMSIN vars: DELTAZ, DZRAT, DZMAX.
!

subroutine SetZmHeights(Nz, ZmHeights)
  implicit none

  real, dimension(0:Nz) :: ZmHeights
  real :: DeltaZ, DzRat, DzMax
  integer :: Nz, iz

  ! all the simulations are using the same z scheme so for now just
  ! use the settings
  ! WARNING: this needs to change is the z scheme in the simulations changes

  DeltaZ = 300.0
  DzRat = 1.065
  DzMax = 1000.0

  ZmHeights(0) = 0.0
  do iz = 1, Nz
    ZmHeights(iz) = ZmHeights(iz - 1) + DeltaZ
    DeltaZ = DeltaZ * DzRat
    if (DeltaZ .ge. DzMax) then
      DeltaZ = DzMax
    endif
  enddo

  return
end subroutine

end module gdata_utils
