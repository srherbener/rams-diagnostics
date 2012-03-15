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
  integer, parameter :: InUnit=8
  integer, parameter :: OutUnit=10
  integer, parameter :: BinRecFactor=4
  
  type GradsDataDescription
    character (len=MaxString) :: DataFile
    real :: UndefVal
    real :: Xstart, Xinc
    real :: Ystart, Yinc
    real, dimension(1:MaxCoords) :: Zlevels
    character (len=MaxString) :: Tstart, Tinc
    character (len=MaxString), dimension(1:MaxVars) :: VarNames
    integer, dimension(1:MaxVars) :: VarNz
    integer :: Nx
    integer :: Ny
    integer :: Nz
    integer :: Nt
    integer :: Nvars
  end type GradsDataDescription

  type GradsVar
    ! Fnum, Vnum define where the variable description lives within
    !   the set of input GRADS control files
    !
    !   Fnum - File Number
    !   Vnum - Variable Number
    character (len=MaxString) :: Vname
    integer :: Nx
    integer :: Ny
    integer :: Nz
    integer :: Nt
    real, dimension(:), allocatable :: Xcoords
    real, dimension(:), allocatable :: Ycoords
    real, dimension(:), allocatable :: Zcoords
    real, dimension(:), allocatable :: Tcoords
    integer :: Fnum
    integer :: Vnum
    character (len=MaxString) :: Tstart, Tinc
    character (len=MaxString) :: DataFile
    real, dimension(:,:,:,:), allocatable :: Vdata
  end type GradsVar

  type GradsOutDescription
    character (len=MaxString) :: CtlFile
    character (len=MaxString) :: DataFile
    character (len=MaxString) :: Title
    character (len=MaxString) :: VarName
    real :: UndefVal
    real :: Xstart, Xinc
    real :: Ystart, Yinc
    real, dimension(1:MaxCoords) :: Zlevels
    character (len=MaxString) :: Tstart, Tinc
    integer :: Nx
    integer :: Ny
    integer :: Nz
    integer :: Nt
  end type GradsOutDescription

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
! ReadGradsCtlFile()
!
! This routine will read in the GRADS control files
!
! Args
!   1. Name of grads control file
!   2. Grads data description
!

subroutine  ReadGradsCtlFile(GradsCtlFile, GdataDescrip, MaxCoords)
  implicit none

  integer, parameter :: MaxLine = 1000
  integer, parameter :: MaxFields = 100

  character (len=*) :: GradsCtlFile
  type (GradsDataDescription) GdataDescrip
  integer :: MaxCoords

  integer :: InError
  character (len=MaxLine) :: InLine
  character (len=MaxLine), dimension(1:MaxFields) :: InFields
  integer :: i
  integer :: Nfields

  character (len=MaxString) :: GfileDir

  open(unit=InUnit, file=GradsCtlFile, form='formatted', &
       status='old', action='read', iostat=InError)
  if (InError .ne. 0) then
    write (0,*) 'ERROR: cannot open GRADS control file: ', trim(GradsCtlFile)
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
      else if ((InFields(1) .eq. 'xdef') .or. (InFields(1) .eq. 'XDEF')) then
        read (InFields(2), '(i)') GdataDescrip%Nx 
        read (InFields(4), '(f)') GdataDescrip%Xstart
        read (InFields(5), '(f)') GdataDescrip%Xinc
      else if ((InFields(1) .eq. 'ydef') .or. (InFields(1) .eq. 'YDEF')) then
        read (InFields(2), '(i)') GdataDescrip%Ny 
        read (InFields(4), '(f)') GdataDescrip%Ystart
        read (InFields(5), '(f)') GdataDescrip%Yinc
      else if ((InFields(1) .eq. 'zdef') .or. (InFields(1) .eq. 'ZDEF')) then
        read (InFields(2), '(i)') GdataDescrip%Nz
        do i = 1, GdataDescrip%Nz
          read(InFields(i+3), '(f)') GdataDescrip%Zlevels(i)
        enddo
      else if ((InFields(1) .eq. 'tdef') .or. (InFields(1) .eq. 'TDEF')) then
        read (InFields(2), '(i)') GdataDescrip%Nt
        GdataDescrip%Tstart = InFields(4)
        GdataDescrip%Tinc   = InFields(5)
      else if ((InFields(1) .eq. 'vars') .or. (InFields(1) .eq. 'VARS')) then
        read (InFields(2), '(i)') GdataDescrip%Nvars
        if (GdataDescrip%Nvars .gt. MaxVars) then
          write (0,*) 'ERROR: Maximum number of variables (', MaxVars, &
                      ') exceeded in GRADS control file specs'
          stop
        endif
        call GetVarInfo(InUnit, GdataDescrip%VarNames, GdataDescrip%VarNz, GdataDescrip%Nvars, MaxLine, MaxFields)
      endif
    endif
  enddo

  close(unit=InUnit, status='keep')

  if ((GdataDescrip%Nx .gt. MaxCoords) .or. (GdataDescrip%Ny .gt. MaxCoords) .or. &
      (GdataDescrip%Nz .gt. MaxCoords) .or. (GdataDescrip%Nt .gt. MaxCoords)) then
    write (0,*) 'ERROR: Maximum number of coordinates (', MaxCoords, ') exceeded in GRADS control file: ', trim(GradsCtlFile)
    stop
  end if

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

subroutine InitGvarFromGdescrip(Nf, GdataDescrip, Gvar, Vname)
  implicit none

  integer :: Nf
  type (GradsDataDescription), dimension(1:Nf) :: GdataDescrip
  type (GradsVar) :: Gvar
  character (len=*) :: Vname

  ! i -> file number, j -> var number
  integer :: i, j
  integer :: Ierror

  Gvar%Fnum = 0
  Gvar%Vnum = 0

  do i = 1, Nf
    do j =1, GdataDescrip(i)%Nvars
      if (GdataDescrip(i)%VarNames(j) .eq. Vname) then
        Gvar%Fnum = i
        Gvar%Vnum = j
      endif
    enddo
  enddo

  ! If we found the variable, then fill in the dimensions, coordinate data, and allocate
  ! the space for the data.
  if (Gvar%Fnum .gt. 0) then
    Gvar%Vname = Vname
    Gvar%Tstart = GdataDescrip(Gvar%Fnum)%Tstart
    Gvar%Tinc = GdataDescrip(Gvar%Fnum)%Tinc
    Gvar%DataFile = GdataDescrip(Gvar%Fnum)%DataFile

    Gvar%Nx = GdataDescrip(Gvar%Fnum)%Nx;
    Gvar%Ny = GdataDescrip(Gvar%Fnum)%Ny;
    Gvar%Nz = GdataDescrip(Gvar%Fnum)%Nz;
    Gvar%Nt = GdataDescrip(Gvar%Fnum)%Nt;

    allocate(Gvar%Xcoords(1:Gvar%Nx), Gvar%Ycoords(1:Gvar%Ny), &
             Gvar%Zcoords(1:Gvar%Nz), Gvar%Tcoords(1:Gvar%Nt), &
             Gvar%Vdata(1:Gvar%Nt,1:Gvar%Nz,1:Gvar%Nx,1:Gvar%Ny), stat=Ierror)
    if (Ierror .ne. 0) then
      write (*,*) 'ERROR: Data array memory allocation failed'
      stop
    endif

    ! fill in the coordinate arrays
    do i = 1, Gvar%Nx
      Gvar%Xcoords(i) = GdataDescrip(Gvar%Fnum)%Xstart + (real(i - 1) * GdataDescrip(Gvar%Fnum)%Xinc)
    enddo
    do i = 1, Gvar%Ny
      Gvar%Ycoords(i) = GdataDescrip(Gvar%Fnum)%Ystart + (real(i - 1) * GdataDescrip(Gvar%Fnum)%Yinc)
    enddo
    do i = 1, Gvar%Nz
      Gvar%Zcoords(i) = GdataDescrip(Gvar%Fnum)%Zlevels(i)
    enddo
    do i = 1, Gvar%Nt
      Gvar%Tcoords(i) = real(i)
    enddo
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
! All of the data is filled into the GRADS variable except
! for the actual 3D or 2D data field which is deferred so that all
! of the GRADS variables can be checked for consistency before trying
! to read in the data field (can be time consuming).

subroutine InitGvarFromGvar(GvarSrc, GvarDest, Vname)
  implicit none

  type (GradsVar) :: GvarSrc, GvarDest
  character (len=*) :: Vname

  integer :: i
  integer :: Ierror

  GvarDest%Fnum = 0
  GvarDest%Vnum = 0

  GvarDest%Vname = Vname
  GvarDest%Tstart = GvarSrc%Tstart
  GvarDest%Tinc = GvarSrc%Tinc
  GvarDest%DataFile = '<NONE>'

  GvarDest%Nx = GvarSrc%Nx;
  GvarDest%Ny = GvarSrc%Ny;
  GvarDest%Nz = GvarSrc%Nz;
  GvarDest%Nt = GvarSrc%Nt;

  allocate(GvarDest%Xcoords(1:GvarDest%Nx), GvarDest%Ycoords(1:GvarDest%Ny), &
           GvarDest%Zcoords(1:GvarDest%Nz), GvarDest%Tcoords(1:GvarDest%Nt), &
           GvarDest%Vdata(1:GvarDest%Nt,1:GvarDest%Nz,1:GvarDest%Nx,1:GvarDest%Ny), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  endif

  ! fill in the coordinate arrays
  do i = 1, GvarDest%Nx
    GvarDest%Xcoords(i) = GvarSrc%Xcoords(i)
  enddo
  do i = 1, GvarDest%Ny
    GvarDest%Ycoords(i) = GvarSrc%Ycoords(i)
  enddo
  do i = 1, GvarDest%Nz
    GvarDest%Zcoords(i) = GvarSrc%Zcoords(i)
  enddo
  do i = 1, GvarDest%Nt
    GvarDest%Tcoords(i) = GvarSrc%Tcoords(i)
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

subroutine ReadGradsData(GdataDescrip, Gvar)
  implicit none

  type (GradsDataDescription) :: GdataDescrip
  type (GradsVar) :: Gvar

  integer RecLen, Ierror
  integer ix,iy,iz,it
  integer RecNum
  integer NumRecs

  ! Each record is a horizontal slice so it has Nx * Ny elements in it
  ! of BinRecFactor size
  RecLen = Gvar%Nx * Gvar%Ny * BinRecFactor

  open (unit=InUnit, file=GdataDescrip%DataFile, form='unformatted', access='direct', &
        recl=RecLen, status='old', action='read', iostat=Ierror)
  if (Ierror .ne. 0) then
    write (0,*) 'ERROR: Cannot open GRADS data file: ', trim(GdataDescrip%DataFile)
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

  NumRecs = Gvar%Nz * Gvar%Nt * GdataDescrip%Nvars

  write (0,*) 'Reading GRADS data file: ', trim(GdataDescrip%DataFile)
  write (0,*) '  Var: ', trim(Gvar%Vname)
  write (0,*) '  X points:          ', Gvar%Nx
  write (0,*) '  Y points:          ', Gvar%Ny
  write (0,*) '  Z points:          ', Gvar%Nz
  write (0,*) '  T points:          ', Gvar%Nt
  write (0,*) '  Number of Vars:    ', GdataDescrip%Nvars
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
      RecNum = (it-1)*(Gvar%Nz*GdataDescrip%Nvars) + (Gvar%Nz*(Gvar%Vnum-1)) + iz
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

!******************************************************************************
! BuildGoutDescrip()
!
! This routine will fill in the GoutDescrip record for the AzAvg data. This
! record is used for creating the GRADS control file.
!

subroutine BuildGoutDescrip(Nx, Ny, Nz, Nt, AzAvg, OfileBase, UndefVal, VarName, &
          Xstart, Xinc, Ystart, Yinc, Zlevels, Tstart, Tinc, GoutDescrip, Diagnostic)

  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: AzAvg
  character (len=*) :: OfileBase, VarName
  real :: UndefVal, Xstart, Xinc, Ystart, Yinc
  real, dimension(1:Nz) :: Zlevels
  character (len=*) :: Tstart, Tinc
  type (GradsOutDescription) :: GoutDescrip
  character (len=*) :: Diagnostic

  integer :: iz

  GoutDescrip%CtlFile = trim(OfileBase) // '.ctl'
  GoutDescrip%DataFile = trim(OfileBase) // '.gra'
  GoutDescrip%Title = 'Diagnostic: ' // trim(Diagnostic)
  GoutDescrip%UndefVal = UndefVal
  GoutDescrip%Nx = Nx
  GoutDescrip%Ny = Ny
  GoutDescrip%Nz = Nz
  GoutDescrip%Nt = Nt
  GoutDescrip%Xstart = Xstart
  GoutDescrip%Xinc = Xinc
  GoutDescrip%Ystart = Ystart
  GoutDescrip%Yinc = Yinc
  GoutDescrip%VarName = trim(VarName)
  do iz = 1, Nz
    GoutDescrip%Zlevels(iz) = Zlevels(iz)
  end do
  GoutDescrip%Tstart = Tstart
  GoutDescrip%Tinc = Tinc

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
  implicit none

  type (GradsOutDescription) :: GoutDescrip
  real, dimension(1:GoutDescrip%Nx,1:GoutDescrip%Ny,1:GoutDescrip%Nz,1:GoutDescrip%Nt) :: AzAvg
  integer :: Nx, Ny, Nz, Nt

  integer :: OutRecLen
  integer :: RecNum
  integer :: Ierror
  integer :: ix
  integer :: iy
  integer :: iz
  integer :: it

  write (0,*) 'Writing out result in GRADS format:'
  write (0,*) '  Control file: ', trim(GoutDescrip%CtlFile)
  write (0,*) '  Data file:    ', trim(GoutDescrip%DataFile)
  write (0,*) '  Total number of data points: ', &
              GoutDescrip%Nx * GoutDescrip%Ny * GoutDescrip%Nz * GoutDescrip%Nt
  write (0,*) ''

  ! Control (data description) file
  open (OutUnit, file=GoutDescrip%CtlFile, form='formatted', action='write', status='replace', iostat=Ierror)
  if (Ierror .ne. 0) then
    write (0,*) 'ERROR: cannot open output GRADS control file for writing: ', trim(GoutDescrip%CtlFile)
    stop
  end if

  write (OutUnit, '(a,a)')          'DSET ^./', trim(GoutDescrip%DataFile)
  write (OutUnit, '(a,a)')          'TITLE ', trim(GoutDescrip%Title)
  write (OutUnit, '(a,g)')          'UNDEF ', GoutDescrip%UndefVal
  write (OutUnit, '(a)')            'OPTIONS  little_endian'
  write (OutUnit, '(a,i,a,g,g)')    'XDEF ', GoutDescrip%Nx, ' LINEAR ', GoutDescrip%Xstart, GoutDescrip%Xinc
  write (OutUnit, '(a,i,a,g,g)')    'YDEF ', GoutDescrip%Ny, ' LINEAR ', GoutDescrip%Ystart, GoutDescrip%Yinc
  write (OutUnit, '(a,i,a,100g)')   'ZDEF ', GoutDescrip%Nz, ' LEVELS ', &
                                    (GoutDescrip%Zlevels(iz), iz = 1, GoutDescrip%Nz)
  write (OutUnit, '(a,i,a,a,5x,a)') 'TDEF ', GoutDescrip%Nt, ' LINEAR ', &
                                    trim(GoutDescrip%Tstart), trim(GoutDescrip%Tinc)
  write (OutUnit, '(a)')            'VARS 1 '
  write (OutUnit, '(a,i,a)')        trim(GoutDescrip%VarName) , GoutDescrip%Nz, ' 99 Diagnostic'
  write (OutUnit, '(a)')            'ENDVARS'

  close (OutUnit, status='keep')

  ! Data file
  ! dimensions from fastest changing to slowest changing are: x, y, z, t
  ! One record is a single horizontal slice -> (nx * ny) data points
  OutRecLen = GoutDescrip%Nx * GoutDescrip%Ny * BinRecFactor
  open (OutUnit, file=GoutDescrip%DataFile, form='unformatted', access='direct', &
        recl=OutRecLen, action='write', status='replace', iostat=Ierror)
  if (Ierror .ne. 0) then
    write (0,*) 'ERROR: cannot open output GRADS data file for writing: ', trim(GoutDescrip%DataFile)
    stop
  end if

  RecNum = 1
  do it = 1, GoutDescrip%Nt
    do iz = 1, GoutDescrip%Nz
      if (modulo(RecNum,100) .eq. 0) then
        write (0,*) '  Writing: ', trim(GoutDescrip%VarName), ': ', RecNum
      end if
      write (OutUnit, rec=RecNum) ((AzAvg(ix,iy,iz,it), ix = 1, GoutDescrip%Nx), iy = 1, GoutDescrip%Ny)
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
