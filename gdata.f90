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

!*************************************************************************
! CheckDataDescrip_CI()
!
! This routine will check the consistency of the data description obtained
! from the GRADS control files. The number of x, y, z, t points should match
! in all files. Also, the total number of variables are calculated and returned
! in Nvars.
!

subroutine CheckDataDescrip_CI(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, &
           DensLoc, TempLoc, VarLoc, VarName)
  use GfileTypes
  implicit none

  type (GradsDataDescription), dimension(*) :: GdataDescrip
  integer Nx, Ny, Nz, Nt, Nfiles
  integer Nvars
  type (GradsVarLocation) :: DensLoc, TempLoc, VarLoc
  character (len=*) :: VarName

  ! i -> file number, j -> var number
  integer i, j
  logical BadData

  DensLoc%Fnum = 0
  DensLoc%Vnum = 0
  TempLoc%Fnum = 0
  TempLoc%Vnum = 0
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
      if (GdataDescrip(i)%VarNames(j) .eq. 'dn0') then
        DensLoc%Fnum = i
        DensLoc%Vnum = j
      end if
      if (GdataDescrip(i)%VarNames(j) .eq. 'tempc') then
        TempLoc%Fnum = i
        TempLoc%Vnum = j
      end if
      if (GdataDescrip(i)%VarNames(j) .eq. VarName) then
        VarLoc%Fnum = i
        VarLoc%Vnum = j
      end if
    end do
  end do

  ! Check to see if you got all vars (dn0, var)
  if (DensLoc%Fnum .eq. 0) then
    write (*,*) 'ERROR: cannot find grid var "dn0" in the GRADS data files'
    BadData = .true.
  end if
  if (TempLoc%Fnum .eq. 0) then
    write (*,*) 'ERROR: cannot find grid var "tempc" in the GRADS data files'
    BadData = .true.
  end if
  if (VarLoc%Fnum .eq. 0) then
    write (*,*) 'ERROR: cannot find grid var "', trim(VarName), '" in the GRADS data files'
    BadData = .true.
  end if

  if (BadData) then
    stop
  end if
end subroutine

!*************************************************************************
! CheckDataDescrip_VS()
!
! This routine will check the consistency of the data description obtained
! from the GRADS control files. The number of x, y, z, t points should match
! in all files. Also, the total number of variables are calculated and returned
! in Nvars.
!

subroutine CheckDataDescrip_VS(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, &
           VarLoc, VarName)
  use GfileTypes
  implicit none

  type (GradsDataDescription), dimension(*) :: GdataDescrip
  integer Nx, Ny, Nz, Nt, Nfiles
  integer Nvars
  type (GradsVarLocation) :: VarLoc
  character (len=*) :: VarName

  ! i -> file number, j -> var number
  integer i, j
  logical BadData

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
      if (GdataDescrip(i)%VarNames(j) .eq. VarName) then
        VarLoc%Fnum = i
        VarLoc%Vnum = j
      end if
    end do
  end do

  ! Check to see if you got all vars (dn0, var)
  if (VarLoc%Fnum .eq. 0) then
    write (*,*) 'ERROR: cannot find grid var "', trim(VarName), '" in the GRADS data files'
    BadData = .true.
  end if

  if (BadData) then
    stop
  end if
end subroutine

!*************************************************************************
! CheckDataDescrip_SW()
!
! This routine will check the consistency of the data description obtained
! from the GRADS control files. The number of x, y, z, t points should match
! in all files. Also, the total number of variables are calculated and returned
! in Nvars.
!

subroutine CheckDataDescrip_SW(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, &
           Uloc, Vloc)
  use GfileTypes
  implicit none

  type (GradsDataDescription), dimension(*) :: GdataDescrip
  integer Nx, Ny, Nz, Nt, Nfiles
  integer Nvars
  type (GradsVarLocation) :: Uloc, Vloc

  ! i -> file number, j -> var number
  integer i, j
  logical BadData

  Uloc%Fnum = 0
  Uloc%Vnum = 0
  Vloc%Fnum = 0
  Vloc%Vnum = 0

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
    end do
  end do

  ! Check to see if you got all vars (dn0, var)
  if (Uloc%Fnum .eq. 0) then
    write (*,*) 'ERROR: cannot find grid var "u" in the GRADS data files'
    BadData = .true.
  end if
  if (Vloc%Fnum .eq. 0) then
    write (*,*) 'ERROR: cannot find grid var "v" in the GRADS data files'
    BadData = .true.
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

!******************************************************************************
! BuildGoutDescrip()
!
! This routine will fill in the GoutDescrip record for the AzAvg data. This
! record is used for creating the GRADS control file.
!

subroutine BuildGoutDescrip(Nx, Ny, Nz, Nt, AzAvg, OfileBase, UndefVal, VarName, &
          Xstart, Xinc, Ystart, Yinc, Zcoords, Tstart, Tinc, GoutDescrip, Diagnostic)

  use GfileTypes
  implicit none

  integer :: Nx, Ny, Nz, Nt
  real, dimension(1:Nx,1:Ny,1:Nz,1:Nt) :: AzAvg
  character (len=*) :: OfileBase, VarName
  real :: UndefVal, Xstart, Xinc, Ystart, Yinc
  real, dimension(1:Nz) :: Zcoords
  character (len=*) :: Tstart, Tinc
  type (GradsOutDescription) :: GoutDescrip
  character (len=*) :: Diagnostic

  integer :: iz

  GoutDescrip%CtlFile = trim(OfileBase) // '.ctl'
  GoutDescrip%DataFile = trim(OfileBase) // '.gra'
  GoutDescrip%Title = 'Diagnostic: ' // trim(Diagnostic)
  GoutDescrip%UndefVal = UndefVal
  GoutDescrip%nx = Nx
  GoutDescrip%ny = Ny
  GoutDescrip%nz = Nz
  GoutDescrip%nt = Nt
  GoutDescrip%Xstart = Xstart
  GoutDescrip%Xinc = Xinc
  GoutDescrip%Ystart = Ystart
  GoutDescrip%Yinc = Yinc
  GoutDescrip%VarName = trim(VarName) // '_' // trim(Diagnostic)
  do iz = 1, Nz
    GoutDescrip%Zcoords(iz) = Zcoords(iz)
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
