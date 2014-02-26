!***************************************************************
! Program to do averaging over the domain and yield a single
! data point at each time step
!
! This program will read in HDF5 data from a RAMS simulation, 
! perform an averaging function, and output the single point time
! series in HDF5 format
!

program tsavg
  use rhdf5_utils
  use diag_utils
  implicit none

  integer, parameter :: LargeString  = 512
  integer, parameter :: MediumString = 256
  integer, parameter :: LittleString = 128
  integer, parameter :: MaxArgFields = 20

  real, parameter :: UndefVal = -999.0

  character (len=MediumString) :: InDir
  character (len=MediumString) :: InSuffix
  character (len=MediumString) :: OutFile
  character (len=MediumString) :: FilterFile
  character (len=LittleString) :: AvgFunc
  character (len=MediumString), dimension(MaxArgFields) :: ArgList
  integer :: Nfields
  character (len=LittleString) :: VarDim, XvarDim, YvarDim
  logical :: UseFilter

  ! Data arrays
  ! Dims: x, y, z, t
  type (Rhdf5Var) :: InXcoords, InYcoords, InZcoords, InTcoords
  type (Rhdf5Var) :: OutXcoords, OutYcoords, OutZcoords, OrigDimSize
  type (Rhdf5Var) :: U, V, AzWind, Speed10m, Dens, Var, Filter, TserAvg, RadMaxWind, Rad34ktWind, Rad50ktWind, Rad64ktWind
  type (Rhdf5Var) :: AzSlp
  type (Rhdf5var) :: PrecipRate, Lwp, Ltss, Theta
  type (Rhdf5var) :: Xvar, Yvar
  character (len=MediumString) :: Ufile, Vfile, AzWindFile, Speed10mFile, DensFile, VarFile, InCoordFile
  character (len=MediumString) :: AzSlpFile
  character (len=MediumString) :: PrecipRateFile, LwpFile, LtssFile, ThetaFile
  character (len=MediumString) :: XvarFile, YvarFile
  character (len=LittleString) :: VarFprefix
  character (len=LittleString) :: VelInType
  character (len=LittleString) :: rh5f_facc
  real :: PrecipRateLimit, Zbot, Ztop
  integer :: Kbot, Ktop
  
  integer :: Xnbins, Ynbins
  real :: Xbstart, Xbinc, Ybstart, Ybinc
  real, dimension(:), allocatable :: Xbins, Ybins
  real :: HfracThresh

  real, dimension(:,:), allocatable :: HorizSpeed

  integer :: rh5f_azwind, rh5f_azslp, rh5f_u, rh5f_v, rh5f_speed10m, rh5f_dens, rh5f_var, rh5f_filter, rh5f_out
  integer :: rh5f_pcprate, rh5f_lwp, rh5f_ltss, rh5f_theta
  integer :: rh5f_xvar, rh5f_yvar

  integer :: id, ib, ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt
  integer :: XvarNz, YvarNz
  real :: DeltaX, DeltaY
  real, dimension(:), allocatable :: InXcoordsKm, InYcoordsKm
  logical :: BadDims

  ! Get the command line arguments
  call GetMyArgs(InDir, InSuffix, OutFile, AvgFunc, FilterFile, UseFilter)
  if ((AvgFunc(1:4) .eq. 'min:') .or. (AvgFunc(1:4) .eq. 'max:') .or. (AvgFunc(1:4) .eq. 'hda:') .or. (AvgFunc(1:10) .eq. 'turb_mmts:')) then
    call String2List(AvgFunc, ':', ArgList, MaxArgFields, Nfields, 'avg spec')
    if (Nfields .eq. 4) then
      ! got the right amount of fields
      !   field    value
      !    1       'min', 'max', or 'hda' 
      !    2       name of variable inside the REVU file
      !    3       prefix for the REVU file name
      !    4       dimensionality of variable
      AvgFunc    = trim(ArgList(1))
      Var%vname  = trim(ArgList(2))
      VarFprefix = trim(ArgList(3))
      VarFile    = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      VarDim     = trim(ArgList(4))
    else
      write (*,*) 'ERROR: average function requires four fields: <avg_func>:<var>:<file>:<dim>'
      stop
    endif
  endif

  if (AvgFunc(1:5) .eq. 'hist:') then
    call String2List(AvgFunc, ':', ArgList, MaxArgFields, Nfields, 'hist spec')
    if (Nfields .eq. 7) then
      ! got the right amount of fields
      !   field    value
      !    1       'hist'
      !    2       name of variable inside the REVU file
      !    3       prefix for the REVU file name
      !    4       dimensionality of variable
      !    5       number of bins
      !    6       bin start
      !    7       bin increment
      AvgFunc    = trim(ArgList(1))
      Var%vname  = trim(ArgList(2))
      VarFprefix = trim(ArgList(3))
      VarFile    = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      VarDim     = trim(ArgList(4))
      read(ArgList(5), '(i)') Xnbins
      read(ArgList(6), '(f)') Xbstart
      read(ArgList(7), '(f)') Xbinc
    else
      write (*,*) 'ERROR: average function hist requires seven fields: hist:<var>:<file>:<dim>:<num_bins>:<bin_start>:<bin_inc>'
      stop
    endif
  endif

  if (AvgFunc(1:6) .eq. 'hfrac:') then
    call String2List(AvgFunc, ':', ArgList, MaxArgFields, Nfields, 'hfrac spec')
    if (Nfields .eq. 5) then
      ! got the right amount of fields
      !   field    value
      !    1       'hist'
      !    2       name of variable inside the REVU file
      !    3       prefix for the REVU file name
      !    4       dimensionality of variable
      !    5       threshold for determining numerator count
      AvgFunc    = trim(ArgList(1))
      Var%vname  = trim(ArgList(2))
      VarFprefix = trim(ArgList(3))
      VarFile    = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      VarDim     = trim(ArgList(4))
      read(ArgList(5), '(f)') HfracThresh
    else
      write (*,*) 'ERROR: average function hfrac requires five fields: hist:<var>:<file>:<dim>:<threshold>'
      stop
    endif
  endif

  if (AvgFunc(1:4) .eq. 'pop:') then
    call String2List(AvgFunc, ':', ArgList, MaxArgFields, Nfields, 'pop spec')
    if (Nfields .eq. 14) then
      ! got the right amount of fields
      !   field    value
      !    1       'pop'
      !    2       name of precip rate variable inside the REVU file
      !    3       prefix for the REVU file name containing the precip rate
      !    4       precip rate threshold to deteriming if rainging or not
      !    5       name of lwp variable inside the REVU file
      !    6       prefix for the REVU file name containing the lwp
      !    7       lwp: number of bins
      !    8       lwp: bin start
      !    9       lwp: bin increment
      !   10       name of ltss variable inside the REVU file
      !   11       prefix for the REVU file name containing the ltss
      !   12       ltss: number of bins
      !   13       ltss: bin start
      !   14       ltss: bin increment
      AvgFunc           = trim(ArgList(1))
      PrecipRate%vname  = trim(ArgList(2))
      VarFprefix        = trim(ArgList(3))
      PrecipRateFile    = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      read(ArgList(4), '(f)') PrecipRateLimit
      Lwp%vname         = trim(ArgList(5))
      VarFprefix        = trim(ArgList(6))
      LwpFile           = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      read(ArgList(7), '(i)') Xnbins
      read(ArgList(8), '(f)') Xbstart
      read(ArgList(9), '(f)') Xbinc
      Ltss%vname        = trim(ArgList(10))
      VarFprefix        = trim(ArgList(11))
      LtssFile          = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      read(ArgList(12), '(i)') Ynbins
      read(ArgList(13), '(f)') Ybstart
      read(ArgList(14), '(f)') Ybinc
    else
      write (*,*) 'ERROR: average function pop requires fourteen fields:'
      write (*,*) 'ERROR:   pop:<pcp_var>:<pcp_file>:<pcp_limit>:<lwp_var>:<lwp_file>:<lwp_nbins>:<lwp_bstart>:<lwp_binc>:'
      write (*,*) 'ERROR:       <ltss_var>:<ltss_file>:<ltss_nbins>:<ltss_bstart>:<ltss_binc>'
      stop
    endif
  endif

  if (AvgFunc(1:7) .eq. 'hist2d:') then
    call String2List(AvgFunc, ':', ArgList, MaxArgFields, Nfields, 'hist2d spec')
    if (Nfields .eq. 13) then
      ! got the right amount of fields
      !   field    value
      !    1       'hist2d'
      !    2       name of x variable inside the REVU file
      !    3       prefix for the REVU file name containing the x variable
      !    4       dimensionalitiy of x variable (2d or 3d)
      !    5       x: number of bins
      !    6       x: bin start
      !    7       x: bin increment
      !    8       name of y variable inside the REVU file
      !    9       prefix for the REVU file name containing the y variable
      !   10       dimensionalitiy of y variable (2d or 3d)
      !   11       y: number of bins
      !   12       y: bin start
      !   13       y: bin increment
      AvgFunc     = trim(ArgList(1))
      Xvar%vname  = trim(ArgList(2))
      VarFprefix  = trim(ArgList(3))
      XvarFile    = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      XvarDim     = trim(ArgList(4))
      read(ArgList(5), '(i)') Xnbins
      read(ArgList(6), '(f)') Xbstart
      read(ArgList(7), '(f)') Xbinc
      Yvar%vname  = trim(ArgList(8))
      VarFprefix  = trim(ArgList(9))
      YvarFile    = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      YvarDim     = trim(ArgList(10))
      read(ArgList(11), '(i)') Ynbins
      read(ArgList(12), '(f)') Ybstart
      read(ArgList(13), '(f)') Ybinc
    else
      write (*,*) 'ERROR: average function hist2d requires thirteen fields:'
      write (*,*) 'ERROR:   hist2d:<x_var>:<x_file>:<x_dim>:<x_nbins><x_bstart><x_binc>:<y_var>:<y_file>:<y_dim>:<y_nbins>:<y_bstart>:<y_binc>:'
      stop
    endif
  endif

  if (AvgFunc(1:9) .eq. 'turb_cov:') then
    call String2List(AvgFunc, ':', ArgList, MaxArgFields, Nfields, 'turb_cov spec')
    if (Nfields .eq. 7) then
      ! got the right amount of fields
      !   field    value
      !    1       'turb_cov'
      !    2       name of x variable inside the REVU file
      !    3       prefix for the REVU file name containing the x variable
      !    4       dimensionalitiy of x variable (2d or 3d)
      !    5       name of y variable inside the REVU file
      !    6       prefix for the REVU file name containing the y variable
      !    7       dimensionalitiy of y variable (2d or 3d)
      AvgFunc     = trim(ArgList(1))

      Xvar%vname  = trim(ArgList(2))
      VarFprefix  = trim(ArgList(3))
      XvarFile    = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      XvarDim     = trim(ArgList(4))

      Yvar%vname  = trim(ArgList(5))
      VarFprefix  = trim(ArgList(6))
      YvarFile    = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      YvarDim     = trim(ArgList(7))

    else
      write (*,*) 'ERROR: average function turb_cov requires seven fields:'
      write (*,*) 'ERROR:   turb_cov:<x_var>:<x_file>:<x_dim>:<y_var>:<y_file>:<y_dim>'
      stop
    endif
  endif

  if (AvgFunc(1:5) .eq. 'ltss:') then
    call String2List(AvgFunc, ':', ArgList, MaxArgFields, Nfields, 'pop spec')
    if (Nfields .eq. 5) then
      ! got the right amount of fields
      !   field    value
      !    1       'ltss'
      !    2       name of theta variable inside the REVU file
      !    3       prefix for the REVU file name containing theta
      !    4       Z (m) of bottom level for taking theta diff
      !    5       Z (m) of top level for taking theta diff
      AvgFunc           = trim(ArgList(1))
      Theta%vname     = trim(ArgList(2))
      VarFprefix        = trim(ArgList(3))
      ThetaFile    = trim(InDir) // '/' //trim(VarFprefix) // trim(InSuffix)
      read(ArgList(4), '(f)') Zbot
      read(ArgList(5), '(f)') Ztop
    else
      write (*,*) 'ERROR: average function ltss requires five fields: ltss:<theta_var>:<theta_file>:<z_bot>:<z_top>'
      stop
    endif
  endif

  if ((AvgFunc(1:9) .eq. 'horiz_ke:') .or. (AvgFunc(1:11) .eq. 'max_azwind:')) then
    call String2List(AvgFunc, ':', ArgList, MaxArgFields, Nfields, 'horiz_ke spec')
    if (Nfields .eq. 2) then
      ! got the right amount of fields
      !   field    value
      !    1       'horiz_ke'
      !    2       type of input
      AvgFunc    = trim(ArgList(1))
      VelInType  = trim(ArgList(2))

      if ((VelInType .ne. 'uv') .and. (VelInType .ne. 's10')) then
        write (*,*) 'ERROR: <in_type> for average functions horiz_ke and max_azwind must be one of: "uv" or "s10"'
        stop
      endif
    else
      write (*,*) 'ERROR: average function ltss requires five fields: ltss:<theta_var>:<theta_file>:<z_bot>:<z_top>'
      stop
    endif
  endif

  write (*,*) 'Time seris of average for RAMS data:'
  write (*,*) '  Input directory: ', trim(InDir)
  write (*,*) '  Input file suffix: ', trim(InSuffix)
  write (*,*) '  Output file:  ', trim(OutFile)
  write (*,*) '  Averaging function: ', trim(AvgFunc)
  if ((AvgFunc .eq. 'horiz_ke') .or. (AvgFunc .eq. 'max_azwind')) then
    write (*,*) '    Input Type: ', trim(VelInType)
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'turb_mmts')) then
    write (*,*) '    Variable name: ', trim(Var%vname)
    write (*,*) '    File name: ', trim(VarFile)
    write (*,*) '    Dimensionality: ', trim(VarDim)
  else if (AvgFunc .eq. 'hist') then
    write (*,*) '    Variable name: ', trim(Var%vname)
    write (*,*) '    File name: ', trim(VarFile)
    write (*,*) '    Dimensionality: ', trim(VarDim)
    write (*,*) '    Binning specs:'
    write (*,*) '      Number of bins: ', Xnbins
    write (*,*) '      Bins start at: ', Xbstart
    write (*,*) '      Delta between bins: ', Xbinc
  else if (AvgFunc .eq. 'hfrac') then
    write (*,*) '    Variable name: ', trim(Var%vname)
    write (*,*) '    File name: ', trim(VarFile)
    write (*,*) '    Dimensionality: ', trim(VarDim)
    write (*,*) '    Threshold for fraction calculation: ', HfracThresh 
  else if (AvgFunc .eq. 'pop') then
    write (*,*) '    Precip rate variable name: ', trim(PrecipRate%vname)
    write (*,*) '    Precip rate File name: ', trim(PrecipRateFile)
    write (*,*) '    Liquid water path variable name: ', trim(Lwp%vname)
    write (*,*) '    Liquid water path file name: ', trim(LwpFile)
    write (*,*) '    Lower troposhperic static stability variable name: ', trim(Ltss%vname)
    write (*,*) '    Lower troposhperic static stability file name: ', trim(LtssFile)
    write (*,*) '    Precip rate threshold: ', PrecipRateLimit
    write (*,*) '    Liquid water path binning specs:'
    write (*,*) '      Number of bins: ', Xnbins
    write (*,*) '      Bins start at: ', Xbstart
    write (*,*) '      Delta between bins: ', Xbinc
    write (*,*) '    Lower tropospheric static stability binning specs:'
    write (*,*) '      Number of bins: ', Ynbins
    write (*,*) '      Bins start at: ', Ybstart
    write (*,*) '      Delta between bins: ', Ybinc
  else if (AvgFunc .eq. 'hist2d') then
    write (*,*) '    X variable name: ', trim(Xvar%vname)
    write (*,*) '    X file name: ', trim(XvarFile)
    write (*,*) '    X variable dimensionality: ', trim(XvarDim)
    write (*,*) '    X variable binning specs:'
    write (*,*) '      Number of bins: ', Xnbins
    write (*,*) '      Bins start at: ', Xbstart
    write (*,*) '      Delta between bins: ', Xbinc
    write (*,*) '    Y variable name: ', trim(Yvar%vname)
    write (*,*) '    Y file name: ', trim(YvarFile)
    write (*,*) '    Y variable dimensionality: ', trim(YvarDim)
    write (*,*) '    Y variable binning specs:'
    write (*,*) '      Number of bins: ', Ynbins
    write (*,*) '      Bins start at: ', Ybstart
    write (*,*) '      Delta between bins: ', Ybinc
  else if (AvgFunc .eq. 'turb_cov') then
    write (*,*) '    X variable name: ', trim(Xvar%vname)
    write (*,*) '    X file name: ', trim(XvarFile)
    write (*,*) '    X variable dimensionality: ', trim(XvarDim)
    write (*,*) '    Y variable name: ', trim(Yvar%vname)
    write (*,*) '    Y file name: ', trim(YvarFile)
    write (*,*) '    Y variable dimensionality: ', trim(YvarDim)
  else if (AvgFunc .eq. 'ltss') then
    write (*,*) '    Theta variable name: ', trim(Theta%vname)
    write (*,*) '    Theta file name: ', trim(ThetaFile)
    write (*,*) '    Z at bottom: ', Zbot
    write (*,*) '    Z at top: ', Ztop
  endif
  write (*,*) '  Filter file: ', trim(FilterFile)
  write (*,*) '    Using filter: ', UseFilter
  write (*,*) ''
  flush(6)

  ! set up file and variable names
  if (VelInType .eq. 'uv') then
    AzWindFile = trim(InDir) // '/speed_t' // trim(InSuffix)
    AzWind%vname = 'speed_t'
  else if (VelInType .eq. 's10') then
    AzWindFile = trim(InDir) // '/speed10m' // trim(InSuffix)
    AzWind%vname = 'speed10m'
  endif

  AzSlpFile = trim(InDir) // '/sea_press' // trim(InSuffix)
  AzSlp%vname = 'sea_press'

  ! FilterFile is set by command line arguments
  Filter%vname = 'filter'

  DensFile = trim(InDir) // '/dn0' // trim(InSuffix)
  Dens%vname = 'dn0'

  Ufile = trim(InDir) // '/u' // trim(InSuffix)
  U%vname = 'u'

  Vfile = trim(InDir) // '/v' // trim(InSuffix)
  V%vname = 'v'

  Speed10mFile = trim(InDir) // '/speed10m' // trim(InSuffix)
  Speed10m%vname = 'speed10m'

  ! Check that the dimensions are consistent between the variables needed for
  ! the selected averaging function.
  !
  ! There is no associated filter with max_azwind and min_azslp since a filter has already been
  ! applied by the azavg program (which created the azwind data). All other functions
  ! need the filter data.
  !
  ! Expect 3D vars to be: (x,y,z,t)
  !        2D vars to be: (x,y,t)
  !

  ! Set Filter%ndims to one here. This will result with Filter%ndims being equal to zero (since
  ! we need to chop of the time dimension to read the filter time step by time step) when we
  ! are not using a filter. When using a filter, Filter%ndims will get set by what's contained
  ! in FilterFile.
  if ((AvgFunc .eq. 'max_azwind') .or. (AvgFunc .eq. 'min_azslp')) then
    ! max_azwind, min_azslp must not use filter (azavg already applied a filter)
    if (UseFilter) then
      write (*,*) 'ERROR: cannot use a filter with function: max_azwind, min_azslp'
      stop
    endif
    
    ! Read in the filter and use it to check against all the other variables
    if (AvgFunc .eq. 'max_azwind') then
      call rhdf5_read_init(AzWindFile, AzWind)

      Nx = AzWind%dims(1)
      Ny = AzWind%dims(2)
      Nz = AzWind%dims(3)
      Nt = AzWind%dims(4)
    else
      ! min_azslp
      call rhdf5_read_init(AzSlpFile, AzSlp)

      Nx = AzSlp%dims(1)
      Ny = AzSlp%dims(2)
      Nz = AzSlp%dims(3)
      Nt = AzSlp%dims(4)
    endif
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
           (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist') .or. &
           (AvgFunc .eq. 'hfrac') .or. (AvgFunc .eq. 'turb_mmts')) then
    call rhdf5_read_init(VarFile, Var)

    if (VarDim .eq. '2d') then
      Nx = Var%dims(1)
      Ny = Var%dims(2)
      Nz = 1
      Nt = Var%dims(3)
    else
      Nx = Var%dims(1)
      Ny = Var%dims(2)
      Nz = Var%dims(3)
      Nt = Var%dims(4)
    endif

    ! filter is optional
    if (UseFilter) then
      call rhdf5_read_init(FilterFile, Filter)
      if (.not.(DimsMatch(Filter, Var))) then
        write (*,*) 'ERROR: dimensions of filter do not match dimensions of input variable: ', trim(Var%vname)
        stop
      endif
    endif
  else if (AvgFunc .eq. 'pop') then
    ! check that the horzontal dimensions of precip rate and lwp match
    call rhdf5_read_init(PrecipRateFile, PrecipRate)
    call rhdf5_read_init(LwpFile, Lwp)

    if (.not.(DimsMatch(PrecipRate, Lwp))) then
      write (*,*) 'ERROR: dimensions of precip rate variable do not match dimensions of liquid water path variable'
      stop
    endif

    ! get info for the LTSS data - this will be 1D (single number for each time step)
    call rhdf5_read_init(LtssFile, Ltss)

    ! record dims - 2d data
    Nx = PrecipRate%dims(1)
    Ny = PrecipRate%dims(2)
    Nz = 1
    Nt = PrecipRate%dims(3)

    ! filter is optional
    if (UseFilter) then
      call rhdf5_read_init(FilterFile, Filter)
      if (.not.(DimsMatch(Filter, PrecipRate))) then
        write (*,*) 'ERROR: dimensions of filter do not match dimensions of precip rate variable'
        stop
      endif
    endif
  else if ((AvgFunc .eq. 'hist2d') .or. (AvgFunc .eq. 'turb_cov')) then
    ! check that the horzontal dimensions of x and y match
    call rhdf5_read_init(XvarFile, Xvar)
    call rhdf5_read_init(YvarFile, Yvar)

    if (.not.(DimsMatch(Xvar, Yvar))) then
      write (*,*) 'ERROR: dimensions of x and y variables do not match'
      stop
    endif

    ! Record dims: if we got to here, then x and y have matching x, y and t dimensions
    ! so get these from the x variable.
    !
    ! Want to allow for x and y to be 2d or 3d. To do this record the z dimensions
    ! separately for x and y. If either of x and y are 3d, then make the output 3d
    ! (do histograms for each level).
    Nx = Xvar%dims(1)
    Ny = Xvar%dims(2)
    Nz = 1                       ! Assume both x and y are 2d, if either is 3d then
                                 ! set Nz to the 3d size. If both are 3d, there is a
                                 ! check that x and y have matching z dimensions so it
                                 ! doesn't matter if Nz gets set from x or y.
    if (XvarDim .eq. '2d') then
      XvarNz = 1
      Nt = Xvar%dims(3)
    else
      XvarNz = Xvar%dims(3)
      Nt = Xvar%dims(4)

      Nz = XvarNz
    endif
    if (YvarDim .eq. '2d') then
      YvarNz = 1
    else
      YvarNz = Yvar%dims(3)

      Nz = YvarNz
    endif

    ! make sure if x and y are both 3d, that their z dimensions match
    if ((XvarDim .eq. '3d') .and. (YvarDim .eq. '3d')) then
      if (XvarNz .ne. YvarNz) then
        write (*,*) 'ERROR: x and y variables are both 3d, but their z dimensions do not match'
        stop
      endif
    endif

    ! make sure if doing turb_cov, that x and y z dimensions match
    if (AvgFunc .eq. 'turb_cov') then
      if (XvarNz .ne. YvarNz) then
        write (*,*) 'ERROR: x and y variables need to have their z dimensions match for "turb_cov" function'
        stop
      endif
    endif

    ! filter is optional
    if (UseFilter) then
      call rhdf5_read_init(FilterFile, Filter)
      if (.not.(DimsMatch(Filter, Xvar))) then
        write (*,*) 'ERROR: dimensions of filter do not match dimensions of x and y variables'
        stop
      endif
    endif
  else if (AvgFunc .eq. 'ltss') then
    ! get dimensions from theta file
    call rhdf5_read_init(ThetaFile, Theta)

    ! record dims - 3d data
    Nx = Theta%dims(1)
    Ny = Theta%dims(2)
    Nz = Theta%dims(3)
    Nt = Theta%dims(4)

    ! filter is optional
    if (UseFilter) then
      call rhdf5_read_init(FilterFile, Filter)
      if (.not.(DimsMatch(Filter, Theta))) then
        write (*,*) 'ERROR: dimensions of filter do not match dimensions of theta variable'
        stop
      endif
    endif
  else if ((AvgFunc .eq. 'horiz_ke') .or. (AvgFunc .eq. 'storm_int')) then
    ! horiz_ke and storm_int require a filter
    if (.not. UseFilter) then
      write (*,*) 'ERROR: must use a filter with functions: horiz_ke and storm_int'
      stop
    endif
    
    ! Read in the filter and use it to check against all the other variables
    call rhdf5_read_init(FilterFile, Filter)
  
    Nx = Filter%dims(1)
    Ny = Filter%dims(2)
    Nz = Filter%dims(3)
    Nt = Filter%dims(4)

    BadDims = .false.
    if (AvgFunc .eq. 'horiz_ke') then
      call rhdf5_read_init(DensFile, Dens)
      BadDims = BadDims .or. (.not.(DimsMatch(Filter, Dens)))
  
      if (VelInType .eq. 'uv') then
        call rhdf5_read_init(Ufile, U)
        BadDims = BadDims .or. (.not.(DimsMatch(Filter, U)))
  
        call rhdf5_read_init(Vfile, V)
        BadDims = BadDims .or. (.not.(DimsMatch(Filter, V)))
      else if (VelInType .eq. 's10') then
        ! speed10m is a 2D variable
        call rhdf5_read_init(Speed10mFile, Speed10m)
        BadDims = BadDims .or. (.not.(DimsMatch(Filter, Speed10m)))

        Nz = 1
      endif

      if (BadDims) then
        write (*,*) 'ERROR: dimensions of filter, dn0, u and v do not match'
        stop
      endif
    else if (AvgFunc .eq. 'storm_int') then
      ! speed10m is a 2D variable
      call rhdf5_read_init(Speed10mFile, Speed10m)
      BadDims = BadDims .or. (.not.(DimsMatch(Filter, Speed10m)))
      Nz = 1
  
      if (BadDims) then
        write (*,*) 'ERROR: dimensions of filter and speed10m do not match'
        stop
      endif
    endif
  endif

  ! Set up the dimensions for reading in the input field data, one time step per read. 
  ! In other words, remove the time dimension from the input dimensions since we will 
  ! be incrementing through every time step in a loop. The time dimension is always the
  ! last dimension so what this boils down to is to decrement number of dimensions by one.
  write (*,*) 'Input variable information:'
  if (AvgFunc .eq. 'max_azwind') then
    AzWind%ndims = AzWind%ndims - 1
    write (*,*) '  Number of dimensions: ', AzWind%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, AzWind%ndims
      write (*,*), '    ', trim(AzWind%dimnames(id)), ': ', AzWind%dims(id)
    enddo
  else if (AvgFunc .eq. 'min_azslp') then
    AzSlp%ndims = AzSlp%ndims - 1
    write (*,*) '  Number of dimensions: ', AzSlp%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, AzSlp%ndims
      write (*,*), '    ', trim(AzSlp%dimnames(id)), ': ', AzSlp%dims(id)
    enddo
  else if (AvgFunc .eq. 'horiz_ke') then
    Dens%ndims = Dens%ndims - 1
    if (VelInType .eq. 'uv') then
      U%ndims = U%ndims - 1
      V%ndims = V%ndims - 1
    else if (VelInType .eq. 's10') then
      Speed10m%ndims = Speed10m%ndims - 1
    endif
    write (*,*) '  Number of dimensions: ', Dens%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Dens%ndims
      write (*,*), '    ', trim(Dens%dimnames(id)), ': ', Dens%dims(id)
    enddo
  else if (AvgFunc .eq. 'storm_int') then
    Speed10m%ndims = Speed10m%ndims - 1
    write (*,*) '  Number of dimensions: ', Speed10m%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Speed10m%ndims
      write (*,*), '    ', trim(Speed10m%dimnames(id)), ': ', Speed10m%dims(id)
    enddo
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
           (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist') .or. &
           (AvgFunc .eq. 'hfrac') .or. (AvgFunc .eq. 'turb_mmts')) then
    Var%ndims = Var%ndims - 1
    write (*,*) '  Number of dimensions: ', Var%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Var%ndims
      write (*,*), '    ', trim(Var%dimnames(id)), ': ', Var%dims(id)
    enddo
  else if (AvgFunc .eq. 'pop') then
    PrecipRate%ndims = PrecipRate%ndims - 1
    Lwp%ndims = Lwp%ndims - 1
    Ltss%ndims = Ltss%ndims - 1
    write (*,*) '  Number of dimensions: ', PrecipRate%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, PrecipRate%ndims
      write (*,*), '    ', trim(PrecipRate%dimnames(id)), ': ', PrecipRate%dims(id)
    enddo
  else if ((AvgFunc .eq. 'hist2d') .or. (AvgFunc .eq. 'turb_cov')) then
    Xvar%ndims = Xvar%ndims - 1
    Yvar%ndims = Yvar%ndims - 1
    write (*,*) '  X variable: '
    write (*,*) '    Number of dimensions: ', Xvar%ndims
    write (*,*) '    Dimension sizes:'
    do id = 1, Xvar%ndims
      write (*,*), '      ', trim(Xvar%dimnames(id)), ': ', Xvar%dims(id)
    enddo
    write (*,*) '  Y variable: '
    write (*,*) '    Number of dimensions: ', Yvar%ndims
    write (*,*) '    Dimension sizes:'
    do id = 1, Yvar%ndims
      write (*,*), '      ', trim(Yvar%dimnames(id)), ': ', Yvar%dims(id)
    enddo
  else if (AvgFunc .eq. 'ltss') then
    Theta%ndims = Theta%ndims - 1
    write (*,*) '  Number of dimensions: ', Theta%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Theta%ndims
      write (*,*), '    ', trim(Theta%dimnames(id)), ': ', Theta%dims(id)
    enddo
  endif
  write (*,*) ''

  if (UseFilter) then
    Filter%ndims = Filter%ndims - 1

    write (*,*) 'Filter variable information:'
    write (*,*) '  Number of dimensions: ', Filter%ndims
    write (*,*) '  Dimension sizes:'
    do id = 1, Filter%ndims
      write (*,*), '    ', trim(Filter%dimnames(id)), ': ', Filter%dims(id)
    enddo
    write (*,*) ''
  endif

  ! Set up the dimensions for the output and allocate the output data array. Always
  ! set up as if the output were 3D. This is done so that the output file can
  ! be read into GRADS which expects 3D variables. Always have (x,y,z) for the
  ! dimension names, but set the sizes of the dimensions according to the averaging
  ! function asked for.
  TserAvg%vname = trim(AvgFunc)
  if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
      (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist') .or. &
      (AvgFunc .eq. 'hfrac') .or. (AvgFunc .eq. 'turb_mmts')) then
    TserAvg%vname = trim(TserAvg%vname) // '_' // trim(VarFprefix)
  endif
  TserAvg%descrip = 'time series averaged ' // trim(AvgFunc) 
  TserAvg%ndims = 3 
  TserAvg%dimnames(1) = 'x' 
  TserAvg%dimnames(2) = 'y' 
  TserAvg%dimnames(3) = 'z' 

  if (AvgFunc .eq. 'max_azwind') then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = 'm/s'
  else if (AvgFunc .eq. 'min_azslp') then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = 'mb'
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max')) then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = Var%units
  else if ((AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hfrac')) then
    ! y has size 2, one for the sum and the other for the count
    ! all z levels
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 2
    TserAvg%dims(3) = Nz
    TserAvg%units = Var%units
  else if ((AvgFunc .eq. 'turb_cov') .or. (AvgFunc .eq. 'turb_mmts')) then
    ! y has size 4
    !   For turb_cov
    !     y(1) - sum of mean values of x
    !     y(2) - sum of mean values of y
    !     y(3) - sum of x'y' products
    !     y(4) - number of points
    !   For turb_mmts
    !     y(1) - sum of mean values of x
    !     y(2) - sum of x'x' products
    !     y(3) - sum of x'x'x' products
    !     y(4) - number of points
    ! all z levels
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 4
    TserAvg%dims(3) = Nz
    TserAvg%units = Var%units
  else if (AvgFunc .eq. 'hist') then
    ! put bin values in the x dimension
    TserAvg%dims(1) = Xnbins
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = Nz
    TserAvg%units = Var%units
  else if (AvgFunc .eq. 'pop') then
    ! put LWP bin values in the x dimension
    ! put LTSS bin values in the y dimension
    ! put Nr and Nt results in z dimension
    TserAvg%dims(1) = Xnbins
    TserAvg%dims(2) = Ynbins
    TserAvg%dims(3) = 2
    TserAvg%units = PrecipRate%units // ':' // Lwp%units
  else if (AvgFunc .eq. 'hist2d') then
    TserAvg%dims(1) = Xnbins
    TserAvg%dims(2) = Ynbins
    TserAvg%dims(3) = Nz
    TserAvg%units = Xvar%units // ':' // Lwp%units
  else if (AvgFunc .eq. 'ltss') then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = Theta%units
  else if (AvgFunc .eq. 'horiz_ke') then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = 'Joules'
  else if (AvgFunc .eq. 'storm_int') then
    ! single point result
    TserAvg%dims(1) = 1
    TserAvg%dims(2) = 1
    TserAvg%dims(3) = 1
    TserAvg%units = 'int'
  endif

  allocate(TserAvg%vdata(TserAvg%dims(1)*TserAvg%dims(2)*TserAvg%dims(3)))

  ! If doing max_azwind, set up the output var for the radius of max wind
  if (AvgFunc .eq. 'max_azwind') then
    RadMaxWind%vname = 'rmw'
    RadMaxWind%descrip = 'time series averaged ' // trim(AvgFunc) 
    RadMaxWind%units = 'km'
    RadMaxWind%ndims = 3 
    RadMaxWind%dimnames(1) = 'x' 
    RadMaxWind%dimnames(2) = 'y' 
    RadMaxWind%dimnames(3) = 'z' 
    RadMaxWind%dims(1) = 1
    RadMaxWind%dims(2) = 1
    RadMaxWind%dims(3) = 1
    allocate(RadMaxWind%vdata(RadMaxWind%dims(1)*RadMaxWind%dims(2)*RadMaxWind%dims(3)))

    Rad34ktWind%vname = 'r34kt'
    Rad34ktWind%descrip = 'time series averaged ' // trim(AvgFunc) 
    Rad34ktWind%units = 'km'
    Rad34ktWind%ndims = 3 
    Rad34ktWind%dimnames(1) = 'x' 
    Rad34ktWind%dimnames(2) = 'y' 
    Rad34ktWind%dimnames(3) = 'z' 
    Rad34ktWind%dims(1) = 1
    Rad34ktWind%dims(2) = 1
    Rad34ktWind%dims(3) = 1
    allocate(Rad34ktWind%vdata(Rad34ktWind%dims(1)*Rad34ktWind%dims(2)*Rad34ktWind%dims(3)))

    Rad50ktWind%vname = 'r50kt'
    Rad50ktWind%descrip = 'time series averaged ' // trim(AvgFunc) 
    Rad50ktWind%units = 'km'
    Rad50ktWind%ndims = 3 
    Rad50ktWind%dimnames(1) = 'x' 
    Rad50ktWind%dimnames(2) = 'y' 
    Rad50ktWind%dimnames(3) = 'z' 
    Rad50ktWind%dims(1) = 1
    Rad50ktWind%dims(2) = 1
    Rad50ktWind%dims(3) = 1
    allocate(Rad50ktWind%vdata(Rad50ktWind%dims(1)*Rad50ktWind%dims(2)*Rad50ktWind%dims(3)))

    Rad64ktWind%vname = 'r64kt'
    Rad64ktWind%descrip = 'time series averaged ' // trim(AvgFunc) 
    Rad64ktWind%units = 'km'
    Rad64ktWind%ndims = 3 
    Rad64ktWind%dimnames(1) = 'x' 
    Rad64ktWind%dimnames(2) = 'y' 
    Rad64ktWind%dimnames(3) = 'z' 
    Rad64ktWind%dims(1) = 1
    Rad64ktWind%dims(2) = 1
    Rad64ktWind%dims(3) = 1
    allocate(Rad64ktWind%vdata(Rad64ktWind%dims(1)*Rad64ktWind%dims(2)*Rad64ktWind%dims(3)))
  endif

  ! Report the dimensions
  write (*,*) 'Output variable information:'
  write (*,*) '  Name: ', trim(TserAvg%vname)
  write (*,*) '  Units: ', trim(TserAvg%units)
  write (*,*) '  Description: ', trim(TserAvg%descrip)
  write (*,*) '  Number of dimensions: ', TserAvg%ndims
  write (*,*) '  Dimension sizes:'
  do id = 1, TserAvg%ndims
    write (*,*), '    ', trim(TserAvg%dimnames(id)), ': ', TserAvg%dims(id)
  enddo
  write (*,*) ''
  write (*,*) '  Number of time steps: ', Nt
  write (*,*) ''
  flush(6)

  ! Read in the input coordinates
  if (AvgFunc .eq. 'max_azwind') then
    InCoordFile = trim(AzWindFile)
  else if (AvgFunc .eq. 'min_azslp') then
    InCoordFile = trim(AzSlpFile)
  else if (AvgFunc .eq. 'horiz_ke') then
    InCoordFile = trim(DensFile)
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
           (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist') .or. &
           (AvgFunc .eq. 'hfrac') .or. (AvgFunc .eq. 'turb_mmts')) then
    InCoordFile = trim(VarFile)
  else if (AvgFunc .eq. 'pop') then
    InCoordFile = trim(PrecipRateFile)
  else if ((AvgFunc .eq. 'hist2d') .or. (AvgFunc .eq. 'turb_cov')) then
    ! Want to end up using a 3d var if one of x and y is 3d, that is use 2d only
    ! if both x and y are 2d. Set InCoordFile to x, and switch to y only if y is 3d.
    ! This gives you:
    !     x   y    InCoordFile
    !    2d  2d       2d  (x)
    !    2d  3d       3d  (y)
    !    3d  2d       3d  (x)
    !    3d  3d       3d  (y)
    InCoordFile = trim(XvarFile)
    if (YvarDim .eq. '3d') then
      InCoordFile = trim(YvarFile)
    endif
  else if (AvgFunc .eq. 'ltss') then
    InCoordFile = trim(ThetaFile)
  else if (AvgFunc .eq. 'storm_int') then
    InCoordFile = trim(Speed10mFile)
  endif

  InXcoords%vname = 'x_coords'
  call rhdf5_read_init(InCoordFile, InXcoords)
  call rhdf5_read(InCoordFile, InXcoords)

  InYcoords%vname = 'y_coords'
  call rhdf5_read_init(InCoordFile, InYcoords)
  call rhdf5_read(InCoordFile, InYcoords)

  if (AvgFunc .eq. 'storm_int') then
    ! fake it for the storm intensity metric
    InZcoords%vname = 'z_coords'
    InZcoords%ndims = 1
    InZcoords%dims(1) = 1
    InZcoords%dimnames(1) = 'z'
    InZcoords%units = 'meter'
    InZcoords%descrip = 'sigma-z'
    allocate(InZcoords%vdata(1))
    InZcoords%vdata(1) = 10.0
  else
    InZcoords%vname = 'z_coords'
    call rhdf5_read_init(InCoordFile, InZcoords)
    call rhdf5_read(InCoordFile, InZcoords)
  endif

  InTcoords%vname = 't_coords'
  call rhdf5_read_init(InCoordFile, InTcoords)
  call rhdf5_read(InCoordFile, InTcoords)

  ! Need to get delta for x and y (in meters) and the z heights (in meters)
  ! for DoHorizKe to be able to calculate volume * density -> mass
  allocate(InXcoordsKm(Nx))
  allocate(InYcoordsKm(Ny))
  if ((AvgFunc .eq. 'max_azwind') .or. (AvgFunc .eq. 'min_azslp')) then
    ! x,y coords are in meters, convert to km
    do ix = 1, Nx
      InXcoordsKm(ix) = InXcoords%vdata(ix) / 1000.0
    enddo
    do iy = 1, Ny
      InYcoordsKm(iy) = InYcoords%vdata(iy) / 1000.0
    enddo
  else
    ! x,y coords are in degrees lon,lat respectively
    call ConvertGridCoords(Nx, Ny, Nz, InXcoords%vdata, InYcoords%vdata, InXcoordsKm, InYcoordsKm)
  endif

  DeltaX = InXcoordsKm(2) - InXcoordsKm(1)
  DeltaY = InYcoordsKm(2) - InYcoordsKm(1)

  write (*,*) 'Horizontal grid info:'
  write (*,*) '  X range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', InXcoords%vdata(1), InXcoords%vdata(Nx), InXcoordsKm(1), InXcoordsKm(Nx)
  write (*,*) '  Y range (min lat, max lat) --> (min y, max y): '
  write (*,*) '    ', InYcoords%vdata(1), InYcoords%vdata(Ny), InYcoordsKm(1), InYcoordsKm(Ny)
  write (*,*) ''
  write (*,*) 'Vertical grid info:'
  do iz = 1, Nz
    write (*,*) '  ', iz, ' --> ', InZcoords%vdata(iz)
  end do
  write (*,*) ''
  flush(6)

  ! If doing histogram or pop, calculate the bin values
  if ((AvgFunc .eq. 'hist') .or. (AvgFunc .eq. 'pop') .or. (AvgFunc .eq. 'hist2d')) then
    allocate(Xbins(Xnbins))
    Xbins(1) = Xbstart
    do ib = 2, Xnbins
      Xbins(ib) = Xbins(ib-1) + Xbinc
    enddo

    if ((AvgFunc .eq. 'pop') .or. (AvgFunc .eq. 'hist2d')) then
      allocate(Ybins(Ynbins))
      Ybins(1) = Ybstart
      do ib = 2, Ynbins
        Ybins(ib) = Ybins(ib-1) + Ybinc
      enddo
    endif
  endif

  ! if doing ltss, find the indices associated with Zbot and Ztop
  if (AvgFunc .eq. 'ltss') then
    Kbot = Nt
    do iz = Nz,1,-1
      if (InZcoords%vdata(iz) .ge. Zbot) then
        Kbot = iz
      endif
    enddo

    Ktop = 1
    do iz = 1,Nz
      if (InZcoords%vdata(iz) .le. Ztop) then
        Ktop = iz
      endif
    enddo

    write (*,*) 'Height indices for LTSS calculation:'
    write (*,*) '   Zbot: ', Zbot
    write (*,*) '     Index for Zbot: ', Kbot
    write (*,*) '   Ztop: ', Ztop
    write (*,*) '     Index for Ztop: ', Ktop
    write (*,*) ''
  endif

  ! Prepare the output coordinates
  !
  ! Create dummy coordinates (for GRADS sake) and write out the
  ! time series as a 4D var, (x,y,z,t) regardless of how many
  ! true dimensions exist in the output.
  OutXcoords%vname = 'x_coords'
  OutXcoords%ndims = 1
  OutXcoords%dimnames(1) = 'x'
  OutXcoords%units = 'degrees_east'
  OutXcoords%descrip = 'longitude'
  if ((AvgFunc .eq. 'hist') .or. (AvgFunc .eq. 'pop') .or. (AvgFunc .eq. 'hist2d')) then
    OutXcoords%dims(1) = Xnbins
    allocate(OutXcoords%vdata(Xnbins))
    do ib = 1, Xnbins
      OutXcoords%vdata(ib) = Xbins(ib)
    enddo
  else
    OutXcoords%dims(1) = 1
    allocate(OutXcoords%vdata(1))
    OutXcoords%vdata(1) = 1.0
  endif
  
  OutYcoords%vname = 'y_coords'
  OutYcoords%ndims = 1
  OutYcoords%dimnames(1) = 'y'
  OutYcoords%units = 'degrees_north'
  OutYcoords%descrip = 'latitude'
  if ((AvgFunc .eq. 'pop') .or. (AvgFunc .eq. 'hist2d')) then
    OutYcoords%dims(1) = Ynbins
    allocate(OutYcoords%vdata(Ynbins))
    do ib = 1, Ynbins
      OutYcoords%vdata(ib) = Ybins(ib)
    enddo
  elseif ((AvgFunc .eq. 'hfrac') .or. (AvgFunc .eq. 'hda')) then
    OutYcoords%dims(1) = 2
    allocate(OutYcoords%vdata(2))
    OutYcoords%vdata(1) = 1.0
    OutYcoords%vdata(2) = 2.0
  elseif ((AvgFunc .eq. 'turb_cov') .or. (AvgFunc .eq. 'turb_mmts')) then
    OutYcoords%dims(1) = 4
    allocate(OutYcoords%vdata(4))
    OutYcoords%vdata(1) = 1.0
    OutYcoords%vdata(2) = 2.0
    OutYcoords%vdata(3) = 3.0
    OutYcoords%vdata(4) = 4.0
  else
    OutYcoords%dims(1) = 1
    allocate(OutYcoords%vdata(1))
    OutYcoords%vdata(1) = 1.0
  endif
  
  OutZcoords%vname = 'z_coords'
  OutZcoords%ndims = 1
  OutZcoords%dimnames(1) = 'z'
  OutZcoords%units = 'meter'
  OutZcoords%descrip = 'sigma-z'
  if ((AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist') .or. (AvgFunc .eq. 'hist2d') .or.  &
      (AvgFunc .eq. 'hfrac') .or. (AvgFunc .eq. 'turb_cov') .or. (AvgFunc .eq. 'turb_mmts')) then
    OutZcoords%dims(1) = Nz
    allocate(OutZcoords%vdata(Nz))
    do iz = 1, Nz
      OutZcoords%vdata(iz) = InZcoords%vdata(iz)
    enddo
  elseif (AvgFunc .eq. 'pop') then
    ! need a size of 2 for the z dimension
    !   index 1 --> number of grids that are raining
    !   index 2 --> number of total grids (selected per bin on lwp)
    OutZcoords%dims(1) = 2
    allocate(OutZcoords%vdata(2))
    ! dummy values
    OutZcoords%vdata(1) = 1.0
    OutZcoords%vdata(2) = 2.0
  else
    OutZcoords%dims(1) = 1
    allocate(OutZcoords%vdata(1))
    OutZcoords%vdata(1) = 1.0
  endif

  ! Perform the averaging function.
  rh5f_facc = 'W'
  call rhdf5_open_file(OutFile, rh5f_facc, 1, rh5f_out)

  rh5f_facc = 'R'
  if (UseFilter) then
    call rhdf5_open_file(FilterFile, rh5f_facc, 1, rh5f_filter)
  endif
  if (AvgFunc .eq. 'max_azwind') then
    call rhdf5_open_file(AzWindFile, rh5f_facc, 1, rh5f_azwind)
  else if (AvgFunc .eq. 'min_azslp') then
    call rhdf5_open_file(AzSlpFile, rh5f_facc, 1, rh5f_azslp)
  else if (AvgFunc .eq. 'horiz_ke') then
    call rhdf5_open_file(DensFile, rh5f_facc, 1, rh5f_dens)
    if (VelInType .eq. 'uv') then
      call rhdf5_open_file(Ufile, rh5f_facc, 1, rh5f_u)
      call rhdf5_open_file(Vfile, rh5f_facc, 1, rh5f_v)
    else if (VelInType .eq. 's10') then
      call rhdf5_open_file(Speed10mFile, rh5f_facc, 1, rh5f_speed10m)
    endif
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
           (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist') .or. &
           (AvgFunc .eq. 'hfrac') .or. (AvgFunc .eq. 'turb_mmts')) then
    call rhdf5_open_file(VarFile, rh5f_facc, 1, rh5f_var)
  else if (AvgFunc .eq. 'pop') then
    call rhdf5_open_file(PrecipRateFile, rh5f_facc, 1, rh5f_pcprate)
    call rhdf5_open_file(LwpFile, rh5f_facc, 1, rh5f_lwp)
    call rhdf5_open_file(LtssFile, rh5f_facc, 1, rh5f_ltss)
  else if ((AvgFunc .eq. 'hist2d') .or. (AvgFunc .eq. 'turb_cov')) then
    call rhdf5_open_file(XvarFile, rh5f_facc, 1, rh5f_xvar)
    call rhdf5_open_file(YvarFile, rh5f_facc, 1, rh5f_yvar)
  else if (AvgFunc .eq. 'ltss') then
    call rhdf5_open_file(ThetaFile, rh5f_facc, 1, rh5f_theta)
  else if (AvgFunc .eq. 'storm_int') then
    call rhdf5_open_file(Speed10mFile, rh5f_facc, 1, rh5f_speed10m)
  endif

  do it = 1, Nt
    ! if using a filter read in the data
    if (UseFilter) then
      call rhdf5_read_variable(rh5f_filter, Filter%vname, Filter%ndims, it, Filter%dims, rdata=Filter%vdata)
    endif

    ! do the averaging function
    if (AvgFunc .eq. 'max_azwind') then
      call rhdf5_read_variable(rh5f_azwind, AzWind%vname, AzWind%ndims, it, AzWind%dims, rdata=AzWind%vdata)
      call DoMaxAzWind(Nx, Ny, Nz, AzWind%vdata, InXcoordsKm, UndefVal, TserAvg%vdata(1), RadMaxWind%vdata(1), &
                      Rad34ktWind%vdata(1), Rad50ktWind%vdata(1), Rad64ktWind%vdata(1))
      deallocate(AzWind%vdata)
    else if (AvgFunc .eq. 'min_azslp') then
      call rhdf5_read_variable(rh5f_azslp, AzSlp%vname, AzSlp%ndims, it, AzSlp%dims, rdata=AzSlp%vdata)
      call DoMinAzSlp(Nx, Ny, Nz, AzSlp%vdata, UndefVal, TserAvg%vdata(1))
      deallocate(AzSlp%vdata)
    else if (AvgFunc .eq. 'horiz_ke') then
      call rhdf5_read_variable(rh5f_dens, Dens%vname, Dens%ndims, it, Dens%dims, rdata=Dens%vdata)

      if (VelInType .eq. 'uv') then
        call rhdf5_read_variable(rh5f_u, U%vname, U%ndims, it, U%dims, rdata=U%vdata)
        call rhdf5_read_variable(rh5f_v, V%vname, V%ndims, it, V%dims, rdata=V%vdata)

        allocate(HorizSpeed(Nx,Ny))
        call UvToSpeed(Nx, Ny, Nz, U%vdata, V%vdata, HorizSpeed)
        call DoHorizKe(Nx, Ny, Nz, Dens%vdata, HorizSpeed, Filter%vdata, DeltaX, DeltaY, InZcoords%vdata, TserAvg%vdata(1))

        deallocate(U%vdata)
        deallocate(V%vdata)
        deallocate(HorizSpeed)
      else if (VelInType .eq. 's10') then
        call rhdf5_read_variable(rh5f_speed10m, Speed10m%vname, Speed10m%ndims, it, Speed10m%dims, rdata=Speed10m%vdata)
        call DoHorizKe(Nx, Ny, Nz, Dens%vdata, Speed10m%vdata, Filter%vdata, DeltaX, DeltaY, InZcoords%vdata, TserAvg%vdata(1))
        deallocate(Speed10m%vdata)
      endif

      deallocate(Dens%vdata)
    else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
             (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist') .or. &
             (AvgFunc .eq. 'hfrac') .or. (AvgFunc .eq. 'turb_mmts')) then
      call rhdf5_read_variable(rh5f_var, Var%vname, Var%ndims, it, Var%dims, rdata=Var%vdata)

      if (AvgFunc .eq. 'min') then
        call DoMin(Nx, Ny, Nz, Filter%dims(3), Var%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata(1))
      else if (AvgFunc .eq. 'max') then
        call DoMax(Nx, Ny, Nz, Filter%dims(3), Var%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata(1))
      else if (AvgFunc .eq. 'hda') then
        call DoHda(Nx, Ny, Nz, Filter%dims(3), Var%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata)
      else if (AvgFunc .eq. 'hist') then
        call DoHist(Nx, Ny, Nz, Filter%dims(3), Xnbins, Var%vdata, Filter%vdata, UseFilter, UndefVal, Xbins, TserAvg%vdata)
      else if (AvgFunc .eq. 'hfrac') then
        call DoHfrac(Nx, Ny, Nz, Filter%dims(3), HfracThresh, Var%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata)
      else if (AvgFunc .eq. 'turb_mmts') then
        call DoTurbMmts(Nx, Ny, Nz, Filter%dims(3), Var%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata)
      endif

      deallocate(Var%vdata)
    else if (AvgFunc .eq. 'pop') then
      call rhdf5_read_variable(rh5f_pcprate, PrecipRate%vname, PrecipRate%ndims, it, PrecipRate%dims, rdata=PrecipRate%vdata)
      call rhdf5_read_variable(rh5f_lwp, Lwp%vname, Lwp%ndims, it, Lwp%dims, rdata=Lwp%vdata)
      call rhdf5_read_variable(rh5f_ltss, Ltss%vname, Ltss%ndims, it, Ltss%dims, rdata=Ltss%vdata)

      call DoPop(Nx, Ny, Nz, Filter%dims(3), Xnbins, Ynbins, PrecipRate%vdata, Lwp%vdata, Ltss%vdata(1), Filter%vdata, UseFilter, UndefVal, PrecipRateLimit, Xbins, Ybins, TserAvg%vdata)

      deallocate(PrecipRate%vdata)
      deallocate(Lwp%vdata)
      deallocate(Ltss%vdata)
    else if ((AvgFunc .eq. 'hist2d') .or. (AvgFunc .eq. 'turb_cov')) then
      call rhdf5_read_variable(rh5f_xvar, Xvar%vname, Xvar%ndims, it, Xvar%dims, rdata=Xvar%vdata)
      call rhdf5_read_variable(rh5f_yvar, Yvar%vname, Yvar%ndims, it, Yvar%dims, rdata=Yvar%vdata)

      if (AvgFunc .eq. 'hist2d') then
        call DoHist2d(Nx, Ny, Nz, XvarNz, YvarNz, Filter%dims(3), Xnbins, Ynbins, Xvar%vdata, Yvar%vdata, Filter%vdata, UseFilter, UndefVal, Xbins, Ybins, TserAvg%vdata)
      else if (AvgFunc .eq. 'turb_cov') then
        call DoTurbCov(Nx, Ny, Nz, Filter%dims(3), Xvar%vdata, Yvar%vdata, Filter%vdata, UseFilter, UndefVal, TserAvg%vdata)
      endif

      deallocate(Xvar%vdata)
      deallocate(Yvar%vdata)
    else if (AvgFunc .eq. 'ltss') then
      call rhdf5_read_variable(rh5f_theta, Theta%vname, Theta%ndims, it, Theta%dims, rdata=Theta%vdata)

      call DoLtss(Nx, Ny, Nz, Filter%dims(3), Theta%vdata, Filter%vdata, UseFilter, UndefVal, Kbot, Ktop, TserAvg%vdata(1))

      deallocate(Theta%vdata)
    else if (AvgFunc .eq. 'storm_int') then
      call rhdf5_read_variable(rh5f_speed10m, Speed10m%vname, Speed10m%ndims, it, Speed10m%dims, rdata=Speed10m%vdata)
      call DoStormInt(Nx, Ny, Nz, Speed10m%vdata, Filter%vdata, TserAvg%vdata(1))
      deallocate(Speed10m%vdata)
    endif

    ! if using a filter, deallocate the space for the next time around
    if (UseFilter) then
      deallocate(Filter%vdata)
    endif

    ! write out the averaged data
    call rhdf5_write_variable(rh5f_out, TserAvg%vname, TserAvg%ndims, it, TserAvg%dims, &
       TserAvg%units, TserAvg%descrip, TserAvg%dimnames, rdata=TserAvg%vdata)
    if (AvgFunc .eq. 'max_azwind') then
      call rhdf5_write_variable(rh5f_out, RadMaxWind%vname, RadMaxWind%ndims, it, RadMaxWind%dims, &
         RadMaxWind%units, RadMaxWind%descrip, RadMaxWind%dimnames, rdata=RadMaxWind%vdata)
      call rhdf5_write_variable(rh5f_out, Rad34ktWind%vname, Rad34ktWind%ndims, it, Rad34ktWind%dims, &
         Rad34ktWind%units, Rad34ktWind%descrip, Rad34ktWind%dimnames, rdata=Rad34ktWind%vdata)
      call rhdf5_write_variable(rh5f_out, Rad50ktWind%vname, Rad50ktWind%ndims, it, Rad50ktWind%dims, &
         Rad50ktWind%units, Rad50ktWind%descrip, Rad50ktWind%dimnames, rdata=Rad50ktWind%vdata)
      call rhdf5_write_variable(rh5f_out, Rad64ktWind%vname, Rad64ktWind%ndims, it, Rad64ktWind%dims, &
         Rad64ktWind%units, Rad64ktWind%descrip, Rad64ktWind%dimnames, rdata=Rad64ktWind%vdata)
    endif

    ! print a message for the user on longer jobs so that it can be
    ! seen that progress is being made
    if (modulo(it,100) .eq. 0) then
      write (*,*) 'Working: Number of time steps processed so far: ', it
    endif
  enddo

  rh5f_facc = 'W'
  call rhdf5_close_file(rh5f_out)

  rh5f_facc = 'R'
  if (UseFilter) then
    call rhdf5_close_file(rh5f_filter)
  endif
  if (AvgFunc .eq. 'max_azwind') then
    call rhdf5_close_file(rh5f_azwind)
  else if (AvgFunc .eq. 'max_azslp') then
    call rhdf5_close_file(rh5f_azslp)
  else if (AvgFunc .eq. 'horiz_ke') then
    call rhdf5_close_file(rh5f_dens)
    if (VelInType .eq. 'uv') then
      call rhdf5_close_file(rh5f_u)
      call rhdf5_close_file(rh5f_v)
    else if (VelInType .eq. 's10') then
      call rhdf5_close_file(rh5f_speed10m)
    endif
  else if ((AvgFunc .eq. 'min') .or. (AvgFunc .eq. 'max') .or. &
           (AvgFunc .eq. 'hda') .or. (AvgFunc .eq. 'hist') .or. &
           (AvgFunc .eq. 'hfrac') .or. (AvgFunc .eq. 'turb_mmts')) then
    call rhdf5_close_file(rh5f_var)
  else if (AvgFunc .eq. 'pop') then
    call rhdf5_close_file(rh5f_pcprate)
    call rhdf5_close_file(rh5f_lwp)
    call rhdf5_close_file(rh5f_ltss)
  else if ((AvgFunc .eq. 'hist2d') .or. (AvgFunc .eq. 'turb_cov')) then
    call rhdf5_close_file(rh5f_xvar)
    call rhdf5_close_file(rh5f_yvar)
  else if (AvgFunc .eq. 'ltss') then
    call rhdf5_close_file(rh5f_theta)
  else if (AvgFunc .eq. 'storm_int') then
    call rhdf5_close_file(rh5f_speed10m)
  endif
 
  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''

  ! Finish off output file
  ! write out the coordinate data
  call rhdf5_write(OutFile, OutXcoords, 1)
  call rhdf5_write(OutFile, OutYcoords, 1)
  call rhdf5_write(OutFile, OutZcoords, 1)
  call rhdf5_write(OutFile, InTcoords, 1)

  ! If doing hist function, write out the input dimension sizes
  ! for downstream analyses. Eg. if you want to do fractional
  ! area calculations then the counts in the histogram do not
  ! tell you how many total points are in the domain.
  if ((AvgFunc .eq. 'hist') .or. (AvgFunc .eq. 'hist2d')) then
    ! common settings
    OrigDimSize%ndims = 1
    OrigDimSize%dims(1) = 1
    OrigDimSize%units = 'number'
    allocate(OrigDimSize%vdata(1))

    ! X
    OrigDimSize%vname = 'Nx'
    OrigDimSize%dimnames(1) = 'x'
    OrigDimSize%descrip = 'number of domain x points'
    OrigDimSize%vdata(1) = float(Nx)
    call rhdf5_write(OutFile, OrigDimSize, 1)

    ! Y
    OrigDimSize%vname = 'Ny'
    OrigDimSize%dimnames(1) = 'y'
    OrigDimSize%descrip = 'number of domain y points'
    OrigDimSize%vdata(1) = float(Ny)
    call rhdf5_write(OutFile, OrigDimSize, 1)

    ! Z
    OrigDimSize%vname = 'Nz'
    OrigDimSize%dimnames(1) = 'z'
    OrigDimSize%descrip = 'number of domain z points'
    OrigDimSize%vdata(1) = float(Nz)
    call rhdf5_write(OutFile, OrigDimSize, 1)

    ! T
    OrigDimSize%vname = 'Nt'
    OrigDimSize%dimnames(1) = 't'
    OrigDimSize%descrip = 'number of domain t points'
    OrigDimSize%vdata(1) = float(Nt)
    call rhdf5_write(OutFile, OrigDimSize, 1)

    deallocate(OrigDimSize%vdata)
  endif

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(OutFile, OutXcoords, 'x')
  call rhdf5_set_dimension(OutFile, OutYcoords, 'y')
  call rhdf5_set_dimension(OutFile, OutZcoords, 'z')
  call rhdf5_set_dimension(OutFile, InTcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(OutFile, TserAvg)
  if (AvgFunc .eq. 'max_azwind') then
    call rhdf5_attach_dimensions(OutFile, RadMaxWind)
    call rhdf5_attach_dimensions(OutFile, Rad34ktWind)
    call rhdf5_attach_dimensions(OutFile, Rad50ktWind)
    call rhdf5_attach_dimensions(OutFile, Rad64ktWind)
  endif
  
  stop

contains
!**********************************************************************
! SUBROUTINES
!**********************************************************************

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the command line arguments
!
subroutine GetMyArgs(InDir, InSuffix, OutFile, AvgFunc, FilterFile, UseFilter)
  implicit none

  character (len=*) :: InDir, InSuffix, OutFile, AvgFunc, FilterFile
  logical :: UseFilter

  integer :: iargc
  character (len=128) :: arg
  integer :: Nitems

  logical :: BadArgs

  if (iargc() .ne. 5) then
    write (*,*) 'ERROR: must supply exactly 5 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: tsavg <in_dir> <in_suffix> <out_file> <avg_function> <filter_file>'
    write (*,*) '        <in_dir>: directory where input files live'
    write (*,*) '        <in_suffix>: suffix on input file names'
    write (*,*) '        <out_file>: name of output file, HDF5 format'
    write (*,*) '        <avg_function>: averaging function to use on input data'
    write (*,*) '            horiz_ke:<in_type> -> total kinetic energy form horizontal winds'
    write (*,*) '              <in_type>: "uv" calculate from lowest level of u and v fields,'
    write (*,*) '                         "s10" calculate from 10m wind speed field,'
    write (*,*) '            storm_int -> storm intensity metric from horizontal wind speeds'
    write (*,*) '            max_azwind:<in_type> -> max value of azimuthially averaged wind'
    write (*,*) '              <in_type>: "uv" calculate from lowest level of u and v fields,'
    write (*,*) '                         "s10" calculate from 10m wind speed field,'
    write (*,*) '            min_azslp -> min value of azimuthially averaged sea-level pressure'
    write (*,*) '            min:<var>:<file>:<dim> -> domain minimum for <var>'
    write (*,*) '              <var>: revu var name inside the file'
    write (*,*) '              <file>: prefix for the revu file'
    write (*,*) '              <dim>: dimensionality of variable, either "2d" or "3d"'
    write (*,*) '            max:<var>:<file>:<dim> -> domain maximum for <var>'
    write (*,*) '              <var>: revu var name inside the file'
    write (*,*) '              <file>: prefix for the revu file'
    write (*,*) '              <dim>: dimensionality of variable, either "2d" or "3d"'
    write (*,*) '            hda:<var>:<file>:<dim> -> horizontal domain average at each z level for <var>'
    write (*,*) '              <var>: revu var name inside the file'
    write (*,*) '              <file>: prefix for the revu file'
    write (*,*) '              <dim>: dimensionality of variable, either "2d" or "3d"'
    write (*,*) '            hist:<var>:<file>:<dim>:<num_bins>:<bin_start>:<bin_inc>'
    write (*,*) '              <var>,<file>,<dim> same as for hda'
    write (*,*) '              <num_bins>: number of bins'
    write (*,*) '              <bin_start>: starting value for bins'
    write (*,*) '              <bin_inc>: delta between bins'
    write (*,*) '            hfrac:<var>:<file>:<dim>:<thresh>'
    write (*,*) '              <var>,<file>,<dim> same as for hda'
    write (*,*) '              <threshold>: if var > threshold, that hoizontal grid cell gets a 1, otherwise a zero.'
    write (*,*) '                           Then the fraction for that level = Number of cells with 1s divided by total Number of cells'
    write (*,*) '            pop:<pcp_var>:<pcp_file>:<pcp_limit>:<lwp_var>:<lwp_file>:<lwp_nbins>:<lwp_bstart>:<lwp_binc>:'
    write (*,*) '                <ltss_var>:<ltss_file>:<ltss_nbins>:<ltss_bstart>:<ltss_binc>'
    write (*,*) '              precip rate variable (2d):'
    write (*,*) '                <pcp_var>: revu var name inside the file'
    write (*,*) '                <pcp_file>: prefix for the revu file'
    write (*,*) '                <pcp_limit>: threshold to determine if raining'
    write (*,*) '                  precip rate >= threshold --> raining'
    write (*,*) '                  precip rate <  threshold --> not raining'
    write (*,*) '              liquid water path variable (2d):'
    write (*,*) '                <lwp_var>: revu var name inside the file'
    write (*,*) '                <lwp_file>: prefix for the revu file'
    write (*,*) '                <lwp_nbins>: number of bins'
    write (*,*) '                <lwp_bstart>: starting value for bins'
    write (*,*) '                <lwp_binc>: delta between bins'
    write (*,*) '              lower tropospheric static stabilty variable (1d):'
    write (*,*) '                <ltss_var>: revu var name inside the file'
    write (*,*) '                <ltss_file>: prefix for the revu file'
    write (*,*) '                <ltss_nbins>: number of bins'
    write (*,*) '                <ltss_bstart>: starting value for bins'
    write (*,*) '                <ltss_binc>: delta between bins'
    write (*,*) '            ltss:<theta_var>:<theta_file>:<k_bot>:<k_top>'
    write (*,*) '                <theta_var>: revu var name inside the file'
    write (*,*) '                <theta_file>: prefix for the revu file'
    write (*,*) '                <z_bot>: height (Z) for bottom'
    write (*,*) '                <z_top>: height (Z) for top'
    write (*,*) '            hist2d:<x_var>:<x_file>:<x_dim>:<x_nbins><x_bstart><x_binc>:<y_var>:<y_file>:<y_dim>:<y_nbins>:<y_bstart>:<y_binc>:'
    write (*,*) '              for both X and Y bins:'
    write (*,*) '                <*_var>: revu var name inside the file'
    write (*,*) '                <*_file>: prefix for the revu file'
    write (*,*) '                <*_dim>: dimensionality of variable, either "2d" or "3d"'
    write (*,*) '                <*_nbins>: number of bins'
    write (*,*) '                <*_bstart>: starting value for bins'
    write (*,*) '                <*_binc>: delta between bins'
    write (*,*) '            turb_cov:<x_var>:<x_file>:<x_dim>:<y_var>:<y_file>:<y_dim>'
    write (*,*) '                <[xy]_var>: revu var name inside the file'
    write (*,*) '                <[xy]_file>: prefix for the revu file'
    write (*,*) '                <[xy]_dim>: dimensionality of variable, either "2d" or "3d"'
    write (*,*) '            turb_mmts:<x_var>:<x_file>:<x_dim>'
    write (*,*) '                <x_var>: revu var name inside the file'
    write (*,*) '                <x_file>: prefix for the revu file'
    write (*,*) '                <x_dim>: dimensionality of variable, either "2d" or "3d"'
    write (*,*) '        <filter_file>: file containing the filter mask'
    stop
  end if

  call getarg(1, InDir)
  call getarg(2, InSuffix)
  call getarg(3, OutFile)
  call getarg(4, AvgFunc)
  call getarg(5, FilterFile)

  if ((FilterFile .eq. 'none') .or. (FilterFile .eq. 'NONE')) then
    UseFilter = .false.
  else
    UseFilter = .true.
  endif

  BadArgs = .false.

  if ((AvgFunc(1:9)  .ne. 'horiz_ke:')    .and. &
      (AvgFunc(1:4)  .ne. 'min:')         .and. &
      (AvgFunc(1:4)  .ne. 'max:')         .and. &
      (AvgFunc(1:4)  .ne. 'hda:')         .and. &
      (AvgFunc(1:5)  .ne. 'hist:')        .and. &
      (AvgFunc(1:9)  .ne. 'turb_cov:')    .and. &
      (AvgFunc(1:10) .ne. 'turb_mmts:')   .and. &
      (AvgFunc(1:6)  .ne. 'hfrac:')       .and. &
      (AvgFunc(1:4)  .ne. 'pop:')         .and. &
      (AvgFunc(1:7)  .ne. 'hist2d:')      .and. &
      (AvgFunc(1:5)  .ne. 'ltss:')        .and. &
      (AvgFunc(1:10) .ne. 'storm_int:')   .and. &
      (AvgFunc(1:9)  .ne. 'min_azslp')   .and. &
      (AvgFunc(1:11) .ne. 'max_azwind:')) then
    write (*,*) 'ERROR: <avg_function> must be one of:'
    write (*,*) '          horiz_ke:<in_type>'
    write (*,*) '          min:<var>:<file>:<dim>'
    write (*,*) '          max:<var>:<file>:<dim>'
    write (*,*) '          hda:<var>:<file>:<dim>'
    write (*,*) '          hist:<var>:<file>:<dim>:<num_bins>:<bin_start>:<bin_inc>'
    write (*,*) '          turb_cov:<x_var>:<x_file>:<x_dim>:<y_var>:<y_file>:<y_dim>'
    write (*,*) '          turb_mmts:<x_var>:<x_file>:<x_dim>'
    write (*,*) '          hfrac:<var>:<file>:<dim>:<threshold>'
    write (*,*) '          pop:<pcp_var>:<pcp_file>:<pcp_limit>:<lwp_var>:<lwp_file>:<lwp_nbins>:<lwp_bstart>:<lwp_binc>:'
    write (*,*) '              <ltss_var>:<ltss_file>:<ltss_nbins>:<ltss_bstart>:<ltss_binc>'
    write (*,*) '          hist2d:<x_var>:<x_file>:<x_dim>:<x_nbins><x_bstart><x_binc>:<y_var>:<y_file>:<y_dim>:<y_nbins>:<y_bstart>:<y_binc>:'
    write (*,*) '          ltss:<theta_var>:<theta_file>:<k_bot>:<k_top>'
    write (*,*) '          storm_int'
    write (*,*) '          min_azslp'
    write (*,*) '          max_azwind:<in_type>'
    write (*,*) ''
    BadArgs = .true.
  end if

  if (BadArgs) then
    stop
  end if

  return
end subroutine GetMyArgs

!************************************************************************************
! DoMaxAzWind()
!
! This subroutine will find the maximum wind speed in AzWind and copy that
! to TserAvg. The intent of this metric is to mimic the Saffir-Simpson scale
! which uses the average 10m tangential wind speed that has been sustained
! over 10 minutes.
!
! Typically, don't have enough RAMS files to do 10m time averaging so as a
! proxy use the azimuthally averaged tangetial wind speed.
!
! Keep track of the index into the radius values that produced the maximum
! wind speed. Pass this back to the caller so that the caller can look
! up the radius value. Radius is the first dimension, x.
!
! Also keep track of the greatest radial index where the wind speed was
! >= 34 kt, 50 kt and 64 kt. AzWind is in m/s so the corresponding values
! are:
!   34 kt -> 17.5 m/s
!   50 kt -> 25.7 m/s
!   64 kt -> 32.9 m/s
!
subroutine DoMaxAzWind(Nx, Ny, Nz, AzWind, Radius, UndefVal, AzWindMax, RadMax, Rad34kt, Rad50kt, Rad64kt)
  implicit none

  ! 1 kt = 0.514444 m/s
  real, parameter :: w_34kt_in_ms = 34.0 * 0.514444
  real, parameter :: w_50kt_in_ms = 50.0 * 0.514444
  real, parameter :: w_64kt_in_ms = 64.0 * 0.514444

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: AzWind
  real, dimension(Nx) :: Radius
  real :: UndefVal
  real :: AzWindMax
  real :: RadMax, Rad34kt, Rad50kt, Rad64kt

  integer :: ix,iy,iz
  integer :: i_rmw, i_34kt, i_50kt, i_64kt

  ! dimension order is: x,y,z

  ! RAMS places the first z level below the surface so use the second level
  ! to approximate the 10m winds.

  AzWindMax = 0.0
  i_rmw = 0
  i_34kt = 0
  i_50kt = 0
  i_64kt = 0
  iy = 1 ! dummy dimension in azmuthally averaged wind speed
  if (Nz .eq. 1) then
    iz = 1
  else
    iz = 2
  endif
  do ix = 1, Nx
    ! Max wind
    if ((AzWind(ix,iy,iz) .gt. AzWindMax) .and. (anint(AzWind(ix,iy,iz)) .ne. UndefVal)) then
      AzWindMax = AzWind(ix,iy,iz)
      i_rmw = ix
    endif
    ! Last index where wind >= 34 kt
    if ((AzWind(ix,iy,iz) .ge. w_34kt_in_ms) .and. (anint(AzWind(ix,iy,iz)) .ne. UndefVal)) then
      i_34kt = ix
    endif
    ! Last index where wind >= 50 kt
    if ((AzWind(ix,iy,iz) .ge. w_50kt_in_ms) .and. (anint(AzWind(ix,iy,iz)) .ne. UndefVal)) then
      i_50kt = ix
    endif
    ! Last index where wind >= 64 kt
    if ((AzWind(ix,iy,iz) .ge. w_64kt_in_ms) .and. (anint(AzWind(ix,iy,iz)) .ne. UndefVal)) then
      i_64kt = ix
    endif
  end do

  RadMax = Radius(i_rmw)
  Rad34kt = InterpRadius(Nx, Ny, Nz, AzWind, Radius, w_34kt_in_ms, i_34kt, UndefVal)
  Rad50kt = InterpRadius(Nx, Ny, Nz, AzWind, Radius, w_50kt_in_ms, i_50kt, UndefVal)
  Rad64kt = InterpRadius(Nx, Ny, Nz, AzWind, Radius, w_64kt_in_ms, i_64kt, UndefVal)

  return
end subroutine DoMaxAzWind

!************************************************************************************
! DoMinAzSlp()
!
! This subroutine will find the min sea-level pressure in AzSlp and copy that
! to TserAvg. The intent of this metric is to mimic the usage of mininum central
! pressure for an intensity measurement.
!
!
subroutine DoMinAzSlp(Nx, Ny, Nz, AzSlp, UndefVal, AzSlpMin)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: AzSlp
  real :: UndefVal
  real :: AzSlpMin

  integer :: ix,iy,iz

  ! dimension order is: x,y,z

  ! RAMS places the first z level below the surface so use the second level
  ! to approximate the 10m winds.

  AzSlpMin = 1000000
  iy = 1 ! dummy dimension in azmuthally averaged sea level pressure
  iz = 1 ! dummy dimension in azmuthally averaged sea level pressure
  do ix = 1, Nx
    ! Min SLP
    if ((AzSlp(ix,iy,iz) .lt. AzSlpMin) .and. (anint(AzSlp(ix,iy,iz)) .ne. UndefVal)) then
      AzSlpMin = AzSlp(ix,iy,iz)
    endif
  end do

  return
end subroutine DoMinAzSlp

!**************************************************************************************
! InterpRadius()
!
! This function will interpolate a radial value (x) given the cooresponding wind speed.
!
real function InterpRadius(Nx, Ny, Nz, AzWind, Radius, W, Ileft, UndefVal)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: AzWind
  real, dimension(Nx) :: Radius
  real :: W, UndefVal
  integer :: Ileft

  integer :: ix1, ix2, iy, iz
  real :: W1, W2, R1, R2

  iy = 1 ! dummy dimension in azmuthally averaged wind speed
  if (Nz .eq. 1) then
    iz = 1
  else
    iz = 2
  endif

  if (Ileft .eq. 0) then
    ! Didn't find any entry in AzWind >= W -> return UndefVal
    InterpRadius = UndefVal
  else if (Ileft .eq. Nx) then
    ! The last entry in AzWind is >= W -> can't interpolate so just return the last entry
    InterpRadius = AzWind(Nx,iy,iz)
  else
    ! W has fallen between two entries in AzWind -> linearly interpolate radius value
    ! Ileft points to the last entry in AzWind that was >= W (left side of interval)
    ix1 = Ileft
    ix2 = Ileft + 1

    R1 = Radius(ix1)
    R2 = Radius(ix2)
    W1 = AzWind(ix1,iy,iz)
    W2 = AzWind(ix2,iy,iz)

    InterpRadius = R1 + (((W - W1)/(W2 - W1)) * (R2 - R1))
  endif

end function InterpRadius

!**************************************************************************************
! DoHorizKe()
!
! This routine will calculate the total kinetic energy over the given cylindrical
! volume. Do not want average since we want the size of the storm reflected in
! this diagnostic.

subroutine DoHorizKe(Nx, Ny, Nz, Dens, Speed, Filter, DeltaX, DeltaY, Zcoords, TserAvg)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: Dens, Filter
  real, dimension(Nx,Ny) :: Speed
  real :: DeltaX, DeltaY
  real, dimension(Nz) :: Zcoords
  real :: TserAvg

  integer ix,iy,iz
  integer NumPoints
  real SumKe, CurrKe, LevThickness

  ! KE is 1/2 * m * v^2
  !   - calculate this a every point inside the defined cylindrical volume
  !   - get m by density * volume
  !   - v^2 is based on horiz velocity so is equal to Speed^2
  !
  ! Zcoords are technically the center points of the levels in the RAMS simulation. Since we don't have
  ! the level definition from the RAMS runs here, just use the difference from the i+1st z coord minus the
  ! ith z coord to approzimate the ith level thickness. This will be close enough for the measurement.
  !
  ! The integrated kinetic engery measurement done on hurricanes is taken over the domain on the 10m
  ! level. Just use the first model level above the surface (z == 2) for this measurement which will
  ! be close enough.

  SumKe = 0.0
  NumPoints = 0

  iz = 2

  do iy = 2, Ny-1
    do ix = 2, Nx-1
      if (anint(Filter(ix,iy,iz)) .eq. 1.0) then
        if (iz .eq. Nz) then
          ! Use the level below for this case (since no level above)
          LevThickness = Zcoords(iz) - Zcoords(iz-1)
        else
          LevThickness = Zcoords(iz+1) - Zcoords(iz)
        end if
        CurrKe = 0.5 * DeltaX * DeltaY * LevThickness * Dens(ix,iy,iz) * Speed(ix,iy)**2.
        SumKe = SumKe + CurrKe
        NumPoints = NumPoints + 1
      endif
    enddo
  enddo

  TserAvg = SumKe;
  
  return
end subroutine DoHorizKe

!**************************************************************************************
! DoHda()
!
! This routine will do horizontal domain average for all z levels.
!
subroutine DoHda(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomStats)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(2,Nz) :: DomStats ! (1,z) <-- domain sum, (2,z) <-- domain count

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint
  integer :: NumPoints
  real    :: DomSum

  do iz = 1, Nz
    DomSum = 0.0
    NumPoints = 0

    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          DomSum = DomSum + Var(ix,iy,iz)
          NumPoints = NumPoints + 1
        endif
      enddo
    enddo

    DomStats(1,iz) = DomSum
    DomStats(2,iz) = float(NumPoints)
  enddo

  return
end subroutine DoHda

!**************************************************************************************
! DoHfrac()
!
! This routine will do horizontal domain fraction calculations
!
subroutine DoHfrac(Nx, Ny, Nz, FilterNz, Threshold, Var, Filter, UseFilter, UndefVal, DomStats)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal, Threshold
  real, dimension(2,Nz) :: DomStats ! (1,x) <-- fraction count, (2,x) <-- total count

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint
  integer :: NumPoints, NumFracPoints

  do iz = 1, Nz
    NumPoints = 0
    NumFracPoints = 0

    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          if (Var(ix,iy,iz) .gt. Threshold) then
            NumFracPoints = NumFracPoints + 1
          endif
          NumPoints = NumPoints + 1
        endif
      enddo
    enddo

    DomStats(1,iz) = float(NumFracPoints)
    DomStats(2,iz) = float(NumPoints)
  enddo

  return
end subroutine DoHfrac

!**************************************************************************************
! DoMin()
!
! This routine will return the domain minimum value. Values equal to UndefVal will
! be ignored.
!
subroutine DoMin(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomMin)
  implicit none

  real, parameter :: BigPosNum = 10.0e50

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real :: DomMin

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint

  DomMin = BigPosNum
  do iz = 1, Nz
    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          if (Var(ix,iy,iz) .lt. DomMin) then
            DomMin = Var(ix,iy,iz)
          endif
        endif
      enddo
    enddo
  enddo

  if (DomMin .eq. BigPosNum) then
    ! all entries were UndefVal
    DomMin = UndefVal
  endif

  return
end subroutine DoMin

!**************************************************************************************
! DoMax()
!
! This routine will return the domain maximum value. Values equal to UndefVal will
! be ignored.
!
subroutine DoMax(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomMax)
  implicit none

  real, parameter :: BigNegNum = -10.0e50

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real :: DomMax

  integer :: ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint

  DomMax = BigNegNum
  do iz = 1, Nz
    ! RAMS reserves the first and last x and y values for lateral
    ! boundaries. These only contain valid field data under certain
    ! circumstances such as cyclic boundary cases. Most of the time
    ! we want these to be excluded so for now always exclude them
    ! (shouldn't hurt results with cyclic boundaries where the
    ! boundary values could have been included).

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          if (Var(ix,iy,iz) .gt. DomMax) then
            DomMax = Var(ix,iy,iz)
          endif
        endif
      enddo
    enddo
  enddo

  if (DomMax .eq. BigNegNum) then
    ! all entries were UndefVal
    DomMax = UndefVal
  endif

  return
end subroutine DoMax

!**************************************************************************************
! DoHist()
!
! This routine will do histogram binning over all of the domain.
!

subroutine DoHist(Nx, Ny, Nz, FilterNz, Nb, Var, Filter, UseFilter, UndefVal, Bins, Counts)
  implicit none

  integer :: Nx, Ny, Nz, Nb, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(Nb,Nz) :: Counts
  real, dimension(Nb) :: Bins

  integer :: ib, ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint

  ! Emulate matlab histc command, ie the values in Bins are treated as the edges
  ! of the bin ranges. Do not count any values from Var that fall outside the
  ! bin range (< Bins(1) or > Bins(Nb)). The count is incremented for bin ib when:
  !
  !   Bins(ib) <= Var(ix,iy,iz) < Bins(ib+1)
  !
  ! The count in Bin(Nb), last bin, is incremented when:
  !
  !   Var(ix,iy,iz) == Bins(ib)
  !
  ! Build the histogram (counts)
  do iz = 1, Nz
    do ib = 1, Nb
      Counts(ib,iz) = 0.0
    enddo

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = 2, Ny-1
      do ix = 2, Nx-1
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          ! Check all bins except the last.
          !
          ! Exiting out of the loop when finding the bin will help a lot when the
          ! distribution of values is biased toward smaller values. After exiting
          ! out of the loop you can either just check the last bin (which will be wasted)
          ! or put in a logical variable and check that variable saying you can skip the
          ! check of the last bin. Since you would have to check the logical variable and
          ! the last bin every time you might as well just check the last bin instead.
          do ib = 1, Nb-1
             if ((Bins(ib) .le. Var(ix,iy,iz)) .and. (Var(ix,iy,iz) .lt. Bins(ib+1))) then
                Counts(ib,iz) = Counts(ib,iz) + 1.0
                exit
             endif
          enddo
          ! check the last bin
          if (Bins(Nb) .eq. Var(ix,iy,iz)) then
            Counts(Nb,iz) = Counts(Nb,iz) + 1.0
          endif
        endif
      enddo
    enddo
  enddo

  return
end subroutine DoHist

!**************************************************************************************
! DoPop()
!
! This routine will create counts for generating a probability of precip statistic.
! The data is binned on LWP, and each bin gets two counts: one for the number of grid cells that
! it is raining in that bin and the other for the total number of cells in that bin. A cell counts
! as raining when its precip rate is >= to the precip limit
!

subroutine DoPop(Nx, Ny, Nz, FilterNz, Nb_lwp, Nb_ltss, PrecipRate, Lwp, Ltss, Filter, UseFilter, UndefVal, PrecipRateLimit, Bins_lwp, Bins_ltss, Counts)
  implicit none

  integer :: Nx, Ny, Nz, Nb_lwp, Nb_ltss, FilterNz
  real, dimension(Nx,Ny,Nz) :: PrecipRate, Lwp
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal, PrecipRateLimit, Ltss
  real, dimension(Nb_lwp) :: Bins_lwp
  real, dimension(Nb_ltss) :: Bins_ltss
  real, dimension(Nb_lwp,Nb_ltss,2) :: Counts

  integer :: ib_lwp, ib_ltss, ix, iy, iz
  integer :: filter_z
  logical :: SelectPoint

  ! Emulate matlab histc command, ie the values in Bins are treated as the edges
  ! of the bin ranges. Do not count any values from Lwp that fall outside the
  ! bin range (< Bins(1) or > Bins(Nb)). The count is incremented for bin ib when:
  !
  !   Bins(ib) <= Lwp(ix,iy,iz) < Bins(ib+1)
  !
  ! The count in Bin(Nb), last bin, is incremented when:
  !
  !   Lwp(ix,iy,iz) == Bins(ib)
  !
  ! Same for the Ltts bins
  !
  do ib_lwp = 1, Nb_lwp
    do ib_ltss = 1, Nb_ltss
      Counts(ib_lwp,ib_ltss,1) = 0.0
      Counts(ib_lwp,ib_ltss,2) = 0.0
    enddo
  enddo

  ! Figure out the LTSS bin
  ! Then loop through the LWP and figure out the counts
  ib_ltss = FindBin(Nb_ltss, Bins_ltss, Ltss)

  ! Only if an LTSS bin has been selected
  if (ib_ltss .ne. -1) then
    do iz = 1, Nz
      if (Nz .eq. 1) then
        ! 2D var, use the z = 2 level (first model level above the surface)
        ! since typically have the 3D filter z = 1 level all zeros (don't want
        ! to include below surface model level in analysis).
        filter_z = 2
      else
        filter_z = iz
      endif
  
      do iy = 2, Ny-1
        do ix = 2, Nx-1
          SelectPoint = anint(Lwp(ix,iy,iz)) .ne. UndefVal
          if (UseFilter) then
            SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
          endif
          if (SelectPoint) then
            ib_lwp = FindBin(Nb_lwp, Bins_lwp, Lwp(ix,iy,iz))
            if (ib_lwp .ne. -1) then
              Counts(ib_lwp,ib_ltss,1) = Counts(ib_lwp,ib_ltss,1) + 1.0
              if (PrecipRate(ix,iy,iz) .ge. PrecipRateLimit) then
                Counts(ib_lwp,ib_ltss,2) = Counts(ib_lwp,ib_ltss,2) + 1.0
              endif
            endif
          endif
        enddo
      enddo
    enddo
  endif

  return
end subroutine DoPop

!**************************************************************************************
! DoHist2d()
!
! This routine will create a 2D histogram.
!
subroutine DoHist2d(Nx, Ny, Nz, XvarNz, YvarNz, FilterNz, Xnbins, Ynbins, Xvar, Yvar, Filter, UseFilter, UndefVal, Xbins, Ybins, Counts)
  implicit none

  integer :: Nx, Ny, Nz, XvarNz, YvarNz, FilterNz, Xnbins, Ynbins
  real, dimension(Nx,Ny,XvarNz) :: Xvar
  real, dimension(Nx,Ny,YvarNz) :: Yvar
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(Xnbins) :: Xbins
  real, dimension(Ynbins) :: Ybins
  real, dimension(Xnbins,Ynbins,Nz) :: Counts

  integer :: ix, iy, iz, iz_x, iz_y, iz_filter
  integer :: ib_x, ib_y

  logical :: SelectPoint

  ! go level by level
  do iz = 1, Nz
    ! zero out the counts
    do ib_x = 1, Xnbins
      do ib_y = 1, Ynbins
        Counts(ib_x,ib_y,iz) = 0.0
      enddo
    enddo

    ! Since Xvar and Yvar can have differing z dimensions, check to make sure we are not
    ! going past the max z coordinate.
    if (iz .gt. XvarNz) then
      iz_x = XvarNz
    else
      iz_x = iz
    endif
    if (iz .gt. YvarNz) then
      iz_y = YvarNz
    else
      iz_y = iz
    endif
   
    ! If both x and y are 2D vars, use the z = 2 level (first model level above the surface)
    ! since typically have the 3D filter z = 1 level all zeros (don't want
    ! to include below surface model level in analysis).
    if (Nz .eq. 1) then
      iz_filter = 2
    else
      iz_filter = iz
    endif

    ! Walk through each grid cell on this level. Simultaneously bin the Xvar and Yvar values.
    ! Increment the count for this grid cell only if the Xvar and Yvar values got placed into
    ! an Xbin and Ybin.
    do iy = 2, Ny-1
      do ix = 2, Nx-1
        ! don't select if either x or y is undefined
        SelectPoint = ((anint(Xvar(ix,iy,iz_x)) .ne. UndefVal) .and. &
                       (anint(Yvar(ix,iy,iz_y)) .ne. UndefVal))

        ! if using a filter, check the filter value
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,iz_filter)) .eq. 1.0)
        endif

        ! get the bins and only count if both x and y bins were found
        ib_x = FindBin(Xnbins, Xbins, Xvar(ix,iy,iz_x))
        ib_y = FindBin(Ynbins, Ybins, Yvar(ix,iy,iz_y))
        SelectPoint = SelectPoint .and. (ib_x .ne. -1) .and. (ib_y .ne. -1)
        
        if (SelectPoint) then
          Counts(ib_x,ib_y,iz) = Counts(ib_x,ib_y,iz) + 1.0
        endif
      enddo
    enddo
  enddo
end subroutine DoHist2d

!**************************************************************************************
! DoLtss()
!
! This routine will calculate the lower tropsheric static stability (LTSS) from the
! theta field. 
!
! Use an average of differences. Tried difference of averages and the traces are noisy.

subroutine DoLtss(Nx, Ny, Nz, FilterNz, Theta, Filter, UseFilter, UndefVal, Kbot, Ktop, Ltss)
  implicit none

  integer :: Nx, Ny, Nz, Nb, FilterNz
  real, dimension(Nx,Ny,Nz) :: Theta
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal, Ltss
  integer :: Kbot, Ktop

  integer :: ib, ix, iy
  integer :: filter_zbot, filter_ztop
  logical :: SelectPoint
  integer :: Npts

  ! Figure out which levels to use in the filter data
  if (Nz .eq. 1) then
    ! 2D var, use the z = 2 level (first model level above the surface)
    ! since typically have the 3D filter z = 1 level all zeros (don't want
    ! to include below surface model level in analysis).
    filter_zbot = 2
    filter_ztop = 2
  else
    filter_zbot = Kbot
    filter_ztop = Ktop
  endif

  ! Calculate theta difference at each point and then take the average of
  ! these differences.
  Ltss = 0.0
  Npts = 0
  do iy = 2, Ny-1
    do ix = 2, Nx-1
      SelectPoint = (anint(Theta(ix,iy,Kbot)) .ne. UndefVal) .and. (anint(Theta(ix,iy,Ktop)) .ne. UndefVal)
      if (UseFilter) then
        SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_zbot)) .eq. 1.0) &
                                  .and. (anint(Filter(ix,iy,filter_ztop)) .eq. 1.0)
      endif

      ! Sum up the theta differences
      if (SelectPoint) then
        Ltss = Ltss + (Theta(ix,iy,Ktop) - Theta(ix,iy,Kbot))
        Npts = Npts + 1
      endif
    enddo
  enddo

  if (Npts .eq. 0) then
    Ltss = UndefVal
  else
    Ltss = Ltss / float(Npts)
  endif

  return
end subroutine DoLtss

!**************************************************************************************
! DoStormInt()
!
! This routine will calculate a storm intensity metric based on the horizontal wind
! speeds. The metric is based on the hurricane categories from the Saffir-Simpson
! and their associated wind speeds.
!
!  convert the surface wind data into weighted coverage data
!  run through each time point and each x,y location of the sfc wind data
!  weight each count by the Saffir-Simpson category number
!
!   Wind speed   Weight
!    <  33 m/s -> 0  (not hurricane force)
!    33-42 m/s -> 1  (Category 1 hurricane)
!    43-49 m/s -> 2  (etc.)
!    50-58 m/s -> 3
!    59-69 m/s -> 4
!    >= 70 m/s -> 5
!

subroutine DoStormInt(Nx, Ny, Nz, Speed10m, Filter, TserAvg)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: Filter
  real, dimension(Nx,Ny) :: Speed10m
  real :: TserAvg

  integer ix,iy,iz
  integer nCat0, nCat1, nCat2, nCat3, nCat4, nCat5, NumPoints
  real Wspeed, SiMetric

  nCat0 = 0
  nCat1 = 0
  nCat2 = 0
  nCat3 = 0
  nCat4 = 0
  nCat5 = 0

  iz = 2 ! use next to bottom layer in the filter

  do iy = 2, Ny-1
    do ix = 2, Nx-1
      if (anint(Filter(ix,iy,iz)) .eq. 1.0) then
        ! Count up the number of grid points with wind speeds fitting each of the
        ! Saffir-Simpson categories. Then form the metric by weighting each category
        ! count.
        Wspeed = Speed10m(ix,iy)
        if (Wspeed .ge. 70.0) then
          nCat5 = nCat5 + 1
        else if (Wspeed .ge. 59.0) then
          nCat4 = nCat4 + 1
        else if (Wspeed .ge. 50.0) then
          nCat3 = nCat3 + 1
        else if (Wspeed .ge. 43.0) then
          nCat2 = nCat2 + 1
        else if (Wspeed .ge. 33.0) then
          nCat1 = nCat1 + 1
        else
          nCat0 = nCat0 + 1
        endif
      endif
    enddo
  enddo

  !Linear weighting
  NumPoints = nCat0 + nCat1 + nCat2 + nCat3 + nCat4 + nCat5
  SiMetric = float(nCat1) + (float(nCat2)*2.0) + (float(nCat3)*3.0) + (float(nCat4)*4.0) + (float(nCat5)*5.0)

  if (NumPoints .eq. 0) then
    TserAvg = 0.0
  else
    TserAvg = SiMetric / float(NumPoints)
  end if
  
  return
end subroutine DoStormInt

!**************************************************************************************
! DoTurbCov()
!
! This routine will do summing for the calculation of the mean, variance and skew
! (first three moments) of the input variable.
!
subroutine DoTurbCov(Nx, Ny, Nz, FilterNz, Xvar, Yvar, Filter, UseFilter, UndefVal, DomStats)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Xvar, Yvar
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(4,Nz) :: DomStats ! (1,Nz) <-- sum of mean values from running mean of Xvar
                                    ! (2,Nz) <-- sum of mean values from running mean of Yvar
                                    ! (3,Nz) <-- sum of (x-xmean)*(y-ymean)
                                    ! (4,Nz) <-- number of points

  integer :: ix, iy, iz
  integer :: Xstart, Xend, Ystart, Yend
  integer :: filter_z
  logical :: SelectPoint
  integer :: NumPoints
  real    :: DomS1
  real    :: DomS2
  real    :: DomS3
  real    :: Xmean, Ymean

  ! RAMS reserves the first and last x and y values for lateral
  ! boundaries. These only contain valid field data under certain
  ! circumstances such as cyclic boundary cases. Most of the time
  ! we want these to be excluded so for now always exclude them
  ! (shouldn't hurt results with cyclic boundaries where the
  ! boundary values could have been included).
  Xstart = 2
  Xend = Nx - 1
  Ystart = 2
  Yend = Ny - 1

  do iz = 1, Nz
    DomS1 = 0.0
    DomS2 = 0.0
    DomS3 = 0.0
    NumPoints = 0

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = Ystart, Yend
      do ix = Xstart, Xend
        SelectPoint = ((anint(Xvar(ix,iy,iz)) .ne. UndefVal) .and. (anint(Yvar(ix,iy,iz)) .ne. UndefVal))
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          Xmean = Calc2dMean(Nx, Ny, Nz, Xvar, ix, iy, iz, Xstart, Xend, Ystart, Yend)
          Ymean = Calc2dMean(Nx, Ny, Nz, Yvar, ix, iy, iz, Xstart, Xend, Ystart, Yend)
          DomS1 = DomS1 + Xmean
          DomS2 = DomS2 + Ymean
          DomS3 = DomS3 + ((Xvar(ix,iy,iz)-Xmean)*(Yvar(ix,iy,iz)-Ymean))
          NumPoints = NumPoints + 1
        endif
      enddo
    enddo

    DomStats(1,iz) = DomS1
    DomStats(2,iz) = DomS2
    DomStats(3,iz) = DomS3
    DomStats(4,iz) = float(NumPoints)
  enddo

  return
end subroutine DoTurbCov

!**************************************************************************************
! DoTurbMmts()
!
! This routine will do summing for the calculation of the mean, variance and skew
! (first three moments) of the input variable.
!
subroutine DoTurbMmts(Nx, Ny, Nz, FilterNz, Var, Filter, UseFilter, UndefVal, DomStats)
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx,Ny,Nz) :: Var
  real, dimension(Nx,Ny,FilterNz) :: Filter
  logical :: UseFilter
  real :: UndefVal
  real, dimension(4,Nz) :: DomStats ! (1,Nz) <-- sum of mean values from running mean
                                    ! (2,Nz) <-- sum of (x-xmean)**2
                                    ! (3,Nz) <-- sum of (x-xmean)**3
                                    ! (4,Nz) <-- number of points

  integer :: ix, iy, iz
  integer :: Xstart, Xend, Ystart, Yend
  integer :: filter_z
  logical :: SelectPoint
  integer :: NumPoints
  real    :: DomS1
  real    :: DomS2
  real    :: DomS3
  real    :: Vmean

  ! RAMS reserves the first and last x and y values for lateral
  ! boundaries. These only contain valid field data under certain
  ! circumstances such as cyclic boundary cases. Most of the time
  ! we want these to be excluded so for now always exclude them
  ! (shouldn't hurt results with cyclic boundaries where the
  ! boundary values could have been included).
  Xstart = 2
  Xend = Nx - 1
  Ystart = 2
  Yend = Ny - 1

  do iz = 1, Nz
    DomS1 = 0.0
    DomS2 = 0.0
    DomS3 = 0.0
    NumPoints = 0

    if (Nz .eq. 1) then
      ! 2D var, use the z = 2 level (first model level above the surface)
      ! since typically have the 3D filter z = 1 level all zeros (don't want
      ! to include below surface model level in analysis).
      filter_z = 2
    else
      filter_z = iz
    endif

    do iy = Ystart, Yend
      do ix = Xstart, Xend
        SelectPoint = anint(Var(ix,iy,iz)) .ne. UndefVal
        if (UseFilter) then
          SelectPoint = SelectPoint .and. (anint(Filter(ix,iy,filter_z)) .eq. 1.0)
        endif
        if (SelectPoint) then
          Vmean = Calc2dMean(Nx, Ny, Nz, Var, ix, iy, iz, Xstart, Xend, Ystart, Yend)
          DomS1 = DomS1 + Vmean
          DomS2 = DomS2 + ((Var(ix,iy,iz)-Vmean)**2)
          DomS3 = DomS3 + ((Var(ix,iy,iz)-Vmean)**3)
          NumPoints = NumPoints + 1
        endif
      enddo
    enddo

    DomStats(1,iz) = DomS1
    DomStats(2,iz) = DomS2
    DomStats(3,iz) = DomS3
    DomStats(4,iz) = float(NumPoints)
  enddo

  return
end subroutine DoTurbMmts

!*****************************************************************************
! Calc2dMean()
!
! This function will find the mean for the given ix, iy, iz point. The mean
! is taken from the surrounding N points about (ix,iy) on the level given
! by iz. When called for each (x,y) point on a level, this effective does
! a 2D running average over that level.
!
! 
!
real function Calc2dMean(Nx, Ny, Nz, Var, ix, iy, iz, xs, xe, ys, ye)
  integer :: Nx, Ny, Nz, ix, iy, iz, xs, xe, ys, ye
  real, dimension(Nx,Ny,Nz) :: Var

  integer :: i, j, Fi, Fj
  integer :: Fsize
  real :: Npts

  ! use a 5x5 box for computing the average
  Fsize = 2
  Npts = 25.0

  Calc2dMean = 0.0
  do Fj = iy-Fsize, iy+Fsize
    do Fi = ix-Fsize, ix+Fsize
      ! If the box extends beyond the borders of Var, then repeat
      ! the value of Var from the border.
      if (Fi .lt. xs) then
         i = xs
      else if (Fi .gt. xe) then
         i = xe
      else
         i = Fi
      endif
      if (Fj .lt. ys) then
         j = ys
      else if (Fj .gt. ye) then
         j = ye
      else
         j = Fj
      endif

      Calc2dMean = Calc2dMean + Var(i,j,iz)
    enddo
  enddo
  
  Calc2dMean = Calc2dMean / Npts

end function Calc2dMean

!*****************************************************************************
! FindBin()
!
! This function will find the bin that a given data value belongs to. Emulate
! the matlab binning where the bin values are treated like edges. Val belongs
! to a bin if:
!
!   Bins(ib) <= Val < Bins(i+1)
!
! except for the last bin where Val belongs to it if:
!
!   Bins(ib) == Val
!
integer function FindBin(Nb, Bins, Val)
  implicit none

  integer :: Nb
  real, dimension(Nb) :: Bins
  real :: Val

  integer :: ib

  ! if Val doesn't fall into any bins, FindBin will remain -1
  FindBin = -1
  do ib = 1, Nb
    if (ib .lt. Nb) then
      if ((Bins(ib) .le. Val) .and. (Val .lt. Bins(ib+1))) then
        FindBin = ib
        exit
      endif
    else
      !last bin
      if (Bins(ib) .eq. Val) then
        FindBin = ib
      endif
    endif
  enddo

  return
end function FindBin

!************************************************************************
! UvToSpeed()
!
! This function will convert the U, V fields into speed (vector magnitude)
! on the lowest above ground model level.

subroutine UvToSpeed(Nx, Ny, Nz, U, V, HorizSpeed)
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx,Ny,Nz) :: U, V
  real, dimension(Nx,Ny) :: HorizSpeed

  integer :: ix, iy, iz

  ! RAMS has first model level underground, so use z = 2 (1st model level above ground)
  iz = 2
  do iy = 1, Ny
    do ix = 1, Nx
      HorizSpeed(ix,iy) = sqrt(U(ix,iy,iz)**2. + V(ix,iy,iz)**2.)
    enddo
  enddo

  return
end subroutine UvToSpeed

!!! !*****************************************************************************
!!! ! DoCloud()
!!! !
!!! ! This subroutine will perform the cloud droplet total mass time series
!!! ! averaging. Can select between supercooled or warm rain droplets.
!!! !
!!! 
!!! subroutine DoCloud(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, Dens, CintLiq, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   logical :: DoSc
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
!!!   type (GradsVar) :: Cloud, TempC, Dens, CintLiq, TserAvg
!!!   integer, dimension(1:Cloud%Nt) :: StmIx, StmIy
!!!   real, dimension(1:Cloud%Nx) :: Xcoords
!!!   real, dimension(1:Cloud%Ny) :: Ycoords
!!!   real, dimension(1:Cloud%Nz) :: Zcoords
!!! 
!!!   real, dimension(0:Cloud%Nz) :: ZmHeights
!!!   integer :: ix, iy, iz, it, NumPoints
!!! 
!!!   call SetZmHeights (Cloud%Nz, ZmHeights)
!!! 
!!!   ! Convert the cloud mixing ratio to mass for each grid point. 
!!!   !
!!!   ! Mixing ratios in GRADS files are g/kg
!!!   ! Density is in kg/m**3
!!!   ! Heights are in m
!!!   ! So, express the mass value in g using the formula
!!!   !   (mix ratio) * (density) * (layer thickness) * (layer horiz area)
!!!   ! The layer thickness for layer k is: ZmHeights(k) - ZmHeights(k-1)
!!! 
!!!   do it = 1, Cloud%Nt
!!!     ! Sum up the cloud droplet mass over the specified radial band. Only include the
!!!     ! grid points where tempc is 0 or less (supercooled)
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, Cloud%Nz
!!!       do ix = 1, Cloud%Nx
!!!         do iy = 1, Cloud%Ny
!!!           if ((InsideCylVol(Cloud%Nx, Cloud%Ny, Cloud%Nz, Cloud%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords))  .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
!!!             if (DoSc) then
!!!               if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
!!!                 TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + (Cloud%Vdata(it,iz,ix,iy) * Dens%Vdata(it,iz,ix,iy) * (ZmHeights(iz)-ZmHeights(iz-1)))
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             else
!!!               if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
!!!                 TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + (Cloud%Vdata(it,iz,ix,iy) * Dens%Vdata(it,iz,ix,iy) * (ZmHeights(iz)-ZmHeights(iz-1)))
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     ! At this point TserAvg%Vdata(it,1,1,1) holds g/m**2, multiply by grid cell horizontal area. Note this assumes
!!!     ! each grid cell has the same horizontal area.
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) * DeltaX * DeltaY
!!!     if (NumPoints .eq. 0) then
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       write (*,*) 'ScCloud: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     endif
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoWup()
!!! !
!!! ! This subroutine will do the average vertical velocity in regions of significant
!!! ! updrafts.
!!! !
!!! 
!!! subroutine DoWup(DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, W, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: W, TserAvg
!!!   integer, dimension(1:W%Nt) :: StmIx, StmIy
!!!   real, dimension(1:W%Nx) :: Xcoords
!!!   real, dimension(1:W%Ny) :: Ycoords
!!!   real, dimension(1:W%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it, NumPoints
!!! 
!!!   do it = 1, W%Nt
!!!     ! Average w over regions where significant updrafts occur
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, W%Nz
!!!       do ix = 1, W%Nx
!!!         do iy = 1, W%Ny
!!!           if (InsideCylVol(W%Nx, W%Ny, W%Nz, W%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             if (W%Vdata(it,iz,ix,iy) .ge. Wthreshold) then
!!!               TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + W%Vdata(it,iz,ix,iy)
!!!               NumPoints = NumPoints + 1
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) / float(NumPoints)
!!!       write (*,*) 'Wup: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!   end do
!!! 
!!! end subroutine
!!! 
!!! 
!!! !************************************************************************************
!!! ! DoCloudDiam()
!!! !
!!! ! This routine will calculate an averaged cloud droplet diameter. Can select between
!!! ! warm rain droplets or supercooled droplets.
!!! 
!!! subroutine DoCloudDiam(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, DoSc, Cloud, TempC, CloudDiam, CintLiq, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   logical :: DoSc
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
!!!   type (GradsVar) :: Cloud, TempC, CloudDiam, CintLiq, TserAvg
!!!   integer, dimension(1:CloudDiam%Nt) :: StmIx, StmIy
!!!   real, dimension(1:CloudDiam%Nx) :: Xcoords
!!!   real, dimension(1:CloudDiam%Ny) :: Ycoords
!!!   real, dimension(1:CloudDiam%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it, NumPoints
!!!   real SumQ, SumQD
!!!   real MaxQ, Climit, SumD
!!! 
!!!   do it = 1, CloudDiam%Nt
!!!     ! Calculate a mass-weighted mean diameter for supercooled cloud droplets.
!!!     ! 
!!!     !    Mean diameter (TserAvg value) = Sum(cloud * cloud_d) / Sum(cloud)
!!!     !    where cloud and cloud_d are only included in the sum when tempc is <= 0
!!!     !
!!!     ! Find the max Q (mass) value and use it to filter data - select data points if
!!!     ! the Q is within 20% of the max Q (.2 to 1.0). Then form the average of the diameters
!!!     ! of the selected points.
!!!     SumQD = 0.0
!!!     SumQ = 0.0
!!!     SumD = 0.0
!!!     NumPoints = 0
!!!     do iz = 1, CloudDiam%Nz
!!!       do ix = 1, CloudDiam%Nx
!!!         do iy = 1, CloudDiam%Ny
!!!           if ((InsideCylVol(CloudDiam%Nx, CloudDiam%Ny, CloudDiam%Nz, CloudDiam%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
!!!             if (DoSc) then
!!!               if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
!!!                  SumQD = SumQD + (Cloud%Vdata(it,iz,ix,iy) * CloudDiam%Vdata(it,iz,ix,iy))
!!!                  SumQ = SumQ + Cloud%Vdata(it,iz,ix,iy)
!!!                  SumD = SumD + CloudDiam%Vdata(it,iz,ix,iy)
!!!                  NumPoints = NumPoints + 1
!!!               end if
!!!             else
!!!               if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
!!!                  SumQD = SumQD + (Cloud%Vdata(it,iz,ix,iy) * CloudDiam%Vdata(it,iz,ix,iy))
!!!                  SumQ = SumQ + Cloud%Vdata(it,iz,ix,iy)
!!!                  SumD = SumD + CloudDiam%Vdata(it,iz,ix,iy)
!!!                  NumPoints = NumPoints + 1
!!!               end if
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (SumQ .eq. 0.0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = SumQD / SumQ
!!!       write (*,*) 'ScCloudDiam: Ts:', it, ', NumPoints: ', NumPoints
!!!     end if
!!! 
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoCloudConc()
!!! !
!!! ! This subroutine will do the average cloud droplet concentration
!!! 
!!! subroutine DoCloudConc(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, DoSc, TempC, CloudConc, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: CloudConc, TempC, TserAvg
!!!   integer, dimension(1:CloudConc%Nt) :: StmIx, StmIy
!!!   real, dimension(1:CloudConc%Nx) :: Xcoords
!!!   real, dimension(1:CloudConc%Ny) :: Ycoords
!!!   real, dimension(1:CloudConc%Nz) :: Zcoords
!!!   logical :: DoSc
!!! 
!!!   integer ix,iy,iz,it
!!!   integer NumPoints
!!!   real SumCloudConc
!!! 
!!!   ! Calculate the average cloud droplet concentration near the eyewall region.
!!!   ! 
!!!   ! Call the cloud base to be between 1000m and 1200m. The z level corresponding to
!!!   ! that is z = 4 which is 1138m. Cover from surface to the cloud base level. This
!!!   ! results in using z = 1 to z = 4.
!!! 
!!!   do it = 1, CloudConc%Nt
!!!     SumCloudConc = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, CloudConc%Nz
!!!       do ix = 1, CloudConc%Nx
!!!         do iy = 1, CloudConc%Ny
!!!           if (InsideCylVol(CloudConc%Nx, CloudConc%Ny, CloudConc%Nz, CloudConc%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             if (DoSc) then
!!!               if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
!!!                 SumCloudConc = SumCloudConc + CloudConc%Vdata(it,iz,ix,iy)
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             else
!!!               if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
!!!                 SumCloudConc = SumCloudConc + CloudConc%Vdata(it,iz,ix,iy)
!!!                 NumPoints = NumPoints + 1
!!!               end if
!!!             end if
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = SumCloudConc / float(NumPoints)
!!!       write (*,*) 'CloudConc: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoPrecipR()
!!! !
!!! ! This subroutine will do the average cloud droplet concentration near the eyewall
!!! ! region.
!!! 
!!! subroutine DoPrecipR(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CilThresh, PrecipR, CintLiq, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, CilThresh
!!!   ! PrecipR is 2D var
!!!   type (GradsVar) :: PrecipR, CintLiq, TserAvg
!!!   integer, dimension(1:PrecipR%Nt) :: StmIx, StmIy
!!!   real, dimension(1:PrecipR%Nx) :: Xcoords
!!!   real, dimension(1:PrecipR%Ny) :: Ycoords
!!!   real, dimension(1:PrecipR%Nz) :: Zcoords
!!! 
!!!   integer :: ix,iy,iz,it
!!!   integer :: NumPoints
!!!   real :: SumPrecip
!!! 
!!!   do it = 1, PrecipR%Nt
!!!     SumPrecip = 0.0
!!!     NumPoints = 0
!!! 
!!!     do ix = 1, PrecipR%Nx
!!!       do iy = 1, PrecipR%Ny
!!!         if (InsideCylVol(PrecipR%Nx, PrecipR%Ny, PrecipR%Nz, PrecipR%Nt, ix, iy, 1, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords) .and. (CintLiq%Vdata(it,1,ix,iy) .ge. CilThresh)) then
!!!           SumPrecip = SumPrecip + PrecipR%Vdata(it,1,ix,iy)
!!!           NumPoints = NumPoints + 1
!!!         end if
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) '  WARNING: no data points selected for this time step'
!!!     else
!!!       ! At this point TserAvg%Vdata(it,1,1,1) holds mm/hr, multiply by grid cell horizontal area.
!!!       ! Note this assumes each grid cell has the same horizontal area. What this does
!!!       ! is convert mm/hr to kg/hr when assuming that the density of water is 1000kg/m**3.
!!!       !   mm/hr * m**2 * 1000 kg/m**3 * 0.001 m/mm -> kg/hr
!!!       TserAvg%Vdata(it,1,1,1) = (SumPrecip / float(NumPoints)) * DeltaX * DeltaY
!!!       write (*,*) 'Precip: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!     write (*,*) ''
!!!     flush(6)
!!!   end do
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoCcnConc()
!!! !
!!! ! This subroutine will do the average CCN concentration.
!!! !
!!! 
!!! subroutine DoCcnConc(DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, CcnConc, TserAvg)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, Wthreshold, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: CcnConc, TserAvg
!!!   integer, dimension(1:CcnConc%Nt) :: StmIx, StmIy
!!!   real, dimension(1:CcnConc%Nx) :: Xcoords
!!!   real, dimension(1:CcnConc%Ny) :: Ycoords
!!!   real, dimension(1:CcnConc%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it, NumPoints
!!! 
!!!   do it = 1, CcnConc%Nt
!!!     ! Average w over regions where significant updrafts occur
!!! 
!!!     TserAvg%Vdata(it,1,1,1) = 0.0
!!!     NumPoints = 0
!!! 
!!!     do iz = 1, CcnConc%Nz
!!!       do ix = 1, CcnConc%Nx
!!!         do iy = 1, CcnConc%Ny
!!!           if (InsideCylVol(CcnConc%Nx, CcnConc%Ny, CcnConc%Nz, CcnConc%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) + CcnConc%Vdata(it,iz,ix,iy)
!!!             NumPoints = NumPoints + 1
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!! 
!!!     if (NumPoints .eq. 0) then
!!!       TserAvg%Vdata(it,1,1,1) = 0.0
!!!       write (*,*) 'WARNING: no data points selected for time step: ', it
!!!     else
!!!       TserAvg%Vdata(it,1,1,1) = TserAvg%Vdata(it,1,1,1) / float(NumPoints)
!!!       write (*,*) 'Wup: Timestep:', it, ', Number of points selected: ', NumPoints
!!!     end if
!!!   end do
!!! 
!!! end subroutine
!!! 
!!! !************************************************************************************
!!! ! DoTestCvs()
!!! !
!!! ! This subroutine will perform a test on the cylindrical volume selection routine.
!!! ! Just runs through the entire grid and outputs a '1' when selection occurs otherwise
!!! ! outputs a '0'. Then view the result in grads and see if selection is correct.
!!! 
!!! subroutine DoTestCvs(DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords, TestSelect)
!!!   use gdata_utils
!!!   use azavg_utils
!!!   implicit none
!!! 
!!!   real :: DeltaX, DeltaY, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ
!!!   type (GradsVar) :: TestSelect
!!!   integer, dimension(1:TestSelect%Nt) :: StmIx, StmIy
!!!   real, dimension(1:TestSelect%Nx) :: Xcoords
!!!   real, dimension(1:TestSelect%Ny) :: Ycoords
!!!   real, dimension(1:TestSelect%Nz) :: Zcoords
!!! 
!!!   integer ix,iy,iz,it
!!!   integer NumPoints
!!! 
!!!   write (*,*) 'Testing cylindrical volume selection:'
!!!   do it = 1, TestSelect%Nt
!!!     NumPoints = 0
!!!     do iz = 1, TestSelect%Nz
!!!       do ix = 1, TestSelect%Nx
!!!         do iy = 1, TestSelect%Ny
!!!           if (InsideCylVol(TestSelect%Nx, TestSelect%Ny, TestSelect%Nz, TestSelect%Nt, ix, iy, iz, it, MinR, MaxR, MinPhi, MaxPhi, MinZ, MaxZ, StmIx, StmIy, Xcoords, Ycoords, Zcoords)) then
!!!             TestSelect%Vdata(it,iz,ix,iy) = 1.0
!!!             NumPoints = NumPoints + 1
!!!           else
!!!             TestSelect%Vdata(it,iz,ix,iy) = 0.0
!!!           end if
!!!         end do
!!!       end do
!!!     end do
!!!     ! mark the storm center
!!!     TestSelect%Vdata(it,iz,StmIx(it),StmIy(it)) = 2.0
!!!     write (*,*) '  Timestep, Number of points selected: ', it, NumPoints
!!!   end do
!!! end subroutine
!!! 
!!! !******************************************************************************
!!! ! ConvertTinc()
!!! !
!!! ! This function will convert the time increment spec'd in the GRAD control
!!! ! file into a number of seconds.
!!! !
!!! 
!!! real function ConvertTinc(Tinc)
!!!   implicit none
!!! 
!!!   character (len=*) :: Tinc
!!! 
!!!   character (len=128) :: Tval, Tunits, Tfmt
!!!   integer :: i, Uloc, Tlen, Itval
!!!   
!!!   ! Walk through the Tinc string. Concatenate the numbers onto Tval and
!!!   ! the alaph characters onto Tunits. This algorithm assumes that GRADS
!!!   ! will use a format like <numeric_value><units> for Tinc where <units>
!!!   ! is a string like 'hr' or 'mn'.
!!!   Uloc = 0
!!!   Tlen = len_trim(Tinc)
!!!   do i = 1, Tlen
!!!     if ((Tinc(i:i) .ge. 'A') .and. (Tinc(i:i) .le. 'z')) then
!!!       if (Uloc .eq. 0) then
!!!         ! the first alpha character is the beginning of the spec for units
!!!         Uloc = i
!!!       end if
!!!     end if
!!!   end do
!!! 
!!!   write (Tfmt, '(a,i2.2,a,i2.2,a)') '(a', (Uloc-1), 'a' , ((Tlen-Uloc)+1), ')'
!!!   read (Tinc, Tfmt) Tval, Tunits
!!! 
!!!   read(Tval, '(i)') Itval
!!!   if (Tunits .eq. 'hr') then
!!!     ConvertTinc = float(Itval) * 3600.0
!!!   else
!!!     if (Tunits .eq. 'mn') then
!!!       ConvertTinc = float(Itval) * 60.0
!!!     else
!!!       ConvertTinc = float(Itval)
!!!     end if
!!!   end if
!!!   return
!!! end function
!!! 
!!! !******************************************************************************
!!! ! CalcRates()
!!! !
!!! ! This routine will calculate time derivatives of the input data using
!!! ! a centered difference method.
!!! !
!!! 
!!! subroutine CalcRates(DeltaT, TserAvg, Rates)
!!!   use gdata_utils
!!!   implicit none
!!! 
!!!   type (GradsVar) :: TserAvg, Rates
!!!   real :: DeltaT
!!! 
!!!   real :: f1, f2
!!!   integer :: ix, iy, iz, it
!!! 
!!!   ! use a centered difference, uses points at t-1, t and t+1
!!!   do ix = 1, TserAvg%Nx
!!!     do iy = 1, TserAvg%Ny
!!!       do iz = 1, TserAvg%Nz
!!!         do it = 2, TserAvg%Nt-1
!!!           f1 = (TserAvg%Vdata(it,1,1,1) + TserAvg%Vdata(it-1,1,1,1)) / 2.0
!!!           f2 = (TserAvg%Vdata(it+1,1,1,1) + TserAvg%Vdata(it,1,1,1)) / 2.0
!!! 
!!!           Rates%Vdata(it-1,1,1,1) = (f2 - f1) / DeltaT
!!!         end do
!!!       end do
!!!     end do
!!!   end do
!!! end subroutine

end program tsavg
