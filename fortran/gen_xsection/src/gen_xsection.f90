!***************************************************************
! Program to do extract a cross section out of 2D and 3D fields
!

program gen_xsection
  use rhdf5_utils
  use diag_utils

  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10
  integer, parameter :: MaxCoords=1000
  real, parameter :: UndefVal=-999.0
  integer, parameter :: MaxArgFields = 20

  type FileSpec
     character (len=RHDF5_MAX_STRING) :: fname
     character (len=RHDF5_MAX_STRING) :: vname
  endtype

  type ArgList
    logical :: DoHorizVel
    logical :: DoTangential
    real :: Lon1
    real :: Lon2
    real :: Lat1
    real :: Lat2
    type (FileSpec) :: Output
    type (FileSpec) :: U
    type (FileSpec) :: V
    type (FileSpec) :: Var
  endtype

  type (ArgList) :: Args

  ! Data arrays: need one for w (vertical velocity), press (pressure)
  ! and the var we are doing the averaging on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number

  type (Rhdf5Var) :: U, V, Var, Xcoords, Ycoords, Zcoords, Tcoords
  type (Rhdf5Var) :: Lcoords, Xsection, Dcoords
  real, dimension(:,:,:), allocatable :: Vt, Vr
  character (len=RHDF5_MAX_STRING) :: rh5f_facc
  integer :: rh5f_u, rh5f_v, rh5f_var, rh5f_out
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm
  integer, dimension(:), allocatable :: Xindices, Yindices
  real :: LineX1, LineY1

  integer :: i

  integer :: ix, iy, iz, it
  integer :: Nx, Ny, Nz, Nt
  integer :: il, Nl

  real :: VarVal

  ! Get the command line arguments
  call GetMyArgs(Args)

  write (*,'("Extracting cross section for RAMS data:")')
  write (*,'("  Xsection line (lon,lat): (",f5.1,",",f5.1,"), (",f5.1,",",f5.1,")")') Args%Lon1, Args%Lat1, Args%Lon2, Args%Lat2
  write (*,'("  Input file specs:")')
  if (Args%DoHorizVel) then
    write (*,'("    U: ",a," (",a,")")') trim(Args%U%fname), trim(Args%U%vname)
    write (*,'("    V: ",a," (",a,")")') trim(Args%V%fname), trim(Args%V%vname)
  else
    write (*,'("    Var: ",a," (",a,")")') trim(Args%Var%fname), trim(Args%Var%vname)
  endif
  write (*,'("    Output: ",a," (",a,")")') trim(Args%Output%fname), trim(Args%Output%vname)

  write (*,*) ''

  ! Read the variable information from the HDF5 files and check for consistency.
  !
  ! Always use filter file
  !
  ! Get variable and coordinate descriptions from the filter file
  ! Get the z coordinate description from the variable file since this
  ! can change between 2D and 3D variables.
  Xcoords%vname = 'x_coords'
  Ycoords%vname = 'y_coords'
  Zcoords%vname = 'z_coords'
  Tcoords%vname = 't_coords'

  if (Args%DoHorizVel) then
    U%vname = trim(Args%U%vname)
    call rhdf5_read_init(Args%U%fname, U)

    call rhdf5_read_init(Args%U%fname, Xcoords)
    call rhdf5_read_init(Args%U%fname, Ycoords)
    call rhdf5_read_init(Args%U%fname, Zcoords)
    call rhdf5_read_init(Args%U%fname, Tcoords)

    call rhdf5_read(Args%U%fname, Xcoords)
    call rhdf5_read(Args%U%fname, Ycoords)
    call rhdf5_read(Args%U%fname, Zcoords)
    call rhdf5_read(Args%U%fname, Tcoords)

    V%vname = trim(Args%V%vname)
    call rhdf5_read_init(Args%V%fname, V)

    ! Initialize the elements in Var
    Var%ndims = U%ndims
    Var%dims = U%dims
    Var%dimnames = U%dimnames
    Var%units = 'm/s'
    if (Args%DoTangential) then
      Var%vname = 'speed_t'
      Var%descrip = 'tangential wind speed'
    else
      Var%vname = 'speed_r'
      Var%descrip = 'radial wind speed'
    endif
  else
    Var%vname = trim(Args%Var%vname)
    call rhdf5_read_init(Args%Var%fname, Var)

    call rhdf5_read_init(Args%Var%fname, Xcoords)
    call rhdf5_read_init(Args%Var%fname, Ycoords)
    call rhdf5_read_init(Args%Var%fname, Zcoords)
    call rhdf5_read_init(Args%Var%fname, Tcoords)

    call rhdf5_read(Args%Var%fname, Xcoords)
    call rhdf5_read(Args%Var%fname, Ycoords)
    call rhdf5_read(Args%Var%fname, Zcoords)
    call rhdf5_read(Args%Var%fname, Tcoords)
  endif

  ! check that the variable dimensions (size and coordinate values) match up, if this
  ! isn't true, then the subsequent anlysis will be bogus
  if (Args%DoHorizVel) then
    if (.not. DimsMatch(U,V)) then
      write (*,*) 'ERROR: dimensions of u and v do not match'
      stop
    endif
  endif
 
  ! Var is either 3D (x,y,z,t) or 2D (x,y,t)
  if (Var%ndims .eq. 4) then
    ! (x,y,z,t)
    Nx = Var%dims(1)
    Ny = Var%dims(2)
    Nz = Var%dims(3)
    Nt = Var%dims(4)
  elseif (Var%ndims .eq. 3) then
    ! (x,y,t)
    Nx = Var%dims(1)
    Ny = Var%dims(2)
    Nz = 1
    Nt = Var%dims(3)
  else
    write (*,*) 'ERROR: input variable must possess either 2D or 3D spatial dimensions'
    stop
  endif

  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:               ', Nx
  write (*,*) '  Number of y (latitude) points:                ', Ny
  write (*,*) '  Number of z (vertical level) points (var):    ', Nz
  write (*,*) '  Number of t (time) points:                    ', Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''

  write (*,*) ''

  ! Convert the GRADS grid coordinates from longitude, latitude to flat plane (x and y).
  ! XcoordsKm, YcoordsKm are in units of km
  allocate(XcoordsKm(Nx))
  allocate(YcoordsKm(Ny))
  call ConvertGridCoords(Nx, Ny, Xcoords%vdata, Ycoords%vdata, XcoordsKm, YcoordsKm)

  write (*,*) 'Horzontal Grid Coordinate Info:'
  write (*,*) '  X Range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', Xcoords%vdata(1), Xcoords%vdata(Nx), XcoordsKm(1), XcoordsKm(Nx)
  write (*,*) '  Y Range (min lat, max lat) --> (min y, max y): '
  write (*,*) '    ', Ycoords%vdata(1), Ycoords%vdata(Ny), YcoordsKm(1), YcoordsKm(Ny)
  write (*,*) ''

  ! Form list of indices that represent the cross section line.
  ! CalcLineIndices will also calculate the corresponding distances for the coordinate
  ! values, and will allocate memory for Lcoords%vdata.
  ! 
  call CalcLineIndices(Nx, Ny, Xcoords%vdata, Ycoords%vdata, Args%Lon1, Args%Lat1, &
                       Args%Lon2, Args%Lat2, XcoordsKm, YcoordsKm, Nl, Xindices, Yindices, Lcoords%vdata)

  ! Finish of variable specs for the line coordinates
  Lcoords%vname = 'x_coords'
  Lcoords%ndims = 1
  Lcoords%dims(1) = Nl
  Lcoords%dimnames(1) = 'x'
  Lcoords%units = 'degrees_east'
  Lcoords%descrip  = 'linear distance'

  ! Open files for read and write. The variables in the HDF5 files include the
  ! time dimension which is always the last dimension. We want to read and write
  ! one time step at at time which requires the variables to not indlude the time
  ! dimension so we need to remove the time dimension from the variables we are
  ! using. This turns out to be easy - since time is always the last dimension, we
  ! simply decrement the number of dimensions by one to remove time.
  !
  ! Chop off the time dimensions of the variables that have time dimensions
  !
  ! set up the input variable data
  rh5f_facc = 'R'
  if (Args%DoHorizVel) then
    call rhdf5_open_file(Args%U%fname, rh5f_facc, 0, rh5f_u)
    call rhdf5_open_file(Args%V%fname, rh5f_facc, 0, rh5f_v)
    U%ndims = U%ndims - 1
    V%ndims = V%ndims - 1
  else
    call rhdf5_open_file(Args%Var%fname, rh5f_facc, 0, rh5f_var)
  endif

  Var%ndims = Var%ndims - 1

  ! set up the output variable
  rh5f_facc = 'W'
  call rhdf5_open_file(Args%Output%fname, rh5f_facc, 1, rh5f_out)
  write (*,*) 'Writing HDF5 output: ', trim(Args%Output%fname)
  write (*,*) ''

  Xsection%vname = trim(Args%Output%vname)
  Xsection%ndims = 3
  Xsection%dims(1) = Nl
  Xsection%dims(2) = 1
  Xsection%dims(3) = Nz
  Xsection%dimnames(1) = 'x'
  Xsection%dimnames(2) = 'y'
  Xsection%dimnames(3) = 'z'
  Xsection%units = Var%units 
  Xsection%descrip = trim(Args%Output%vname) // ' cross section'
  allocate(Xsection%vdata(Xsection%dims(1)*Xsection%dims(2)*Xsection%dims(3)))

  ! Do the selection - one time step at a time
  LineX1 = XcoordsKm(Xindices(1))  ! for tangential/radial wind speed
  LineY1 = YcoordsKm(Yindices(1))
  do it = 1, Nt
    if (Args%DoHorizVel) then
      call rhdf5_read_variable(rh5f_u, U%vname, U%ndims, it, U%dims, rdata=U%vdata)
      call rhdf5_read_variable(rh5f_v, V%vname, V%ndims, it, V%dims, rdata=V%vdata)

      ! Calculate tangential or radial wind speed relative to the starting point of
      ! the cross section line.
      allocate(Vt(Nx,Ny,Nz))
      allocate(Vr(Nx,Ny,Nz))
      call ConvertHorizVelocity(Nx, Ny, Nz, U%vdata, V%vdata, LineX1, LineY1, &
                                XcoordsKm, YcoordsKm, Vt, Vr)

      ! Free up variable memory
      deallocate(U%vdata)
      deallocate(V%vdata)
    else
      call rhdf5_read_variable(rh5f_var, Var%vname, Var%ndims, it, Var%dims, rdata=Var%vdata)
    endif

    ! Selection
    do il = 1, Nl
       ix = Xindices(il)
       iy = Yindices(il)
       do iz = 1, Nz
         if (Args%DoHorizVel) then
           if (Args%DoTangential) then
             VarVal = Vt(ix,iy,iz)
           else
             VarVal = Vr(ix,iy,iz)
           endif
         else
           VarVal = MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=Var%vdata)
         endif
         call MultiDimAssign(Nl, 1, Nz, il, 1, iz, VarVal, Var3d=Xsection%vdata)
       enddo
    enddo

    call rhdf5_write_variable(rh5f_out, Xsection%vname, Xsection%ndims, it, Xsection%dims, &
      Xsection%units, Xsection%descrip, Xsection%dimnames, rdata=Xsection%vdata)

    ! Free up variable memory
    if (Args%DoHorizVel) then
      deallocate(Vt)
      deallocate(Vr)
    else
      deallocate(Var%vdata)
    endif
    
    ! Write out status to screen every 10 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,10) .eq. 0) then
      write (*,*) 'Working: Timestep, Time: ', it, Tcoords%vdata(it)
      flush(6)
    endif
  enddo
  write (*,*) ''

  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''
  flush(6)

  ! close the files
  if (Args%DoHorizVel) then
    call rhdf5_close_file(rh5f_u)
    call rhdf5_close_file(rh5f_v)
  else
    call rhdf5_close_file(rh5f_var)
  endif
  call rhdf5_close_file(rh5f_out)

  ! Write out the coordinate data
  ! need to create a dummy y coordinate to keep grads happy
  Dcoords%vname = 'y_coords'
  Dcoords%ndims = 1
  Dcoords%dims(1) = 1
  Dcoords%dimnames(1) = 'y'
  Dcoords%units = 'degrees_north'
  Dcoords%descrip = 'dummy coordinates'
  allocate (Dcoords%vdata(1))
  Dcoords%vdata(1) = 1.0
  
  call rhdf5_write(Args%Output%fname, Lcoords, 1)
  call rhdf5_write(Args%Output%fname, Dcoords, 1)
  call rhdf5_write(Args%Output%fname, Zcoords, 1)
  call rhdf5_write(Args%Output%fname, Tcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(Args%Output%fname, Lcoords, 'x')
  call rhdf5_set_dimension(Args%Output%fname, Dcoords, 'y')
  call rhdf5_set_dimension(Args%Output%fname, Zcoords, 'z')
  call rhdf5_set_dimension(Args%Output%fname, Tcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(Args%Output%fname, Xsection)

  stop

Contains

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the command line arguments
!

subroutine GetMyArgs(Args)

  use getoptions

  implicit none

  type (ArgList) :: Args

  character :: optval
  logical :: BadArgs
  character (len=MediumString), dimension(MaxArgFields) :: ArgList
  integer :: Nfields

  ! default values
  Args%DoHorizVel = .false.
  Args%DoTangential = .false.

  ! initialization
  Args%Lon1 = -999.0
  Args%Lon2 = -999.0
  Args%Lat1 = -999.0
  Args%Lat2 = -999.0

  Args%Output%fname = 'none'
  Args%Output%vname = 'none'

  Args%U%fname = 'none'
  Args%U%vname = 'none'

  Args%V%fname = 'none'
  Args%V%vname = 'none'

  Args%Var%fname = 'none'
  Args%Var%vname = 'none'

  ! loop through all of the command line tokens
  ! optarg is a character string variable that the getoptions module supplies
  !   optarg gets set to the argument value for an option that uses an argument
  ! getopt returns single character:
  !    '>' finished
  !    '!' unknown option
  !    '.' command line argument (no option)
  BadArgs = .false.
  do
    optval = getopt('t:')

    select case (optval)
      case ('>')  ! finished
        exit

      case ('!')  ! unrecognized argument
        write(*,*) 'ERROR: unknown option: ', trim(optarg)
        write(*,*) ''
        BadArgs = .true.

      case ('t')
        Args%DoHorizVel = .true.
        select case (optarg)
          case ('h_tan')
            Args%DoTangential = .true.

          case ('h_rad')
            Args%DoTangential = .false.

          case default
            write(*,*) 'ERROR: must use one of "h_tan" or "h_rad" for argement to -t option: ', trim(optarg)
            write(*,*) ''
            BadArgs = .true.
        endselect

      case ('.')
        ! file spec ->   <file_type>:<file_name>:<variable_name>
        ! line spec ->   line:lon1:lat2:lon2:lat2
        call String2List(optarg, ':', ArgList, MaxArgFields, Nfields, 'file spec') 
        select case (ArgList(1))
          ! line spec
          case ('line')
            if (Nfields .eq. 5) then
              read(ArgList(2), '(f)') Args%Lon1
              read(ArgList(3), '(f)') Args%Lat1
              read(ArgList(4), '(f)') Args%Lon2
              read(ArgList(5), '(f)') Args%Lat2
            else
              write(*,*) 'ERROR: must use line:lon1:lat1:lon2:lat2 for line spec: ', trim(optarg)
              write(*,*) ''
              BadArgs = .true.
            endif

          ! file specs
          case ('var')
            if (Nfields .eq. 3) then
              Args%Var%fname = trim(ArgList(2))
              Args%Var%vname = trim(ArgList(3))
            else
              write(*,*) 'ERROR: must use <file_type>:<file_name>:<variable_name> for var file spec: ', trim(optarg)
              write(*,*) ''
              BadArgs = .true.
            endif

          case ('u')
            if (Nfields .eq. 3) then
              Args%U%fname = trim(ArgList(2))
              Args%U%vname = trim(ArgList(3))
            else
              write(*,*) 'ERROR: must use <file_type>:<file_name>:<variable_name> for u file spec: ', trim(optarg)
              write(*,*) ''
              BadArgs = .true.
            endif

          case ('v')
            if (Nfields .eq. 3) then
              Args%V%fname = trim(ArgList(2))
              Args%V%vname = trim(ArgList(3))
            else
              write(*,*) 'ERROR: must use <file_type>:<file_name>:<variable_name> for v file spec: ', trim(optarg)
              write(*,*) ''
              BadArgs = .true.
            endif

          case ('out')
            if (Nfields .eq. 3) then
              Args%Output%fname = trim(ArgList(2))
              Args%Output%vname = trim(ArgList(3))
            else
              write(*,*) 'ERROR: must use <file_type>:<file_name>:<variable_name> for out file spec: ', trim(optarg)
              write(*,*) ''
              BadArgs = .true.
            endif

          case default
            write(*,*) 'ERROR: unknown <file_type>: ', trim(ArgList(1))
            write(*,*) ''
            BadArgs = .true.
        endselect
    endselect
  enddo

  ! check for required files
  if (trim(Args%Output%fname) .eq. 'none') then
    write(*,*) 'ERROR: must always specify out file'
    write(*,*) ''
    BadArgs = .true.
  endif

  if (Args%DoHorizVel) then
    if ((trim(Args%U%fname) .eq. 'none') .or. (trim(Args%V%fname) .eq. 'none')) then
      write(*,*) 'ERROR: must specify u and v files when using -t option'
      write(*,*) ''
      BadArgs = .true.
    endif
  else
    if (trim(Args%Var%fname) .eq. 'none') then
      write(*,*) 'ERROR: must specify var file when not using -t option'
      write(*,*) ''
      BadArgs = .true.
    endif
  endif

  if ((Args%Lon1 .lt. -200.0) .or. (Args%Lat1 .lt. -200.0) .or. &
      (Args%Lon2 .lt. -200.0) .or. (Args%Lat2 .lt. -200.0)) then
    write(*,*) 'ERROR: must specify line for cross section'
    write(*,*) ''
    BadArgs = .true.
  endif

  if (BadArgs) then
    write (*,*) 'USAGE: gen_xsection [-t <input_type>] <line_spec> <list_of_files>'
    write (*,*) ''
    write (*,*) '   -t: for input data processing, valid values of <input_type>:'
    write (*,*) '     "h_tan": do tangential speed of horizontal winds'
    write (*,*) '     "h_rad": do radial speed of horizontal winds'
    write (*,*) '     both of these values require u and v files'
    write (*,*) ''
    write (*,*) '   <line_spec>: horizontal line for basis of cross section'
    write (*,*) '     <line_spec> := line:lon1:lat1:lon2:lat2'
    write (*,*) '       lon1,lat1: first endpoint of line, longitude and latitude'
    write (*,*) '       lon2,lat2: second endpoint of line, longitude and latitude'
    write (*,*) ''
    write (*,*) '   <list_of_files>: supply specs for input and output files'
    write (*,*) '     <list_of_files> := <file_spec> ...'
    write (*,*) '       <file_spec> := <file_type>:<file_name>:<variable_name>'
    write (*,*) '          <file_type> is one of:'
    write (*,*) '             "var": variable to average, required except when -t is used'
    write (*,*) '             "u": zonal component of momentum, required only if -t used'
    write (*,*) '             "v": meridional component of momentum, required only if -t used'
    write (*,*) '             "out": output file, required'
    write (*,*) '          <file_name> is the complete path to the file'
    write (*,*) '          <variable_name> is the name of the variable inside the file'
    write (*,*) ''
    stop
  endif

  return
end subroutine GetMyArgs


!**********************************************************************
! CalcLineIndices()
!
! This program will calculate the 2D (horizontal domain) array indices
! that correspond to the given line (endpoints) chopped up into Nl
! uniform pieces.
!
subroutine CalcLineIndices(Nx, Ny, X, Y, X1, Y1, X2, Y2, Xkm, Ykm, Nl, Xindices, Yindices, Lcoords)
  implicit none

  integer :: Nx, Ny, Nl
  real, dimension(Nx) :: X, Xkm
  real, dimension(Ny) :: Y, Ykm
  real :: X1, Y1, X2, Y2
  integer, dimension(:), allocatable :: Xindices, Yindices
  real, dimension(:), allocatable :: Lcoords

  integer :: ix, iy, il
  integer :: ix1, ix2, iy1, iy2
  real :: xind, yind, slope, intercept
  real :: xstart, xincr

  ! Figure out the indices corresponding to the endpoints.
  ix1 = FindNearestIndex(Nx, X, X1)
  ix2 = FindNearestIndex(Nx, X, X2)
  iy1 = FindNearestIndex(Ny, Y, Y1)
  iy2 = FindNearestIndex(Ny, Y, Y2)

  ! Calulate the number of points in the line
  ! MATLAB requires the coordinates of the line (for plotting) to be
  ! strictly increasing or strictly decreasing. The 0.9 factor in the
  ! following formula is an attempt to prevent a point in the line
  ! from being repeated as a result of rounding in the algorithm.
  Nl = nint(0.9 * sqrt(float(iabs(ix2-ix1)**2 + iabs(iy2-iy1)**2)))
  
  ! Form a line equation between the two endpoint, using the index
  ! values. Then chop up the x range into Nl evenly spaced values
  ! and apply the line equation to get the y values. Then round
  ! the results to get each x,y pairs of indices.
  slope = float(iy2-iy1) / float(ix2-ix1)
  intercept = float(iy1) - (slope * float(ix1))
  xstart = float(ix1)
  xincr = float(ix2-ix1) / float(Nl-1)

  allocate(Xindices(Nl))
  allocate(Yindices(Nl))
  allocate(Lcoords(Nl))
  do il = 1, Nl
    xind = xstart + (float(il-1) * xincr)
    yind = (slope * xind) + intercept

    Xindices(il) = nint(xind)
    Yindices(il) = nint(yind)

    Lcoords(il) = sqrt((Xkm(Xindices(il))-Xkm(Xindices(1)))**2 + (Ykm(Yindices(il))-Ykm(Yindices(1)))**2)
  enddo

  return
end subroutine CalcLineIndices


!**********************************************************************
! FindNearestIndex()
!
! This routine will return the index of the 1D array element that most
! closely matches the given value.
!
integer function FindNearestIndex(Na, Array, Val)
  implicit none

  integer :: Na
  real, dimension(Na) :: Array
  real :: Val

  integer :: i
  real :: Diff, MinDiff

  FindNearestIndex = 1
  MinDiff = abs(Array(1) - Val)
  do i = 2, Na
    Diff = abs(Array(i) - Val)
    if (Diff .lt. MinDiff) then
      MinDiff = Diff
      FindNearestIndex = i
    endif
  enddo

  return
end function FindNearestIndex

end program gen_xsection
