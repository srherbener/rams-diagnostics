!***************************************************************
! Program to generate storm relative winds (Vt, Vr) from
! zonal, meridional winds (u,v).
!

program gen_vt_vr
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
    logical :: SubtractSmotion
    logical :: ReverseConv
    type (FileSpec) :: Vt
    type (FileSpec) :: Vr
    type (FileSpec) :: StormCenter
    type (FileSpec) :: U
    type (FileSpec) :: V
  endtype

  type (ArgList) :: Args
  integer :: Uv2Tr

  ! Data arrays: need one for w (vertical velocity), press (pressure)
  ! and the var we are doing the averaging on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number

  type (Rhdf5Var) :: U, V, Ustorm, Vstorm, Vt, Vr, Xcoords, Ycoords, Zcoords, Tcoords
  type (Rhdf5Var) :: StormXindex, StormYindex
  character (len=RHDF5_MAX_STRING) :: rh5f_facc
  integer :: rh5f_center, rh5f_u, rh5f_v, rh5f_vt, rh5f_vr
  integer :: i_storm, j_storm
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm
  real :: StormX, StormY

  integer :: i
  integer :: AvarNelems

  integer :: ix, iy, it
  integer :: Nx, Ny, Nz, Nt

  ! Get the command line arguments
  call GetMyArgs(Args)

  if (Args%ReverseConv) then
    ! do (vt,vr) to (u,v) conversion
    Uv2Tr = 0
  else
    ! do (u,v) to (vt,vr) conversion
    Uv2Tr = 1
  endif

  write (*,'("Calculating storm relative winds:")')
  write (*,'("  Subtract storm motion from input winds: ",l)') Args%SubtractSmotion
  write (*,'("  Reverse conversion, (vt,vr) to (u,v): ",l)') Args%ReverseConv
  write (*,'("  Input file specs:")')
  write (*,'("    U: ",a," (",a,")")') trim(Args%U%fname), trim(Args%U%vname)
  write (*,'("    V: ",a," (",a,")")') trim(Args%V%fname), trim(Args%V%vname)
  write (*,'("    Vt: ",a," (",a,")")') trim(Args%Vt%fname), trim(Args%Vt%vname)
  write (*,'("    Vr: ",a," (",a,")")') trim(Args%Vr%fname), trim(Args%Vr%vname)
  write (*,'("    Storm center: ",a," (",a,")")') trim(Args%StormCenter%fname), trim(Args%StormCenter%vname)
  write (*,*) ''

  ! Read the variable information from the HDF5 files and check for consistency.
  StormXindex%vname = 'press_cent_x_index'
  StormYindex%vname = 'press_cent_y_index'
  Xcoords%vname = 'x_coords'
  Ycoords%vname = 'y_coords'
  Zcoords%vname = 'z_coords'
  Tcoords%vname = 't_coords'

  call rhdf5_read_init(Args%StormCenter%fname, StormXindex)
  call rhdf5_read_init(Args%StormCenter%fname, StormYindex)

  call rhdf5_read(Args%StormCenter%fname, StormXindex)
  call rhdf5_read(Args%StormCenter%fname, StormYindex)

  U%vname = trim(Args%U%vname)
  V%vname = trim(Args%V%vname)
  call rhdf5_read_init(Args%U%fname, U)
  call rhdf5_read_init(Args%V%fname, V)
  call rhdf5_read_init(Args%U%fname, Xcoords)
  call rhdf5_read_init(Args%U%fname, Ycoords)
  call rhdf5_read_init(Args%U%fname, Zcoords)
  call rhdf5_read_init(Args%U%fname, Tcoords)

  ! Read in coordinate values, only for copying to output files
  call rhdf5_read(Args%U%fname, Xcoords)
  call rhdf5_read(Args%U%fname, Ycoords)
  call rhdf5_read(Args%U%fname, Zcoords)
  call rhdf5_read(Args%U%fname, Tcoords)

  ! If subtracting storm motion, initialize the storm u and v vars
  if (Args%SubtractSmotion) then
    Ustorm%vname = 'storm_speed_x'
    call rhdf5_read_init(Args%StormCenter%fname, Ustorm)

    Vstorm%vname = 'storm_speed_y'
    call rhdf5_read_init(Args%StormCenter%fname, Vstorm)
  endif

  ! check that the variable dimensions (size and coordinate values) match up, if this
  ! isn't true, then the subsequent anlysis will be bogus
  if (.not. DimsMatch(U, V)) then
    write (*,*) 'ERROR: dimensions of u and v do not match'
    stop
  endif
 
  Nx = U%dims(1)
  Ny = U%dims(2)
  Nz = U%dims(3)
  Nt = U%dims(4)

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

  ! Open files for read and write. The variables in the HDF5 files include the
  ! time dimension which is always the last dimension. We want to read and write
  ! one time step at at time which requires the variables to not indlude the time
  ! dimension so we need to remove the time dimension from the variables we are
  ! using. This turns out to be easy - since time is always the last dimension, we
  ! simply decrement the number of dimensions by one to remove time.
  !
  ! Chop off the time dimensions of the variables that have time dimensions
  !
  rh5f_facc = 'R'
  call rhdf5_open_file(Args%StormCenter%fname, rh5f_facc, 0, rh5f_center)

  ! set up the input variable data
  rh5f_facc = 'R'
  call rhdf5_open_file(Args%U%fname, rh5f_facc, 0, rh5f_u)
  call rhdf5_open_file(Args%V%fname, rh5f_facc, 0, rh5f_v)
  U%ndims = U%ndims - 1
  V%ndims = V%ndims - 1

  ! set up the output variable
  rh5f_facc = 'W'
  call rhdf5_open_file(Args%Vt%fname, rh5f_facc, 1, rh5f_vt)
  call rhdf5_open_file(Args%Vr%fname, rh5f_facc, 1, rh5f_vr)
  write (*,*) 'Writing: ', trim(Args%Vt%fname)
  write (*,*) 'Writing: ', trim(Args%Vr%fname)
  write (*,*) ''

  Vt%vname = trim(Args%Vt%vname)
  Vt%ndims = 3
  Vt%dims(1) = Nx
  Vt%dims(2) = Ny
  Vt%dims(3) = Nz
  Vt%dimnames(1) = 'x'
  Vt%dimnames(2) = 'y'
  Vt%dimnames(3) = 'z'
  Vt%units = U%units 
  Vt%descrip = 'Storm Vt'
  allocate(Vt%vdata(Vt%dims(1)*Vt%dims(2)*Vt%dims(3)))

  Vr%vname = trim(Args%Vr%vname)
  Vr%ndims = 3
  Vr%dims(1) = Nx
  Vr%dims(2) = Ny
  Vr%dims(3) = Nz
  Vr%dimnames(1) = 'x'
  Vr%dimnames(2) = 'y'
  Vr%dimnames(3) = 'z'
  Vr%units = U%units 
  Vr%descrip = 'Storm Vr'
  allocate(Vr%vdata(Vr%dims(1)*Vr%dims(2)*Vr%dims(3)))

  ! If subtracting out storm motion, read in the storm u and v
  ! Read in entire dataset (time step number == 0, Argument 4 to rhdf5_read_variable). Because
  ! we are doing this outside the time step loop, don't alter the number of dimensions on [UV]storm.
  if (Args%SubtractSmotion) then
    call rhdf5_read_variable(rh5f_center, Ustorm%vname, Ustorm%ndims, 0, Ustorm%dims, rdata=Ustorm%vdata)
    call rhdf5_read_variable(rh5f_center, Vstorm%vname, Vstorm%ndims, 0, Vstorm%dims, rdata=Vstorm%vdata)
  endif

  ! Do the averaging - one time step at a time
  do it = 1, Nt
    ! find storm center in km
    i_storm = nint(StormXindex%vdata(it))
    j_storm = nint(StormYindex%vdata(it))
    StormX = XcoordsKm(i_storm)
    StormY = YcoordsKm(j_storm)

    call rhdf5_read_variable(rh5f_u, U%vname, U%ndims, it, U%dims, rdata=U%vdata)
    call rhdf5_read_variable(rh5f_v, V%vname, V%ndims, it, V%dims, rdata=V%vdata)

    ! If subtracting storm motion, adjust u and v
    if (Args%SubtractSmotion) then
      call TranslateField(Nx, Ny, Nz, U%vdata, Ustorm%vdata(it))
      call TranslateField(Nx, Ny, Nz, V%vdata, Vstorm%vdata(it))
    endif

    ! convert u,v to tangential,radial (last arg to ConvertHorizVelocity indicates
    ! direction of conversion: 1 -> (u,v) to (vt,vr), otherwise (vt,vr) to (u,v)
    call ConvertHorizVelocity(Nx, Ny, Nz, U%vdata, V%vdata, StormX, StormY, &
      XcoordsKm, YcoordsKm, Vt%vdata, Vr%vdata, Uv2Tr)

    ! Free up variable memory
    deallocate(U%vdata)
    deallocate(V%vdata)

    ! write out the results to the output files
    call rhdf5_write_variable(rh5f_vt, Vt%vname, Vt%ndims, it, Vt%dims, &
      Vt%units, Vt%descrip, Vt%dimnames, rdata=Vt%vdata)
    call rhdf5_write_variable(rh5f_vr, Vr%vname, Vr%ndims, it, Vr%dims, &
      Vr%units, Vr%descrip, Vr%dimnames, rdata=Vr%vdata)

    ! Write out status to screen every 100 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,50) .eq. 0) then
      write (*,*) 'Working: Timestep, Time: ', it, Tcoords%vdata(it)

      write (*,'(a,4f15.4)') '   Storm center: lon, lat, x, y: ', Xcoords%vdata(i_storm), Ycoords%vdata(j_storm), StormX, StormY
      if (Args%SubtractSmotion) then
        write (*,'(a,4f15.4)') '   Storm motion: u, v: ', Ustorm%vdata(it), Vstorm%vdata(it)
      endif
      write (*,*) ''
      flush(6)
    endif
  enddo

  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''
  flush(6)

  ! close the files
  call rhdf5_close_file(rh5f_center)
  call rhdf5_close_file(rh5f_u)
  call rhdf5_close_file(rh5f_v)
  call rhdf5_close_file(rh5f_vt)
  call rhdf5_close_file(rh5f_vr)

  ! write out the coordinate data
  call rhdf5_write(Args%Vt%fname, Xcoords, 1)
  call rhdf5_write(Args%Vt%fname, Ycoords, 1)
  call rhdf5_write(Args%Vt%fname, Zcoords, 1)
  call rhdf5_write(Args%Vt%fname, Tcoords, 1)

  ! write out the coordinate data
  call rhdf5_write(Args%Vr%fname, Xcoords, 1)
  call rhdf5_write(Args%Vr%fname, Ycoords, 1)
  call rhdf5_write(Args%Vr%fname, Zcoords, 1)
  call rhdf5_write(Args%Vr%fname, Tcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(Args%Vt%fname, Xcoords, 'x')
  call rhdf5_set_dimension(Args%Vt%fname, Ycoords, 'y')
  call rhdf5_set_dimension(Args%Vt%fname, Zcoords, 'z')
  call rhdf5_set_dimension(Args%Vt%fname, Tcoords, 't')

  call rhdf5_set_dimension(Args%Vr%fname, Xcoords, 'x')
  call rhdf5_set_dimension(Args%Vr%fname, Ycoords, 'y')
  call rhdf5_set_dimension(Args%Vr%fname, Zcoords, 'z')
  call rhdf5_set_dimension(Args%Vr%fname, Tcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(Args%Vt%fname, Vt)
  call rhdf5_attach_dimensions(Args%Vr%fname, Vr)

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
  Args%SubtractSmotion = .false.
  Args%ReverseConv     = .false.

  ! initialization
  Args%Vt%fname = 'none'
  Args%Vr%vname = 'none'

  Args%StormCenter%fname = 'none'
  Args%StormCenter%vname = 'none'

  Args%U%fname = 'none'
  Args%U%vname = 'none'

  Args%V%fname = 'none'
  Args%V%vname = 'none'

  ! loop through all of the command line tokens
  ! optarg is a character string variable that the getoptions module supplies
  !   optarg gets set to the argument value for an option that uses an argument
  ! getopt returns single character:
  !    '>' finished
  !    '!' unknown option
  !    '.' command line argument (no option)
  BadArgs = .false.
  do
    optval = getopt('mr')

    select case (optval)
      case ('>')  ! finished
        exit

      case ('!')  ! unrecognized argument
        write(*,*) 'ERROR: unknown option: ', trim(optarg)
        write(*,*) ''
        BadArgs = .true.

      case ('m')
        Args%SubtractSmotion = .true.

      case ('r')
        Args%ReverseConv = .true.

      case ('.')
        ! file spec ->   <file_type>:<file_name>:<variable_name>
        call String2List(optarg, ':', ArgList, MaxArgFields, Nfields, 'file spec') 
        if (Nfields .eq. 3) then
          select case (ArgList(1))
            case ('storm')
              Args%StormCenter%fname = trim(ArgList(2))
              Args%StormCenter%vname = trim(ArgList(3))

            case ('u')
              Args%U%fname = trim(ArgList(2))
              Args%U%vname = trim(ArgList(3))

            case ('v')
              Args%V%fname = trim(ArgList(2))
              Args%V%vname = trim(ArgList(3))

            case ('vt')
              Args%Vt%fname = trim(ArgList(2))
              Args%Vt%vname = trim(ArgList(3))

            case ('vr')
              Args%Vr%fname = trim(ArgList(2))
              Args%Vr%vname = trim(ArgList(3))

            case default
              write(*,*) 'ERROR: unknown <file_type>: ', trim(ArgList(1))
              write(*,*) ''
              BadArgs = .true.
          endselect
        else
          write(*,*) 'ERROR: must use <file_type>:<file_name>:<variable_name> for file specs: ', trim(optarg)
          write(*,*) ''
          BadArgs = .true.
        endif
    endselect
  enddo

  ! check for required files
  if ((trim(Args%U%fname) .eq. 'none') .or. &
      (trim(Args%V%fname) .eq. 'none') .or. &
      (trim(Args%Vt%fname) .eq. 'none') .or. &
      (trim(Args%Vr%fname) .eq. 'none') .or. &
      (trim(Args%StormCenter%fname) .eq. 'none')) then
    write(*,*) 'ERROR: must specify u, v, vt, vr, and storm center files'
    write(*,*) ''
    BadArgs = .true.
  endif

  if (BadArgs) then
    write (*,*) 'USAGE: gen_vt_vr [-m] <list_of_files>'
    write (*,*) ''
    write (*,*) '   -m: subtract storm motion from input winds'
    write (*,*) '   -r: do reverse conversion (vt,vr) to (u,v)'
    write (*,*) ''
    write (*,*) '   <list_of_files>: supply specs for input and output files'
    write (*,*) '     <list_of_files> := <file_spec> ...'
    write (*,*) '       <file_spec> := <file_type>:<file_name>:<variable_name>'
    write (*,*) '          <file_type> is one of:'
    write (*,*) '             "storm": storm track and motion, required'
    write (*,*) '             "u": zonal component of momentum, required'
    write (*,*) '             "v": meridional component of momentum, required'
    write (*,*) '             "vt": output file for Vt, required'
    write (*,*) '             "vr": output file for Vr, required'
    write (*,*) '          <file_name> is the complete path to the file'
    write (*,*) '          <variable_name> is the name of the variable inside the file'
    write (*,*) ''
    stop
  endif

  return
end subroutine GetMyArgs

end program gen_vt_vr
