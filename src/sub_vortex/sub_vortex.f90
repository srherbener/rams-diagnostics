!***************************************************************
! Program to subtract out a TC vortex from the mean wind flow.
!

program sub_vortex
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
    type (FileSpec) :: FieldVt
    type (FileSpec) :: FieldVr
    type (FileSpec) :: StormVt
    type (FileSpec) :: StormVr
    type (FileSpec) :: StormCenter
    type (FileSpec) :: OutVt
    type (FileSpec) :: OutVr
  endtype

  type (ArgList) :: Args

  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number

  type (Rhdf5Var) :: FieldVt, FieldVr,StormVt, StormVr
  type (Rhdf5Var) :: OutVt, OutVr
  type (Rhdf5Var) :: Radius, StormXindex, StormYindex
  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords, Rcoords
  character (len=RHDF5_MAX_STRING) :: rh5f_facc
  integer :: rh5f_fvt, rh5f_fvr, rh5f_svt, rh5f_svr, rh5f_sctr
  integer :: rh5f_ovt, rh5f_ovr
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm, RcoordsKm
  real :: StormX, StormY

  integer :: ir, it, isx, isy
  integer :: Nx, Ny, Nz, Nt, Nr

  ! Get the command line arguments
  call GetMyArgs(Args)

  write (*,'("Calculating azimuthal average for RAMS data:")')
  write (*,'("  Input file specs:")')
  write (*,'("    Field Vt: ",a," (",a,")")') trim(Args%FieldVt%fname), trim(Args%FieldVt%vname)
  write (*,'("    Field Vr: ",a," (",a,")")') trim(Args%FieldVr%fname), trim(Args%FieldVr%vname)
  write (*,'("    Storm Vt: ",a," (",a,")")') trim(Args%StormVt%fname), trim(Args%StormVt%vname)
  write (*,'("    Storm Vr: ",a," (",a,")")') trim(Args%StormVr%fname), trim(Args%StormVr%vname)
  write (*,'("    Storm Center: ",a," (",a,")")') trim(Args%StormCenter%fname), trim(Args%StormCenter%vname)
  write (*,'("    Output Vt: ",a," (",a,")")') trim(Args%OutVt%fname), trim(Args%OutVt%vname)
  write (*,'("    Output Vr: ",a," (",a,")")') trim(Args%OutVr%fname), trim(Args%OutVr%vname)

  write (*,*) ''

  ! Read the variable information from the HDF5 files and check for consistency.

  ! Vt, Vr info for the entire field and averaged values for the storm
  FieldVt%vname = Args%FieldVt%vname
  FieldVr%vname = Args%FieldVr%vname
  StormVt%vname = Args%StormVt%vname
  StormVr%vname = Args%StormVr%vname
  call rhdf5_read_init(Args%FieldVt%fname, FieldVt)
  call rhdf5_read_init(Args%FieldVr%fname, FieldVr)
  call rhdf5_read_init(Args%StormVt%fname, StormVt)
  call rhdf5_read_init(Args%StormVr%fname, StormVr)

  ! Radius info
  Radius%vname = 'radius'
  call rhdf5_read_init(Args%StormCenter%fname, Radius)

  ! check that the variable dimensions (size and coordinate values) match up, if this
  ! isn't true, then the subsequent anlysis will be bogus
  !
  ! FieldV[tr] and Radius are (x,y,z,t)
  ! StormV[tr] are (r,z,t)
  !
  if (.not. (DimsMatch(FieldVt, FieldVr) .and. DimsMatch(FieldVt, Radius))) then
    write (*,*) 'ERROR: dimensions of field_v[tr], and radius do not match'
    stop
  endif
  if (.not. (DimsMatch(StormVt, StormVr))) then
    write (*,*) 'ERROR: dimensions of storm_vt and storm_vt do not match'
    stop
  endif
 
  Nx = FieldVt%dims(1)
  Ny = FieldVt%dims(2)
  Nz = FieldVt%dims(3)
  Nt = FieldVt%dims(4)

  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:         ', Nx
  write (*,*) '  Number of y (latitude) points:          ', Ny
  write (*,*) '  Number of z (vertical level) points:    ', Nz
  write (*,*) '  Number of t (time) points:              ', Nt
  write (*,*) ''
  write (*,*) '  Number of r (radial bands) points:      ', Nr
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''

  ! Read in vector data
  ! Coordinate values (from the field Vt file).
  Xcoords%vname = 'x_coords'
  Ycoords%vname = 'y_coords'
  Zcoords%vname = 'z_coords'
  Tcoords%vname = 't_coords'
  call rhdf5_read_init(Args%FieldVt%fname, Xcoords)
  call rhdf5_read_init(Args%FieldVt%fname, Ycoords)
  call rhdf5_read_init(Args%FieldVt%fname, Zcoords)
  call rhdf5_read_init(Args%FieldVt%fname, Tcoords)
  call rhdf5_read(Args%FieldVt%fname, Xcoords)
  call rhdf5_read(Args%FieldVt%fname, Ycoords)
  call rhdf5_read(Args%FieldVt%fname, Zcoords)
  call rhdf5_read(Args%FieldVt%fname, Tcoords)

  ! Coordinate data (from the storm Vt file).
  Rcoords%vname = 'x_coords'
  call rhdf5_read_init(Args%StormVt%fname, Rcoords)
  call rhdf5_read(Args%StormVt%fname, Rcoords)
  Nr = StormVt%dims(1)

  ! Read in the (x,y) indices for the storm center
  StormXindex%vname = 'press_cent_x_index'
  StormYindex%vname = 'press_cent_y_index'
  call rhdf5_read_init(Args%StormCenter%fname, StormXindex)
  call rhdf5_read_init(Args%StormCenter%fname, StormYindex)
  call rhdf5_read(Args%StormCenter%fname, StormXindex)
  call rhdf5_read(Args%StormCenter%fname, StormYindex)

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

  ! Radius values are in km, and Rcoords values are in m
  ! Convert the Rcoords to km (instead of Radius to m) since Rcoords is
  ! much less memory.
  allocate(RcoordsKm(Nr))
  do ir = 1, Nr
    RcoordsKm(ir) = Rcoords%vdata(ir) * 1.0e-3
  enddo

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
  call rhdf5_open_file(Args%FieldVt%fname, rh5f_facc, 0, rh5f_fvt)
  call rhdf5_open_file(Args%FieldVr%fname, rh5f_facc, 0, rh5f_fvr)
  call rhdf5_open_file(Args%StormVt%fname, rh5f_facc, 0, rh5f_svt)
  call rhdf5_open_file(Args%StormVr%fname, rh5f_facc, 0, rh5f_svr)
  call rhdf5_open_file(Args%StormCenter%fname, rh5f_facc, 0, rh5f_sctr)

  FieldVt%ndims = FieldVt%ndims - 1
  FieldVr%ndims = FieldVr%ndims - 1
  StormVt%ndims = StormVt%ndims - 1
  StormVr%ndims = StormVr%ndims - 1
  Radius%ndims = Radius%ndims - 1

  ! set up the output variables
  rh5f_facc = 'W'
  call rhdf5_open_file(Args%OutVt%fname, rh5f_facc, 1, rh5f_ovt)
  call rhdf5_open_file(Args%OutVr%fname, rh5f_facc, 1, rh5f_ovr)
  write (*,*) 'Writing HDF5 output: ', trim(Args%OutVt%fname)
  write (*,*) 'Writing HDF5 output: ', trim(Args%OutVr%fname)
  write (*,*) ''

  ! NumBins will be one when not doing histograms. Otherwise it will be
  ! the number of bins the user specified.

  OutVt%vname = trim(Args%OutVt%vname)
  OutVt%ndims = 3
  OutVt%dims(1) = Nx
  OutVt%dims(2) = Ny
  OutVt%dims(3) = Nz
  OutVt%dimnames(1) = 'x'
  OutVt%dimnames(2) = 'y'
  OutVt%dimnames(3) = 'z'
  OutVt%units = FieldVt%units 
  OutVt%descrip = 'no vortex ' // trim(Args%FieldVt%vname) 
  allocate(OutVt%vdata(OutVt%dims(1)*OutVt%dims(2)*OutVt%dims(3)))

  OutVr%vname = trim(Args%OutVr%vname)
  OutVr%ndims = 3
  OutVr%dims(1) = Nx
  OutVr%dims(2) = Ny
  OutVr%dims(3) = Nz
  OutVr%dimnames(1) = 'x'
  OutVr%dimnames(2) = 'y'
  OutVr%dimnames(3) = 'z'
  OutVr%units = FieldVr%units 
  OutVr%descrip = 'no vortex ' // trim(Args%FieldVr%vname) 
  allocate(OutVr%vdata(OutVr%dims(1)*OutVr%dims(2)*OutVr%dims(3)))

  ! Do the averaging - one time step at a time
  do it = 1, Nt
    ! find storm center in km
    isx = nint(StormXindex%vdata(it))
    isy = nint(StormYindex%vdata(it))
    StormX = XcoordsKm(isx)
    StormY = YcoordsKm(isy)

    ! Read in the radius and variable data
    call rhdf5_read_variable(rh5f_sctr, Radius%vname, Radius%ndims, it, Radius%dims, rdata=Radius%vdata)
    call rhdf5_read_variable(rh5f_fvt, FieldVt%vname, FieldVt%ndims, it, FieldVt%dims, rdata=FieldVt%vdata)
    call rhdf5_read_variable(rh5f_fvr, FieldVr%vname, FieldVr%ndims, it, FieldVr%dims, rdata=FieldVr%vdata)
    call rhdf5_read_variable(rh5f_svt, StormVt%vname, StormVt%ndims, it, StormVt%dims, rdata=StormVt%vdata)
    call rhdf5_read_variable(rh5f_svr, StormVr%vname, StormVr%ndims, it, StormVr%dims, rdata=StormVr%vdata)

    call SubtractVortex(Nx,Ny,Nz,Nr,FieldVt%vdata,FieldVr%vdata,StormVt%vdata,StormVr%vdata, &
                        Radius%vdata,RcoordsKm,OutVt%vdata,OutVr%vdata)

    call rhdf5_write_variable(rh5f_ovt, OutVt%vname, OutVt%ndims, it, OutVt%dims, &
      OutVt%units, OutVt%descrip, OutVt%dimnames, rdata=OutVt%vdata)
    call rhdf5_write_variable(rh5f_ovr, OutVr%vname, OutVr%ndims, it, OutVr%dims, &
      OutVr%units, OutVr%descrip, OutVr%dimnames, rdata=OutVr%vdata)

    ! Free up variable memory
    deallocate(Radius%vdata)
    deallocate(FieldVt%vdata)
    deallocate(FieldVr%vdata)
    deallocate(StormVt%vdata)
    deallocate(StormVr%vdata)
    
    ! Write out status to screen every 100 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,50) .eq. 0) then
      write (*,*) 'Working: Timestep, Time: ', it, Tcoords%vdata(it)

      write (*,'(a,4f15.4)') '   Storm center: lon, lat, x, y: ', Xcoords%vdata(isx), Ycoords%vdata(isy), StormX, StormY
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
  call rhdf5_close_file(rh5f_fvt)
  call rhdf5_close_file(rh5f_fvr)
  call rhdf5_close_file(rh5f_svt)
  call rhdf5_close_file(rh5f_svr)
  call rhdf5_close_file(rh5f_sctr)

  call rhdf5_close_file(rh5f_ovt)
  call rhdf5_close_file(rh5f_ovr)

  ! write out the coordinate data
  call rhdf5_write(Args%OutVt%fname, Xcoords, 1)
  call rhdf5_write(Args%OutVt%fname, Ycoords, 1)
  call rhdf5_write(Args%OutVt%fname, Zcoords, 1)
  call rhdf5_write(Args%OutVt%fname, Tcoords, 1)

  call rhdf5_write(Args%OutVr%fname, Xcoords, 1)
  call rhdf5_write(Args%OutVr%fname, Ycoords, 1)
  call rhdf5_write(Args%OutVr%fname, Zcoords, 1)
  call rhdf5_write(Args%OutVr%fname, Tcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(Args%OutVt%fname, Xcoords, 'x')
  call rhdf5_set_dimension(Args%OutVt%fname, Ycoords, 'y')
  call rhdf5_set_dimension(Args%OutVt%fname, Zcoords, 'z')
  call rhdf5_set_dimension(Args%OutVt%fname, Tcoords, 't')

  call rhdf5_set_dimension(Args%OutVr%fname, Xcoords, 'x')
  call rhdf5_set_dimension(Args%OutVr%fname, Ycoords, 'y')
  call rhdf5_set_dimension(Args%OutVr%fname, Zcoords, 'z')
  call rhdf5_set_dimension(Args%OutVr%fname, Tcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(Args%OutVt%fname, OutVt)
  call rhdf5_attach_dimensions(Args%OutVr%fname, OutVr)

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

  ! initialization
  Args%FieldVt%fname = 'none'
  Args%FieldVt%vname = 'none'

  Args%FieldVr%fname = 'none'
  Args%FieldVr%vname = 'none'

  Args%StormVt%fname = 'none'
  Args%StormVt%vname = 'none'

  Args%StormVr%fname = 'none'
  Args%StormVr%vname = 'none'

  Args%StormCenter%fname = 'none'
  Args%StormCenter%vname = 'none'

  Args%OutVt%fname = 'none'
  Args%OutVt%vname = 'none'

  Args%OutVr%fname = 'none'
  Args%OutVr%vname = 'none'

  ! loop through all of the command line tokens
  ! optarg is a character string variable that the getoptions module supplies
  !   optarg gets set to the argument value for an option that uses an argument
  ! getopt returns single character:
  !    '>' finished
  !    '!' unknown option
  !    '.' command line argument (no option)
  BadArgs = .false.
  do
    optval = getopt('')

    select case (optval)
      case ('>')  ! finished
        exit

      case ('!')  ! unrecognized argument
        write(*,*) 'ERROR: unknown option: ', trim(optarg)
        write(*,*) ''
        BadArgs = .true.

      case ('.')
        ! file spec ->   <file_type>:<file_name>:<variable_name>
        call String2List(optarg, ':', ArgList, MaxArgFields, Nfields, 'file spec') 
        if (Nfields .eq. 3) then
          select case (ArgList(1))
            case ('field_vt')
              Args%FieldVt%fname = trim(ArgList(2))
              Args%FieldVt%vname = trim(ArgList(3))

            case ('field_vr')
              Args%FieldVr%fname = trim(ArgList(2))
              Args%FieldVr%vname = trim(ArgList(3))

            case ('storm_vt')
              Args%StormVt%fname = trim(ArgList(2))
              Args%StormVt%vname = trim(ArgList(3))

            case ('storm_vr')
              Args%StormVr%fname = trim(ArgList(2))
              Args%StormVr%vname = trim(ArgList(3))

            case ('storm_center')
              Args%StormCenter%fname = trim(ArgList(2))
              Args%StormCenter%vname = trim(ArgList(3))

            case ('out_vt')
              Args%OutVt%fname = trim(ArgList(2))
              Args%OutVt%vname = trim(ArgList(3))

            case ('out_vr')
              Args%OutVr%fname = trim(ArgList(2))
              Args%OutVr%vname = trim(ArgList(3))

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
  if ((trim(Args%FieldVt%fname) .eq. 'none') .or. &
      (trim(Args%FieldVr%fname) .eq. 'none') .or. &
      (trim(Args%StormVt%fname) .eq. 'none') .or. &
      (trim(Args%StormVr%fname) .eq. 'none') .or. &
      (trim(Args%StormCenter%fname) .eq. 'none') .or. &
      (trim(Args%OutVt%fname) .eq. 'none') .or. &
      (trim(Args%OutVr%fname) .eq. 'none')) then
    write(*,*) 'ERROR: must always specify field_v[tr], storm_v[tr], storm_center, and out_v[tr], and files'
    write(*,*) ''
    BadArgs = .true.
  endif

  if (BadArgs) then
    write (*,*) 'USAGE: sub_vortex <list_of_files>'
    write (*,*) ''
    write (*,*) '   <list_of_files>: supply specs for input and output files'
    write (*,*) '     <list_of_files> := <file_spec> ...'
    write (*,*) '       <file_spec> := <file_type>:<file_name>:<variable_name>'
    write (*,*) '          <file_type> is one of:'
    write (*,*) '             "field_vt": tangential velocity, required'
    write (*,*) '             "field_vr": radial velocity, required'
    write (*,*) '             "storm_vt": vortex average tangential velocity, required'
    write (*,*) '             "storm_vr": vortex average radial velocity, required'
    write (*,*) '             "storm_center": storm track and motion, required'
    write (*,*) '             "out_vt": output tangential velocity, required'
    write (*,*) '             "out_vr": output radial velocity, required'
    write (*,*) '          <file_name> is the complete path to the file'
    write (*,*) '          <variable_name> is the name of the variable inside the file'
    write (*,*) ''
    stop
  endif

  return
end subroutine GetMyArgs

!****************************************************************************
! In a storm relative sense, this routine will subtract out the average 
! (axisymmetric) storm (Vt,Vr) winds from the field (Vt,Vr) winds.
!
subroutine SubtractVortex(Nx,Ny,Nz,Nr,FieldVt,FieldVr,StormVt,StormVr,Radius,Rbins,OutVt,OutVr)
  implicit none

  integer :: Nx, Ny, Nz, Nr
  real, dimension(Nx,Ny) :: Radius
  real, dimension(Nx,Ny,Nz) :: FieldVt, FieldVr, OutVt, OutVr
  real, dimension(Nr,Nz) :: StormVt, StormVr
  real, dimension(Nr) :: Rbins

  integer :: ix, iy, iz, ir

  do iz = 1, Nz
    do iy = 1, Ny
      do ix = 1, Nx
        ! Find the index into the StormV[tr] arrays that represents
        ! the average Vt,Vr values of the vortex. These are the numbers
        ! that need to be subtracted out of the field Vt,Vr
        ir = FindBin(Nr,Rbins,Radius(ix,iy))

        if (ir .ne. -1) then
          ! Found the vortex Vt,Vr values
          OutVt(ix,iy,iz) = FieldVt(ix,iy,iz) - StormVt(ir,iz)
          OutVr(ix,iy,iz) = FieldVr(ix,iy,iz) - StormVr(ir,iz)
        else
          ! Copy field value 
          OutVt(ix,iy,iz) = FieldVt(ix,iy,iz)
          OutVr(ix,iy,iz) = FieldVr(ix,iy,iz)
        endif
      enddo
    enddo
  enddo

  return
end subroutine SubtractVortex

end program sub_vortex
