!***************************************************************
! Program to generate advection terms around the periphery of
! a volume in 3D space.
!

program gen_adv
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
    real :: Zmin
    real :: Zmax
    type (FileSpec) :: Output
    type (FileSpec) :: Filter
    type (FileSpec) :: U
    type (FileSpec) :: V
    type (FileSpec) :: W
    type (FileSpec) :: Var
  endtype

  type (ArgList) :: Args

  integer :: Nx
  integer :: Ny
  integer :: Nz
  integer :: Nt
  integer :: FilterNz

  real, dimension(:), allocatable :: XcoordsKm
  real, dimension(:), allocatable :: YcoordsKm

  integer :: ix
  integer :: iy
  integer :: iz
  integer :: it

  integer :: Ktop
  integer :: Kbot

  ! variables
  type (Rhdf5Var) :: U, V, W, Var
  type (Rhdf5Var) :: Filter, Advect
  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords, FilterZcoords
  type (Rhdf5Var) :: AdvectXcoords, AdvectYcoords, AdvectZcoords
  character (len=RHDF5_MAX_STRING) :: rh5f_facc
  integer :: rh5f_filter, rh5f_u, rh5f_v, rh5f_w, rh5f_var, rh5f_adv

  ! Get the command line arguments
  call GetMyArgs(Args)

  write (*,'("Generating advection terms:")')
  write (*,'("    Filter: ",a," (",a,")")') trim(Args%Filter%fname), trim(Args%Filter%vname)
  write (*,'("    U: ",a," (",a,")")') trim(Args%U%fname), trim(Args%U%vname)
  write (*,'("    V: ",a," (",a,")")') trim(Args%V%fname), trim(Args%V%vname)
  write (*,'("    W: ",a," (",a,")")') trim(Args%W%fname), trim(Args%W%vname)
  write (*,'("    Var: ",a," (",a,")")') trim(Args%Var%fname), trim(Args%Var%vname)
  write (*,'("    Output: ",a," (",a,")")') trim(Args%Output%fname), trim(Args%Output%vname)
  write (*,*) ''

  ! Read the variable information from the HDF5 files and check for consistency.
  Filter%vname = trim(Args%Filter%vname)
  U%vname = trim(Args%U%vname)
  V%vname = trim(Args%V%vname)
  W%vname = trim(Args%W%vname)
  Var%vname = trim(Args%Var%vname)
  Advect%vname = trim(Args%Output%vname)

  Xcoords%vname = 'x_coords'
  Ycoords%vname = 'y_coords'
  Zcoords%vname = 'z_coords'
  FilterZcoords%vname = 'z_coords'
  Tcoords%vname = 't_coords'

  call rhdf5_read_init(Args%Filter%fname, Filter)
  call rhdf5_read_init(Args%U%fname, U)
  call rhdf5_read_init(Args%V%fname, V)
  call rhdf5_read_init(Args%W%fname, W)
  call rhdf5_read_init(Args%Var%fname, Var)

  call rhdf5_read_init(Args%Var%fname, Xcoords)
  call rhdf5_read_init(Args%Var%fname, Ycoords)
  call rhdf5_read_init(Args%Var%fname, Zcoords)
  call rhdf5_read_init(Args%Filter%fname, FilterZcoords)
  call rhdf5_read_init(Args%Var%fname, Tcoords)

  ! check that the variable dimensions (size and coordinate values) match up, if this
  ! isn't true, then the subsequent anlysis will be bogus
  if (.not. (DimsMatch(Filter, U) .and. DimsMatch(Filter, V) .and. &
             DimsMatch(Filter, W) .and. DimsMatch(Filter, Var))) then
    write (*,*) 'ERROR: dimensions of filter, u, v, w, and var do not match'
    stop
  endif
 
  ! Avar needs to be 3D
  if (Var%ndims .eq. 4) then
    ! (x,y,z,t)
    Nx = Var%dims(1)
    Ny = Var%dims(2)
    Nz = Var%dims(3)
    Nt = Var%dims(4)
  else
    write (*,*) 'ERROR: <var_to_average> must be a 3D field'
    stop
  endif

  FilterNz = Filter%dims(3)

  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:               ', Nx
  write (*,*) '  Number of y (latitude) points:                ', Ny
  write (*,*) '  Number of z (vertical level) points (var):    ', Nz
  write (*,*) '  Number of z (vertical level) points (filter): ', FilterNz
  write (*,*) '  Number of t (time) points:                    ', Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''

  write (*,*) ''

  ! Read in coordinate data
  call rhdf5_read(Args%Var%Fname, Xcoords)
  call rhdf5_read(Args%Var%Fname, Ycoords)
  call rhdf5_read(Args%Var%Fname, Zcoords)
  call rhdf5_read(Args%Var%Fname, Tcoords)

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

  ! Set the z range indicies for the analysis region
  !  Default is cued by receiving a negative value in either of zmin or zmax
  !  Default is k = 2 and k = Nz - 1 for bottom and top
  !  If zmin is spec'd, then use first k level that is >= to zmin for the bottom
  !  If zmax is spec'd, then use last k level that is <= to zmax for the top
  !  When searching for Kbot and Ktop, make sure that range is restricted to: 2 through Nz-1
  Kbot = 2
  Ktop = Nz - 1

  if (Args%Zmin .ge. 0.0) then
    do iz = 2, Nz-1
      if (Zcoords%vdata(iz)/1000.0 .ge. Args%Zmin) then
        exit
      endif
    enddo
    Kbot = iz
  endif
  
  if (Args%Zmax .ge. 0.0) then
    do iz = Nz-1, 2, -1
      if (Zcoords%vdata(iz)/1000.0 .le. Args%Zmax) then
        exit
      endif
    enddo
    Ktop = iz
  endif

  write (*,*) 'Height range for analysis region:'
  write (*,*) '  Bottom: height (km), k: ', Zcoords%vdata(Kbot)/1000.0, Kbot 
  write (*,*) '  Top: height (km), k: ', Zcoords%vdata(Ktop)/1000.0, Ktop 
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
  call rhdf5_open_file(Args%Filter%fname, rh5f_facc, 0, rh5f_filter)
  call rhdf5_open_file(Args%U%fname, rh5f_facc, 0, rh5f_u)
  call rhdf5_open_file(Args%V%fname, rh5f_facc, 0, rh5f_v)
  call rhdf5_open_file(Args%W%fname, rh5f_facc, 0, rh5f_w)
  call rhdf5_open_file(Args%Var%fname, rh5f_facc, 0, rh5f_var)

  Filter%ndims = Filter%ndims - 1
  U%ndims = U%ndims - 1
  V%ndims = V%ndims - 1
  W%ndims = W%ndims - 1
  Var%ndims = Var%ndims - 1

  ! set up the output variable
  rh5f_facc = 'W'
  call rhdf5_open_file(Args%Output%fname, rh5f_facc, 1, rh5f_adv)
  write (*,*) 'Writing HDF5 output: ', trim(Args%Output%fname)
  write (*,*) ''

  ! The results of the advection calculation will be the individual sums and
  ! counts associated with each of the six faces of the analysis region defined
  ! by filter.
  !
  !  x(1)  - west face
  !  x(2)  - east face
  !  x(3)  - south face
  !  x(4)  - north face
  !  x(5)  - bottom face
  !  x(6)  - top face
  !
  !  y(1)  - sum of advective terms
  !  y(2)  - count of grids cells used in sum of advective terms
  !
  !  z     - height
  !  t     - time

  Advect%vname = trim(Args%Output%vname)
  Advect%ndims = 3
  Advect%dims(1) = 6
  Advect%dims(2) = 2
  Advect%dims(3) = 1
  Advect%dimnames(1) = 'x'
  Advect%dimnames(2) = 'y'
  Advect%dimnames(3) = 'z'
  Advect%units = Var%units 
  Advect%descrip = 'advection of ' // trim(Args%Output%vname) 
  allocate(Advect%vdata(Advect%dims(1)*Advect%dims(2)*Advect%dims(3)))

  ! output coordinates
  AdvectXcoords%vname = 'x_coords'
  AdvectXcoords%ndims = 1
  AdvectXcoords%dims(1) = 6
  AdvectXcoords%dimnames(1) = 'x'
  AdvectXcoords%units = 'degrees_east'
  AdvectXcoords%descrip = 'dummy coordinates'
  allocate (AdvectXcoords%vdata(6))
  do ix = 1, 6
    AdvectXcoords%vdata(ix) = ix
  enddo

  AdvectYcoords%vname = 'y_coords'
  AdvectYcoords%ndims = 1
  AdvectYcoords%dims(1) = 2
  AdvectYcoords%dimnames(1) = 'y'
  AdvectYcoords%units = 'degrees_north'
  AdvectYcoords%descrip = 'dummy coordinates'
  allocate (AdvectYcoords%vdata(2))
  do iy = 1, 2
    AdvectYcoords%vdata(iy) = iy
  enddo

  AdvectZcoords%vname = 'z_coords'
  AdvectZcoords%ndims = 1
  AdvectZcoords%dims(1) = 1
  AdvectZcoords%dimnames(1) = 'z'
  AdvectZcoords%units = 'meter'
  AdvectZcoords%descrip = 'dummy coordinates'
  allocate (AdvectZcoords%vdata(2))
  do iz = 1, 1
    AdvectZcoords%vdata(iz) = iz
  enddo

  ! Do the averaging - one time step at a time
  do it = 1, Nt
    ! Read in the filter, u, v, w, and variable data
    call rhdf5_read_variable(rh5f_filter, Filter%vname, Filter%ndims, it, Filter%dims, rdata=Filter%vdata)
    call rhdf5_read_variable(rh5f_u, U%vname, U%ndims, it, U%dims, rdata=U%vdata)
    call rhdf5_read_variable(rh5f_v, V%vname, V%ndims, it, V%dims, rdata=V%vdata)
    call rhdf5_read_variable(rh5f_w, W%vname, W%ndims, it, W%dims, rdata=W%vdata)
    call rhdf5_read_variable(rh5f_var, Var%vname, Var%ndims, it, Var%dims, rdata=Var%vdata)

    ! sum up the advection terms and write out the results to the output file
    ! pass in x, y and z coordinate values in units of meters
    call SumAdvectionTerms(Nx, Ny, Nz, FilterNz, XcoordsKm*1000.0, YcoordsKm*1000.0, Zcoords%vdata, &
                           Filter%vdata, U%vdata, V%vdata, W%vdata, Var%vdata, Advect%vdata, Kbot, Ktop);

    call rhdf5_write_variable(rh5f_adv, Advect%vname, Advect%ndims, it, Advect%dims, &
      Advect%units, Advect%descrip, Advect%dimnames, rdata=Advect%vdata)

    ! Free up variable memory
    deallocate(Filter%vdata)
    deallocate(U%vdata)
    deallocate(V%vdata)
    deallocate(W%vdata)
    deallocate(Var%vdata)
    
    ! Write out status to screen every 100 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,50) .eq. 0) then
      write (*,*) 'Working: Timestep, Time: ', it, Tcoords%vdata(it)
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
  call rhdf5_close_file(rh5f_filter)
  call rhdf5_close_file(rh5f_u)
  call rhdf5_close_file(rh5f_v)
  call rhdf5_close_file(rh5f_w)
  call rhdf5_close_file(rh5f_var)

  call rhdf5_close_file(rh5f_adv)

  ! write out the coordinate data
  call rhdf5_write(Args%Output%fname, AdvectXcoords, 1)
  call rhdf5_write(Args%Output%fname, AdvectYcoords, 1)
  call rhdf5_write(Args%Output%fname, AdvectZcoords, 1)
  call rhdf5_write(Args%Output%fname, Tcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(Args%Output%fname, AdvectXcoords, 'x')
  call rhdf5_set_dimension(Args%Output%fname, AdvectYcoords, 'y')
  call rhdf5_set_dimension(Args%Output%fname, Zcoords, 'z')
  call rhdf5_set_dimension(Args%Output%fname, Tcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(Args%Output%fname, Advect)

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
  Args%Zmin = -1.0
  Args%Zmax = -1.0

  ! initialization
  Args%Output%fname = 'none'
  Args%Output%vname = 'none'

  Args%Filter%fname = 'none'
  Args%Filter%vname = 'none'

  Args%U%fname = 'none'
  Args%U%vname = 'none'

  Args%V%fname = 'none'
  Args%V%vname = 'none'

  Args%W%fname = 'none'
  Args%W%vname = 'none'

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
    optval = getopt('z:')

    select case (optval)
      case ('>')  ! finished
        exit

      case ('!')  ! unrecognized argument
        write(*,*) 'ERROR: unknown option: ', trim(optarg)
        write(*,*) ''
        BadArgs = .true.

      case ('z')
         ! height range spec -> Zmin:Zmax
         call String2List(optarg, ':', ArgList, MaxArgFields, Nfields, 'height spec')
         read(ArgList(1), '(f)') Args%Zmin
         read(ArgList(2), '(f)') Args%Zmax

      case ('.')
        ! file spec ->   <file_type>:<file_name>:<variable_name>
        call String2List(optarg, ':', ArgList, MaxArgFields, Nfields, 'file spec') 
        if (Nfields .eq. 3) then
          select case (ArgList(1))
            case ('filter')
              Args%Filter%fname = trim(ArgList(2))
              Args%Filter%vname = trim(ArgList(3))

            case ('var')
              Args%Var%fname = trim(ArgList(2))
              Args%Var%vname = trim(ArgList(3))

            case ('u')
              Args%U%fname = trim(ArgList(2))
              Args%U%vname = trim(ArgList(3))

            case ('v')
              Args%V%fname = trim(ArgList(2))
              Args%V%vname = trim(ArgList(3))

            case ('w')
              Args%W%fname = trim(ArgList(2))
              Args%W%vname = trim(ArgList(3))

            case ('out')
              Args%Output%fname = trim(ArgList(2))
              Args%Output%vname = trim(ArgList(3))

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
  if ((trim(Args%Filter%fname) .eq. 'none') .or. (trim(Args%Output%fname) .eq. 'none') .or. &
      (trim(Args%U%fname) .eq. 'none') .or. (trim(Args%V%fname) .eq. 'none') .or. &
      (trim(Args%W%fname) .eq. 'none') .or. (trim(Args%Var%fname) .eq. 'none')) then
    write(*,*) 'ERROR: must always specify filter, u, v, w, var, and out files'
    write(*,*) ''
    BadArgs = .true.
  endif

  if (BadArgs) then
    write (*,*) 'USAGE: gen_adv [-z <height_range>] <list_of_files>'
    write (*,*) '   -z <height_range>: define bottom and top of analysis region'
    write (*,*) '     <height_range>: zmin:zmax in km'
    write (*,*) '        default is k = 2 and k = Nz - 1 levels'
    write (*,*) ''
    write (*,*) '   <list_of_files>: supply specs for input and output files'
    write (*,*) '     <list_of_files> := <file_spec> ...'
    write (*,*) '       <file_spec> := <file_type>:<file_name>:<variable_name>'
    write (*,*) '          <file_type> is one of:'
    write (*,*) '             "filter": filter mask, required'
    write (*,*) '             "var": variable to average, required'
    write (*,*) '             "u": zonal component of momentum, required'
    write (*,*) '             "v": meridional component of momentum, required'
    write (*,*) '             "w": vertical component of momentum, required'
    write (*,*) '             "out": output file, required'
    write (*,*) '          <file_name> is the complete path to the file'
    write (*,*) '          <variable_name> is the name of the variable inside the file'
    write (*,*) ''
    stop
  endif

  return
end subroutine GetMyArgs

!****************************************************************************
! SumAdvectionTerms
!
! This routine will calculate advection terms for each grid cell that is touching
! one of the 6 faces of the volume defined by Filter. For each of the 6 faces, the
! advection terms will be counted and summed. These results will then be output into
! Advect.
!
! Entries in Advect
!
!  x(1)  - west face
!  x(2)  - east face
!  x(3)  - south face
!  x(4)  - north face
!  x(5)  - bottom face
!  x(6)  - top face
!
!  y(1)  - sum of advective terms
!  y(2)  - count of grids cells used in sum of advective terms
!
! Ktop and Kbot have been restricted to the range 2 through Nz-1 so no need to check this.

subroutine SumAdvectionTerms(Nx, Ny, Nz, FilterNz, X, Y, Z, Filter, U, V, W, Var, Advect, Kbot, Ktop);
  implicit none

  integer :: Nx, Ny, Nz, FilterNz
  real, dimension(Nx) :: X
  real, dimension(Ny) :: Y
  real, dimension(Nz) :: Z
  real, dimension(Nx,Ny,FilterNz) :: Filter
  real, dimension(Nx,Ny,Nz) :: U, V, W, Var
  real, dimension(6,2,1) :: Advect
  integer :: Kbot, Ktop

  integer :: ix, iy, iz, filter_iz
  logical :: OnWest, OnEast, OnSouth, OnNorth, OnBot, OnTop
  real :: Var1, Var2
  real :: DeltaX, DeltaY, DeltaZ, DeltaVar

  Advect = 0.0

  do iz = Kbot, Ktop
    if (FilterNz .eq. 1) then
      filter_iz = 1
    else
      filter_iz = iz
    endif

    ! The active region of the domain is from 2 to Nx-1, and 2 to Ny-1.
    !
    ! Allow the horizontal edges (x,y directions) to vary, but assume that Kbot and Ktop
    ! remain constant across the entire top and bottom faces.
    !
    ! In x direction, indices 1 and Nx are boundaries which may have 1's in the filter.
    ! So if sitting on index 2, you want this to be called the West edge if the
    ! filter has a 1 regardless  of the value in index 1 so test for this case.  Same for
    ! index Nx-1 and identidification of the East edge.
    !
    ! In the x direction the RAMS staggered grid looks like:
    !             T  U  T  U  T  U  T  U  T  U  T  U
    ! index:      1  1  2  2  3  3  4  4  5  5  6  6
    !
    !   T represents Var
    !
    ! In this example Nx = 6. The domain goes from 2 through 5, and indices 1 and 6 are
    ! the boundaries. Note that the U value at index 6 is extraneous - it's only included
    ! so that the arrays for T and U are consistent size. For the east edge you want the
    ! value of U from the cell to the left, U(ix-1), and for the west edge you want the
    ! value of U from the current cell, U(ix). Then you want to average the surrounding T
    ! values onto the location of U.
    !
    ! If the edge of the filter region is the same as the edge of the active domain, then
    ! zero-gradient boundary conditions will force the advection terms to zero. Don't want
    ! this on lateral boundaries so the solution is to calculate the advection between the
    ! cell you are on and the adjacent cell that is toward the center of the domain
    ! (away from the boundary). This can be automatically handled by constraining the
    ! x and y indicies to the range: 3 to N-2. Then the code below that does the advection
    ! calculation will use the correct edge.
    ! 
    ! If Kbot is 2, we would actually want the zero advection values, so it's okay to allow
    ! Kbot to be 2. It doesn't make sense to place Ktop at Nz-1 since this would likely be
    ! in the sponge layer.

    do iy = 3, Ny-2
      do ix = 3, Nx-2
        ! check if on west or east edges
        if (ix .eq. 3) then
          OnWest = (anint(Filter(ix,iy,filter_iz)) .eq. 1.0)
          OnEast = .false.
        elseif (ix .eq. Nx-2) then
          OnWest = .false.
          OnEast = (anint(Filter(ix,iy,filter_iz)) .eq. 1.0)
        else
          OnWest = ((anint(Filter(ix,iy,filter_iz)) .eq. 1.0) .and. &
                    (anint(Filter(ix-1,iy,filter_iz)) .eq. 0.0))
          OnEast = ((anint(Filter(ix,iy,filter_iz)) .eq. 1.0) .and. &
                    (anint(Filter(ix+1,iy,filter_iz)) .eq. 0.0))
        endif

        ! check if on south or north edges
        if (iy .eq. 3) then
          OnSouth = (anint(Filter(ix,iy,filter_iz)) .eq. 1.0)
          OnNorth = .false.
        elseif (iy .eq. Ny-2) then
          OnSouth = .false.
          OnNorth = (anint(Filter(ix,iy,filter_iz)) .eq. 1.0)
        else
          OnSouth = ((anint(Filter(ix,iy,filter_iz)) .eq. 1.0) .and. &
                    (anint(Filter(ix,iy-1,filter_iz)) .eq. 0.0))
          OnNorth = ((anint(Filter(ix,iy,filter_iz)) .eq. 1.0) .and. &
                    (anint(Filter(ix,iy+1,filter_iz)) .eq. 0.0))
        endif

        ! check if on bottom or top edges
        OnBot = ((iz .eq. Kbot) .and. (anint(Filter(ix,iy,filter_iz)) .eq. 1.0))
        OnTop = ((iz .eq. Ktop) .and. (anint(Filter(ix,iy,filter_iz)) .eq. 1.0))

        ! If on any of the 6 faces, calculate the grid cell advection term
        ! add to the running sum and count.
        ! If at a point on the west face (x-axis) and:
        !   u > 0, dM/dx > 0 makes their product > 0, but there is a lesser amount
        !   of M entering the box. Ie, a positive u*dM/dx represents decreasing
        !   M so add a negative sign to the term: -u*dM/dx.
        !
        !   But if this were taking place on the east face, it would mean a lesser
        !   amount leaving the box. Ie, a postitive u*dM/dx represents increasing
        !   M so do not add teh negative sign.

        ! X - axis
        if (OnWest) then
          DeltaX = X(ix) - X(ix-1)
          DeltaVar = Var(ix,iy,iz) - Var(ix-1,iy,iz)

          Advect(1,1,1) = Advect(1,1,1) - (U(ix-1,iy,iz) * (DeltaVar/DeltaX))
          Advect(1,2,1) = Advect(1,2,1) + 1.0
        endif

        if (OnEast) then
          DeltaX = X(ix+1) - X(ix)
          DeltaVar = Var(ix+1,iy,iz) - Var(ix,iy,iz)

          Advect(2,1,1) = Advect(2,1,1) + (U(ix,iy,iz) * (DeltaVar/DeltaX))
          Advect(2,2,1) = Advect(2,2,1) + 1.0
        endif

        ! Y - axis
        if (OnSouth) then
          DeltaY = Y(iy) - Y(iy-1)
          DeltaVar = Var(ix,iy,iz) - Var(ix,iy-1,iz)

          Advect(3,1,1) = Advect(3,1,1) - (V(ix,iy-1,iz) * (DeltaVar/DeltaY))
          Advect(3,2,1) = Advect(3,2,1) + 1.0
        endif

        if (OnNorth) then
          DeltaY = Y(iy+1) - Y(iy)
          DeltaVar = Var(ix,iy+1,iz) - Var(ix,iy,iz)

          Advect(4,1,1) = Advect(4,1,1) + (V(ix,iy,iz) * (DeltaVar/DeltaY))
          Advect(4,2,1) = Advect(4,2,1) + 1.0
        endif

        ! Z - axis
        if (OnBot) then
          DeltaZ = Z(iz) - Z(iz-1)
          DeltaVar = Var(ix,iy,iz) - Var(ix,iy,iz-1)

          Advect(5,1,1) = Advect(5,1,1) - (W(ix,iy,iz-1) * (DeltaVar/DeltaZ))
          Advect(5,2,1) = Advect(5,2,1) + 1.0
        endif

        if (OnTop) then
          DeltaZ = Z(iz+1) - Z(iz)
          DeltaVar = Var(ix,iy,iz) - Var(ix,iy,iz+1)

          Advect(6,1,1) = Advect(6,1,1) + (W(ix,iy,iz) * (DeltaVar/DeltaZ))
          Advect(6,2,1) = Advect(6,2,1) + 1.0
        endif

      enddo
    enddo
  enddo

  return
end subroutine SumAdvectionTerms

end program gen_adv
