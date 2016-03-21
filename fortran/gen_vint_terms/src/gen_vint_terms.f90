!***************************************************************
! Program to calculate volume integration terms on a 3D field.
!

program gen_vint_terms
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
    type (FileSpec) :: Output
    type (FileSpec) :: Var
  endtype

  type (ArgList) :: Args

  integer :: Nx
  integer :: Ny
  integer :: Nz
  integer :: Nt

  real, dimension(:), allocatable :: XcoordsKm
  real, dimension(:), allocatable :: YcoordsKm

  integer :: it

  ! variables
  type (Rhdf5Var) :: Var, VintTerms
  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords
  character (len=RHDF5_MAX_STRING) :: rh5f_facc
  integer :: rh5f_var, rh5f_vint_terms

  ! Get the command line arguments
  call GetMyArgs(Args)

  write (*,'("Generating volume integration terms:")')
  write (*,'("    Var: ",a," (",a,")")') trim(Args%Var%fname), trim(Args%Var%vname)
  write (*,'("    Output: ",a," (",a,")")') trim(Args%Output%fname), trim(Args%Output%vname)
  write (*,*) ''

  ! Read the variable information from the HDF5 files and check for consistency.
  Var%vname = trim(Args%Var%vname)
  VintTerms%vname = trim(Args%Output%vname)

  Xcoords%vname = 'x_coords'
  Ycoords%vname = 'y_coords'
  Zcoords%vname = 'z_coords'
  Tcoords%vname = 't_coords'

  call rhdf5_read_init(Args%Var%fname, Var)

  call rhdf5_read_init(Args%Var%fname, Xcoords)
  call rhdf5_read_init(Args%Var%fname, Ycoords)
  call rhdf5_read_init(Args%Var%fname, Zcoords)
  call rhdf5_read_init(Args%Var%fname, Tcoords)

  ! Var needs to be 3D
  if (Var%ndims .eq. 4) then
    ! (x,y,z,t)
    Nx = Var%dims(1)
    Ny = Var%dims(2)
    Nz = Var%dims(3)
    Nt = Var%dims(4)
  else
    write (*,*) 'ERROR: input variable must be a 3D field'
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
  call rhdf5_open_file(Args%Var%fname, rh5f_facc, 0, rh5f_var)
  Var%ndims = Var%ndims - 1

  ! set up the output variable
  rh5f_facc = 'W'
  call rhdf5_open_file(Args%Output%fname, rh5f_facc, 1, rh5f_vint_terms)
  write (*,*) 'Writing HDF5 output: ', trim(Args%Output%fname)
  write (*,*) ''

  VintTerms%vname = trim(Args%Output%vname)
  VintTerms%ndims = 3
  VintTerms%dims(1) = Nx
  VintTerms%dims(2) = Ny
  VintTerms%dims(3) = Nz
  VintTerms%dimnames(1) = 'x'
  VintTerms%dimnames(2) = 'y'
  VintTerms%dimnames(3) = 'z'
  VintTerms%units = Var%units 
  VintTerms%descrip = 'vertical integration of ' // trim(Args%Output%vname) 
  allocate(VintTerms%vdata(VintTerms%dims(1)*VintTerms%dims(2)*VintTerms%dims(3)))

  ! Do the calculations one time step at a time
  do it = 1, Nt
    ! Read in the filter, u, v, w, and variable data
    call rhdf5_read_variable(rh5f_var, Var%vname, Var%ndims, it, Var%dims, rdata=Var%vdata)

    ! calculate the tems: quantity X grid volume for each grid cell
    ! pass in x, y and z coordinate values in units of meters
    call CalcVintTerms(Nx, Ny, Nz, XcoordsKm*1000.0, YcoordsKm*1000.0, Zcoords%vdata, Var%vdata, VintTerms%vdata);

    call rhdf5_write_variable(rh5f_vint_terms, VintTerms%vname, VintTerms%ndims, it, VintTerms%dims, &
      VintTerms%units, VintTerms%descrip, VintTerms%dimnames, rdata=VintTerms%vdata)

    ! Free up variable memory
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
  call rhdf5_close_file(rh5f_var)

  call rhdf5_close_file(rh5f_vint_terms)

  ! write out the coordinate data
  call rhdf5_write(Args%Output%fname, Xcoords, 1)
  call rhdf5_write(Args%Output%fname, Ycoords, 1)
  call rhdf5_write(Args%Output%fname, Zcoords, 1)
  call rhdf5_write(Args%Output%fname, Tcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(Args%Output%fname, Xcoords, 'x')
  call rhdf5_set_dimension(Args%Output%fname, Ycoords, 'y')
  call rhdf5_set_dimension(Args%Output%fname, Zcoords, 'z')
  call rhdf5_set_dimension(Args%Output%fname, Tcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(Args%Output%fname, VintTerms)

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
  Args%Output%fname = 'none'
  Args%Output%vname = 'none'

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
            case ('var')
              Args%Var%fname = trim(ArgList(2))
              Args%Var%vname = trim(ArgList(3))

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
  if ((trim(Args%Output%fname) .eq. 'none') .or. (trim(Args%Var%fname) .eq. 'none')) then
    write(*,*) 'ERROR: must always specify var and out files'
    write(*,*) ''
    BadArgs = .true.
  endif

  if (BadArgs) then
    write (*,*) 'USAGE: gen_vint_terms <list_of_files>'
    write (*,*) '   <list_of_files>: supply specs for input and output files'
    write (*,*) '     <list_of_files> := <file_spec> ...'
    write (*,*) '       <file_spec> := <file_type>:<file_name>:<variable_name>'
    write (*,*) '          <file_type> is one of:'
    write (*,*) '             "var": variable to average, required'
    write (*,*) '             "out": output file, required'
    write (*,*) '          <file_name> is the complete path to the file'
    write (*,*) '          <variable_name> is the name of the variable inside the file'
    write (*,*) ''
    stop
  endif

  return
end subroutine GetMyArgs

!************************************************************************
! CalcVintTerms
!
! This routine will calculate the terms that are used in a volume integration.
!  

subroutine CalcVintTerms(Nx, Ny, Nz, Xcoords, Ycoords, Zcoords, Var, VintTerms);
  implicit none

  integer :: Nx, Ny, Nz
  real, dimension(Nx) :: Xcoords
  real, dimension(Ny) :: Ycoords
  real, dimension(Nz) :: Zcoords
  real, dimension(Nx,Ny,Nz) :: Var, VintTerms

  integer :: ix, iy, iz
  real :: dX, dY, dZ
  real :: GridVol

  do iz = 1, Nz
    ! Use the Nz-1 interval for delta z when at the top
    ! of the domain. This is the case for 99% of RAMS configurations
    if (iz .eq. Nz) then
      dZ = Zcoords(Nz) - Zcoords(Nz-1)
    else
      dZ = Zcoords(iz+1) - Zcoords(iz)
    endif

    do iy = 1, Ny
      ! Use the Ny-1 interval for delta y when at the north edge of the
      ! of the domain.
      if (iy .eq. Ny) then
        dY = Ycoords(Ny) - Ycoords(Ny-1)
      else
        dY = Ycoords(iy+1) - Ycoords(iy)
      endif

      do ix = 1, Nx
        ! Use the Nx-1 interval for delta x when at the east edge of the
        ! of the domain.
        if (ix .eq. Nx) then
          dX = Xcoords(Nx) - Xcoords(Nx-1)
        else
          dX = Xcoords(ix+1) - Xcoords(ix)
        endif

        VintTerms(ix,iy,iz) = Var(ix,iy,iz) * dX * dY * dZ

      enddo
    enddo
  enddo
  
  return
end subroutine CalcVintTerms
  

end program gen_vint_terms
