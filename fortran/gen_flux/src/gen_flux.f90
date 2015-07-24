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

program gen_flux
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

  character (len=MediumString) :: InDir
  character (len=MediumString) :: InSuffix
  character (len=MediumString) :: OutFile
  character (len=LittleString) :: VarName
  character (len=MediumString) :: FluxSpec

  character (len=MediumString), dimension(MaxArgFields) :: ArgList
  integer :: Nfields
  character (len=LittleString) :: Ftype
  character (len=MediumString) :: Units, Descrip

  type (Rhdf5Var) :: MixRatio, Velocity, Density
  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords
  type (Rhdf5Var) :: Flux
  character (len=RHDF5_MAX_STRING) :: MRfile, Vfile, Dfile
  character (len=RHDF5_MAX_STRING) :: rh5f_facc
  integer :: rh5f_d, rh5f_v, rh5f_mr, rh5f_flux

  integer :: Nx, Ny, Nz, Nt
  integer :: i, it

  ! Get the command line arguments
  call GetMyArgs(InDir, InSuffix, OutFile, VarName, FluxSpec)

  if (FluxSpec(1:5) .eq. 'mass:') then
    call String2List(FluxSpec, ':', ArgList, MaxArgFields, Nfields, 'flux spec')
    if (Nfields .eq. 7) then
      ! got the right amount of fields
      !   field    value
      !    1       'mass'
      !    2       density file prefix
      !    3       density variable name
      !    4       velocity file prefix
      !    5       velocity variable name
      !    6       mixing ratio file prefix
      !    7       mixing ratio variable name
      Ftype = 'mass'

      Dfile  = trim(InDir) // '/' // trim(ArgList(2)) // trim(InSuffix)
      Vfile  = trim(InDir) // '/' // trim(ArgList(4)) // trim(InSuffix)
      MRfile = trim(InDir) // '/' // trim(ArgList(6)) // trim(InSuffix)

      Density%vname  = trim(ArgList(3))
      Velocity%vname = trim(ArgList(5))
      MixRatio%vname = trim(ArgList(7))
    else
      write (*,*) 'ERROR: flux function "mass" requires seven fields -> mass:<dens_fprefix>:<dens_vname>:<velocity_fprefix>:<velocity_vname>:<mix_ratio_fprefix>:<mix_ratio_vname>'
      stop
    endif
  endif

  write (*,*) 'Calculating flux for RAMS data:'
  write (*,*) '  Input directory: ', trim(InDir)
  write (*,*) '  Input file name suffix: ', trim(InSuffix)
  write (*,*) '  Output file name:  ', trim(OutFile)
  write (*,*) '  Output variable name:  ', trim(VarName)
  write (*,*) ''

  ! Read the variable information from the HDF5 files and check for consistency.
  Xcoords%vname = 'x_coords'
  Ycoords%vname = 'y_coords'
  Zcoords%vname = 'z_coords'
  Tcoords%vname = 't_coords'

  if (Ftype .eq. 'mass') then
    ! read in descriptions of variables
    call rhdf5_read_init(Dfile, Density)
    call rhdf5_read_init(Vfile, Velocity)
    call rhdf5_read_init(MRfile, MixRatio)

    Nx = Density%dims(1)
    Ny = Density%dims(2)
    Nz = Density%dims(3)
    Nt = Density%dims(4)
    
    Units   = MixRatio%units
    Descrip = MixRatio%descrip

    ! get coordinates from the density file
    call rhdf5_read_init(Dfile, Xcoords)
    call rhdf5_read_init(Dfile, Ycoords)
    call rhdf5_read_init(Dfile, Zcoords)
    call rhdf5_read_init(Dfile, Tcoords)

    call rhdf5_read(Dfile, Xcoords)
    call rhdf5_read(Dfile, Ycoords)
    call rhdf5_read(Dfile, Zcoords)
    call rhdf5_read(Dfile, Tcoords)
    
    ! check that the dimensions all match up between the variables
    if (.not. (DimsMatch(Density, Velocity) .and. DimsMatch(Density, MixRatio))) then
      write (*,*) 'ERROR: dimensions of density, velocity, and mixing ratio do not match'
      stop
    endif

    ! Set up the input files. We are going to process one time step at a time
    ! so set up the input variables as if they had no time dimension (3D).
    rh5f_facc = 'R'
    call rhdf5_open_file(Dfile,  rh5f_facc, 0, rh5f_d)
    call rhdf5_open_file(Vfile,  rh5f_facc, 0, rh5f_v)
    call rhdf5_open_file(MRfile, rh5f_facc, 0, rh5f_mr)
    Density%ndims = 3
    Velocity%ndims = 3
    MixRatio%ndims = 3
    allocate(Density%vdata(Nx*Ny*Nz))
    allocate(Velocity%vdata(Nx*Ny*Nz))
    allocate(MixRatio%vdata(Nx*Ny*Nz))
  endif
 
  ! write out dimension info
  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:       ', Nx
  write (*,*) '  Number of y (latitude) points:        ', Ny
  write (*,*) '  Number of z (vertical level) points:  ', Nz
  write (*,*) '  Number of t (time) points:            ', Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''

  ! set up output var, only set up 3D field, the time dimension
  ! will get formed as we go through the loop
  Flux%vname = trim(VarName)
  Flux%ndims = 3
  Flux%dims(1) = Nx
  Flux%dims(2) = Ny
  Flux%dims(3) = Nz
  Flux%dimnames(1) = 'x'
  Flux%dimnames(2) = 'y'
  Flux%dimnames(3) = 'z'
  Flux%units = Units
  Flux%descrip = 'flux of ' // trim(Descrip) 
  allocate(Flux%vdata(Nx*Ny*Nz))
  
  rh5f_facc = 'W'
  call rhdf5_open_file(OutFile, rh5f_facc, 1, rh5f_flux)
  write (*,*) 'Writing HDF5 output: ', trim(OutFile)
  write (*,*) ''

  ! Do the flux calculation, one time step at at time so that we
  ! can work on a huge file (lots of time steps) while only requiring
  ! the memory for one copy of the 3D field.
  do it = 1, Nt
    if (Ftype .eq. 'mass') then
      ! read in density, velocity and mixing ratio and multiply them together to get the flux
      call rhdf5_read_variable(rh5f_d, Density%vname, Density%ndims, it, Density%dims, rdata=Density%vdata)
      call rhdf5_read_variable(rh5f_v, Velocity%vname, Velocity%ndims, it, Velocity%dims, rdata=Velocity%vdata)
      call rhdf5_read_variable(rh5f_mr, MixRatio%vname, MixRatio%ndims, it, MixRatio%dims, rdata=MixRatio%vdata)

      ! flux calculation - all %vdata arrays hold the 3D fields, but are organized
      ! as 1D arrays
      do i = 1, Nx*Ny*Nz
        Flux%vdata(i) = Density%vdata(i) * Velocity%vdata(i) * MixRatio%vdata(i)
      enddo
    endif

    ! Write out the flux
    call rhdf5_write_variable(rh5f_flux, Flux%vname, Flux%ndims, it, Flux%dims, &
      Flux%units, Flux%descrip, Flux%dimnames, rdata=Flux%vdata)

    ! Write out status to screen every 100 timesteps so that the user can see that a long
    ! running job is progressing okay.
    if (modulo(it,100) .eq. 0) then
      write (*,*) 'Working: Timestep, Time: ', it, Tcoords%vdata(it)
      flush(6)
    endif
  enddo

  ! 'it' will be one beyond its loop limit (Nt) so subtract one
  ! from 'it' when reporting how many times steps were processed
  write (*,*) 'Finished: Total number of time steps processed: ', it-1
  write (*,*) ''
  flush(6)

  ! close the files
  if (Ftype .eq. 'mass') then
    call rhdf5_close_file(rh5f_d)
    call rhdf5_close_file(rh5f_v)
    call rhdf5_close_file(rh5f_mr)
  endif
  call rhdf5_close_file(rh5f_flux)

  ! write out the coordinate data
  call rhdf5_write(OutFile, Xcoords, 1)
  call rhdf5_write(OutFile, Ycoords, 1)
  call rhdf5_write(OutFile, Zcoords, 1)
  call rhdf5_write(OutFile, Tcoords, 1)

  ! set up four (x,y,z,t) dimensions for use by GRADS
  call rhdf5_set_dimension(OutFile, Xcoords, 'x')
  call rhdf5_set_dimension(OutFile, Ycoords, 'y')
  call rhdf5_set_dimension(OutFile, Zcoords, 'z')
  call rhdf5_set_dimension(OutFile, Tcoords, 't')

  ! attach the dimension specs to the output variable
  call rhdf5_attach_dimensions(OutFile, Flux)

  stop

Contains

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the command line arguments
!

subroutine GetMyArgs(InDir, InSuffix, OutFile, VarName, FluxSpec)

  implicit none

  character (len=*) :: InDir, InSuffix, OutFile, VarName, FluxSpec

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 5) then
    write (*,*) 'ERROR: must supply exactly 5 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: gen_flux <in_dir> <in_suffix> <out_file> <var_name> <flux_spec>'
    write (*,*) '         <in_dir>: directory where input files live'
    write (*,*) '         <in_suffix>: suffix on input file names'
    write (*,*) '         <out_file>: name of output file, HDF5 format'
    write (*,*) '         <var_name> is name of output variable'
    write (*,*) '         <flux_spec>: specification of flux calculation'
    write (*,*) '            <flux_spec>: <type>:<varlist>'
    write (*,*) '               <type> is one of: "mass"'
    write (*,*) '               <varlist> specifies the required input variables'
    write (*,*) '                  colon separated pairs of <fprefix>:<vname>'
    write (*,*) '                     <fprefix>: prefix of REVU file'
    write (*,*) '                     <vname>: name of variable in the REVU file'
    write (*,*) ''
    stop
  end if

  call getarg(1, InDir)
  call getarg(2, InSuffix)
  call getarg(3, OutFile)
  call getarg(4, VarName)
  call getarg(5, FluxSpec)

  if (FluxSpec(1:5) .ne. 'mass:') then
    write (*,*) 'ERROR: must use one of the following for the flux type:'
    write (*,*) '         mass'
    stop
  endif

  return
end subroutine GetMyArgs

end program gen_flux
