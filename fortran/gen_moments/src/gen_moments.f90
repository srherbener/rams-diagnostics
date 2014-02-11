!***************************************************************
! Program to calculate terms (sums) used for generating
! statistical moments.
!
! Args
!   1. directory containing input files
!   2. suffix to tag onto the end of input file names
!   3. output file name 
!   4. selection of averaging function
!
! Output
!   For each time step:
!     1. number of points per level
!     2. moment terms per level
!
! This program is going to break up the formulas for moments into the individual
! summation terms so that they can be assembled as desired by a post processing
! step. Doing it this way allows for the combination of multiple time steps
! into one moment calculation. That is, the terms and number of point counts
! can simply be added together across multiple time steps and then the terms
! can be combined into the moment calculation.
!
! At this point, only up to the fourth moment is supported. The first moment
! is the mean, and subsequent moments use the definition:
!
!   Mr = sum(Xi - Xbar)^r / N
!
! where r is the degree of the moment (2, 3 or 4). Note that covariances can
! also be calculated by specifying 2 different variables to this program.
!
! For example, if doing the 3rd moment of the vertical velocity, w, you would
! specify w three times to gen_moments. Then the terms of
!
!   (Wi - Wbar)^3
!
! are actually calculated using (Wi - Wbar) * (Wi - Wbar) * (Wi - Wbar). Doing
! it the algorithm this way allows for covariance when you specify 2 different
! variables to gen_moments. Say you specify x and y, then you get terms:
! 
!   (Xi - Xbar) * (Yi - Ybar)
!
! that lead to covariance.
!
! In order to accomodate combining multiple time steps, the means are never
! calculate in this program. Instead all of the summation terms are calculated
! and output. These terms are obtained by multiplying out the expressions
! for sum(Xi - Xbar)^r. Using the notation, X1 for the first variable,
! X2 for the second, X3 for the third and X4 for the fourth: The list of terms
! for the first four moments are:
!
! 1st moment
!    Sum(X1)
!
!    where Sum(X1) means summing up all of the elements of X1(i,j)
!    across the horizontal domain
!
! 2nd moment
!    Sum(X1), Sum(X2)
!    Sum(X1 * X2)
!
!    where Sum(X1 * X2) means summing up X1(i,j) * X2(i,j)
!    across the horizontal domain
!
! 3rd moment
!   Sum(X1), Sum(X2), Sum(X3)
!   Sum(X1 * X2), Sum(X1 * X3), Sum(X2 * X3)
!   Sum(X1 * X2 * X3)
!
! 4th moment
!   Sum(X1), Sum(X2), Sum(X3), Sum(X4)
!   Sum(X1 * X2), Sum(X1 * X3), Sum(X1, X4), Sum(X2 * X3), Sum(X2 * X4), Sum(X3 * X4)
!   Sum(X1 * X2 * X3), Sum(X1 * X2 * X4), Sum(X1 * X3 * X4), Sum(X2 * X3 * X4)
!   Sum(X1 * X2 * X3 * X4)
!
! For all of these terms the associate N (number of points) is also saved
!
! Then a postprocessor, can finish the moment calculation using:
!
! 1st moment
!    M = Sum(X1) / N
!
! 2nd moment
!   M = ( Sum(X1 * X2) - (Sum(X1) * Sum(X2) / N) ) / N
!
! 3rd moment
!  M = ( Sum(X1 * X2 * X3) - ((Sum(X1)*Sum(X2*X3) + Sum(X2)*Sum(X1*X3 + Sum(X3)*Sum(X1*X2)) / N) + ( 2 * (Sum(X1*X2*X3) / N^2)) / N
!
! Output is organized as 3D array (Nterms, Order of Terms, Nz)
!  The single variable terms (X1, X2, ...) are in the first column
!  The double variable terms (X1*X2, X1*X3, X2*X3, ...) are in the second column
!  The triple variable terms are in the third column
!  The quadruple variable terms are in the fourth column
!
! The loop below that computes these terms and fills in the output variable array will order terms by varying
! the first index (denoting the variable number) the slowest. Take the exmaple of 4 input variables. This
! will create 4 single variable terms, 6 double variable terms, 4 triple variable terms and 1 quadruple
! variable term. The output array (sums, 6 rows, 4 columns) will be filled in as (just one level shown):
!
!      Col->  1        2        3            4
!   Row
!    |
!    v
!
!    1       X1      X1*X2   X1*X2*X3    X1*X2*X3*X4
!    2       X2      X1*X3   X1*X2*X4              0
!    3       X3      X1*X4   X1*X3*X4              0
!    4       X4      X2*X3   X2*X3*X4              0
!    5        0      X2*X4          0              0
!    6        0      X3*X4          0              0
!

program gen_moments
  use rhdf5_utils
  use diag_utils
  implicit none

  integer, parameter :: LargeString  = 512
  integer, parameter :: MediumString = 256
  integer, parameter :: LittleString = 128
  integer, parameter :: MaxVars = 4
  integer, parameter :: MaxTerms = 6

  real, parameter :: UndefVal = -999.0

  type VarSpec
    character (len=LittleString) :: Vname
    character (len=LittleString) :: Vfprefix
  end type VarSpec

  character (len=LargeString) :: InDir, InSuffix, OutFile, FilterFile
  type (VarSpec), dimension(MaxVars) :: VarList
  integer :: Nvars, Nterms
  logical :: UseFilter, SelectPoint

  logical :: BadDims, ZdimMatches

  integer :: id, ivar, ix, iy, iz, it, filter_z
  integer :: iv1, iv2, iv3, iv4
  integer :: iterm1, iterm2, iterm3, iterm4
  real :: v1, v2, v3, v4
  real :: s_temp
  integer :: Nx, Ny, Nz, Nt

  character (len=LargeString), dimension(MaxVars) :: InFiles, InVarNames
  type (Rhdf5Var), dimension(MaxVars) :: InVars
  integer, dimension(MaxVars) :: InFileIds

  character (len=RHDF5_MAX_STRING) :: FileAcc

  type (Rhdf5Var) :: Xcoords, Ycoords, Zcoords, Tcoords
  real, dimension(:), allocatable :: XcoordsKm, YcoordsKm

  real :: DeltaX, DeltaY

  integer :: OutFileId, FilterFileId
  type (Rhdf5Var) :: OutVar, NumPoints, Filter
  character (len=RHDF5_MAX_STRING) :: OutVarName, OutVarUnits, OutVarDescrip

  ! Get the command line arguments
  call GetMyArgs(LargeString, MaxVars, InDir, InSuffix, OutFile, FilterFile, UseFilter, VarList, Nvars, Nterms)
  do ivar = 1, Nvars
    InFiles(ivar) = trim(InDir) // '/' // trim(VarList(ivar)%Vfprefix) // trim(InSuffix)
    InVarNames(ivar) = trim(VarList(ivar)%Vname)
  enddo

  write (*,*) 'Generating moments:'
  write (*,*) '  Number of variables: ', Nvars
  write (*,*) '  Max number of terms: ', Nterms
  write (*,*) '  Input variables:'
  do ivar = 1, Nvars
    write (*,*) '    File: ', trim(InFiles(ivar))
    write (*,*) '    Variable: ', trim(InVarNames(ivar))
  enddo
  write (*,*) '  Output file:  ', trim(OutFile)
  if (UseFilter) then
    write (*,*) '  Filter file:  ', trim(FilterFile)
  else
    write (*,*) '  No filtering'
  endif
  write (*,*) ''
  flush(6)

  ! Since we will be processing one time step at a time, the effective number of dimensions
  ! on each of the input is decreased by one. We want to eliminate the time dimension which
  ! is always the last one. This works out conveniently since all we need to do is decrement
  ! the number of dimensions by one.

  BadDims = .false.
  do ivar = 1, Nvars
    InVars(ivar)%vname = trim(InVarNames(ivar))
    call rhdf5_read_init(InFiles(ivar), InVars(ivar))

    if (ivar .eq. 1) then
      ! If first variable, record the dimensions
      Nx = InVars(ivar)%dims(1)
      Ny = InVars(ivar)%dims(2)
      if (InVars(ivar)%ndims .eq. 3) then
        ! 2D field
        Nz = 1
        Nt = InVars(ivar)%dims(3)
      else
        ! 3D field
        Nz = InVars(ivar)%dims(3)
        Nt = InVars(ivar)%dims(4)
      endif
    else
      ! If beyond the first variable, check to make sure its dimensions match
      ! The DimsMatch() function only checks horzontal and time dimensions. We also want to
      ! force the z dimensions to match as well.
      if (Nz .eq. 1) then
        ! 2D vars
        ZdimMatches = InVars(ivar)%ndims .eq. 3
      else
        ! 3D vars
        ZdimMatches = Nz .eq. InVars(ivar)%dims(3)
      endif

      if (DimsMatch(InVars(1), InVars(ivar)) .and. ZdimMatches) then
      else
        write (*,*) 'ERROR: dimensions of variable do not match the first variable: ', trim(InVars(ivar)%vname)
        BadDims = .true.
      endif
    endif
  enddo

  ! check dimensions of filter if being used
  if (UseFilter) then
    Filter%vname = 'filter'
    call rhdf5_read_init(FilterFile, Filter)

   if (.not. DimsMatch(InVars(1), Filter)) then
     write (*,*) 'ERROR: dimensions of filter do not match the variables: '
     BadDims = .true.
   endif
  endif

  if (BadDims) then
    stop
  endif

  ! Prepare for reading
  ! Create names for ouput variable
  do ivar = 1, Nvars
    InVars(ivar)%ndims = InVars(ivar)%ndims - 1
 
    if (ivar .eq. 1) then
      OutVarName = trim(InVars(ivar)%vname)
      OutVarUnits = '(' // trim(InVars(ivar)%units) // ')'
      OutVarDescrip = 'moment: (' // trim(InVars(ivar)%vname) // ')'
    else
      OutVarName = trim(OutVarName) // '-' // trim(InVars(ivar)%vname)
      OutVarUnits = trim(OutVarUnits) // '(' // trim(InVars(ivar)%units) // ')'
      OutVarDescrip = trim(OutVarDescrip) // '(' // trim(InVars(ivar)%vname) // ')'
    endif
  enddo

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

  ! Set the output dimensions and coordinates to those of the selected input var
  call SetOutCoords(InFiles(1), Xcoords, Ycoords, Zcoords, Tcoords)

  ! Convert lat (x coords) and lon (y coords) to distances in km
  allocate(XcoordsKm(Nx))
  allocate(YcoordsKm(Ny))
  call ConvertGridCoords(Nx, Ny, Nz, Xcoords%vdata, Ycoords%vdata, XcoordsKm, YcoordsKm)

  DeltaX = (XcoordsKm(2) - XcoordsKm(1)) * 1000.0
  DeltaY = (YcoordsKm(2) - YcoordsKm(1)) * 1000.0

  write (*,*) 'Horizontal grid info:'
  write (*,*) '  X range (min lon, max lon) --> (min x, max x): '
  write (*,*) '    ', Xcoords%vdata(1), Xcoords%vdata(Nx), XcoordsKm(1), XcoordsKm(Nx)
  write (*,*) '  Y range (min lat, max lat) --> (min y, max y): '
  write (*,*) '    ', Ycoords%vdata(1), Ycoords%vdata(Ny), YcoordsKm(1), YcoordsKm(Ny)
  write (*,*) ''
  write (*,*) '  Delta x of domain:     ', DeltaX
  write (*,*) '  Delta y of domain:     ', DeltaY
  write (*,*) ''
  write (*,*) 'Vertical grid info:'
  do iz = 1, Nz
    write (*,*) '  ', iz, ' --> ', Zcoords%vdata(iz)
  end do
  write (*,*) ''
  flush(6)

  ! Stats
  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', Nx
  write (*,*) '  Number of y (latitude) points:           ', Ny
  write (*,*) '  Number of z (vertical level) points:     ', Nz
  write (*,*) '  Number of t (time) points:               ', Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''
  write (*,*) '  Grid delta x: ', DeltaX
  write (*,*) '  Grid delta y: ', DeltaY
  write (*,*) ''
  flush(6)

  ! Prepare output variables
  NumPoints%vname = 'num_points'
  NumPoints%units = 'count'
  NumPoints%descrip = 'number of points selected'
  NumPoints%ndims = 1
  NumPoints%dims(1) = Nz
  NumPoints%dimnames(1) = 'z'
  allocate(NumPoints%vdata(Nz))

  OutVar%vname = trim(OutVarName)
  OutVar%units = trim(OutVarUnits)
  OutVar%descrip = trim(OutVarDescrip)
  OutVar%ndims = 3
  OutVar%dims(1) = Nterms
  OutVar%dims(2) = Nvars
  OutVar%dims(3) = Nz
  OutVar%dimnames(1) = 'x'
  OutVar%dimnames(2) = 'y'
  OutVar%dimnames(3) = 'z'
  allocate(OutVar%vdata(Nterms*Nvars*Nz))

  ! Open the input files and the output file
  FileAcc = 'R'
  do ivar = 1, Nvars
    call rhdf5_open_file(InFiles(ivar), FileAcc, 0, InFileIds(ivar))
    write (*,*) 'Reading HDF5 file: ', trim(InFiles(ivar))
  enddo
  write (*,*) ''

  if (UseFilter) then
    call rhdf5_open_file(FilterFile, FileAcc, 0, FilterFileId)
    write (*,*) 'Reading HDF5 file: ', trim(FilterFile)
  endif

  FileAcc = 'W'
  call rhdf5_open_file(OutFile, FileAcc, 1, OutFileId)
  write (*,*) 'Writing HDF5 file: ', trim(OutFile)
  write (*,*) ''

  ! The necessary input files have been opened, data buffers allocated,
  ! plus the time dimensions have been "removed" from the descriptions
  ! of the input variables.
  !

  do it = 1, Nt
    ! zero out the terms (sums) an number of points counters
    do iz = 1, Nz
      NumPoints%vdata(iz) = 0.0
      do ivar = 1, Nvars
        do iterm1 = 1, Nterms
          call MultiDimAssign(Nterms, Nvars, Nz, iterm1, ivar, iz, 0.0, Var3d=OutVar%vdata)
        enddo
      enddo
    enddo

    ! Read in the input variables and filter
    do ivar = 1, Nvars
      call rhdf5_read_variable(InFileIds(ivar), InVars(ivar)%vname, InVars(ivar)%ndims, it, InVars(ivar)%dims, rdata=InVars(ivar)%vdata)
    enddo

    if (UseFilter) then
      call rhdf5_read_variable(FilterFileId, Filter%vname, Filter%ndims, it, Filter%dims, rdata=Filter%vdata)
    endif

    ! Form the sums
    ! Output is 3D array (iterm, ivar, iz)
    !   The ivar dimension represent the "order" of the term stored in that column. The order
    !   is 1 for sum(Xi) terms, the order is 2 for sum(Xi * Xj) terms, etc.
    do iz = 1, Nz
      if (Nz .eq. 1) then
        ! Use z = 2 for the filter (first model level above zero)
        filter_z = 2
      else
        filter_z = iz
      endif

      do iy = 2, Ny-1
        do ix = 2, Nx-1
          SelectPoint = .true.
          if (UseFilter) then
            SelectPoint = SelectPoint .and. (anint(MultiDimLookup(Nx, Ny, Nz, ix, iy, filter_z, Var3d=Filter%vdata)) .eq. 1.0)
          endif

          if (SelectPoint) then
            NumPoints%vdata(iz) = NumPoints%vdata(iz) + 1.0

            iterm1 = 0
            iterm2 = 0
            iterm3 = 0
            iterm4 = 0

            do iv1 = 1, Nvars
              iterm1 = iterm1 + 1
              v1 = MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=InVars(1)%vdata)
              s_temp = MultiDimLookup(Nterms, Nvars, Nz, iterm1, 1, iz, Var3d=OutVar%vdata) + v1
              call MultiDimAssign(Nterms, Nvars, Nz, iterm1, 1, iz, s_temp, Var3d=OutVar%vdata)

              do iv2 = iv1+1, Nvars
                iterm2 = iterm2 + 1
                v2 = MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=InVars(2)%vdata)
                s_temp = MultiDimLookup(Nterms, Nvars, Nz, iterm2, 2, iz, Var3d=OutVar%vdata) + (v1 * v2)
                call MultiDimAssign(Nterms, Nvars, Nz, iterm2, 2, iz, s_temp, Var3d=OutVar%vdata)

                do iv3 = iv2+1, Nvars
                  iterm3 = iterm3 + 1
                  v3 = MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=InVars(3)%vdata)
                  s_temp = MultiDimLookup(Nterms, Nvars, Nz, iterm3, 3, iz, Var3d=OutVar%vdata) + (v1 * v2 * v3)
                  call MultiDimAssign(Nterms, Nvars, Nz, iterm3, 3, iz, s_temp, Var3d=OutVar%vdata)

                  do iv4 = iv3+1, Nvars
                    iterm4 = iterm4 + 1
                    v4 = MultiDimLookup(Nx, Ny, Nz, ix, iy, iz, Var3d=InVars(4)%vdata)
                    s_temp = MultiDimLookup(Nterms, Nvars, Nz, iterm4, 4, iz, Var3d=OutVar%vdata) + (v1 * v2 * v3 * v4)
                    call MultiDimAssign(Nterms, Nvars, Nz, iterm4, 4, iz, s_temp, Var3d=OutVar%vdata)

                  enddo
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

    enddo

    ! Write out the moment data
    call rhdf5_write_variable(OutFileId, NumPoints%vname, NumPoints%ndims, it, NumPoints%dims, &
       NumPoints%units, NumPoints%descrip, NumPoints%dimnames, rdata=NumPoints%vdata)
    call rhdf5_write_variable(OutFileId, OutVar%vname, OutVar%ndims, it, OutVar%dims, &
       OutVar%units, OutVar%descrip, OutVar%dimnames, rdata=OutVar%vdata)

    ! Free up input var space
    do ivar = 1, Nvars
      deallocate(InVars(ivar)%vdata)
    enddo

    if (UseFilter) then
      deallocate(Filter%vdata)
    endif
  enddo

  ! Write out coordinates and attach dimensions

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
  call rhdf5_attach_dimensions(OutFile, NumPoints)
  call rhdf5_attach_dimensions(OutFile, OutVar)

  ! cleanup
  call rhdf5_close_file(OutFileId)
  do ivar = 1, Nvars
    call rhdf5_close_file(InFileIds(ivar))
  enddo

  stop

contains
!**********************************************************************
! Subroutines go here. Want these contained within the main program
! so that the interfaces to these subroutine are 'explicit' (ie, like
! ANSI C external declarations). This is needed to get the passing
! of allocatable arrays working properly.
!**********************************************************************

!**********************************************************************
! GetMyArgs()
!

 subroutine GetMyArgs(MaxStr, MaxVars, InDir, InSuffix, OutFile, FilterFile, UseFilter, VarList, NumVars, Nterms)
  implicit none

  integer, parameter :: MaxFields = 15

  integer :: MaxStr, MaxVars, NumVars, Nterms
  character (len=MaxStr) :: InDir, InSuffix, OutFile, FilterFile
  type (VarSpec), dimension(MaxVars) :: VarList
  logical :: UseFilter

  integer :: iargc, i, Nargs, Nfld
  character (len=MaxStr) :: Arg
  character (len=MaxStr), dimension(MaxFields) :: Fields

  logical :: BadArgs

  BadArgs = .false.
  InDir = ''
  InSuffix = ''
  OutFile = ''
  NumVars = 0

  ! parse args
  i = 1
  Nargs = iargc()

  do while (i .le. Nargs)
    call getarg(i, Arg)

    if (i .eq. 1) then 
      InDir = Arg
      i = i + 1
    else if (i .eq. 2) then 
      InSuffix = Arg
      i = i + 1
    else if (i .eq. 3) then 
      OutFile = Arg
      i = i + 1
    else if (i .eq. 4) then 
      FilterFile = Arg
      i = i + 1
    else
      call String2List(Arg, ':', Fields, MaxFields, Nfld, 'var spec')
      i = i + 1

      NumVars = NumVars + 1
      if (NumVars .gt. MaxVars) then
        write (*,'(a,i1,a)') 'ERROR: Exceeded maximum number of times <var_spec> can be used (', MaxVars, ')'
        BadArgs = .true.
      else
        VarList(NumVars)%Vname = Fields(1)
        VarList(NumVars)%Vfprefix = Fields(2)
      endif
    endif
  enddo

  if (FilterFile .eq. 'none') then
    UseFilter = .false.
  else
    UseFilter = .true.
  endif
  
  if (NumVars .lt. 1) then
    write (*,*) 'ERROR: must specify at least one <var_spec>'
    BadArgs = .true.
  endif
  
  if (BadArgs) then
    write (*,*) 'USAGE: gen_moments <in_dir> <in_suffix> <out_file> <filter_file> <ts_start> <ts_end> <var_spec> [<var_spec>...]'
    write (*,*) '        <in_dir>: directory containing input files'
    write (*,*) '        <in_suffix>: suffix attached to all input files'
    write (*,*) '        <out_file>: output file name'
    write (*,*) '        <filter_file>: filter file name'
    write (*,*) '        <var_spec> = <var_name>:<file_prefix>'
    write (*,*) '          <var_name>: name of variable inside input file'
    write (*,*) '          <file_prefix>: leading name of input HDF5 file'
    write (*,*) '              input file name gets built by joining: <file_prefix><in_suffix>'
    write (*,*) ''

    stop
  end if

  ! The number of variables determines the number of terms
  !     1 variable  --> 1 term
  !     2 variables --> 2 terms
  !     3 variables --> 3 terms
  !     4 variables --> 6 terms
  if (NumVars .eq. 4) then
    Nterms = 6
  else
   Nterms = NumVars
  endif

  return
end subroutine GetMyArgs

end program gen_moments
