!***************************************************************
! Program to create histograms
!
! This program will read in GRADS data from a RAMS simulation, and
! generate histograms of a given variable.
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
!   3. variable used for histgrams
!
! Output
!   The output will be two files which make a GRADS data set. One
!   file is the control file and the other is the binary data.
!

program main
  use gdata_utils
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarName

  type (GradsControlFiles) :: GctlFiles

  integer :: Nx, Ny, Nz, Nt, Nbins
  real :: BinStart, BinSize, DataMin, DataMax
  logical :: DoLogScale, DoLinScale, DoFaScale

  ! Data arrays: need one for dn0 (density) and the var we are doing the
  ! column integration on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  type (GradsVar) :: Gvar, Hist
  real :: Xstart, Xinc, Ystart, Yinc

  integer :: i, Nhoriz

  integer :: ix, iy, iz, it, ib
  logical :: SuperCooled, WarmRain

  ! Get the command line argument
  call GetMyArgs(Infiles, OfileBase, VarName, DoLogScale, DoLinScale, DoFaScale, DataMin, DataMax, &
                     Nbins, BinStart, BinSize)
  call String2List(Infiles, ':', GctlFiles%Fnames, MaxFiles, GctlFiles%Nfiles, 'input files')

  write (*,*) 'Building histograms for GRADS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, GctlFiles%Nfiles
    write (*,*) '  ', i, ': ', trim(GctlFiles%Fnames(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  GRADS variable that is being analyzed: ', trim(VarName)
  write (*,*) '  DoLogScale: ', DoLogScale
  write (*,*) '  DoLinScale: ', DoLinScale
  write (*,*) '  DoFaScale: ', DoFaScale
  write (*,*) '  DataMin: ', DataMin
  write (*,*) '  DataMax: ', DataMax
  write (*,*) '  Number of bins: ', Nbins
  write (*,*) '  Bin start (lower end of first bin): ', BinStart
  write (*,*) '  Bin size: ', BinSize
  write (*,*) ''
  flush(6)


  ! Read the GRADS data description files and collect the information about the data
  call ReadGradsCtlFiles(GctlFiles)

  call InitGvarFromGdescrip(GctlFiles, Gvar, VarName)

  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', Gvar%Nx
  write (*,*) '  Number of y (latitude) points:           ', Gvar%Ny
  write (*,*) '  Number of z (vertical level) points:     ', Gvar%Nz
  write (*,*) '  Number of t (time) points:               ', Gvar%Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Gvar%Nx*Gvar%Ny*Gvar%Nz*Gvar%Nt
  write (*,*) ''

  write (*,*) 'Locations of variables in GRADS data (file, var number):'
  write (*,'(a17,a3,a,a2,i3,a1)') trim(VarName), ': (', trim(Gvar%DataFile), ', ', Gvar%Vnum, ')'
  write (*,*) ''
  flush(6)

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(Gvar)

  ! Initialize the output var
  call InitGradsVar(Hist, VarName, Nbins, 1, Gvar%Nz, Gvar%Nt, &
                    BinStart, BinSize, 0.0, 1.0, Gvar%Zcoords, Gvar%Tstart, Gvar%Tinc, &
                    Gvar%UndefVal, '<NONE>', 0, 0)

  ! Generate a histogram based on the entire horizontal domain for each t and z combination
  !     Bin       Range
  !      1      (BinStart + (0 * BinSize)) to just less than (BinStart + (1 * BinSize))
  !      2      (BinStart + (1 * BinSize)) to just less than (BinStart + (2 * BinSize))
  !      ...
  !      n      (BinStart + ((n-1) * BinSize)) to just less than (BinStart + (n * BinSize))
  !
  ! So, the bin number for a value x is given by: ib = int((x - BinStart)/BinSize) + 1
  ! This formula doesn't check to see if x will fit in one of the defined bins meaning that
  ! an index (ib value) can fall outside the defined bins. Check for this and just snap an
  ! ib value < 1 to 1 and an ib value > Nbins to Nbins.
  Nhoriz = Gvar%Nx * Gvar%Ny
  do it = 1, Gvar%Nt
    do iz = 1, Gvar%Nz

      ! zero out the bin counts
      do ib = 1, Nbins
        Hist%Vdata(it,iz,ib,1) = 0.0
      enddo

      ! walk through the entire horizontal domain and count up the occurrances of each bin value
      do ix = 1, Gvar%Nx
        do iy = 1, Gvar%Ny
          ! figure out which bin to place the current var value
          ib = int((Gvar%Vdata(it,iz,ix,iy) - BinStart)/BinSize) + 1
          if (ib .lt. 1) then
            ib = 1
          endif
          if (ib .gt. Nbins) then
            ib = Nbins
          endif

          if ((Gvar%Vdata(it,iz,ix,iy) .ge. DataMin) .and. (Gvar%Vdata(it,iz,ix,iy) .le. DataMax)) then
            Hist%Vdata(it,iz,ib,1) = Hist%Vdata(it,iz,ib,1) + 1.0
          endif
        enddo
      enddo

     ! apply scaling, if linear scaling just leave counts as is
     if (DoLogScale) then
       ! convert the counts to logarithmic scale
       do ib = 1, Nbins
         if (Hist%Vdata(it,iz,ib,1) .gt. 0.0) then
           Hist%Vdata(it,iz,ib,1) = log10(Hist%Vdata(it,iz,ib,1))
         endif
       enddo
     else if (DoFaScale) then
       ! scale counts by fractional area of entire domain
       ! divide counts by total number of horiz grid points (assume uniform sized grid cells)
       do ib = 1, Nbins
         Hist%Vdata(it,iz,ib,1) = Hist%Vdata(it,iz,ib,1) / float(Nhoriz)
       enddo
     endif

    enddo
  enddo

  call WriteGrads(Hist, OfileBase, 'hist')

  stop
end

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the following command line arguments
!   Infiles - input GRADS file
!   OfileBase - output GRADS file, base name for two files
!   VarName - RAMS variable to do the column integration on
!

subroutine GetMyArgs(Infiles, OfileBase, VarName, DoLogScale, DoLinScale, DoFaScale, DataMin, DataMax, &
                     Nbins, BinStart, BinSize)
  implicit none

  character (len=*) :: Infiles, OfileBase, VarName
  integer :: Nbins
  real :: BinStart, BinSize, DataMin, DataMax
  logical :: DoLogScale, DoLinScale, DoFaScale

  integer :: iargc
  character (len=128) :: arg

  logical :: BadArgs

  if (iargc() .ne. 9) then
    write (*,*) 'ERROR: must supply exactly 9 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: column <in_data_files> <out_data_file> <var_to_integrate> <scale> <data_min> <data_max> <number_of_bins> <bin_start> <bin_size>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <var_to_integrate>: name of RAMS variable to do the integration on'
    write (*,*) '        <scale>: scale for counts - "log", "linear", or "frac_area"'
    write (*,*) '        <data_min> <data_max>: toss out data outside this range'
    write (*,*) '        <number_of_bins>: number of bins to build into the histograms'
    write (*,*) '        <bin_start>: starting value of first bin'
    write (*,*) '        <bin_size>: size of bin (bins are uniform sizes)'
    write (*,*) ''
    stop
  end if

  BadArgs    = .false.
  DoLogScale = .false.
  DoLinScale = .false.
  DoFaScale  = .false.

  call getarg(1, Infiles)

  call getarg(2, OfileBase)

  call getarg(3, VarName)

  call getarg(4, arg)
  if (arg .eq. 'log') then
    DoLogScale = .true.
  else if (arg .eq. 'linear') then
    DoLinScale = .true.
  else if (arg .eq. 'frac_area') then
    DoFaScale = .true.
  else
    write (*,*) 'ERROR: must choose one of "log", "linear", or "frac_area" for <scale> argument'
    BadArgs = .true.
  endif

  call getarg(5, arg)
  read (arg, '(f)') DataMin

  call getarg(6, arg)
  read (arg, '(f)') DataMax
  if (DataMax <= DataMin) then
    write (*,*) 'ERROR: <data_max> must be greater than <data_min>'
    BadArgs = 'true'
  endif

  call getarg(7, arg)
  read (arg, '(i)') Nbins

  call getarg(8, arg)
  read (arg, '(f)') BinStart

  call getarg(9, arg)
  read (arg, '(f)') BinSize

  if (BadArgs) then
    stop
  end if

  return
end subroutine

