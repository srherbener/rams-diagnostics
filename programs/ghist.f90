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
  use gdata_mod
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarName

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles
  integer :: Nfiles

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars, Nbins
  type (GradsOutDescription) :: GoutDescrip
  real :: BinStart, BinSize
  logical :: DoLogScale

  ! Data arrays: need one for dn0 (density) and the var we are doing the
  ! column integration on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: Var, Hist
  real :: Xstart, Xinc, Ystart, Yinc
  type (GradsVarLocation) :: VarLoc

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it, ib
  logical :: SuperCooled, WarmRain

  ! Get the command line argument
  call GetMyArgs(Infiles, OfileBase, VarName, DoLogScale, Nbins, BinStart, BinSize)
  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Building histograms for GRADS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  GRADS variable that is being analyzed: ', trim(VarName)
  write (*,*) '  DoLogScale: ', DoLogScale
  write (*,*) '  Number of bins: ', Nbins
  write (*,*) '  Bin start (lower end of first bin): ', BinStart
  write (*,*) '  Bin size: ', BinSize
  write (*,*) ''
  flush(6)


  ! Read the GRADS data description files and collect the information about the data
  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''
  flush(6)

  call CheckDataDescripOneVar(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, VarLoc, VarName)

  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', Nx
  write (*,*) '  Number of y (latitude) points:           ', Ny
  write (*,*) '  Number of z (vertical level) points:     ', Nz
  write (*,*) '  Number of t (time) points:               ', Nt
  write (*,*) '  Total number of grid variables:          ', Nvars
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''

  write (*,*) 'Location of variable in GRADS data (file number, var number):'
  write (*,'(a17,a3,i3,a2,i3,a1)') trim(VarName), ': (', VarLoc%Fnum, ', ', VarLoc%Vnum, ')'
  write (*,*) ''
  flush(6)

  ! Allocate the data arrays and read in the data from the GRADS data files
  allocate (Var(1:Nx,1:Ny,1:Nz,1:Nt), Hist(1:Nbins, 1:1, 1:Nz, 1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(GdataDescrip, VarName, VarLoc, Var, Nx, Ny, Nz, Nt)

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
  do it = 1, Nt
    do iz = 1, Nz

      ! zero out the bin counts
      do ib = 1, Nbins
        Hist(ib,1,iz,it) = 0.0
      enddo

      ! walk through the entire horizontal domain and count up the occurrances of each bin value
      do ix = 1, Nx
        do iy = 1, Ny
          ! figure out which bin to place the current var value
          ib = int((Var(ix,iy,iz,it) - BinStart)/BinSize) + 1
          if (ib .lt. 1) then
            ib = 1
          endif
          if (ib .gt. Nbins) then
            ib = Nbins
          endif

          Hist(ib,1,iz,it) = Hist(ib,1,iz,it) + 1.0
        enddo
      enddo

     ! convert the counts to logarithmic scale
     if (DoLogScale) then
       do ib = 1, Nbins
         if (Hist(ib,1,iz,it) .gt. 0.0) then
           Hist(ib,1,iz,it) = log10(Hist(ib,1,iz,it))
         endif
       enddo
     endif

    enddo
  enddo

  Xstart = BinStart
  Xinc = BinSize
  Ystart = 0.0
  Yinc = 1.0
  call BuildGoutDescrip(Nbins, 1, Nz, Nt, Hist, OfileBase, GdataDescrip(1)%UndefVal, &
          VarName, Xstart, Xinc, Ystart, Yinc, GdataDescrip(1)%Zcoords, &
          GdataDescrip(1)%Tstart, GdataDescrip(1)%Tinc, GoutDescrip, 'hist')

  call WriteGrads(GoutDescrip, Hist)

  ! Clean up
  deallocate (Var, Hist, stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory de-allocation failed'
    stop
  end if

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

subroutine GetMyArgs(Infiles, OfileBase, VarName, DoLogScale, Nbins, BinStart, BinSize)
  implicit none

  character (len=*) :: Infiles, OfileBase, VarName
  integer :: Nbins
  real :: BinStart, BinSize
  logical :: DoLogScale

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 7) then
    write (*,*) 'ERROR: must supply exactly 7 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: column <in_data_files> <out_data_file> <var_to_integrate> <scale> <number_of_bins> <bin_start> <bin_size>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <var_to_integrate>: name of RAMS variable to do the integration on'
    write (*,*) '        <scale>: scale for counts - "log" or "linear"'
    write (*,*) '        <number_of_bins>: number of bins to build into the histograms'
    write (*,*) '        <bin_start>: starting value of first bin'
    write (*,*) '        <bin_size>: size of bin (bins are uniform sizes)'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)

  call getarg(2, OfileBase)

  call getarg(3, VarName)

  call getarg(4, arg)
  if (arg .eq. 'log') then
    DoLogScale = .true.
  else
    DoLogScale = .false.
  endif

  call getarg(5, arg)
  read (arg, '(i)') Nbins

  call getarg(6, arg)
  read (arg, '(f)') BinStart

  call getarg(7, arg)
  read (arg, '(f)') BinSize

  return
end subroutine

