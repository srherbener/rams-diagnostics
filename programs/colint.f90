!***************************************************************
! Program to do azimuthial averaging
!
! This program will read in GRADS data from a RAMS simulation, and
! generate the column integrated value of a given variable.
!
! Args
!   1. input GRADS file names (control files, colon separated list)
!   2. output GRADS file name (basename, this program will tag
!      on the .ctl and .dat suffixes)
!   3. variable used for column integration
!
! Output
!   The output will be two files which make a GRADS data set. One
!   file is the control file and the other is the binary data.
!

program main
  use GfileTypes
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128
  integer, parameter :: MaxFiles=10

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarToInt

  character (len=MediumString), dimension(1:MaxFiles) :: GradsCtlFiles
  integer :: Nfiles

  type (GradsDataDescription), dimension(1:MaxFiles) :: GdataDescrip
  integer :: Nx, Ny, Nz, Nt, Nvars
  type (GradsOutDescription) :: GoutDescrip

  ! Data arrays: need one for dn0 (density) and the var we are doing the
  ! column integration on
  ! Dims: x, y, z, t
  ! The *Loc vars hold the locations of w, press, var in the GRADS
  ! data files: the first index is the file number, the second index is the
  ! var number
  real, dimension(:,:,:,:), allocatable :: Dens, TempC, ColVar, ColInt
  real, dimension(:), allocatable :: ZmHeights
  real, dimension(1:1) :: DummyZcoords
  real :: Xstart, Xinc, Ystart, Yinc
  type (GradsVarLocation) :: DensLoc, TempLoc, VarLoc

  integer :: i
  integer :: Ierror

  integer :: ix, iy, iz, it
  logical :: SuperCooled, WarmRain

  ! Get the command line argument
  call GetMyArgs(Infiles, OfileBase, VarToInt)
  if (VarToInt .eq. 'liquid_s') then
    SuperCooled = .true.
    VarToInt = 'liquid'
  else
    SuperCooled = .false.
  end if

  if (VarToInt .eq. 'liquid_w') then
    WarmRain = .true.
    VarToInt = 'liquid'
  else
    WarmRain = .false.
  end if

  call String2List(Infiles, ':', GradsCtlFiles, MaxFiles, Nfiles, 'input files')

  write (*,*) 'Calculating column integration for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, Nfiles
    write (*,*) '  ', i, ': ', trim(GradsCtlFiles(i))
  end do
  write (*,*) '  Output file base name:  ', trim(OfileBase)
  write (*,*) '  RAMS variable that is being integrated: ', trim(VarToInt)
  if (SuperCooled) then
    write (*,*) '    Applying tempc filter for supercooled liquid' 
  end if
  if (WarmRain) then
    write (*,*) '    Applying tempc filter for warm rain liquid' 
  end if
  write (*,*) ''


  ! Read the GRADS data description files and collect the information about the data
  do i = 1, Nfiles
    write (*,*) 'Reading GRADS Control File: ', trim(GradsCtlFiles(i))
    call ReadGradsCtlFile(GradsCtlFiles(i), GdataDescrip(i))
  end do
  write (*,*) ''

  call CheckDataDescrip_CI(GdataDescrip, Nfiles, Nx, Ny, Nz, Nt, Nvars, DensLoc, TempLoc, VarLoc, VarToInt)

  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', Nx
  write (*,*) '  Number of y (latitude) points:           ', Ny
  write (*,*) '  Number of z (vertical level) points:     ', Nz
  write (*,*) '  Number of t (time) points:               ', Nt
  write (*,*) '  Total number of grid variables:          ', Nvars
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', Nx*Ny*Nz*Nt
  write (*,*) ''

  write (*,*) 'Locations of variables in GRADS data (file number, var number):'
  write (*,'(a20,i3,a2,i3,a1)') 'dn0: (', DensLoc%Fnum, ', ', DensLoc%Vnum, ')'
  write (*,'(a20,i3,a2,i3,a1)') 'tempc: (', TempLoc%Fnum, ', ', TempLoc%Vnum, ')'
  write (*,'(a17,a3,i3,a2,i3,a1)') trim(VarToInt), ': (', VarLoc%Fnum, ', ', VarLoc%Vnum, ')'
  write (*,*) ''

  ! Allocate the data arrays and read in the data from the GRADS data files
  ! Make sure ZmHeights goes from 0 to n
  allocate (Dens(1:Nx,1:Ny,1:Nz,1:Nt), TempC(1:Nx,1:Ny,1:Nz,1:Nt), &
            ColVar(1:Nx,1:Ny,1:Nz,1:Nt), ZmHeights(0:Nz), &
            ColInt(1:Nx, 1:Ny, 1:1, 1:Nt), stat=Ierror)
  if (Ierror .ne. 0) then
    write (*,*) 'ERROR: Data array memory allocation failed'
    stop
  end if

  ! Read in the data for the vars using the description and location information
  call ReadGradsData(GdataDescrip, 'dn0', DensLoc, Dens, Nx, Ny, Nz, Nt)
  call ReadGradsData(GdataDescrip, 'tempc', TempLoc, TempC, Nx, Ny, Nz, Nt)
  call ReadGradsData(GdataDescrip, VarToInt, VarLoc, ColVar, Nx, Ny, Nz, Nt)

  ! Set the z heights which can be used to determine the layer thicknesses
  call SetZmHeights(Nz, ZmHeights)

  write (*,*) 'Zm Heigths:'
  do iz = 0, Nz
    write (*,*) '  Zm(', iz, '): ', ZmHeights(iz)
  enddo
  write (*,*) ''
 
  ! Generate a column integration at each time step for ColVar
  ! Mixing ratios in GRADS files are g/kg
  ! Density is in kg/m**3
  ! Heights are in m
  ! So, express the column integrated value in g/m**2 using the formula
  !   (mix ratio) * (density) * (layer thickness)
  !   sum up each layer value
  ! The layer thickness for layer k is: ZmHeights(k) - ZmHeights(k-1)
  do it = 1, Nt
    do ix = 1, Nx
      do iy = 1, Ny
        ColInt(ix,iy,1,it) = 0.0
        do iz = 1, Nz
          if (SuperCooled) then
            if (TempC(ix,iy,iz,it) .le. 0.0) then
              ColInt(ix,iy,1,it) = ColInt(ix,iy,1,it) + &
                (ColVar(ix,iy,iz,it) * Dens(ix,iy,iz,it) * (ZmHeights(iz) - ZmHeights(iz-1)))
            end if
          else
            if (WarmRain) then
              if (TempC(ix,iy,iz,it) .gt. 0.0) then
                ColInt(ix,iy,1,it) = ColInt(ix,iy,1,it) + &
                  (ColVar(ix,iy,iz,it) * Dens(ix,iy,iz,it) * (ZmHeights(iz) - ZmHeights(iz-1)))
              end if
            else
              ColInt(ix,iy,1,it) = ColInt(ix,iy,1,it) + &
                (ColVar(ix,iy,iz,it) * Dens(ix,iy,iz,it) * (ZmHeights(iz) - ZmHeights(iz-1)))
            end if
          end if
        enddo
      enddo
    enddo
  enddo

  DummyZcoords(1) = 0.0
  Xstart = GdataDescrip(1)%Xcoords(1)
  Xinc = GdataDescrip(1)%Xcoords(2) - Xstart
  Ystart = GdataDescrip(1)%Ycoords(1)
  Yinc = GdataDescrip(1)%Ycoords(2) - Ystart
  call BuildGoutDescrip(Nx, Ny, 1, Nt, ColInt, OfileBase, GdataDescrip(1)%UndefVal, &
          VarToInt, Xstart, Xinc, Ystart, Yinc, DummyZcoords, &
          GdataDescrip(1)%Tstart, GdataDescrip(1)%Tinc, GoutDescrip, 'colint')

  call WriteGrads(GoutDescrip, ColInt)

  ! Clean up
  deallocate (Dens, TempC, ColVar, ZmHeights, ColInt, stat=Ierror)
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
!   VarToInt - RAMS variable to do the column integration on
!

subroutine GetMyArgs(Infiles, OfileBase, VarToInt)
  implicit none

  character (len=*) :: Infiles, OfileBase, VarToInt

  integer :: iargc
  character (len=128) :: arg

  if (iargc() .ne. 3) then
    write (*,*) 'ERROR: must supply exactly 3 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: column <in_data_files> <out_data_file> <var_to_integrate>'
    write (*,*) '        <in_data_files>: GRADS format, control file, colon separated list'
    write (*,*) '        <out_data_file>: GRADS format, this programe will tag on .ctl, .dat suffixes'
    write (*,*) '        <var_to_integrate>: name of RAMS variable to do the integration on'
    write (*,*) ''
    stop
  end if

  call getarg(1, Infiles)
  call getarg(2, OfileBase)

  call getarg(3, VarToInt)

  return
end subroutine

