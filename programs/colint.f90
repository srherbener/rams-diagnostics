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
  use gdata_utils
  implicit none

  integer, parameter :: LargeString=512
  integer, parameter :: MediumString=256
  integer, parameter :: LittleString=128

  character (len=LargeString) :: Infiles
  character (len=MediumString) :: OfileBase
  character (len=LittleString) :: VarToInt

  type (GradsControlFiles) :: GctlFiles

  ! Data arrays: need one for dn0 (density) and the var we are doing the
  ! column integration on
  ! Dims: t, z, x, y
  type (GradsVar) :: Dens, TempC, ColVar, ColInt
  real, dimension(:), allocatable :: ZmHeights
  real, dimension(1:1) :: DummyZcoords

  integer :: i

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

  call String2List(Infiles, ':', GctlFiles%Fnames, MaxFiles, GctlFiles%Nfiles, 'input files')

  write (*,*) 'Calculating column integration for RAMS data:'
  write (*,*) '  GRADS input control files:'
  do i = 1, GctlFiles%Nfiles
    write (*,*) '  ', i, ': ', trim(GctlFiles%Fnames(i))
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
  call ReadGradsCtlFiles(GctlFiles)

  call InitGvarFromGdescrip(GctlFiles, Dens, 'dn0')
  call InitGvarFromGdescrip(GctlFiles, TempC, 'tempc')
  call InitGvarFromGdescrip(GctlFiles, ColVar, VarToInt)

  ! make sure vars have the same dimensions, third argument to GvarDimsMatch is true if variable
  ! given by second argument is 2D
  if (.not. (GvarDimsMatch(Dens, TempC, .false.) .and. GvarDimsMatch(Dens, ColVar, .false.))) then
    write (*,*) 'ERROR: dimensions of dn0, tempc, and ', trim(VarToInt), ' do not match'
    stop
  endif

  write (*,*) 'Gridded data information:'
  write (*,*) '  Number of x (longitude) points:          ', ColVar%Nx
  write (*,*) '  Number of y (latitude) points:           ', ColVar%Ny
  write (*,*) '  Number of z (vertical level) points:     ', ColVar%Nz
  write (*,*) '  Number of t (time) points:               ', ColVar%Nt
  write (*,*) ''
  write (*,*) '  Number of data values per grid variable: ', ColVar%Nx*ColVar%Ny*ColVar%Nz*ColVar%Nt
  write (*,*) ''

  write (*,*) 'Locations of variables in GRADS data (file, var number):'
  write (*,'(a20,a,a2,i3,a1)') 'dn0: (', trim(Dens%DataFile), ', ', Dens%Vnum, ')'
  write (*,'(a20,a,a2,i3,a1)') 'tempc: (', trim(TempC%DataFile), ', ', TempC%Vnum, ')'
  write (*,'(a17,a3,a,a2,i3,a1)') trim(VarToInt), ': (', trim(ColVar%DataFile), ', ', ColVar%Vnum, ')'
  write (*,*) ''


  ! Read in the data for the vars using the description and location information
  call ReadGradsData(Dens)
  call ReadGradsData(TempC)
  call ReadGradsData(ColVar)

  ! Set the z heights which can be used to determine the layer thicknesses
  allocate(ZmHeights(0:ColVar%Nz))
  call SetZmHeights(ColVar%Nz, ZmHeights)

  write (*,*) 'Zm Heigths:'
  do iz = 0, ColVar%Nz
    write (*,*) '  Zm(', iz, '): ', ZmHeights(iz)
  enddo
  write (*,*) ''
 
  ! Initialize the output GRADS var: ColInt
  DummyZcoords(1) = 0.0
  call InitGradsVar(ColInt, VarToInt, ColVar%Nx, ColVar%Ny, 1, ColVar%Nt, &
                    ColVar%Xstart, ColVar%Xinc, ColVar%Ystart, ColVar%Yinc, DummyZcoords, &
                    ColVar%Tstart, ColVar%Tinc, ColVar%UndefVal, '<NONE>', 0, 0)

  ! Generate a column integration at each time step for ColVar
  ! Mixing ratios in GRADS files are g/kg
  ! Density is in kg/m**3
  ! Heights are in m
  ! So, express the column integrated value in g/m**2 using the formula
  !   (mix ratio) * (density) * (layer thickness)
  !   sum up each layer value
  ! The layer thickness for layer k is: ZmHeights(k) - ZmHeights(k-1)
  do it = 1, ColVar%Nt
    do ix = 1, ColVar%Nx
      do iy = 1, ColVar%Ny
        ColInt%Vdata(it,1,ix,iy) = 0.0
        do iz = 1, ColVar%Nz
          if (SuperCooled) then
            if (TempC%Vdata(it,iz,ix,iy) .le. 0.0) then
              ColInt%Vdata(it,1,ix,iy) = ColInt%Vdata(it,1,ix,iy) + &
                (ColVar%Vdata(it,iz,ix,iy) * Dens%Vdata(it,iz,ix,iy) * (ZmHeights(iz) - ZmHeights(iz-1)))
            end if
          else
            if (WarmRain) then
              if (TempC%Vdata(it,iz,ix,iy) .gt. 0.0) then
                ColInt%Vdata(it,1,ix,iy) = ColInt%Vdata(it,1,ix,iy) + &
                  (ColVar%Vdata(it,iz,ix,iy) * Dens%Vdata(it,iz,ix,iy) * (ZmHeights(iz) - ZmHeights(iz-1)))
              end if
            else
              ColInt%Vdata(it,1,ix,iy) = ColInt%Vdata(it,1,ix,iy) + &
                (ColVar%Vdata(it,iz,ix,iy) * Dens%Vdata(it,iz,ix,iy) * (ZmHeights(iz) - ZmHeights(iz-1)))
            end if
          end if
        enddo
      enddo
    enddo
  enddo

  call WriteGrads(ColInt, OfileBase, 'colint')

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

