!***************************************************************
! Data types for GRADS data routines
!

module GfileTypes
  integer, parameter :: MaxString=256
  integer, parameter :: MaxCoords=1000
  integer, parameter :: MaxVars=50
  integer, parameter :: InUnit=8
  integer, parameter :: OutUnit=10
  integer, parameter :: BinRecFactor=4
  
  type GradsDataDescription
    character (len=MaxString) :: DataFile
    real :: UndefVal
    real, dimension(1:MaxCoords) :: Xcoords
    real, dimension(1:MaxCoords) :: Ycoords
    real, dimension(1:MaxCoords) :: Zcoords
    real, dimension(1:MaxCoords) :: Tcoords
    character (len=MaxString), dimension(1:MaxVars) :: VarNames
    character (len=MaxString) :: Tstart, Tinc
    integer :: nx
    integer :: ny
    integer :: nz
    integer :: nt
    integer :: nvars
  end type GradsDataDescription

  type GradsVarLocation
    integer :: Fnum
    integer :: Vnum
  end type GradsVarLocation

  type GradsOutDescription
    character (len=MaxString) :: CtlFile
    character (len=MaxString) :: DataFile
    character (len=MaxString) :: Title
    character (len=MaxString) :: VarName
    real :: UndefVal
    real :: Xstart, Xinc
    real :: Ystart, Yinc
    real, dimension(1:MaxCoords) :: Zcoords
    character (len=MaxString) :: Tstart, Tinc
    integer :: nx
    integer :: ny
    integer :: nz
    integer :: nt
  end type GradsOutDescription
end module GfileTypes
