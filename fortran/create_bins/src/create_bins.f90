!***************************************************************
! Program to create bins edge values for generating histograms.
!

program create_bins
  implicit none

  integer, parameter :: BTYPE_LINEAR = 1
  integer, parameter :: BTYPE_LOG10  = 2
  integer, parameter :: BTYPE_LN     = 3

  integer :: Nbins, Btype
  real :: Bstart, Binc

  integer :: i
  logical :: DoLinear
  real :: Bval, Bseries

  ! Get the command line arguments
  call GetMyArgs(Nbins, Bstart, Binc, Btype)

  print('(i)'), Nbins
  do i = 1, Nbins
    ! Bseries is either the value for the linear sacle
    ! or the power for the log scales
    if (i .eq. 1) then
      Bseries = Bstart
    else
      Bseries = Bseries + Binc
    endif

    if (Btype .eq. BTYPE_LINEAR) then
      Bval = Bseries
    else if (Btype .eq. BTYPE_LOG10) then
      Bval = 10.0 ** Bseries
    else if (Btype .eq. BTYPE_LN) then
      Bval = exp(Bseries)
    endif

    print('(e)'), Bval
  enddo

  stop

contains
!**********************************************************************
! SUBROUTINES
!**********************************************************************

!**********************************************************************
! GetMyArgs()
!
! This routine will read in the command line arguments
!
subroutine GetMyArgs(Nbins, Bstart, Binc, Btype)
  implicit none

  integer :: Nbins, Btype
  real :: Bstart, Binc

  integer :: iargc
  character (len=128) :: arg

  logical :: BadArgs

  if (iargc() .ne. 4) then
    write (*,*) 'ERROR: must supply exactly 4 arguments'
    write (*,*) ''
    write (*,*) 'USAGE: create_bins <num_bins> <bin_start> <bin_inc> <bin_type>'
    write (*,*) '           <num_bins>: number of bins'
    write (*,*) '           <bin_start>: starting value for bins'
    write (*,*) '           <bin_inc>: delta between bins'
    write (*,*) '           <bin_type>:'
    write (*,*) '               "linear": create even spaced linear bin edges'
    write (*,*) '               "log": create even spaced base 10 logrithmic bin edges'
    write (*,*) '               "ln": create even spaced base e logrithmic bin edges'
    stop
  end if

  BadArgs = .false.

  call getarg(1, arg)
  read(arg, '(i)') Nbins

  call getarg(2, arg)
  read(arg, '(f)') Bstart

  call getarg(3, arg)
  read(arg, '(f)') Binc

  call getarg(4, arg)
  if (arg .eq. 'linear') then
    Btype = BTYPE_LINEAR
  else if (arg .eq. 'log') then
    Btype = BTYPE_LOG10
  else if (arg .eq. 'ln') then
    Btype = BTYPE_LN
  else
    write (*,*) 'ERROR: <bin_type> must be one of: "linear", "log", or "ln"'
    BadArgs = .true.
  endif

  if (BadArgs) then
    stop
  end if

  return
end subroutine GetMyArgs

end program create_bins
