!!****f* source/Grid/Grid_coordTransfm
!!
!! NAME
!!
!!  Grid_coordTransfm
!!
!! SYNOPSIS
!!
!!  call Grid_coordTransfm(real(IN)    :: x,
!!                         real(IN)    :: y,
!!                         real(IN)    :: z,
!!                         real(OUT)   :: xout,
!!                         real(OUT)   :: yout,
!!                         real(OUT)   :: zout,
!!                OPTIONAL,integer(IN) :: geometryIn,
!!                OPTIONAL,integer(IN) :: geometryOut,
!!                OPTIONAL,integer(IN) :: ndim,
!!                OPTIONAL,real(IN)    :: velI,
!!                OPTIONAL,real(IN)    :: velJ,
!!                OPTIONAL,real(IN)    :: velK,
!!                OPTIONAL,real(OUT)   :: velIOut,
!!                OPTIONAL,real(OUT)   :: velJOut,
!!                OPTIONAL,real(OUT)   :: velKOut)
!!
!!
!! DESCRIPTION
!!
!!  Convert between Cartesian and other coordinates.
!!
!!  Supports the FLASH geometries CARTESIAN, CYLINDRICAL, SPHERICAL, and POLAR.
!!
!! ARGUMENTS
!!
!!  x - First input coordinate.
!!  y - Second input coordinate.
!!  z - Third input coordinate.
!!  xout - First Cartesian coordinate.
!!  yout - Second Cartesian coordinate.
!!  zout - Third Cartesian coordinate.
!!  geometryIn - What is the input geometry? (default is the FLASH geometry established by
!!                                            configuration and runtime parameter)
!!  geometryOut - What is the output geometry? (default is the FLASH geometry established by
!!                                            configuration and runtime parameter)
!!  ndim - dimensionality (default N_DIM)
!!  velI,velJ,velK - input velocity components
!!  velIOut,velJOut,velKOut - output velocity components
!!
!!
!!***

subroutine Grid_coordTransfm(x,y,z, xout,yout,zout, geometryIn,geometryOut, ndimArg, &
     velI,velJ,velK,velIOut,velJOut,velKOut)
  implicit none
  real,intent(IN) :: x,y,z
  real,intent(OUT) :: xout,yout,zout
  integer,OPTIONAL,intent(IN) :: geometryIn
  integer,OPTIONAL,intent(IN) :: geometryOut
  integer,OPTIONAL,intent(IN) :: ndimArg
  real,OPTIONAL,intent(IN) :: velI,velJ,velK
  real,OPTIONAL,intent(OUT) :: velIOut,velJOut,velKOut

  if (present(velIOut).AND.present(velI)) velIOut = velI
  if (present(velJOut).AND.present(velJ)) velJOut = velJ
  if (present(velKOut).AND.present(velK)) velKOut = velK
  xout = 0.0
  yout = 0.0
  zout = 0.0

end subroutine Grid_coordTransfm
