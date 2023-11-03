!!****f* source/Grid/Grid_advanceAnisoDiffusion
!!
!! NAME
!!
!!  Grid_advanceAnisoDiffusion
!!
!! SYNOPSIS
!!
!!  call Grid_advanceAnisoDiffusion (integer(IN)  :: iVar,
!!                              integer(IN)         :: iSrc,
!!                              integer(IN),dimension(:) :: iFactorsB(GRID_PDE_DIFFCOEF_NCOMP),
!!                              integer(IN)  :: iFactorA,
!!                              integer(IN)  :: bcTypes(6),
!!                              real(IN)     :: bcValues(2,6),
!!                              real(IN)     :: dt,
!!                              real(IN)     :: chi,
!!                              real(IN)     :: scaleFact,
!!                              real(IN)     :: theta,
!!                              logical(IN)  :: solnIsDelta,
!!                     OPTIONAL,integer(IN)  :: iFactorC,
!!                     OPTIONAL,integer(IN)  :: iFactorD)
!!
!!  DESCRIPTION 
!!
!!      This routine advances a generalized diffusion operator of the form
!!
!!         A*(df/dt) + C*f = div(B*grad(f)) + D ,
!!
!!      where
!!         f = f(x,t) is the  Variable to be diffused (x=1D..3D position);
!!
!!         A,B,C,D are optional given scalar factors/terms that may depend
!!         on position; they are either physcially constant in time, or at
!!         least considered time-independent for the purpose of the operation
!!         implemented here (typically by computing their values from the
!!         solution state reached by the previous time step).
!!
!!      Presently it is used to do heat conduction and multigroup diffusion.
!!
!! ARGUMENTS
!!   iVar           : Variable on which the diffusion operation is performed (e.g., TEMP_VAR)
!!   iFactorA       :| Are factors in the equation with spatial variation.
!!   iFactorsB      :| Factor C,D are optional and are generally used
!!   iFactorC       :| to represent emission/absorption in MGD.
!!   iFactorD       :| iFactorA is not needed (and thus set to constant 1.0) for
!!                   | radiation diffusion.
!!                   | The factor B (whose role is that of diffusion coefficient) is passed to
!!                   | the solver implementation in cell-centered form in UNK; the solver will
!!                   | use averaging of adjacent cells (or perhaps some other mechanism) to
!!                   | derive the required face-cented diffusion coefficent values.
!!   bcTypes        : Presently OUTFLOW, VACUUM, DIRICHLET are supported, with additional
!!                    limited support for OUTSTREAM, in the HYPRE implementation.
!!   bcValues       : Values of iVar,iFactorB on boundary (DIRICHLET).                        
!!   dt             : The time step.
!!   scaleFact      : Factor by which the end solution is scaled (not used).
!!   chi            : useful for constant diffusion problems (not used).
!!   theta          : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   iSrc           : Ignored.
!!   solnIsDelta    : Is the solution only a delta that the caller has to apply to the
!!                    temperature, rather than temperature itself (ignored).
!!
!! SIDE EFFECTS
!!
!!  The iVar component in the solution array UNK is updated on successful completion.
!!
!! NOTES
!!  This implementation
!!    * supports: 
!!           2D, 3D Cartesian in UG (3D untested).
!!           2D Cylindrical in UG (R-Z)
!!
!!    *  uses HYPRE library to solve AX = B
!!
!! SEE ALSO
!! 
!!  Grid_advanceDiffusion
!!  Diffuse_solveScalar
!!  diff_advanceTherm
!!
!!
!!***


subroutine Grid_advanceAnisoDiffusion (iVar, iSrc, iFactorsB, iFactorA, bcTypes, bcValues, &
     dt, chi, scaleFact,theta, solnIsDelta,iFactorC, iFactorD)
  
  
  implicit none
  
  integer, intent(IN) :: iVar
  integer, intent(IN) :: iSrc
  integer, intent(IN) :: iFactorsB(:)
  integer, intent(IN) :: iFactorA
  real, intent(IN)    :: dt 
  real, intent(IN)    :: chi
  real, intent(IN)    :: scaleFact
  real, intent(IN)    :: theta
  logical, intent(IN) :: solnIsDelta
  integer, dimension(6),  intent(IN) :: bcTypes
  real   , dimension(2,6),intent(IN) :: bcValues
  integer, intent(IN), OPTIONAL :: iFactorC
  integer, intent(IN), OPTIONAL :: iFactorD   
  
end subroutine Grid_advanceAnisoDiffusion
