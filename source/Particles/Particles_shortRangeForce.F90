!!****f* source/Particles/Particles_shortRangeForce
!!
!! NAME
!!
!!  Particles_shortRangeForce
!!
!! SYNOPSIS
!!
!!  Particles_shortRangeForce()
!!
!! PARAMETERS
!!
!! ARGUMENTS
!!  
!! DESCRIPTION
!!
!!  Computes short-range forces on particles, ie. forces which couple particles
!!  only to some number of nearest neighbors.
!!  
!!***

!===============================================================================

subroutine Particles_shortRangeForce(particles,p_count,mapType)
 
!-------------------------------------------------------------------------------
#include "Flash.h"

  implicit none

!-------------------------------------------------------------------------------
  integer, intent(IN) :: p_count,mapType
  real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles

!-------------------------------------------------------------------------------

  return

end subroutine Particles_shortRangeForce

!===============================================================================
