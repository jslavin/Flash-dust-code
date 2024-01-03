!!****if* source/Particles/localAPI/pt_advanceDust_AB
!!
!! NAME
!!
!!  pt_advanceDust_AB
!!
!! SYNOPSIS
!!
!!  call pt_advanceDust_AB(real(in)    :: dtold,
!!                      real(in)    :: dtnew,
!!                      real(inout) :: particles(NPART_PROPS,p_count),
!!                      integer(in) :: p_count,
!!                      integer(in) :: ind)
!!
!! DESCRIPTION
!!
!!   Advances particles in time
!!
!! ARGUMENTS
!!
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particles -- particles on which to operate
!!   p_count - the number of particles in the list to advance
!!   ind -- index into pt_typeInfo 
!!
!! DESCRIPTION
!!
!!   Time advancement routine for dust particles. Currently just the Leapfrog
!!   method. When mass evolves as well, may need different technique.
!!
!!***

subroutine pt_advanceDust_AB(dtOld,dtNew,particles,p_count,ind)

  implicit none

#include "Flash.h"
#include "constants.h"  
#include "Particles.h"

  real, intent(in)   :: dtOld, dtNew
  integer,intent(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles

end subroutine pt_advanceDust_AB
