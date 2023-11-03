!!****if* source/physics/sourceTerms/Cool/CoolMain/MeKaL/Cool_computeDt
!!
!! NAME
!!  
!!  Cool_computeDt
!!
!!
!! SYNOPSIS
!! 
!!  Cool_computeDt ( integer(IN) : blockID, 
!!                   
!!                  real,pointer :  solnData(:,:,:,:),   
!!                  real,(INOUT):   dt_check, 
!!                  integer(INOUT): dt_minloc(:) )
!!  
!! DESCRIPTION
!!
!!  Computes the timestep limiter for cooling source term solver.
!!  Version for APEC tabulated cooling function
!!
!!
!! ARGUMENTS
!!
!!  blockID        local block ID
!!  
!!  solnData        the physical, solution data from grid
!!  dt_check        variable to hold timestep constraint
!!  dt_minloc(5)    array to hold limiting zone info:  zone indices
!!
!!***

subroutine Cool_computeDt (blockID, &
                              blkLimits,blkLimitsGC,        &
                              solnData,   &
                              dt_check, dt_minloc )
  
#include "Flash.h"
#include "constants.h"

  implicit none

  integer, intent(IN) :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real,INTENT(INOUT)    :: dt_check
  integer,INTENT(INOUT)    :: dt_minloc(5)
  real, pointer, dimension(:,:,:,:) :: solnData
  
  !==============================================================================

  ! We don't have to do anything here

  return

end subroutine Cool_computeDt

