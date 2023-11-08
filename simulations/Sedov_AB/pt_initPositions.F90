!!
!! NAME
!!
!!   pt_initPositions
!!
!! SYNOPSIS
!!
!!   pt_initPositions(integer, INTENT(in) :: blockId, 
!!                    logical, INTENT(out) :: success)
!!
!! DESCRIPTION
!!   Put grains in ejecta clumps
!!
!!
!! ARGUMENTS
!!
!!   blockId : Id of block in current processor
!!   success: logical argument indicating whether particles initalised 
!!            correctly on this block.
!!
!!
!!***

subroutine pt_initPositions(blockId, success)

use Simulation_data, ONLY : sim_pmass, sim_pdens, xpart, ypart, zpart, &
    vxpart, vypart, vzpart, npart
use Particles_data, ONLY : pt_numLocal, pt_meshMe
use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords

#include "Flash.h"
#include "constants.h"

implicit none

integer, intent(in) :: blockId
logical, intent(out) :: success

integer       :: i, j, k
real          :: xr, yr, zr, vxr, vyr, vzr
integer, save :: p, nleft

logical       :: IsInBlock
real          :: bxLower, bxUpper, byLower, byUpper, bzLower, bzUpper
real, dimension(3) :: blockSize, blockCenter
real :: time

!call Driver_getSimTime(time)
!if(time < sim_tpstart) then
!    success = .true.
!    return
!endif

! Particle slot number (incremented and saved between calls)
p = pt_numLocal

! Location of block faces
call Grid_getBlkPhysicalSize(blockId, blockSize)
call Grid_getBlkCenterCoords(blockId, blockCenter)
bxLower = blockCenter(1) - 0.5*blockSize(1)
bxUpper = blockCenter(1) + 0.5*blockSize(1)
if (NDIM > 1) then
   byLower = blockCenter(2) - 0.5*blockSize(2)
   byUpper = blockCenter(2) + 0.5*blockSize(2)
endif
if (NDIM > 2) then 
   bzLower = blockCenter(3) - 0.5*blockSize(3)
   bzUpper = blockCenter(3) + 0.5*blockSize(3)
endif

!if (pt_meshMe == MASTER_PE) &
!  print *, 'Particles_initPositions:  computing particle positions...'

do k = 1, npart
    xr = xpart(k)
    vxr = vxpart(k)
    yr = 0.
    vyr = 0.
    zr = 0.
    vzr = 0.
    if (NDIM > 1) then
        yr = ypart(k) 
        vyr = vypart(k)
    endif
    if (NDIM > 2) then
        zr = zpart(k)
        vzr = vzpart(k)
    endif
! Check if particle is in this block 
    IsInBlock = (xr >= bxLower) .and. (xr < bxUpper)
    if (NDIM > 1) &
        IsInBlock = IsInBlock.and.((yr >= byLower).and.(yr < byUpper))
    if (NDIM > 2) & 
        IsInBlock = IsInBlock.and.((zr >= bzLower).and.(zr < bzUpper))
! If it is, keep it; otherwise discard it.
    if (IsInBlock) then
        p = p + 1
        call InitSingleParticle(p, xr, yr, zr, vxr, vyr, vzr, sim_pmass, &
            sim_pdens, blockId)
    endif
enddo

pt_numLocal = p
success = .true.

return
end subroutine pt_initPositions

!**********************************************************************
!  Routine:     InitSingleParticle
!  Description: Initialize a single particle with given characteristics.

subroutine InitSingleParticle(p, xpos, ypos, zpos, xvel, yvel, zvel, &
        mass, dens, block)

  use Simulation_data, ONLY : sim_pISM, sim_rhocl
  use Eos_data, ONLY : eos_singleSpeciesA
  use Particles_data, ONLY : pt_maxPerProc, particles,pt_meshMe
  use Driver_interface, ONLY : Driver_abortFlash

#include "Flash.h"
#include "constants.h"

  implicit none

  real, intent(in)    :: xpos, ypos, zpos, xvel, yvel, zvel, mass, dens
  integer, intent(in) :: p, block
  real :: Tclump
  real, parameter :: Rconst = 8.3144598E7

  if (p > pt_maxPerProc) &
    call Driver_abortFlash("InitSingleParticle:  Exceeded max # of particles!")

  particles(BLK_PART_PROP,p)  = real(block)
  particles(PROC_PART_PROP,p) = real(pt_meshMe)
  particles(MASS_PART_PROP,p) = mass
  particles(DENS_PART_PROP,p) = dens
  particles(POSX_PART_PROP,p) = xpos
  particles(VELX_PART_PROP,p) = xvel
  particles(OVLX_PART_PROP,p) = xvel
  particles(ACCX_PART_PROP,p) = 0.
  particles(OACX_PART_PROP,p) = 0.
  particles(POSY_PART_PROP,p) = ypos
  particles(VELY_PART_PROP,p) = yvel
  particles(OVLY_PART_PROP,p) = yvel
  particles(ACCY_PART_PROP,p) = 0.
  particles(OACY_PART_PROP,p) = 0.
  particles(POSZ_PART_PROP,p) = zpos
  particles(VELZ_PART_PROP,p) = zvel
  particles(OVLZ_PART_PROP,p) = zvel
  particles(ACCZ_PART_PROP,p) = 0.
  particles(OACZ_PART_PROP,p) = 0.
  particles(GDEN_PART_PROP,p) = sim_rhocl
  ! This doesn't necessarily hold if we're increasing the pressure 
  ! in the ejecta - see Simulation_initBlock
  Tclump = sim_pISM/sim_rhocl*eos_singleSpeciesA/Rconst
  particles(GTMP_PART_PROP,p) = Tclump
  particles(VREL_PART_PROP,p) = 0.
  ! New values - now storing the gas velocity as well as the relative velocity
  ! Dummy value here
  particles(GVLX_PART_PROP,p) = xvel
  particles(GVLY_PART_PROP,p) = yvel
  particles(GVLZ_PART_PROP,p) = zvel
  particles(DMDT_PART_PROP,p) = 0.
  particles(ODMT_PART_PROP,p) = 0.
  ! Dummy values - should get updated at the next step
  particles(MAGX_PART_PROP,p) = 1.E-7
  particles(MAGY_PART_PROP,p) = 0.
  particles(MAGZ_PART_PROP,p) = 0.
  particles(CHRG_PART_PROP,p) = 0.
return
end subroutine InitSingleParticle
