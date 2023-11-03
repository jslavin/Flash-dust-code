subroutine pt_sputrate(grvel, grrad, nH, Tgas, ntype, dmgdt)
!!     Program currently calculates the combined thermal and inertial 
!!     sputtering (aka thermo-kinetic sputtering rate)
!!     rate (g/s) for the species (from the atomic set H, He, O, C and N).
!!
!!     Parameters:
!!     -----------
!!     grvel            gas-grain relative velocity (cm/s)
!!     grrad            grain radius (cm)
!!     nH               gas (H) density (cm^-3)
!!     Tgas             gas temperature (K)
!!     ntype            grain (target) type (1 => carbonaceous, 2 => silicate)
!!
!!     Returns:
!!     --------
!!     dmgdt            sputtering rate (g/s) 
!!
!!     Extra parameters: (contained in pt_yieldIntgrlData)
!!     -----------------
!!     mproj            mass of projectiles (gas) particles
!!     sputparm         sputtering parameters (see pt_yieldData)
!!     abund            gas phase abundances of projectiles
!!
    use pt_yieldIntgrlData, only : abund, amu, sputparm

    implicit none

#include "Flash.h"
#include "constants.h"

    real, intent(in) :: grvel, grrad, nH, Tgas
    integer, intent(in) :: ntype
    real, intent(out) :: dmgdt
    real :: yieldvel
    integer :: i
    !parameter(ngrid=27)
    !real :: en(ngrid),Y_mean(ngrid),f_M(ngrid),Mprj,fint

    ! Leaving out grain potential for now
    ! Ignore ionization of projectile (gas) particles for now
    ! As a result, don't include He++
    dmgdt = 0.
    ! Note: we're using 4 ions: H, He, O, C
    do i = 1,4
        ! pt_sputYieldVel is the yield averaged over the velocity distribution
        ! for given gas-grain relative velocity and temperature
        call pt_sputYieldVel(ntype,grvel,Tgas,i,yieldvel)
        dmgdt = dmgdt + yieldvel*abund(i)
    enddo
    ! sputparm(ntype,2) is the mass of a sputtered fragment (in amu)
    dmgdt = dmgdt*PI*grrad*grrad*nH*sputparm(ntype,2)*amu
end subroutine pt_sputrate
