subroutine pt_sputrate(grvel, grrad, nH, ntype, dmgdt)
!!     Program currently calculates the non-thermal (or inertial) sputtering 
!!     rate (g/s) for the species (from the atomic set H, He, O, C and N).
!!
!!     Parameters:
!!     -----------
!!     grvel            gas-grain relative velocity (cm/s)
!!     grrad            grain radius (cm)
!!     nH               gas (H) density
!!     ntype            grain (target) type (1 => carbonaceous, 2 => silicate)
!!
!!     Returns:
!!     --------
!!     dmgdt            sputtering rate (g/s) 
!!
!!     Extra parameters: (contained in pt_yieldData)
!!     -----------------
!!     mproj            mass of projectiles (gas) particles
!!     sputparm         sputtering parameters (see pt_yieldData)
!!     abund            gas phase abundances of projectiles
!!
    use pt_yieldData

    implicit none

#include "Flash.h"
#include "constants.h"

    real, intent(in) :: grvel, grrad, nH
    integer, intent(in) :: ntype
    real, intent(out) :: dmgdt
    real :: ergeV,gi,MiMd,Eth,yield,Mprj,Egr
    integer :: i
    !parameter(ngrid=27)
    !real :: en(ngrid),Y_mean(ngrid),f_M(ngrid),Mprj,fint

    ergeV = 1.602176E-12
    ! Leaving out grain potential for now
    ! Ignore ionization of projectile (gas) particles for now
    ! As a result, don't include He++
    dmgdt = 0.
    do i = 1,nion-1
        gi = 4.*mproj(i)*sputparm(ntype,2)/((mproj(i) + sputparm(ntype,2))**2)
        MiMd = mproj(i)/sputparm(ntype,2)
        if(MiMd.le.0.3) then
            Eth = sputparm(ntype,1)/(gi*(1. - gi))
        else 
            Eth = 8.0*sputparm(ntype,1)*MiMd**(1./3.)
        endif
        Mprj = mproj(i)*amu
        Egr = Mprj*grvel**2/(2.*ergeV)
        if(Egr < Eth) cycle
        call pt_sputYield(ntype,Egr,i,yield)
        dmgdt = dmgdt + yield*abund(i)
    enddo
    dmgdt = dmgdt*PI*grrad*grrad*nH*grvel*sputparm(ntype,2)*amu
end subroutine pt_sputrate
