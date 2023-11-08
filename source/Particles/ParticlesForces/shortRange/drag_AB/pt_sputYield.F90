subroutine pt_sputYield(ntype,energy,ion,yield)


!!     Subroutine to calculate the sputtering yield for a given ion ( velocity, 
!!     atomic mass and atomic number specified as input parameters ) as a 
!!     function of target material ( binding energy, atomic mass, atomic
!!     number and sputtering parameter, K, given as input parameters ).
!!     NOTE: routine returns angle averaged yield
!!                 
!!     sptyld  = sputtering yield  Y(E)   ( 'atoms' per incident ion )  output
!!
!!     ntype        = index of grain type
!!     energy       = incident species energy  ( E )             ( eV )  input
!!     ion          = index of projectile ion species 
!!
!!     Data used from module pt_yieldData:
!!     mproj(ion)    = atomic mass of species 'ion' (amu)
!!     Zproj(ion)    = atomic number of species 'ion' 
!!     sputparm(ntype,1) surface binding energy (U) of target (eV)
!!     sputparm(ntype,2) atomic mass of sputtered fragment (amu), mtarg
!!     sputparm(ntype,3) atomic number of sputtered fragment, Ztarg
!!     sputparm(ntype,4) free parameter "K" constant

    use pt_yieldData

    implicit none

#include "constants.h"
#include "Flash.h"

    integer, intent(in) :: ntype,ion
    real, intent(in) :: energy
    real, intent(out) :: yield
    real :: U0,Mi,Md,Zi,Zd,K,mu,Eth,gi,rpovr,asc,zm2m1,alpha,epsln, &
        sqtep,sneps,S,alphai

    ! Connecting to Nozawa et al.'s (2006) notation:
    ! U0 = sputparm(ntype,1) surface binding energy (eV)
    ! Mi = mproj(ion) (in amu) (projectile/ion mass)
    ! Md = sputparm(ntype,2) (") (mean mass of a sputtered atoms)
    ! Zi = Zproj(ion) atomic no. of the incident ion
    ! Zd = sputparm(ntype,3) mean atomic no. of sputtered atoms
    ! K = sputparm(ntype,4) sputtering parameter K
    U0 = sputparm(ntype,1)
    Mi = mproj(ion)
    Md = sputparm(ntype,2)
    Zi = Zproj(ion)
    Zd = sputparm(ntype,3)
    K = sputparm(ntype,4)
    ! Calculation and definition of the appropriate threshold
    ! energy, Eth, for sputtering as a function of ( M1/M2 ). 
    ! for   M1/M2 <= 0.3    Eth = U0/gi*(1 - gi)
    ! for   M1/M2 >  0.3    Eth = 8*U0(Mi/Md)**0.33333
    mu = Md/Mi
    if((Mi/Md) .le. 0.3) then
        gi = 4.0*Mi*Md/(Mi + Md)**2        
        Eth = U0/(gi*(1. - gi))
    else 
        Eth = 8.0*U0/mu**(1./3.)
    endif
    if(energy.gt.Eth) then
        ! Define initial mass independent parameters  ('asc' is in cm)
        ! asc = 0.885*ao/(Z1**0.667 + Z2**0.667)**0.5
        ! ao is the Bohr radius ( 0.529 Angstrom )
        asc = 4.6832E-9/sqrt(Zi**(2./3.) + Zd**(2./3.))
        ! Calculate the dimensionless function alpha.
        ! for   0.5 < M2/M1 < 10    alpha = 0.3 * ( M2/M1 )**0.667
        ! for         M2/M1 < 0.5   alpha = constant ( = 0.2 )
        ! for    5  < M2/M1   correction factor Rp/R is important.
        ! N.B. for the numbers used here M2/M1 is always > 0.5
        ! this one is alpha_sc
        if(mu.le.0.5) then
            alpha = 0.2
        else if(mu.le.1.0) then
            alpha = 0.1/mu + 0.25*(mu - 0.5)**2
        else
            alpha = 0.3*(mu - 0.6)**(2./3.)
        endif
        ! epsilon = Md/(Mi + Md) * asc/(Zi*Zd*e**2) * E
        ! where E is in ergs and e is in electrostatic units (e.s.u.)
        ! e = 4.802 e-10 e.s.u. (6.93E6 = 1.602E-12/(4.802E-10)**2
        epsln = (Md*asc*energy*6.935D+6)/((Mi + Md)*Zi*Zd)
        !                  3.44 * sqrt(eps) * ln( eps + 2.718 )
        ! Sn(eps) = ---------------------------------------------------
        !           1 + 6.35*sqrt(eps) + eps*(-1.708 + 6.882*sqrt(eps))
        sqtep = sqrt(epsln)
        ! si in Nozawa's notation
        sneps = 3.441*sqtep*log(epsln + 2.718)/(1. + 6.35*sqtep &
            + epsln*(6.882*sqtep - 1.708))
        ! S = 4*pi*alpha*Zi*Zd*e**2*Mi/(Mi + Md)*s
        ! nuclear stopping cross section in units of ergs cm**2
        ! 1.43984E-7 = e**2/1.602176E-12
        S = 4.*PI*asc*Zi*Zd*Mi/(Mi + Md)*sneps*1.43984E-7
        ! Sputtering yield calculation. 
        !         3.56      M1           ( Z1 * Z2 )
        ! Y(E) = ------ * ----- * ----------------------- * alpha 
        !        Uo(eV)   M1+M2   (Z1**.67 + Z2**.67)**.5
        !
        !           Rp              (     ( Eth ).67 )     (     ( Eth ) 2 )
        !         * -- * Sn(eps) *  ( 1 - (-----)    )  *  ( 1 - (-----)   )
        !           R               (     (  E  )    )     (     (  E  )   )
        yield = 4.2E14*S/U0*alpha/(K*mu + 1.)*(1. - (Eth/energy)**(2./3.))* &
            (1. - (Eth/energy))**2
        ! Nozawa's expression for the normal incidence yield:
        !               S(E) alpha_i  
        ! Y0 = 4.2E14 * ---  -------- (1 - (Eth/E)**(2/3))(1 - Eth/E)**2
        !               U_0  K*mu + 1 

        ! The angle averaged sputtering yield is twice that for normal
        ! incidence, hence the factor of two in the final yield.
        yield = 2.0 * yield
    else
        yield = 0.
    endif
end subroutine pt_sputYield
