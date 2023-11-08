!!****if* source/Particles/ParticlesForces/shortRange/drag/pt_yieldIntgrlData
!!
!! NAME
!!
!!  pt_yieldIntgrlData
!!
!! SYNOPSIS
!!
!!  use pt_yieldIntgrlData 
!!
!!  DESCRIPTION
!!
!!  Stores the data for calculating the sputtering yield integral for dust.
!!  The yield integral is the integral of the yield over the velocity
!!  distribution function, assumed to be a skewed (or not skewed) maxwellian
!!  
!! PARAMETERS
!!
!!     mproj            atomic mass (A) of gas species 'ion' (amu) 
!!     Zproj            atomic number (Z) of gas species 'ion'
!!     Cproj            charge of gas ion species 'ion'
!!     amu              amu in g
!!     nion             no. of different (projectile) ions included
!!     vtbl             grid of velocties (cm/s) in yield integral table
!!     Ttbl             grid of temperatures (K) in yield integral table
!!     sputparm         array including sputtering parameters:
!!       sputparm(ntype,1) surface binding energy (U) of target (eV)
!!       sputparm(ntype,2) atomic mass of sputtered fragment (amu), mtarg
!!       sputparm(ntype,3) atomic number of sputtered fragment, Ztarg
!!       sputparm(ntype,4) free parameter "K" constant
!!    abund             array of the gas phase abundances of the gas species
!!                      Note: includes both the ionization fraction and 
!!                      elemental abundance abund(i) = Ai * X(ion)
!!                      so n(X) = nH*abund(i) where X is one of H+, He+, etc.
!!    Yield integral tables:
!!    Yv_SiH            silicate grains impacted by H
!!    Yv_SiHe           silicate grains impacted by He
!!    Yv_SiO            silicate grains impacted by O
!!    Yv_SiC            silicate grains impacted by C
!!    Yv_aCH            carbonaceous grains impacted by H
!!    Yv_aCHe           carbonaceous grains impacted by He
!!    Yv_aCO            carbonaceous grains impacted by O
!!    Yv_aCC            carbonaceous grains impacted by C
!!
!!***
module pt_yieldIntgrlData

   real, dimension(6), save :: mproj, zproj, cproj, abund
   real, dimension(2,4), save :: sputparm
   real, save :: amu
   integer, save :: nion
   integer, save :: nvtbl,nTtbl
   real, save :: dvtbl,dTtbl
   real, dimension(61), save :: vtbl, Ttbl
   real, dimension(61,61), save :: Yv_SiH, Yv_SiHe, Yv_SiO, Yv_SiC, &
       Yv_aCH, Yv_aCHe,Yv_aCO,Yv_aCC
end module pt_yieldIntgrlData
