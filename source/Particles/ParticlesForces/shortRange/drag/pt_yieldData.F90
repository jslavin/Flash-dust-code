!!****if* source/Particles/ParticlesForces/shortRange/drag/pt_yieldData
!!
!! NAME
!!
!!  pt_yieldData
!!
!! SYNOPSIS
!!
!!  use pt_yieldData 
!!
!!  DESCRIPTION
!!
!!  Stores the data for calculating the sputtering yield for dust
!!  
!! PARAMETERS
!!
!!     mproj            atomic mass (A) of gas species 'ion' (amu) 
!!     Zproj            atomic number (Z) of gas species 'ion'
!!     Cproj            charge of gas ion species 'ion'
!!     amu              amu in g
!!     nion             no. of different (projectile) ions included
!!     sputparm         array including sputtering parameters:
!!       sputparm(ntype,1) surface binding energy (U) of target (eV)
!!       sputparm(ntype,2) atomic mass of sputtered fragment (amu), mtarg
!!       sputparm(ntype,3) atomic number of sputtered fragment, Ztarg
!!       sputparm(ntype,4) free parameter "K" constant
!!    abund             array of the gas phase abundances of the gas species
!!                      Note: includes both the ionization fraction and 
!!                      elemental abundance abund(i) = Ai * X(ion)
!!                      so n(X) = nH*abund(i) where X is one of H+, He+, etc.
!!
!!***
module pt_yieldData

   real, dimension(6), save :: mproj, zproj, cproj, abund
   real, dimension(2,4), save :: sputparm
   real, save :: amu
   integer, save :: nion
end module pt_yieldData
