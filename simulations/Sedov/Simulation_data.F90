!!****if* source/Simulation/SimulationMain/Sedov/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Sedov
!!  
!! PARAMETERS
!!
!!  sim_Eej             Explosion energy (erg)
!!  sim_Mej             Ejecta mass (g)
!!  sim_vej             Maximum ejecta speed (cm/s)
!!  sim_chi             Density ratio of clumpy ejecta to smooth ejecta
!!  sim_fcl             Volume filling factor of clumpy ejecta within the core
!!  sim_ufac            Factor related to ratio of core radius to envelope
!!                      radius - needed for root finding
!!  sim_eta             ratio of average density to smooth ejecta density in the
!                       ejecta core
!!  sim_ejpl            power law exponent for ejecta envelope
!!  sim_rcl             radius of the clumps in the clumpy ejecta
!!  sim_rhoCSM          pre-shock density in the stellar wind shell at the shock
!!                      radius in 2004
!!  sim_Mprog           Mass of the progenitor star
!!  sim_Msh             Mass of the circumstellar medium (i.e. the stellar wind
!!                      shell) that had been shocked (in 2004)
!!  sim_Rcore           Radius of inner (core) region of ejecta
!!  sim_Rej             Radius of the outer edge of the ejecta
!!  sim_Rb              Shock radius in 2004
!!  sim_Ro              Radius of the outer edge of the stellar wind shell
!!  sim_pISM            Initial ambient (ISM) pressure
!!  sim_rhoISM          Initial ambient density
!!  sim_xctr            Explosion center coordinates
!!  sim_yctr            Explosion center coordinates
!!  sim_zctr            Explosion center coordinates
!!  sim_nsubzones       Number of `sub-zones' in cells for applying 1d profile
!!  sim_mubar           Mean mass per nucleus = 
!!                      (mH + AHe*mHe + ...)/(1 + AHe + ...)
!!                      where mH = H mass, mHe = He mass, AHe = He abundance
!!  sim_nnh             ratio of total number density to H number density
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save :: sim_Eej, sim_Mej, sim_vej, sim_chi, sim_fcl, sim_ufac, sim_eta
  integer, save :: sim_ejpl
  real, save :: sim_rcl, sim_rhoCSM, sim_Mprog, sim_Msh, sim_Rej, sim_Rb
  real, save :: sim_rhosm, sim_rhocl, sim_Ro, sim_Rcore, sim_pISM, sim_rhoISM 
  real, save :: sim_mubar, sim_nnh, sim_gamma, sim_xCenter, sim_yCenter
  real, save :: sim_zCenter, sim_smallX, sim_smallRho, sim_smallP, sim_smallT
  real, save :: sim_smallE, sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin
  real, save :: sim_zMax, pCoreFac, pCorePL
  integer, save :: sim_nSubZones
  character(len=80), save :: sim_dustdata, sim_yielddata, sim_yieldintdata, &
      sim_chargedata, sim_clumpdata, sim_partdata
  logical, save :: sim_usetsput, sim_useisput, sim_usedrag, sim_readparts

  !! *** Variables pertaining to this Simulation *** !!

  integer, parameter                  :: sim_nProfile = 10000
  !! maximum no. of clumps - if we go to 3D we will need more
  integer, parameter                  :: sim_clmax = 500
  !! sim_partmax must be at least sim_ncl*sim_nperclump
  integer, parameter                  :: sim_partmax = 5000
  real   , save                       :: sim_inSubZones, sim_inSubzm1
  real   , save                       :: sim_inszd
  real, dimension(sim_nProfile), save :: sim_rProf, sim_rhoProf, sim_pProf
  real, dimension(sim_nProfile), save :: sim_vProf
  real, dimension(sim_clmax), save  :: xcl,ycl,zcl,clrad,clrho
  real, dimension(sim_partmax), save  :: xpart,ypart,zpart,vxpart,vypart,vzpart
  logical, save :: sim_pactive
  integer, save                       :: sim_ncl, npart  
  real, save                          :: sim_drProf, sim_pExp, sim_vctr

  !! *** variables related to dust *** !!
  integer, save :: sim_nperclump, sim_ptype
  real, save :: sim_pmass, sim_pdens, sim_G0
  !! *** variables related to the B field *** !!
  real, save :: sim_Bx0, sim_By0, sim_Bz0
  logical, save :: sim_killdivb
  !! *** variables related to explosions *** !!
  character(LEN=80) :: sim_exptimesfile
  real, dimension(10), save :: texp,xexp,yexp,zexp
  integer, save :: nexps
  logical, dimension(10), save :: exploded
  !! *** variables related to the stellar wind *** !!
  real, save :: sim_vwind, sim_MLR, sim_rInit, sim_rInject
  real :: sim_SNoutint

  integer, save :: sim_meshMe
  logical, save :: sim_threadBlockList = .false.
  logical, save :: sim_threadWithinBlock = .false.

end module Simulation_data
