!!****if* source/Particles/ParticlesForces/shortRange/drag
!!
!! NAME
!!
!!  Particles_shortRangeForce
!!
!! SYNOPSIS
!!
!!  Particles_shortRangeForce(real,intent(inout) :: particles(:,:),
!!                           integer,intent(in) :: p_count,
!!                           integer,intent(in) :: mapType)
!!
!! DESCRIPTION
!!
!!  Computes short-range forces on particles, ie. forces which couple particles
!!  to the gas.
!!  
!! ARGUMENTS
!!
!!   particles :: the list of particles to be operated on
!!   p_count   :: count of the particles in the list
!!   mapType   :: when mapping grid quantities to particle, method to use
!!
!! NOTES
!!
!!***

!=======================================================================

subroutine Particles_shortRangeForce(particles,p_count,mapType)
  use Driver_interface, only : Driver_getSimTime
  use Grid_interface, ONLY : Grid_mapMeshToParticles
  use Particles_data, ONLY : pt_posAttrib,pt_geometry
  use Simulation_data, ONLY : sim_mubar, sim_nnh, sim_usetsput, sim_useisput, &
      sim_ptype, sim_pactive, sim_smallT
  use pt_chargeData
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Eos_data, ONLY : eos_singleSpeciesA
  use Cool_data, ONLY : cl_logXT, cl_N, cl_logT, cl_dlogTdlogXT, cl_dene, &
      cl_dnedlogT
!-------------------------------------------------------------------------------

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"
!-------------------------------------------------------------------------
  integer, intent(IN) :: p_count,mapType
  real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles
  real,dimension(NPART_PROPS,p_count) :: tmpparts
  integer,parameter :: part_props=NPART_PROPS
  integer :: numAttrib 
  ! numAttrib is the no. of attributes that need to be mapped - for now this
  ! is the velocity components (NDIM) + 1 for the density
  parameter(numAttrib=8)
  integer,dimension(2,numAttrib) :: attrib
  integer :: i, ptype, itbl
  real :: pmass, prad, vxrel, vyrel, vzrel, vrel, gdens, Gfac, gtemp, Tfac, &
      afac, Tconv, dmgdt, dens, hdens, ke2, ec, Z, phi, e, c, kB, Na, &
      vx, vy, vz, phimax, phimin, phang, logchiT, logT, dene, edens, dlgT, &
      gam, r, vp, vr, vtot, vxg, vyg, vzg

  call PhysicalConstants_get("electron charge",e)
  call PhysicalConstants_get("Boltzmann",kB)
  call PhysicalConstants_get("speed of light",c)
  call PhysicalConstants_get("Avogadro",Na)
  ke2 = kB/e**2
  ec = 5.679571E-20 ! e/c*sqrt(4*pi) the sqrt(4*pi) is to make B in gauss
  Tconv = Na*sim_mubar/eos_singleSpeciesA
  Tfac = 6.250305E-16/sim_mubar != 128./(9.*PI)*(kB/sim_mubar)
  attrib(GRID_DS_IND,1) = VELX_VAR
  attrib(GRID_DS_IND,2) = VELY_VAR
  attrib(GRID_DS_IND,3) = VELZ_VAR
  attrib(GRID_DS_IND,4) = DENS_VAR
  attrib(GRID_DS_IND,5) = TEMP_VAR
  attrib(GRID_DS_IND,6) = MAGX_VAR
  attrib(GRID_DS_IND,7) = MAGY_VAR
  attrib(GRID_DS_IND,8) = MAGZ_VAR
  attrib(PART_DS_IND,1) = VELX_PART_PROP
  attrib(PART_DS_IND,2) = VELY_PART_PROP
  attrib(PART_DS_IND,3) = VELZ_PART_PROP
  attrib(PART_DS_IND,4) = GDEN_PART_PROP
  attrib(PART_DS_IND,5) = GTMP_PART_PROP
  attrib(PART_DS_IND,6) = MAGX_PART_PROP
  attrib(PART_DS_IND,7) = MAGY_PART_PROP
  attrib(PART_DS_IND,8) = MAGZ_PART_PROP
  tmpparts = particles
  ! To avoid any problems in mapping for 1D runs
  !if(NDIM < 3) then
  !    tmpparts(POSZ_PART_PROP,:) = 0.
  !    if(NDIM < 2) then
  !        tmpparts(POSY_PART_PROP,:) = 0.
  !    endif
  !endif
  call Grid_mapMeshToParticles(tmpparts, part_props, BLK_PART_PROP, &
      p_count, pt_posAttrib, numAttrib, attrib, mapType)
  do i = 1, p_count
      pmass = particles(MASS_PART_PROP,i)
      gtemp = max(tmpparts(GTMP_PART_PROP,i), sim_smallT)
      logchiT = log10(Tconv*gtemp)
      call ut_hunt(cl_logXT, cl_N, logchiT, itbl)
      logT = cl_logT(itbl) + cl_dlogTdlogXT(itbl)*(logchiT - cl_logXT(itbl))
      dlgT = (logT - cl_logT(itbl))
      dene  = cl_dene(itbl) + cl_dnedlogT(itbl)*dlgT
      gtemp = 10.**logT
      particles(GTMP_PART_PROP,i) = gtemp
      particles(GDEN_PART_PROP,i) = tmpparts(GDEN_PART_PROP,i)
      particles(MAGX_PART_PROP,i) = tmpparts(MAGX_PART_PROP,i)
      particles(MAGY_PART_PROP,i) = tmpparts(MAGY_PART_PROP,i)
      particles(MAGZ_PART_PROP,i) = tmpparts(MAGZ_PART_PROP,i)
      ! Note: DENS_PART_PROP should point to the value of the grain solid
      ! material density (g/cm^3) which should be unchanging (unlike DENS_VAR
      ! which points to the gas density). Could also set a single value,
      ! sim_pdens, and use that, but this way you could set each grain to have
      ! a different grain type if desired.

      prad = (3.*pmass/(4.*PI*particles(DENS_PART_PROP,i)))**(1./3.)
      particles(GVLX_PART_PROP,i) = tmpparts(VELX_PART_PROP,i)
      if(NDIM > 1) then
          particles(GVLY_PART_PROP,i) = tmpparts(VELY_PART_PROP,i)
      else
          particles(GVLY_PART_PROP,i) = 0.
      endif
      ! If using cylindrical symmetry we assume VELZ_VAR holds the
      ! gas angular velocity in phi direction in radians/s (vphi = dphi/dt)
      if(NDIM > 2) then
          particles(GVLZ_PART_PROP,i) = tmpparts(VELZ_PART_PROP,i)
      else
          particles(GVLZ_PART_PROP,i) = 0.
      endif
      gdens = tmpparts(GDEN_PART_PROP,i)
      dens = gdens/sim_mubar
      hdens = dens/sim_nnh
      edens = dene*dens
      if(pt_geometry == CYLINDRICAL) then
          r = particles(POSX_PART_PROP,i)
          vr = particles(VELX_PART_PROP,i) 
          vp =  particles(VELZ_PART_PROP,i)
          ! needed because vzrel is angular speed
          ! here we transform to cartesian coordinates to calculate accelerations
          ! note that vp = dphi/dt, accp = d^2(phi)/dt^2
          phang = particles(POSZ_PART_PROP,i)
          vx = vr*cos(phang) - vp*r*sin(phang)
          vxg = particles(GVLX_PART_PROP,i)*cos(phang) - &
              particles(GVLZ_PART_PROP,i)*r*sin(phang)
          vxrel = vx - vxg
          vy = vr*sin(phang) + vp*r*cos(phang)
          vyg = particles(GVLX_PART_PROP,i)*sin(phang) + &
              particles(GVLZ_PART_PROP,i)*r*cos(phang)
          vyrel = vy - vyg
          vz = particles(VELY_PART_PROP,i)
          vzrel = vz - particles(GVLY_PART_PROP,i)
          vrel = sqrt(vxrel**2 + vyrel**2 + vzrel**2)
          vtot = sqrt(vx**2 + vy**2 + vz**2)
      else
          vx = particles(VELX_PART_PROP,i)
          vy = particles(VELY_PART_PROP,i)
          vz = particles(VELZ_PART_PROP,i)
          vxg = particles(GVLX_PART_PROP,i)
          vyg = particles(GVLY_PART_PROP,i)
          vzg = particles(GVLZ_PART_PROP,i)
          vxrel = vx - vxg
          vyrel = vy - vyg
          vzrel = vz - vzg
          vrel = sqrt(vxrel**2 + vyrel**2 + vzrel**2)
          vtot = sqrt(vx**2 + vy**2 + vz**2)
      endif
      call pt_grain_pot(edens,gtemp,vrel,sim_ptype,phi)
      ! phi = Z*e^2/akT
      ! Added limitation on phi due to field emission - from Draine & Sutin via
      ! McKee et al. (1987) - 5/9/2022
      phimax = 34.8*(prad*1.E6)/(gtemp*1.E-5)
      phimin = -1.16*(prad*1.E6)/(gtemp*1.E-5)
      phi = min(max(phi,phimin),phimax)
      Z = phi*prad*gtemp*ke2
      particles(CHRG_PART_PROP,i) = Z
      ! Now incorporating vrel into Gfac so no danger of divide by 0
      Gfac = sqrt(vrel**2 + Tfac*max(min(gtemp,1.E8),10.))
      ! Note: this now includes the thermal part of the drag (see McKee et al.
      ! 1987) via the second term in Gfac 
      gam = 1./sqrt(1. - (vtot/c)**2)
      afac = -PI*prad*prad*gdens*Gfac/pmass
      ! accelerations now only include the drag - Lorentz force is incorporated
      ! into Particles_advance
      particles(ACCX_PART_PROP,i) = vxrel*afac/gam
      particles(ACCY_PART_PROP,i) = vyrel*afac/gam
      particles(ACCZ_PART_PROP,i) = vzrel*afac/gam
      ! Currently integer values for particle properties are not supported, 
      ! so if we wanted to add a per particle type would need to give it as 
      ! a real value. Instead here we're using a global value sim_ptype.
      ptype = sim_ptype 
      particles(VREL_PART_PROP,i) = vrel
      if(.not.sim_usetsput) gtemp = 1.E3
      if(.not.sim_useisput) vrel = 1.E6
      call pt_sputrate(vrel,prad,hdens,gtemp,ptype,dmgdt)
      particles(DMDT_PART_PROP,i) = -dmgdt
  enddo
  return
end subroutine Particles_shortRangeForce

!=====================================================================
