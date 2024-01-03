!!****if* source/Particles/ParticlesForces/shortRange/drag_AB
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
      sim_usedrag, sim_ptype, sim_pactive, sim_smallT
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
  integer :: i,ptype,itbl
  real :: pmass, prad, vxrel, vyrel, vzrel, vrel, gdens, Gfac, gtemp, Tfac, &
      afac, Tconv, dmgdt, dens, hdens, ke2, ec, EMACCX, EMACCY, EMACCZ, Z, &
      phi, e, c, kB, Na, vx, vy, vr, vz, vp, vtot, r, omx, omy, omz, gam, &
      phimax, phimin, accr, accp, accz, phang, logchiT, logT, dene, dlgT, &
      edens, accx, accy

  call PhysicalConstants_get("electron charge",e)
  call PhysicalConstants_get("Boltzmann",kB)
  call PhysicalConstants_get("speed of light",c)
  call PhysicalConstants_get("Avogadro",Na)
  ke2 = kB/e**2 ! 598.44226
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
  call Grid_mapMeshToParticles(tmpparts,part_props, BLK_PART_PROP, &
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
      if(NDIM == 3) then
          particles(MAGZ_PART_PROP,i) = tmpparts(MAGZ_PART_PROP,i)
      else
          particles(MAGZ_PART_PROP,i) = 0.
      endif
      ! Note: DENS_PART_PROP should point to the value of the grain solid
      ! material density (g/cm^3) which should be unchanging (unlike DENS_VAR
      ! which points to the gas density). Could also set a single value,
      ! sim_pdens, and use that, but this way you could set each grain to have
      ! a different grain type if desired.
      prad = (3.*pmass/(4.*PI*particles(DENS_PART_PROP,i)))**(1./3.)
      vx = particles(VELX_PART_PROP,i)
      vy = particles(VELY_PART_PROP,i)
      vz = particles(VELZ_PART_PROP,i)
      vxrel = vx - tmpparts(VELX_PART_PROP,i)
      particles(GVLX_PART_PROP,i) = tmpparts(VELX_PART_PROP,i)
      if(NDIM > 1) then
          vyrel = vy - tmpparts(VELY_PART_PROP,i)
          particles(GVLY_PART_PROP,i) = tmpparts(VELY_PART_PROP,i)
      else
          vyrel = vy
          particles(GVLY_PART_PROP,i) = 0.
      endif
      ! If using cylindrical symmetry we assume VELZ_VAR holds the
      ! gas angular velocity in phi direction in radians/s (vphi = dphi/dt)
      if(NDIM > 2) then
          vzrel = vz - tmpparts(VELZ_PART_PROP,i)
          particles(GVLZ_PART_PROP,i) = tmpparts(VELZ_PART_PROP,i)
      else
          vzrel = vz
          particles(GVLZ_PART_PROP,i) = 0.
      endif
      gdens = tmpparts(GDEN_PART_PROP,i)
      dens = gdens/sim_mubar
      hdens = dens/sim_nnh
      edens = dene*dens
      if(pt_geometry == CYLINDRICAL) then
          r = particles(POSX_PART_PROP,i)
          ! needed because vzrel is angular speed
          vrel = sqrt(vxrel**2 + vyrel**2 + (vzrel*r)**2)
      else
          vrel = sqrt(vxrel**2 + vyrel**2 + vzrel**2)
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
      omx = ec*Z/pmass*tmpparts(MAGX_PART_PROP,i)
      omy = ec*Z/pmass*tmpparts(MAGY_PART_PROP,i)
      ! Note: in 2D omz = 0 (though not in 2.5D)
      omz = ec*Z/pmass*tmpparts(MAGZ_PART_PROP,i)
      if(pt_geometry == CYLINDRICAL) then
          vr = vx ! assumes cylindrical symmetry
          ! we assume particles(VELZ_PART_PROP) = d(phi)/dt
          vp = vz ! assumes cylindrical symmetry
          ! if VELZ_PART_PROP is angular speed (d(phi)/dt in radians/s) then
          ! Here x is r, y is z and z is phi. In true cyl. coords:
          ! d^2z/dt^2 = omega_phi*v_r - omega_r*(r*dphi/dt)
          EMACCY = omz*vxrel - omx*vzrel*r
          ! d^2r/dt^2 = omega_z*(r*dphi/dt) - omega_phi*v_z + r*(dphi/dt)**2
          EMACCX = omy*vzrel*r - omz*vyrel + r*vp**2
          ! d^2phi/dt^2 = (omega_r*v_z - omega_z*v_r - 2*v_r*dphi/dt)/r
          EMACCZ = omx*vyrel/r - omy*vxrel/r - 2.*vr*vp/r
          vtot = sqrt(vx**2 + vy**2 + (vp*r)**2)
      else
          EMACCY = omz*vxrel - omx*vzrel
          EMACCX = omy*vzrel - omz*vyrel
          EMACCZ = omx*vyrel - omy*vxrel
          vtot = sqrt(vx**2 + vy**2 + vz**2)
      endif
      ! Now incorporating vrel into Gfac so no danger of divide by 0
      Gfac = sqrt(vrel**2 + Tfac*max(min(gtemp,1.E8),10.))
      ! Note: this now includes the thermal part of the drag (see McKee et al.
      ! 1987) via the second term in Gfac 
      gam = 1./sqrt(1. - (vtot/c)**2)
      afac = -PI*prad*prad*gdens*Gfac/pmass
      if(pt_geometry == CYLINDRICAL) then
          if(sim_usedrag) then
              accr = (vxrel*afac + EMACCX)/gam
              accz = (vyrel*afac + EMACCY)/gam
              accp = (vzrel*afac + EMACCZ)/gam
          else
              accr = EMACCX/gam
              accz = EMACCY/gam
              accp = EMACCZ/gam
          endif
          phang = particles(POSZ_PART_PROP,i)
          ! These accelerations are in cartesian coordinates
          particles(ACCX_PART_PROP,i) = (accr - r*vp**2)*cos(phang) - &
              (r*accp + 2.*vr*vp)*sin(phang)
          particles(ACCY_PART_PROP,i) = (accr - r*vp**2)*sin(phang) + &
              (r*accp + 2.*vr*vp)*cos(phang)
          particles(ACCZ_PART_PROP,i) = accz
      else
          accx = (vxrel*afac + EMACCX)/gam
          accy = (vyrel*afac + EMACCY)/gam
          accz = (vzrel*afac + EMACCZ)/gam
          particles(ACCX_PART_PROP,i) = accx
          particles(ACCY_PART_PROP,i) = accy
          particles(ACCZ_PART_PROP,i) = accz
      endif
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
