!!****if* source/Particles/ParticlesMain/active/dust/pt_advanceDust_AB
!!
!! NAME
!!
!!  pt_advanceDust_AB
!!
!! SYNOPSIS
!!
!!  call pt_advanceDust_AB(real(in)    :: dtold,
!!                      real(in)    :: dtnew,
!!                      real(inout) :: particles(NPART_PROPS,p_count),
!!                      integer(in) :: p_count,
!!                      integer(in) :: ind)
!!
!! DESCRIPTION
!!
!!   Advances particles in time
!!
!! ARGUMENTS
!!
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particles -- particles on which to operate
!!   p_count - the number of particles in the list to advance
!!   ind -- index into pt_typeInfo 
!!
!! DESCRIPTION
!!
!!   Time advancement routine for dust particles. Using 2nd order
!!   Adams-Bashforth method. Currently assumes cylindrical symmetry. Mass 
!!   is evolved as well.
!!
!!***

subroutine pt_advanceDust_AB(dtOld,dtNew,particles,p_count,ind)

  use Simulation_data, only : sim_pactive
  use Driver_interface, only : Driver_getSimTime
  use Driver_data, only : dr_dtInit
  use Particles_data, only: useParticles, pt_restart, pt_meshMe, &
       pt_numLocal, pt_maxPerProc, pt_indexCount, pt_indexList, pt_geometry
  use Timers_interface, only : Timers_start, Timers_stop
  use Particles_interface, ONLY : Particles_shortRangeForce

  implicit none

#include "Flash.h"
#include "constants.h"  
#include "Particles.h"

  real, intent(in)   :: dtOld, dtNew
  integer,intent(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles

  integer       :: i
  real          :: wterm, woldterm, mass0, mratio, dtN
  real, parameter :: onethird = 1./3., onesixth = 1./6.
  logical, save :: firstCall = .true., firstActive = .false.
  integer,parameter :: mapType=QUADRATIC
  real, parameter :: minmass = 1.E-20
  logical :: pdiag
  real :: r,vtot,fac,r0,p0,x0,y0,vr0,vp0,vx0,vy0,accx0,accy0,r1,p1,x1,y1, &
      vr1,vp1,vx1,vy1,accx1,accy1,x2,y2,r2,p2,vx2,vy2,vr2,vp2,z0,vz0,z1, &
      vz1,accz0,accz1,z2,vz2
  real, parameter :: vmax = 2.5E9 ! upper limit on grain speed = 15,000 km/s
!-------------------------------------------------------------------------------
  if ((.not. useParticles).or.(.not. sim_pactive)) return

  if (firstCall .and. .not. pt_restart) then
     wterm = 0.5*dtNew
     woldterm = 0.
     firstCall = .false.
  else
     wterm = 0.5*dtNew + onethird*dtOld + onesixth*dtNew**2/dtOld
     woldterm = onesixth*(dtOld**2 - dtNew**2)/dtOld
  endif

  call Timers_start ("particles")

!-------------------------------------------------------------------------------
    
! Compute current components of the acceleration of each particle


  call Timers_start ("particle forces")

  do i = 1, p_count
     particles(ACCX_PART_PROP,i) = 0.
     particles(ACCY_PART_PROP,i) = 0.  
     particles(ACCZ_PART_PROP,i) = 0.
  enddo

  !write(*,*) 'Calling Particles_shortRangeForce from pt_advanceDust'
  call Particles_shortRangeForce(particles,p_count,mapType)
  call Timers_stop ("particle forces")
  !call Driver_getSimTime(t)
  !write(*,'("t =",ES12.5," sim_tpactive =",ES12.5)') t,sim_tpactive
!-------------------------------------------------------------------------------

! Apply the Euler step to update particle positions and velocities

  call Timers_start ("move particles")

  pdiag = .false.
  dtN = dtNew
  if(sim_pactive) then
      if(firstActive) then
          write(*,*) 'firstActive = .true.'
          dtN = dr_dtInit
          wterm = 0.5*dtN
          woldterm = 0.
          firstActive = .false.
          !pdiag = .true.
          !open(95,file='part_init_vals.txt',status='unknown',access='append')
          !write(95,'("# before update")')
      endif
      do i = 1, p_count
          if(pdiag) then
              write(95,'(6(ES12.4))') particles(POSX_PART_PROP,i), &
                  particles(POSY_PART_PROP,i),particles(VELX_PART_PROP,i), &
                  particles(VELY_PART_PROP,i),particles(ACCX_PART_PROP,i), &
                  particles(ACCY_PART_PROP,i)
          endif
          mass0 = particles(MASS_PART_PROP,i)
          particles(MASS_PART_PROP,i) = particles(MASS_PART_PROP,i) + &
              dtN/2.*(3.*particles(DMDT_PART_PROP,i) - & 
              particles(ODMT_PART_PROP,i))
          ! need to guard against divide by zero
          particles(MASS_PART_PROP,i) = max(particles(MASS_PART_PROP,i),minmass)
          particles(ODMT_PART_PROP,i) = particles(DMDT_PART_PROP,i)

          if (pt_geometry == CYLINDRICAL) then
              r0 = particles(OPSX_PART_PROP,i)
              p0 = particles(OPSZ_PART_PROP,i)
              x0 = r0*cos(p0)
              y0 = r0*sin(p0)
              vr0 = particles(OVLX_PART_PROP,i)
              vp0 = particles(OVLZ_PART_PROP,i)
              vx0 = vr0*cos(p0) - r0*vp0*sin(p0)
              vy0 = vr0*sin(p0) + r0*vp0*cos(p0)
              accx0 = particles(OACX_PART_PROP,i)
              accy0 = particles(OACY_PART_PROP,i)
              r1 = particles(POSX_PART_PROP,i)
              p1 = particles(POSZ_PART_PROP,i)
              x1 = r1*cos(p1)
              y1 = r1*sin(p1)
              vr1 = particles(VELX_PART_PROP,i)
              vp1 = particles(VELZ_PART_PROP,i)
              vx1 = vr1*cos(p1) - r1*vp1*sin(p1)
              vy1 = vr1*sin(p1) + r1*vp1*cos(p1)
          else
              x0 = particles(OPSX_PART_PROP,i)
              y0 = particles(OPSY_PART_PROP,i)
              z0 = particles(OPSZ_PART_PROP, i)
              x1 = particles(POSX_PART_PROP, i)
              y1 = particles(POSY_PART_PROP, i)
              z1 = particles(POSZ_PART_PROP, i)
              vx0 = particles(OVLX_PART_PROP,i)
              vx1 = particles(VELX_PART_PROP,i)
              vy0 = particles(OVLY_PART_PROP,i)
              vy1 = particles(VELY_PART_PROP,i)
              vz0 = particles(OVLZ_PART_PROP,i)
              vz1 = particles(VELZ_PART_PROP,i)
          endif
          accx1 = particles(ACCX_PART_PROP,i)
          accy1 = particles(ACCY_PART_PROP,i)
 
          x2 = x1 + dtN/2.*(3.*vx1 - vx0)
          y2 = y1 + dtN/2.*(3.*vy1 - vy0)
          vx2 = vx1 + dtN/2.*(3.*accx1 - accx0)
          vy2 = vy1 + dtN/2.*(3.*accy1 - accy0)
          if (pt_geometry == CYLINDRICAL) then
              r2 = sqrt(x2**2 + y2**2)
              p2 = atan2(y2,x2)
              vr2 = vx2*cos(p2) + vy2*sin(p2)
              vp2 = (vy2*cos(p2) - vx2*sin(p2))/r2
          else
              accz1 = particles(ACCZ_PART_PROP,i)
              z2 = z1 + dtN/2.*(3.*vz1 - vz0)
              vz2 = vz1 + dtN/2.*(3.*accz1 - accz0)
          endif

          if (pt_geometry == CYLINDRICAL) then
              !!! In cylindrical symmetry, X corresponds with r (cylindical) !!!
              particles(OPSX_PART_PROP,i) = particles(POSX_PART_PROP,i) 
              particles(POSX_PART_PROP,i) = r2
              particles(OVLX_PART_PROP,i) = particles(VELX_PART_PROP,i)
              particles(VELX_PART_PROP,i) = vr2
              particles(OACX_PART_PROP,i) = particles(ACCX_PART_PROP,i)
              !!! In cylindrical symmetry, Y corresponds with z !!!
              particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + &
                  dtN/2.*(3.*particles(VELY_PART_PROP,i) - &
                  particles(OVLY_PART_PROP,i))
              particles(OVLY_PART_PROP,i) = particles(VELY_PART_PROP,i)
              particles(VELY_PART_PROP,i) = particles(VELY_PART_PROP,i) + &
                  dtN/2.*(3.*particles(ACCZ_PART_PROP,i) - &
                  particles(OACZ_PART_PROP,i))
              particles(OACY_PART_PROP,i) = particles(ACCY_PART_PROP,i)
              !!! In cylindrical symmetry, Z corresponds with phi !!!
              particles(OPSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) 
              particles(POSZ_PART_PROP,i) = p2
              particles(OVLZ_PART_PROP,i) = particles(VELZ_PART_PROP,i)
              particles(VELZ_PART_PROP,i) = vp2
              particles(OACZ_PART_PROP,i) = particles(ACCZ_PART_PROP,i)
          else
              particles(OPSX_PART_PROP,i) = particles(POSX_PART_PROP,i) 
              particles(POSX_PART_PROP,i) = x2
              particles(OVLX_PART_PROP,i) = particles(VELX_PART_PROP,i)
              particles(VELX_PART_PROP,i) = vx2
              particles(OACX_PART_PROP,i) = particles(ACCX_PART_PROP,i)
              particles(OPSY_PART_PROP,i) = particles(POSY_PART_PROP,i) 
              particles(POSY_PART_PROP,i) = y2
              particles(OVLY_PART_PROP,i) = particles(VELY_PART_PROP,i)
              particles(VELY_PART_PROP,i) = vy2
              particles(OACY_PART_PROP,i) = particles(ACCY_PART_PROP,i)
              particles(OPSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) 
              particles(POSZ_PART_PROP,i) = z2
              particles(OVLZ_PART_PROP,i) = particles(VELZ_PART_PROP,i)
              particles(VELZ_PART_PROP,i) = vz2
              particles(OACZ_PART_PROP,i) = particles(ACCZ_PART_PROP,i)
          endif

          if (pt_geometry == CYLINDRICAL) then
              if(particles(POSZ_PART_PROP,i) < -PI) &
                  particles(POSZ_PART_PROP,i) = 2.*PI +  &
                  particles(POSZ_PART_PROP,i)
              if(particles(POSZ_PART_PROP,i) > PI) &
                  particles(POSZ_PART_PROP,i) = &
                  particles(POSZ_PART_PROP,i) - 2.*PI
              vtot = sqrt(particles(VELX_PART_PROP,i)**2 + &
                  particles(VELY_PART_PROP,i)**2 + &
                  (particles(VELZ_PART_PROP,i)*particles(POSX_PART_PROP,i))**2)
          else
              vtot = sqrt(particles(VELX_PART_PROP,i)**2 + &
                  particles(VELY_PART_PROP,i)**2 + &
                  particles(VELZ_PART_PROP,i)**2)
          endif
          if (vtot > vmax) then
              fac = vmax/vtot
              particles(VELX_PART_PROP,i) = fac*particles(VELX_PART_PROP,i)
              particles(VELY_PART_PROP,i) = fac*particles(VELY_PART_PROP,i)
              particles(VELZ_PART_PROP,i) = fac*particles(VELZ_PART_PROP,i)
          endif
      enddo
      if(pdiag) then
          write(95,'("# after update")')
          do i = 1, p_count
              write(95,'(6(ES12.4))') particles(POSX_PART_PROP,i), &
                  particles(POSY_PART_PROP,i),particles(VELX_PART_PROP,i), &
                  particles(VELY_PART_PROP,i),particles(ACCX_PART_PROP,i), &
                  particles(ACCY_PART_PROP,i)
          enddo
          close(95)
      endif
  endif
  call Timers_stop ("move particles")
!-------------------------------------------------------------------------------
  call Timers_stop ("particles")
  return

end subroutine pt_advanceDust_AB
