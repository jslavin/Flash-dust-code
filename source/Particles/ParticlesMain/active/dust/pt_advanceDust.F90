!!****if* source/Particles/ParticlesMain/active/dust/pt_advanceDust
!!
!! NAME
!!
!!  pt_advanceDust
!!
!! SYNOPSIS
!!
!!  call pt_advanceDust(real(in)    :: dtold,
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

subroutine pt_advanceDust(dtOld,dtNew,particles,p_count,ind)

  use Simulation_data, only : sim_pactive, sim_usedrag
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

  integer       :: i, n, Nsteps
  real          :: mass0, mratio, dtN
  real, parameter :: onethird = 1./3., onesixth = 1./6.
  logical, save :: firstActive = .false.
  integer,parameter :: mapType=QUADRATIC
  real, parameter :: minmass = 1.E-20
  logical :: pdiag = .false., do_STS = .true.
  real :: vtot, fac, vr0, vp0, vx0, vy0, vz0, r1, p1, x1, y1, z1, &
      vx1, vy1, vz1, x2, y2, z2, r2, p2, vx2, vy2, vz2, vr2, vp2, &
      accx, accy, accz, Bp, Br, Bx, By, Bz, Btot, omega, ec, tanom, &
      vgr, vgp, vgx, vgy, vgz, vxp, vyp, vzp
  real :: pmass, charge, tgyro, gfac = 0.01
  real, parameter :: vmax = 2.5E9 ! upper limit on grain speed = 25,000 km/s
!-------------------------------------------------------------------------------
  if (.not. useParticles) return

  call Timers_start ("particles")

!-------------------------------------------------------------------------------
  ec = 5.679571E-20 ! e/c*sqrt(4*pi) the sqrt(4*pi) is to make B in gauss    

  call Timers_start ("particle forces")

  do i = 1, p_count
     particles(ACCX_PART_PROP,i) = 0.
     particles(ACCY_PART_PROP,i) = 0.  
     particles(ACCZ_PART_PROP,i) = 0.
  enddo

  !write(*,*) 'Calling Particles_shortRangeForce from pt_advanceDust'
  ! Compute current components of the acceleration of each particle
  call Particles_shortRangeForce(particles,p_count,mapType)
  call Timers_stop ("particle forces")
  !call Driver_getSimTime(t)
  !write(*,'("t =",ES12.5," sim_tpactive =",ES12.5)') t,sim_tpactive
!-------------------------------------------------------------------------------

! Apply the Euler step to update particle positions and velocities

  call Timers_start ("move particles")

  !if(t > sim_tpactive) then
  pdiag = .false.
  dtN = dtNew
  if(sim_pactive) then
      if(firstActive) then
          write(*,*) 'firstActive = .true.'
          dtN = dr_dtInit
          firstActive = .false.
          !pdiag = .true.
          !open(95,file='part_init_vals.txt',status='unknown',access='append')
          !write(95,'("# before update")')
      endif
      if(do_STS) then
          tgyro = HUGE(1.0)
          do i = 1, p_count
              pmass = particles(MASS_PART_PROP,i)
              Btot = sqrt(particles(MAGX_PART_PROP,i)**2 + &
                particles(MAGY_PART_PROP,i)**2 + particles(MAGZ_PART_PROP,i)**2)
              charge = abs(particles(CHRG_PART_PROP,i))
              if((charge > 0.).and.(Btot > 0.)) tgyro = &
                  min(2.*PI*pmass/(ec*charge*Btot), tgyro)
          enddo
          if(dtNew > gfac*tgyro) then
              Nsteps = int(dtNew/(gfac*tgyro)) + 1
              dtN = dtNew/Nsteps
          else
              Nsteps = 1
              dtN = dtNew
          endif
      else
          Nsteps = 1
          dtN = dtNew
      endif
      do n=1,Nsteps
          if(n > 1) call Particles_shortRangeForce(particles,p_count,mapType)
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
              particles(MASS_PART_PROP,i) = max(particles(MASS_PART_PROP,i), &
                  minmass)
              particles(ODMT_PART_PROP,i) = particles(DMDT_PART_PROP,i)

              if (pt_geometry == CYLINDRICAL) then
                  r1 = particles(POSX_PART_PROP,i)
                  p1 = particles(POSZ_PART_PROP,i)
                  x1 = r1*cos(p1)
                  y1 = r1*sin(p1)
                  z1 = particles(POSY_PART_PROP,i)
                  vr0 = particles(VELX_PART_PROP,i)
                  vp0 = particles(VELZ_PART_PROP,i)
                  vx0 = vr0*cos(p1) - r1*vp0*sin(p1)
                  vy0 = vr0*sin(p1) + r1*vp0*cos(p1)
                  vz0 = particles(VELY_PART_PROP,i)

                  vgr = particles(GVLX_PART_PROP,i)
                  vgz = particles(GVLY_PART_PROP,i)
                  vgp = particles(GVLZ_PART_PROP,i)
                  vgx = vgr*cos(p1) - r1*vgp*sin(p1)
                  vgy = vgr*sin(p1) + r1*vgp*cos(p1)

                  Btot = sqrt(particles(MAGX_PART_PROP,i)**2 + &
                      particles(MAGY_PART_PROP,i)**2 + &
                      particles(MAGZ_PART_PROP,i)**2)
                  omega = particles(CHRG_PART_PROP,i)*ec*Btot/ &
                      particles(MASS_PART_PROP,i)
                  tanom = tan(omega*dtN/2.)
                  Br = particles(MAGX_PART_PROP,i)
                  Bz = particles(MAGY_PART_PROP,i)
                  Bp = particles(MAGZ_PART_PROP,i)
                  Bx = Br*cos(p1) - Bp*sin(p1)
                  By = Br*sin(p1) + Bp*cos(p1)
                  ! This is now just drag in cartesian coords
                  accx = particles(ACCX_PART_PROP,i)
                  accy = particles(ACCY_PART_PROP,i)
                  accz = particles(ACCZ_PART_PROP,i)
                  ! apply drag in operator split fashion
                  if(sim_usedrag) then
                      vx1 = vx0 + accx*dtN
                      vy1 = vy0 + accy*dtN
                      vz1 = vz0 + accz*dtN
                  else
                      vx1 = vx0
                      vy1 = vy0
                      vz1 = vz0
                  endif
                  ! transform to frame of the plasma
                  vx1 = vx1 - vgx
                  vy1 = vy1 - vgy
                  vz1 = vz1 - vgz

                  vxp = vx1 + (vy1*Bz/Btot - vz1*By/Btot)*tanom
                  vyp = vy1 + (vz1*Bx/Btot - vx1*Bz/Btot)*tanom
                  vzp = vz1 + (vx1*By/Btot - vy1*Bx/Btot)*tanom
                  vx2 = vx1 + (vyp*Bz/Btot - vzp*By/Btot)*2.*tanom/ &
                      (1. + tanom**2)
                  vy2 = vy1 + (vzp*Bx/Btot - vxp*Bz/Btot)*2.*tanom/ &
                      (1. + tanom**2)
                  vz2 = vz1 + (vxp*By/Btot - vyp*Bx/Btot)*2.*tanom/ &
                      (1. + tanom**2)
                  ! transform back to the lab frame
                  vx2 = vx2 + vgx
                  vy2 = vy2 + vgy
                  vz2 = vz2 + vgz
                  x2 = x1 + dtN/2.*(vx2 + vx0)
                  y2 = y1 + dtN/2.*(vy2 + vy0)
                  z2 = z1 + dtN/2.*(vz2 + vz0)
                  ! transform back to cylindrical coordinates
                  r2 = sqrt(x2**2 + y2**2)
                  p2 = atan2(y2,x2)
                  vr2 = vx2*cos(p2) + vy2*sin(p2)
                  vp2 = (vy2*cos(p2) - vx2*sin(p2))/r2
                  particles(POSX_PART_PROP,i) = r2
                  particles(VELX_PART_PROP,i) = vr2
                  particles(POSY_PART_PROP,i) = z2
                  particles(VELY_PART_PROP,i) = vz2
                  particles(POSZ_PART_PROP,i) = p2
                  particles(VELZ_PART_PROP,i) = vp2
                  if(particles(POSZ_PART_PROP,i) < -PI) &
                      particles(POSZ_PART_PROP,i) = 2.*PI +  &
                      particles(POSZ_PART_PROP,i)
                  if(particles(POSZ_PART_PROP,i) > PI) &
                      particles(POSZ_PART_PROP,i) = &
                      particles(POSZ_PART_PROP,i) - 2.*PI
                  vtot = sqrt(particles(VELX_PART_PROP,i)**2 + &
                      particles(VELY_PART_PROP,i)**2 + &
                      (particles(VELZ_PART_PROP,i)* &
                      particles(POSX_PART_PROP,i))**2)
              else
                  !!! Cartesian coordinates
                  x1 = particles(POSX_PART_PROP,i)
                  y1 = particles(POSY_PART_PROP,i)
                  z1 = particles(POSZ_PART_PROP,i)
                  vx0 = particles(VELX_PART_PROP,i)
                  vy0 = particles(VELY_PART_PROP,i)
                  vz0 = particles(VELZ_PART_PROP,i)

                  vgx = particles(GVLX_PART_PROP,i)
                  vgy = particles(GVLY_PART_PROP,i)
                  vgz = particles(GVLZ_PART_PROP,i)
                  Btot = sqrt(particles(MAGX_PART_PROP,i)**2 + &
                      particles(MAGY_PART_PROP,i)**2 + &
                      particles(MAGZ_PART_PROP,i)**2)
                  omega = particles(CHRG_PART_PROP,i)*ec*Btot/ &
                      particles(MASS_PART_PROP,i)
                  tanom = tan(omega*dtN/2.)
                  Bx = particles(MAGX_PART_PROP,i)
                  By = particles(MAGY_PART_PROP,i)
                  Bz = particles(MAGZ_PART_PROP,i)
                  ! This is now just drag
                  accx = particles(ACCX_PART_PROP,i)
                  accy = particles(ACCY_PART_PROP,i)
                  accz = particles(ACCZ_PART_PROP,i)
                  ! apply drag in operator split fashion
                  if(sim_usedrag) then
                      vx1 = vx0 + accx*dtN
                      vy1 = vy0 + accy*dtN
                      vz1 = vz0 + accz*dtN
                  else
                      vx1 = vx0
                      vy1 = vy0
                      vz1 = vz0
                  endif
                  ! transform to frame of the plasma
                  vx1 = vx1 - vgx
                  vy1 = vy1 - vgy
                  vz1 = vz1 - vgz

                  vxp = vx1 + (vy1*Bz/Btot - vz1*By/Btot)*tanom
                  vyp = vy1 + (vz1*Bx/Btot - vx1*Bz/Btot)*tanom
                  vzp = vz1 + (vx1*By/Btot - vy1*Bx/Btot)*tanom
                  vx2 = vx1 + (vyp*Bz/Btot - vzp*By/Btot)*2.*tanom/ &
                      (1. + tanom**2)
                  vy2 = vy1 + (vzp*Bx/Btot - vxp*Bz/Btot)*2.*tanom/ &
                      (1. + tanom**2)
                  vz2 = vz1 + (vxp*By/Btot - vyp*Bx/Btot)*2.*tanom/ &
                      (1. + tanom**2)
                  ! transform back to the lab frame
                  vx2 = vx2 + vgx
                  vy2 = vy2 + vgy
                  vz2 = vz2 + vgz
                  x2 = x1 + dtN/2.*(vx2 + vx0)
                  y2 = y1 + dtN/2.*(vy2 + vy0)
                  z2 = z1 + dtN/2.*(vz2 + vz0)

                  particles(POSX_PART_PROP,i) = x2
                  particles(VELX_PART_PROP,i) = vx2
                  particles(POSY_PART_PROP,i) = y2
                  particles(VELY_PART_PROP,i) = vy2
                  particles(POSZ_PART_PROP,i) = z2
                  particles(VELZ_PART_PROP,i) = vz2

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

end subroutine pt_advanceDust
