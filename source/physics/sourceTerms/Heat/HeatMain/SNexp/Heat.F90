!!****if* source/physics/sourceTerms/Heat/HeatMain/Neutrino/Heat
!!
!! NAME
!!  
!!  Heat 
!!
!!
!! SYNOPSIS
!! 
!!  call Heat (integer(IN) :: blkcnt,
!!             integer(IN) :: blklst(blkcnt),
!!             real(IN)    :: dt,
!!             real(IN)    :: time)
!!
!!  
!!  
!! DESCRIPTION
!!
!!   Calculates heating due to a Supernova explosion
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***
!
!==============================================================================
!
subroutine Heat(blkcnt, blklst, dt, time)

#include "Flash.h"
#include "constants.h"

  use Heat_data, ONLY : useHeat, newexp
  use Eos_interface, ONLY : Eos_wrapped
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getCellCoords
  use Simulation_data, ONLY: sim_Rcore, sim_rhosm, sim_Rej, sim_vej, &
      sim_ejpl, sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin, &
      sim_nProfile, sim_drProf, sim_rProf, sim_vProf, sim_pProf, &
      sim_rhoProf, sim_gamma, sim_smallRho, sim_smallP, sim_nSubZones, &
      sim_inSubzm1, sim_inszd, xcl, ycl, zcl, clrad, clrho, sim_ncl, &
      pCoreFac, pCorePL, sim_Bx0, sim_By0, sim_Bz0, &
      sim_killdivb, texp, xexp, yexp, zexp, nexps, exploded, &
      xpart,ypart,zpart,vxpart,vypart,vzpart,npart,sim_pISM,sim_SNoutint, &
      sim_pactive, sim_readparts
  use IO_interface, ONLY : IO_writeCheckpoint
  use Driver_data, ONLY : dr_dtInit, dr_nstep
  use IO_data, ONLY : IO_plotFileIntervalTime
  use IOParticles_data, ONLY : IO_particleFileIntervalTime

  implicit none
  
  integer,intent(IN) :: blkcnt
  integer,dimension(blkcnt),intent(IN)::blklst
  integer :: i, j, k, lb, n, jLo, jHi
  integer  ::  ii, jj, kk
  real     ::  distInv, xDist, yDist, zDist
  real     ::  sumRho, sumP, sumVX, sumVY, sumVZ
  real     ::  vel, diagonal
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac
  real     ::  vx, vy, vz, p, rho, e, ek
  real     ::  dist, rclump, vclump
  integer :: istat, blockId, iexp, Ncl
  real, pointer :: blkPtr(:,:,:,:)
  real,intent(IN) :: dt,time

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  real, dimension(:,:,:,:),pointer :: solnData, facexData, faceyData, facezData

  integer :: sizeX,sizeY,sizeZ,iclump
  integer,dimension(MDIM) :: axis
  logical :: gcell = .true., inclump, restart = .false.
  logical :: test_current_position, success, writeout
  real :: specs(4)
  
  if (.not. useHeat) return
  iexp = 0
  do n=1,nexps
      ! exploded variable keeps track of whether the SN has already gone off
      ! texp is time of the explosion in seconds after start of the simulation
      if((.not. exploded(n)) .and. (time .gt. texp(n))) then
          exploded(n) = .true.
          write(*,*) 'Explosion initiated...'
          iexp = n
          exit
      endif
  enddo
  ! Initiate particles for second explosion
  if(iexp == 2) then
      IO_plotFileIntervalTime = sim_SNoutint
      IO_particleFileIntervalTime = sim_SNoutint
      sim_pactive = .true.
  endif
  if(iexp > 0) then
      writeout = .true. ! write out clump/particle data
      if(sim_readparts) then
          call read_clumps(xcl, ycl, zcl, clrad, clrho)
      else
          call gen_clumps(Ncl,xcl,ycl,zcl,clrad,clrho,writeout)
      endif
      call IO_writeCheckpoint()
      ! Not sure if restart should be false or true here
      !call Particles_init(restart)
  !
  !  Construct the radial samples needed for the initialization.
  !
      diagonal = (sim_xMax-sim_xMin)**2
      diagonal = diagonal + K2D*(sim_yMax-sim_yMin)**2
      diagonal = diagonal + K3D*(sim_zMax-sim_zMin)**2
      diagonal = sqrt(diagonal)

      sim_drProf = diagonal / (sim_nProfile-1)
      do i = 1, sim_nProfile
          sim_rProf(i)   = (i-1) * sim_drProf
          if (sim_rProf(i).le.sim_Rcore) then
              sim_rhoProf(i) = sim_rhosm
              sim_vProf(i) = sim_vej*(sim_rProf(i)/sim_Rej)
              sim_pProf(i) = sim_pISM*pCorefac
          else if(sim_rProf(i).le.sim_Rej) then
              sim_rhoProf(i) = sim_rhosm*(sim_Rcore/sim_rProf(i))**sim_ejpl
              sim_vProf(i) = sim_vej*(sim_rProf(i)/sim_Rej)
              sim_pProf(i) = sim_pISM*pCorefac*(sim_Rcore/sim_rProf(i))**pCorePL
          endif
      enddo
      do lb = 1, blkcnt
         call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
         sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
         allocate(xCoord(sizeX),stat=istat)
         xCoord = 0.0
         sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
         allocate(yCoord(sizeY),stat=istat)
         yCoord = 0.0
         sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
         allocate(zCoord(sizeZ),stat=istat)
         zCoord = 0.0
         blockId = blklst(lb)
         if (NDIM == 3) then
            call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, &
                 zCoord, sizeZ)
         endif
         if (NDIM >= 2) then
            call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, &
                 yCoord, sizeY)
         endif
         call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
         call Grid_getBlkPtr(blockId,solnData)
!#if NFACE_VARS > 0
!         if (sim_killdivb) then
!            call Grid_getBlkPtr(blockID,facexData,FACEX)
!            if (NDIM >= 2) call Grid_getBlkPtr(blockID,faceyData,FACEY)
!            if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
!         endif
!#endif
         do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
            ! Find a real difference between z's if problem is >= 3D
            if (NDIM > 2) then
               if (k .eq. 1) then
                  dzz = zCoord(2) - zCoord(1) 
               else
                  dzz = zCoord(k) - zCoord(k-1) 
               endif
            ! Otherwise this problem is <= 2D, so dzz is meaningless
            else
               dzz = 0.0
            endif
            zz = zCoord(k)
            do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
               if (NDIM > 1) then
                  if (j .eq. 1) then
                     dyy = yCoord(2) - yCoord(1) 
                  else
                     dyy = yCoord(j) - yCoord(j-1) 
                  endif
                  ! Otherwise this problem is <= 1D, so dyy is meaningless
               else
                  dyy = 0.0
               endif
               yy = yCoord(j)
               do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                  xx = xCoord(i)
                  if (i .eq. 1) then
                     dxx = xCoord(2) - xCoord(1) 
                  else
                     dxx = xCoord(i) - xCoord(i-1) 
                  endif
           
                  sumRho = 0.
                  sumP   = 0.
                  sumVX  = 0.
                  sumVY  = 0.
                  sumVZ  = 0.
                 !
                 ! Break the cell into sim_nSubZones^NDIM sub-zones, and look up
                 ! the appropriate quantities along the 1d profile for that
                 ! subzone.  
                 !
                 ! Have the final values for the zone be equal to the average of
                 ! the subzone values.
                 ! 
                  do kk = 0, (sim_nSubZones-1)*K3D
                     zz = zCoord(k) + (kk*sim_inSubzm1-.5)*dzz 
                     zDist = (zz - zexp(iexp)) * K3D
              
                     do jj = 0, (sim_nSubZones-1)*K2D
                        yy = yCoord(j) + (jj*sim_inSubzm1-.5)*dyy
                        yDist = (yy - yexp(iexp)) * K2D
                 
                        do ii = 0, (sim_nSubZones-1)
                           xx = xCoord(i) + (ii*sim_inSubzm1-.5)*dxx
                           xDist = xx - xexp(iexp)
                    
                           dist = sqrt(xDist**2 + yDist**2 + zDist**2)
                           distInv = 1. / max( dist, 1.E-10 )
                           call sim_find(sim_rProf, sim_nProfile, dist, jLo)
                    !
                    !  a point at `dist' is frac-way between jLo and jHi.   We
                    !  do a linear interpolation of the quantities at jLo and
                    !  jHi and sum those.
                    ! 
                           if (jLo .eq. 0) then
                              jLo = 1
                              jHi = 1
                              frac = 0.
                           else if (jLo .eq. sim_nProfile) then
                              jLo = sim_nProfile
                              jHi = sim_nProfile
                              frac = 0.
                           else
                              jHi = jLo + 1
                              frac = (dist - sim_rProf(jLo)) / & 
                                    (sim_rProf(jHi)-sim_rProf(jLo))
                           endif
                    ! 
                    ! Now total these quantities.  Note that v is a radial
                    ! velocity; we multiply by the tangents of the appropriate
                    ! angles to get the projections in the x, y and z
                    ! directions.
                    !
                           if(dist .lt. sim_Rej) then
                               sumP = sumP +  sim_pProf(jLo) + &
                                 frac*(sim_pProf(jHi) - sim_pProf(jLo))
                               sumRho = sumRho + sim_rhoProf(jLo) + &
                                 frac*(sim_rhoProf(jHi)- sim_rhoProf(jLo))
                               vel = sim_vProf(jLo) + frac*(sim_vProf(jHi) - &
                                 sim_vProf(jLo))
                               sumVX  = sumVX  + vel*xDist*distInv
                               sumVY  = sumVY  + vel*yDist*distInv
                               sumVZ  = sumVZ  + vel*zDist*distInv
                           else
                               sumP = sumP +  solnData(PRES_VAR,i,j,k)
                               sumRho = sumRho + solnData(DENS_VAR,i,j,k)
                               sumVX  = sumVX  + solnData(VELX_VAR,i,j,k)
                               sumVY  = sumVY  + solnData(VELY_VAR,i,j,k)
                               sumVZ  = sumVZ  + solnData(VELZ_VAR,i,j,k)
                           endif
                        enddo
                     enddo
                  enddo
           
                  rho = max(sumRho * sim_inszd, sim_smallRho)
                  p   = max(sumP   * sim_inszd, sim_smallP)
                  vx  = sumVX  * sim_inszd
                  vy  = sumVY  * sim_inszd
                  vz  = sumVZ  * sim_inszd
!!*********************************************************************
!! Add clumps
!! The function below will return .true. if there is a clump at the 
!! current coordinates and .false. if there isn't 
!!!! Trying offsetting by dxx/2, dyy/2 to fix issue of density clumps not 
!!!! aligning with the correct positions
                  inclump = test_current_position(xx-dxx/2., yy-dyy/2., &
                      zz-dzz/2., sim_ncl, xcl, ycl, zcl, clrad, clrho, iclump)
                  if(inclump) then
                      rho = clrho(iclump)
                      rclump = sqrt(xcl(iclump)**2 + ycl(iclump)**2 + &
                          zcl(iclump)**2)
                      vclump = sim_vej*(rclump/sim_Rej)
                      vx = vclump*xcl(iclump)/rclump
                      vy = vclump*ycl(iclump)/rclump
                      vz = vclump*zcl(iclump)/rclump
                  endif
!!*********************************************************************
                  ek  = 0.5*(vx*vx + vy*vy + vz*vz)
                  e   = p/(sim_gamma-1.)
                  e   = e/rho + ek
                  e   = max(e, sim_smallP)
           
                  solnData(DENS_VAR,i,j,k) = rho
                  solnData(PRES_VAR,i,j,k) = p
                  solnData(ENER_VAR,i,j,k) = e
                  solnData(GAME_VAR,i,j,k) = sim_gamma
                  solnData(GAMC_VAR,i,j,k) = sim_gamma
                  solnData(VELX_VAR,i,j,k) = vx
                  solnData(VELY_VAR,i,j,k) = vy
                  solnData(VELZ_VAR,i,j,k) = vz
!  Try without altering the magnetic field in the ejecta
!#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
!                     solnData(MAGX_VAR,i,j,k) = sim_Bx0 
!                     solnData(MAGY_VAR,i,j,k) = sim_By0
!                     solnData(MAGZ_VAR,i,j,k) = sim_Bz0
!                     solnData(MAGP_VAR,i,j,k) = dot_product( &
!                           solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
!                           solnData(MAGX_VAR:MAGZ_VAR,i,j,k))/(8.*PI)
!                     solnData(DIVB_VAR,i,j,k) = 0.
!#endif
                enddo ! end of blkLimits IAXIS loop
             enddo ! end of blkLimits JAXIS loop
          enddo ! end of blkLimits KAXIS loop
          call Eos_wrapped(MODE_DENS_PRES,blkLimits,blockID)
          call Grid_releaseBlkPtr(blklst(lb), blkPtr)
!#if NFACE_VARS > 0
!          if (sim_killdivb) then
!             call Grid_releaseBlkPtr(blockID,facexData,FACEX)
!             if (NDIM >= 2) call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
!             if (NDIM == 3) call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
!          endif
!#endif

          deallocate(xCoord)
          if(NDIM > 1) deallocate(yCoord)
          if(NDIM > 2) deallocate(zCoord)
      enddo ! endof blkcnt loop
      call IO_writeCheckpoint()
      newexp = .true.
  endif
  return
end subroutine Heat
