!!****if* source/Simulation/SimulationMain/Sedov/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  call Simulation_initBlock(integer(IN) :: blockId)
!!                       
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the simulation of the Cas A
!!  supernova remnant including fast expanding ejecta inside a stellar wind
!!  shell.  The ejecta are assumed to have a core-envelope structure with a
!!  core that has a constant average density and a sharply falling off, as a
!!  steep power law in radius, envelope.  The core, while having a constant (in
!!  radius) density on average, have a smooth, low density component with dense
!!  clumps sprinkled randomly throughout it.  The aim (eventually) is to study
!!  how dust formed in the clumps evolves as the clumps are destroyed when
!!  overrun by the reverse shock.
!!
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  
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
!!  sim_Rej             Radius of the outer edge of the ejecta
!!  sim_Rb              Shock radius in 2004
!!  sim_Ro              Radius of the outer edge of the stellar wind shell
!!  sim_pISM            Initial ambient (ISM) pressure
!!  sim_rhoISM          Initial ambient density
!!  sim_xctr            Explosion center coordinates
!!  sim_yctr            Explosion center coordinates
!!  sim_zctr            Explosion center coordinates
!!  sim_nsubzones       Number of `sub-zones' in cells for applying 1d profile
!!
!!
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY: sim_xMax, sim_xMin, sim_yMax, sim_yMin, &
  sim_zMax, sim_zMin, sim_nProfile, sim_drProf, sim_rProf, sim_vProf, &
  sim_pProf, sim_rhoProf, sim_gamma, sim_rhoISM, sim_pISM, sim_rhosm, &
  sim_vej, sim_Rej, pCorefac, pCorePL, sim_Rcore, sim_smallX, sim_smallRho, &
  sim_ejpl, sim_smallP, sim_rInit, sim_nSubZones, sim_xCenter, sim_yCenter, &
  sim_zCenter, sim_inSubzm1, sim_inszd, sim_threadBlockList, &
  sim_threadWithinBlock, sim_vwind, sim_MLR, sim_Bx0, sim_By0, sim_Bz0, &
  sim_killdivb, sim_meshMe, sim_smalle, sim_ncl, xcl, ycl, zcl, clrad, clrho, &
  sim_fcl, sim_chi

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
  Grid_releaseBlkPtr, Grid_getCellCoords, Grid_putPointData
  
  implicit none

#include "constants.h"
#include "Flash.h"
  
  integer,intent(IN) ::  blockId
  
  integer  ::  i, j, k, n, jLo, jHi
  integer  ::  ii, jj, kk
  real     ::  distInv, xDist, yDist, zDist
  real     ::  sumRho, sumP, sumVX, sumVY, sumVZ
  real     ::  vel, diagonal
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac
  real     ::  vx, vy, vz, p, rho, e, ek
  real     ::  dist
  logical  ::  validGeom, inclump, test_current_position, add_clumps
  integer :: istat
  real :: magp, rhomean

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ,iclump
  integer,dimension(MDIM) :: axis
  real, dimension(:,:,:,:),pointer :: solnData, facexData, faceyData, facezData

  logical :: gcell = .true.
!!*********************************************************************
!! Variables for clumps
  real :: rho_out, vclump, rclump

  !write(*,*) 'In Simulation_initBlock, blockId =',blockId
  add_clumps = .true.
  !
  !  Construct the radial samples needed for the initialization.
  !
  diagonal = (sim_xMax-sim_xMin)**2
  diagonal = diagonal + K2D*(sim_yMax-sim_yMin)**2
  diagonal = diagonal + K3D*(sim_zMax-sim_zMin)**2
  diagonal = sqrt(diagonal)
  
  sim_drProf = diagonal / (sim_nProfile-1)
  
  rhomean = ((1. - sim_fcl) + sim_chi*sim_fcl)*sim_rhosm
  do i = 1, sim_nProfile
     sim_rProf(i) = (i-1) * sim_drProf
     if (sim_rProf(i).le.sim_Rcore) then
         if(add_clumps) then
             sim_rhoProf(i) = sim_rhosm ! use sim_rhosm with clumps 
         else
             sim_rhoProf(i) = rhomean ! use rhomean without clumps
         endif
         sim_vProf(i) = sim_vej*(sim_rProf(i)/sim_Rej)
         sim_pProf(i) = sim_pISM*pCorefac
     else if(sim_rProf(i).le.sim_Rej) then
         if(add_clumps) then
             sim_rhoProf(i) = sim_rhosm*(sim_Rcore/sim_rProf(i))**sim_ejpl
         else
             sim_rhoProf(i) = rhomean*(sim_Rcore/sim_rProf(i))**sim_ejpl
         endif
         sim_vProf(i) = sim_vej*(sim_rProf(i)/sim_Rej)
         sim_pProf(i) = sim_pISM*pCorefac*(sim_Rcore/sim_rProf(i))**pCorePL
     else
         sim_rhoProf(i) = sim_rhoISM
         sim_vProf(i) = 0.
         sim_pProf(i) = sim_pISM
     endif
  enddo
     
  ! get the coordinate information for the current block from the database

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat); xCoord = 0.0
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat); yCoord = 0.0
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat); zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
  !
  !     For each cell
  !  
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_getBlkPtr(blockId,solnData)
#endif

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     if (NDIM >= 2) call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  !There is no parallel region in Grid_initDomain and so we use the
  !same thread within block code for both multithreading strategies.

  !$omp parallel if (sim_threadBlockList .or. sim_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(blkLimitsGC,xCoord,yCoord,zCoord,blockID,&
  !$omp sim_inSubzm1,sim_nSubZones,sim_rProf,sim_smallRho,sim_smallP,&
  !$omp sim_smallX,sim_pProf,sim_rhoProf,sim_vProf,sim_gamma,sim_inszd,&
  !$omp sim_xCenter,sim_yCenter,sim_zCenter) &
  !$omp private(i,j,k,ii,jj,kk,n,dxx,dyy,dzz,sumRho,sumP,sumVX,sumVY,sumVZ,&
  !$omp xx,yy,zz,xDist,yDist,zDist,dist,distInv,jLo,jHi,frac,vel,axis,&
  !$omp rho,p,vx,vy,vz,ek,e)

#if NDIM == 3
  !$omp do schedule(static)
#endif
  !write(*,*) 'sim_pISM =',sim_pISM,' sim_inszd =',sim_inszd, &
  !    ' sim_smallP =',sim_smallP
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
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

#if NDIM == 2
     !$omp do schedule(static)
#endif
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        ! Find a real difference between y's if problem is >= 2D
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
        
#if NDIM == 1
        !$omp do schedule(static)
#endif
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
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
           ! the appropriate quantities along the 1d profile for that subzone.  
           !
           ! Have the final values for the zone be equal to the average of
           ! the subzone values.
           ! 
           do kk = 0, (sim_nSubZones-1)*K3D
              zz    = zCoord(k) + (kk*sim_inSubzm1-.5)*dzz 
              zDist = (zz - sim_zCenter) * K3D
              
              do jj = 0, (sim_nSubZones-1)*K2D
                 yy    = yCoord(j) + (jj*sim_inSubzm1-.5)*dyy
                 yDist = (yy - sim_yCenter) * K2D
                 
                 do ii = 0, (sim_nSubZones-1)
                    xx    = xCoord(i) + (ii*sim_inSubzm1-.5)*dxx
                    xDist = xx - sim_xCenter
                    
                    dist    = sqrt( xDist**2 + yDist**2 + zDist**2 )
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
                    sumP = sumP +  sim_pProf(jLo) + frac*(sim_pProf(jHi) - &
                            sim_pProf(jLo))
                    
                    sumRho = sumRho + sim_rhoProf(jLo) + &
                        frac*(sim_rhoProf(jHi)- sim_rhoProf(jLo))
                    
                    vel = sim_vProf(jLo) + frac*(sim_vProf(jHi) - &
                        sim_vProf(jLo))
                    
                    sumVX  = sumVX  + vel*xDist*distInv
                    sumVY  = sumVY  + vel*yDist*distInv
                    sumVZ  = sumVZ  + vel*zDist*distInv
                    
                 enddo
              enddo
           enddo
           
           rho = max(sumRho * sim_inszd, sim_smallRho)
           p  = max((sumP   * sim_inszd), sim_smallP)
           vx = sumVX  * sim_inszd
           vy = sumVY  * sim_inszd
           vz = sumVZ  * sim_inszd
!!*********************************************************************
!! Add clumps
!! The function below will return .true. if there is a clump at the 
!! current coordinates and .false. if there isn't 
!!!! Trying offsetting by dxx/2, dyy/2 to fix issue of density clumps not 
!!!! aligning with the correct positions
           if(add_clumps) then
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
           endif
           ek = 0.5*(vx*vx + vy*vy + vz*vz)
           !
           !  assume gamma-law equation of state
           !
           e = p/(sim_gamma-1.)
           e = e/rho + ek
           e = max(e, sim_smallE)
           
           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k


#ifdef FL_NON_PERMANENT_GUARDCELLS
           if (NSPECIES > 0) then
              solnData(SPECIES_BEGIN,i,j,k)=1.0-(NSPECIES-1)*sim_smallX
              solnData(SPECIES_BEGIN+1:SPECIES_END,i,j,k)=sim_smallX
           end if
           solnData(DENS_VAR,i,j,k) = rho
           solnData(PRES_VAR,i,j,k) = p
           solnData(ENER_VAR,i,j,k) = e
           solnData(GAME_VAR,i,j,k) = sim_gamma
           solnData(GAMC_VAR,i,j,k) = sim_gamma
           solnData(VELX_VAR,i,j,k) = vx
           solnData(VELY_VAR,i,j,k) = vy
           solnData(VELZ_VAR,i,j,k) = vz
#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
           solnData(MAGX_VAR,i,j,k) = sim_Bx0/sqrt(4.*PI)
           solnData(MAGY_VAR,i,j,k) = sim_By0/sqrt(4.*PI)
           solnData(MAGZ_VAR,i,j,k) = sim_Bz0/sqrt(4.*PI)
           magp = dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k), &
               solnData(MAGX_VAR:MAGZ_VAR,i,j,k))/2. 
           solnData(MAGP_VAR,i,j,k) = magp
           solnData(DIVB_VAR,i,j,k) = 0.
#endif

#else
           if (NSPECIES > 0) then
              ! putting in the value of the default species
              call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                   axis, 1.0e0-(NSPECIES-1)*sim_smallX)


              !if there is only one species, this loop will not execute
              do n = SPECIES_BEGIN+1, SPECIES_END

                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                      axis, sim_smallX)

              enddo
           end if


           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)    
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, &
               sim_gamma)
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, &
               sim_gamma)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
           call Grid_putPointData(blockId, CENTER, MAGX_VAR, EXTERIOR, axis, &
               sim_Bx0/sqrt(4.*PI))
           call Grid_putPointData(blockId, CENTER, MAGY_VAR, EXTERIOR, axis, &
               sim_By0/sqrt(4.*PI))
           call Grid_putPointData(blockId, CENTER, MAGZ_VAR, EXTERIOR, axis, &
               sim_Bz0/sqrt(4.*PI))
           call Grid_putPointData(blockId, CENTER, MAGP_VAR, EXTERIOR, axis, &
               magp)
           call Grid_putPointData(blockId, CENTER, DIVB_VAR, EXTERIOR, axis, 0.)
#endif
        enddo
#if NDIM == 1
  !$omp end do nowait
#endif
     enddo
#if NDIM == 2
  !$omp end do nowait
#endif
  enddo
#if NDIM == 3
  !$omp end do nowait
#endif
  !$omp end parallel

#if NFACE_VARS > 0
  if (sim_killdivb) then
     facexData(:,:,:,:) = sim_Bx0/sqrt(4.*PI)
     if (NDIM >= 2) faceyData(:,:,:,:) = sim_By0/sqrt(4.*PI)
     if (NDIM == 3) facezData(:,:,:,:) = sim_Bz0/sqrt(4.*PI)

  endif
#endif

#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_releaseBlkPtr(blockID, solnData)
#endif
#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     if (NDIM >= 2) call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  return
end subroutine Simulation_initBlock

!******************************************************************************
!  Routine:     sim_find()
!
!  Description: Given a monotonically increasing table x(N) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

subroutine sim_find(x, N, x0, i)
  implicit none
  integer, intent(IN) :: N
  integer, intent(OUT):: i
  real, intent(IN)    :: x(N), x0

! local variables
  integer  il, ir, im

  if (x0 .lt. x(1)) then
     i = 0
  else if (x0 .gt. x(N)) then
     i = N
  else
     il = 1
     ir = N

     do while (ir-il > 1)
         im = (il + ir) / 2
         if (x(im) .gt. x0) then
            ir = im
         else
            il = im
         endif
     enddo
     i = il

!10   if (ir .eq. il+1) goto 20
!     im = (il + ir) / 2
!     if (x(im) .gt. x0) then
!        ir = im
!     else
!        il = im
!     endif
!     goto 10
!20   i = il
  endif
  return
end subroutine sim_find
