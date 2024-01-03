!!****if* source/Simulation/SimulationMain/Sedov/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sedov Spherical Explosion 
!!  problem.
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  sim_Eej             Explosion energy (erg)
!!  sim_Mej             Ejecta mass (Msun)
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
!! For runs with dust:
!!  sim_
!!
!!***

subroutine Simulation_init()

  use Simulation_data 
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface,           ONLY : Logfile_stamp
  use Eos_data, ONLY : eos_smalle,eos_smallt
  use Driver_interface, ONLY: Driver_getMype
  use Particles_data, ONLY: useParticles


  implicit none

#include "constants.h"
#include "Flash.h"

  real :: Msun, Msw, dfac, u, ufnc, zbrent
  external ufnc
  logical :: threadBlockListBuild, threadWithinBlockBuild  
!!*****************************************************************************
!! Variables for clumps
  real :: Vcore,Vcl,h,xx,yy,rsph,dist,Vclump,rpos_clump,theta_clump
  real :: rp, thetap, phip, rpart, vpart
  real :: cmpc = 3.08568E18
  integer :: NN, Ncl, nin, i, j, np
  logical :: ovlp, writeout
  real, allocatable :: x_try(:), y_try(:), z_try(:)
!!*****************************************************************************

  call Driver_getMype(GLOBAL_COMM, sim_meshMe)

  call RuntimeParameters_get('sim_Eej', sim_Eej)
  call RuntimeParameters_get('sim_Mej', sim_Mej)
  call RuntimeParameters_get('sim_vej', sim_vej)
  call RuntimeParameters_get('sim_Rej', sim_Rej)
  call RuntimeParameters_get('sim_chi', sim_chi)
  call RuntimeParameters_get('sim_fcl', sim_fcl)
  if (useParticles) then
      ! sim_nperclump is number of particles per clump
      call RuntimeParameters_get('sim_nperclump', sim_nperclump)
      call RuntimeParameters_get('sim_pmass', sim_pmass)
      call RuntimeParameters_get('sim_pdens', sim_pdens)
      call RuntimeParameters_get('sim_ptype', sim_ptype)
      call RuntimeParameters_get('sim_dustdata', sim_dustdata)
      call RuntimeParameters_get('sim_yielddata', sim_yielddata)
      call RuntimeParameters_get('sim_yieldintdata', sim_yieldintdata)
      call RuntimeParameters_get('sim_usetsput', sim_usetsput)
      call RuntimeParameters_get('sim_useisput', sim_useisput)
      call RuntimeParameters_get('sim_usedrag', sim_usedrag)
      call RuntimeParameters_get('sim_chargedata', sim_chargedata)
      call RuntimeParameters_get('sim_G0', sim_G0)
      call RuntimeParameters_get('sim_partdata', sim_partdata)
      call RuntimeParameters_get('sim_readparts', sim_readparts)
  endif
  call RuntimeParameters_get('sim_clumpdata', sim_clumpdata)
  call RuntimeParameters_get('sim_ejpl', sim_ejpl)
  call RuntimeParameters_get('sim_rcl', sim_rcl)
  call RuntimeParameters_get('sim_rhoCSM', sim_rhoCSM)
  call RuntimeParameters_get('sim_Rb', sim_Rb)
  !call RuntimeParameters_get('sim_Mprog', sim_Mprog)
  call RuntimeParameters_get('sim_Msh', sim_Msh)
  call RuntimeParameters_get('sim_exptimesfile', sim_exptimesfile)
  call RuntimeParameters_get('sim_pISM', sim_pISM)
  call RuntimeParameters_get('sim_rhoISM', sim_rhoISM)
  call RuntimeParameters_get('sim_vwind', sim_vwind)
  call RuntimeParameters_get('sim_MLR', sim_MLR)
  call RuntimeParameters_get('sim_rInit', sim_rInit)
  call RuntimeParameters_get('sim_rInject', sim_rInject)

  call RuntimeParameters_get('sim_mubar', sim_mubar)
  call RuntimeParameters_get('sim_nnh', sim_nnh)
  call RuntimeParameters_get('sim_nsubzones',sim_nSubZones)
  call RuntimeParameters_get('sim_xctr',sim_xCenter)
  call RuntimeParameters_get('sim_yctr',sim_yCenter)
  call RuntimeParameters_get('sim_zctr',sim_zCenter)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smlrho', sim_smallRho)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallt', sim_smallT)
  call RuntimeParameters_get('smalle', sim_smallE)
  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmax',sim_zMax)
  ! output time interval (s) after the SN has gone off:
  call RuntimeParameters_get('sim_SNoutint',sim_SNoutint)
#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
  call RuntimeParameters_get('killdivb',sim_killdivb)
  call RuntimeParameters_get('Bx0', sim_Bx0)
  call RuntimeParameters_get('By0', sim_By0)
  call RuntimeParameters_get('Bz0', sim_Bz0)
#endif

  if (sim_nSubZones .le. 1) sim_nSubZones = 2
  
  sim_inSubZones = 1./real(sim_nSubZones)
  sim_inSubzm1   = 1./real(sim_nSubZones-1)
  sim_inszd      = sim_inSubZones**NDIM
  eos_smallt = sim_smallT
  eos_smalle = sim_smallE
  !
  ! Calculate values related to density profile:
  !     sim_Rej - outer radius of ejecta and inner edge of stellar wind shell
  !     sim_eta - ratio of mean density to smooth ejecta density in the core
  !             region of the ejecta
  !     sim_ufac - quantity used in calculating u = ratio of radius ejecta
  !             core to ejecta envelope
  !
  ! Model has flat core density made up of a smooth ejecta background with
  ! clumps placed randomly within the core.  Outside of the core is the ejecta
  ! envelope which drops off as r^(-sim_ejpl). I assume that there are no clumps
  ! in the envelope and that it matches up smoothly with the smooth ejecta in
  ! the core. The velocity of ejecta material is assumed to go linearly with
  ! radial distance.  The explosion energy, ejecta mass and maximum ejecta
  ! velocity (v at sim_Rej) constrain the ratio of core radius to envelope
  ! radius. The radius of the inner edge of the stellar wind shell is
  ! constrained by the estimated amount of shocked CSM, sim_Msh

  Msun = 1.989E33
  ! NOTE: Provide values of Mprog and Mej in Solar masses
  !Msw = (sim_Mprog - sim_Mej)
  ! Inner radius of stellar wind shell, set to outer radius of ejecta
  ! Here we're assuming sim_Rb is in cm
  !!! Now using input value for sim_Rej - disconnecting from Cas A !!!
  !sim_Rej = sim_Rb - sim_Msh*Msun/(4.*PI*sim_Rb**2*sim_rhoCSM)

  sim_eta = sim_chi*sim_fcl + 1. - sim_fcl
  ! NOTE: Here we're assuming the sim_Eej is in erg and sim_vej is in cm/s
  sim_ufac = 10./3.*sim_Eej/(sim_Mej*Msun*sim_vej**2)
  u = zbrent(ufnc,0.1,0.99,1.E-5)
  sim_Rcore = u*sim_Rej
  sim_rhosm = Msun*sim_Mej*3./(4.*PI*sim_Rcore**3*(sim_eta + 3.* &
      (1. - u**(sim_ejpl - 3.))/(sim_ejpl - 3.)))
  sim_rhocl = sim_chi*sim_rhosm
  ! Factors to increase the thermal pressure in the ejecta - idea is to
  ! raise the temperature to values >~ 10 K or so
  ! Here this factor is hard coded - should really be a runtime parameter
  pCorefac = 100.
  ! goal is to make clumps and interclump medium hotter
  pCorePL = log10(pCorefac)/log10(sim_Rej/sim_Rcore)
  !!**************************************************************************
  ! Start clump creation
  Vcore = 4./3.*PI*sim_Rcore**3
  Vcl = sim_fcl*Vcore
  ! estimate (high) the no. of clouds needed to achieve the filling factor
  ! assumes cylindrical symmetry so volume per cloud is 2*pi**2*rcl**2*R, where
  ! R is the cylindrical radial distance of the cloud center.  Average value 
  ! for R for randomly placed clouds within Rcore is 4/(3*pi)*Rcore
  Ncl = int(2.*Vcl/(8./3.*PI*sim_rcl**2*sim_Rcore))
  if(sim_meshMe == MASTER_PE) then
      write(*,fmt='("chi =",ES12.5," fcl =",F6.3)') sim_chi, sim_fcl
      write(*,fmt='("Rej =",ES12.5," eta =",F6.3)') sim_Rej, sim_eta
      write(*,fmt='("ufac =",ES12.5," ejpl =",I3)') sim_ufac, sim_ejpl
      write(*,fmt='("Rcore =",ES12.5," u =",F8.5," rhosm =",ES12.5)') &
          sim_Rcore, u, sim_rhosm
      write(*,fmt='("pCorePL =",ES12.5," Vcl =",ES12.5," Ncl =",I4)') &
          pCorePL,Vcl,Ncl
      call pt_yieldInit(trim(sim_dustdata),trim(sim_yielddata), &
        trim(sim_yieldintdata))
      call pt_chargeInit(trim(sim_chargedata))
  endif
  !! Read in explosion data
  exploded(:) = .false.
  open(99,file=trim(sim_exptimesfile),status='old')
  read(99,*) ! First line is a comment
  read(99,fmt='(i3)') nexps ! number of explosions
  !write(*,'("nexps =",i3)') nexps 
  do i=1,nexps
      read(99,*) texp(i),xexp(i),yexp(i),zexp(i)
  enddo
  close(99)
  if(texp(1) == 0.) exploded(1) = .true.
  ! Now not using particles in first explosion
  sim_pactive = .false.
  writeout = .false. ! don't write out clump/particle data yet
  if(sim_readparts) then
      call read_clumps(xcl, ycl, zcl, clrad, clrho)
  else
      call gen_clumps(Ncl,xcl,ycl,zcl,clrad,clrho,writeout)
  endif
  if(useParticles) then
      call pt_yieldInit(trim(sim_dustdata),trim(sim_yielddata), &
        trim(sim_yieldintdata))
      call pt_chargeInit(trim(sim_chargedata))
  endif
end subroutine Simulation_init

real function ufnc(u)
    use Simulation_data, ONLY: sim_ufac,sim_eta,sim_ejpl
    implicit none
    real, intent(in) :: u
    ufnc = (sim_ufac*(sim_eta + 3./(sim_ejpl - 3.)*(1. - u**(sim_ejpl - 3.))) &
        - u**2*(sim_eta + 5./(sim_ejpl - 5.)*(1. - u**(sim_ejpl - 5.))))
    return
end function ufnc 

logical function test_current_position(x, y, z, N, x_pos, y_pos, z_pos, &
        clump_rad, clump_rho, iclump)
    ! return 0 if there is no overlap or the additional density of the region
    ! if there is a clump at the current coordinates.
    ! revised (4/11/2019) - now returns iclump which is equal to the index
    ! of the cloud - note: this is technically a side effect of the function 
    ! and maybe not advisable, but seems fairly benign
    use Simulation_data, ONLY : sim_clmax
    integer :: i,N,iclump
    real :: x,y,z,x_pos(N), y_pos(N), z_pos(N), clump_rad(N), clump_rho(N), &
        distance

    if(N.gt.sim_clmax) then
        write(*,*) 'Problem with value of N =',N
        N = sim_clmax
    endif
    test_current_position = .false.
    do i=1, min(N,sim_clmax)
        distance = sqrt((x - x_pos(i))**2 + (y - y_pos(i))**2 + &
            (z - z_pos(i))**2)
        if (distance <= clump_rad(i)) then
            test_current_position = .true.
            iclump = i
            exit
        endif
    enddo
    return
end function test_current_position

subroutine gen_clumps(Ncl,xcl,ycl,zcl,clrad,clrho,writeout)
    !! Generates clump positions within the smooth ejecta of the core
    !! Also generates particle positions

    use Simulation_data, ONLY: sim_rcl, sim_Rcore, sim_rhocl, sim_nperclump, &
        sim_meshMe,xpart,ypart,zpart,vxpart,vypart,vzpart,sim_ncl,npart, &
        sim_clmax, sim_mubar, sim_fcl, sim_vej, sim_Rej
    use Driver_interface, ONLY: Driver_getMype
    use Particles_data, ONLY: useParticles

    implicit none

#include "constants.h"
#include "Flash.h"
    real, dimension(sim_clmax) :: xcl,ycl,zcl,clrad,clrho
    real :: Vcore,Vcl,h,xx,yy,rsph,dist,Vclump,rpos_clump,theta_clump
    real :: rp, thetap, phip, rpart, vpart
    real :: cmpc = 3.08568E18
    integer :: NN, Ncl, nin, i, j, np
    logical, intent(in) :: writeout
    ! gfortran demands larger rseed for some reason
    !integer, dimension(2) :: rseed
    integer, dimension(33) :: rseed
    logical :: ovlp
    real, allocatable :: x_try(:), y_try(:), z_try(:)
    !data rseed/3,10/

    call Driver_getMype(GLOBAL_COMM, sim_meshMe)
    Vcore = 4./3.*PI*sim_Rcore**3
    Vcl = sim_fcl*Vcore
    ! estimate (high) the no. of clouds needed to achieve the filling factor
    ! assumes cylindrical symmetry so volume per cloud is 2*pi**2*rcl**2*R,
    ! where R is the cylindrical radial distance of the cloud center.  Average
    ! value for R for randomly placed clouds within Rcore is 4/(3*pi)*Rcore
    Ncl = int(2.*Vcl/(8./3.*PI*sim_rcl**2*sim_Rcore))
    allocate(x_try(Ncl))
    allocate(y_try(Ncl))
    allocate(z_try(Ncl))
    nin = 0
    ! gfortran requires 33 seeds for reproducibility
    rseed(1) = 3
    rseed(2) = 10
    call random_seed(PUT=rseed)
    do i=1,Ncl
        call random_number(h)
        ! sim_rcl <= xx <= sim_Rcore - sim_rcl
        !rpos_clump = sim_rcl + (sim_Rcore - 2.*sim_rcl)*sqrt(h)
        ! Try making clumps farther from boundaries (4/25/2022)
        rpos_clump = 1.5*sim_rcl + (sim_Rcore - 2.5*sim_rcl)*sqrt(h)
        ! xx = sim_rcl + (sim_Rcore - 2.*sim_rcl)*h
        call random_number(h)
        ! make sure that distance from z axis is large enough that cloud stays
        ! within the simulation volume
        ! (11/1/2019) changed - multiply sim_r_cl by 2 to prevent problems with
        ! the symmetry axis (also particles leaving the grid)
        ! increase sim_rcl to 2.5 (4/25/2022)
        theta_clump = asin(2.5*sim_rcl/rpos_clump) + (PI - &
            2.*asin(2.5*sim_rcl/rpos_clump))*h
        xx = rpos_clump*sin(theta_clump)
        ! This assumes that we're using full cylindrical volume
        ! -(sim_Rcore - sim_rcl) <= yy <= sim_Rcore - sim_rcl
        !yy = (sim_Rcore - sim_rcl)*(-1. + 2.*h)
        yy = rpos_clump*cos(theta_clump)
        rsph = sqrt(xx**2 + yy**2)
        ! this check should no longer be necessary
        if((rsph + sim_rcl) .lt. sim_Rcore) then
            nin = nin + 1
            x_try(nin) = xx
            y_try(nin) = yy
        endif
    enddo
    !!Allocate the position and radius variables.
    np = 0
    NN = 0
    if((sim_meshMe == MASTER_PE).and.(writeout)) then
        open(unit=30,file='ejecta_clumps.txt',status='unknown',access='append')
        write(30,fmt='("# Clumpy ejecta distribution")')
        write(30,fmt='("# Clump filling factor =",F7.4)') sim_fcl
        write(30,fmt='("#",A7,A10,A10,A10,A10)') 'x','y','z','rcl','rhocl'
        open(unit=31,file='parts_init.txt',status='unknown',access='append')
        write(31,fmt='("# Initial particle data")')
        write(31,fmt='("#",A7,A10,A12,A12)') 'x','y','vx','vy'
    endif
    Vclump = 0.
    do i=1,nin
        xx = x_try(i)
        yy = y_try(i)
        ovlp = .false.
        if(i.gt.1) then
            do j=1,NN
                dist = sqrt((xx - xcl(j))**2 + (yy - ycl(j))**2)
                if(dist .le. 2.*sim_rcl) then
                    ovlp = .true.
                    exit
                endif
            enddo
        endif
        if(.not.ovlp) then
            NN = NN + 1
            xcl(NN) = xx
            ycl(NN) = yy
            zcl(NN) = 0.
            clrad(NN) = sim_rcl
            clrho(NN) = sim_rhocl
            do j=1,sim_nperclump
                call random_number(h)
                if (NDIM == 3) then
                    rp = clrad(NN)*h**(1./3.)
                    call random_number(h)
                    thetap = PI*h
                    call random_number(h)
                    phip = 2.*PI*h
                    np = np + 1
                    xpart(np) = xcl(NN) + rp*sin(thetap)*cos(phip)
                    ypart(np) = ycl(NN) + rp*sin(thetap)*sin(phip)
                    zpart(np) = zcl(NN) + rp*cos(thetap)
                    rpart = sqrt(xpart(np)**2 + ypart(np)**2 + zpart(np)**2)
                    vpart = sim_vej*(rpart/sim_Rej)
                    vxpart(np) = vpart * xpart(np)/rpart
                    vypart(np) = vpart * ypart(np)/rpart
                    vzpart(np) = vpart * zpart(np)/rpart
                else
                    np = np + 1
                    ! factor of 0.95 is so that the particles aren't too 
                    ! close to the edge which, given discretization, could
                    ! lead to them actually not being in a clump
                    ! reduced factor to 0.9
                    rp = 0.9*clrad(NN)*sqrt(h)
                    call random_number(h)
                    thetap = 2.*PI*h
                    xpart(np) = xcl(NN) + rp*cos(thetap)
                    ypart(np) = ycl(NN) + rp*sin(thetap)
                    rpart = sqrt(xpart(np)**2 + ypart(np)**2)
                    vpart = sim_vej*(rpart/sim_Rej)
                    vxpart(np) = vpart * xpart(np)/rpart
                    vypart(np) = vpart * ypart(np)/rpart
                endif
                if((sim_meshMe == MASTER_PE).and.(writeout)) then
                    write(31,fmt='(F10.4,F10.4,ES12.4,ES12.4)') &
                        xpart(np)/cmpc,ypart(np)/cmpc,vxpart(np)/1.E5, &
                        vypart(np)/1.E5
                endif
            enddo
            if((sim_meshMe == MASTER_PE).and.(writeout)) then
                write(30,fmt='(F10.4,F10.4,F10.4,F10.4,F10.2)') xcl(NN)/cmpc, &
                    ycl(NN)/cmpc,zcl(NN)/cmpc,clrad(NN)/cmpc,clrho(NN)/sim_mubar
            endif
            Vclump = Vclump + 2.*PI**2*sim_rcl**2*xx
            if(Vclump .gt. 0.99*Vcl) exit
        endif
    enddo
    sim_ncl = NN
    !write(*,*) 'sim_ncl =',sim_ncl
    npart = np
    if((sim_meshMe == MASTER_PE).and.(writeout)) then
        close(30)
        close(31)
    endif
    deallocate(x_try)
    deallocate(y_try)
    deallocate(z_try)

end subroutine gen_clumps

subroutine read_clumps(xcl, ycl, zcl, clrad, clrho)
    !! Reads in clump positions within the smooth ejecta of the core
    !! Also reads in particle positions

    use Simulation_data, ONLY: sim_nperclump, xpart, ypart, zpart, &
        vxpart, vypart, vzpart, sim_ncl, npart, sim_clmax, sim_mubar, &
        sim_clumpdata, sim_partdata, sim_meshMe

    implicit none

#include "constants.h"
#include "Flash.h"

    real, dimension(sim_clmax) :: xcl, ycl, zcl, clrad, clrho
    real :: cmpc = 3.08568E18
    integer :: i, j, np, NN, stat

    call Driver_getMype(GLOBAL_COMM, sim_meshMe)

    open(unit=30,file=trim(sim_clumpdata),status='old')
    ! First three lines are comments
    read(30,*)
    read(30,*)
    read(30,*)
    open(unit=31,file=trim(sim_partdata),status='old')
    ! First two lines are comments
    read(31,*)
    read(31,*)
    np = 0
    NN = 0
    do i=1,sim_clmax
        read(30,fmt='(F10.4,F10.4,F10.4,F10.4,F10.2)',iostat=stat) xcl(i), &
                ycl(i), zcl(i), clrad(i), clrho(i)
        if(stat == -1) exit
        NN = NN + 1
        xcl(i) = xcl(i)*cmpc
        ycl(i) = ycl(i)*cmpc
        zcl(i) = zcl(i)*cmpc
        clrad(i) = clrad(i)*cmpc
        clrho(i) = clrho(i)*sim_mubar
        do j=1,sim_nperclump
            np = np + 1
            read(31,fmt='(F10.4,F10.4,ES12.4,ES12.4)') &
               xpart(np), ypart(np), vxpart(np), vypart(np)
            xpart(np) = xpart(np)*cmpc
            ypart(np) = ypart(np)*cmpc
            vxpart(np) = vxpart(np)*1.E5
            vypart(np) = vypart(np)*1.E5
        enddo
    enddo
    zpart(:) = 0.
    vzpart(:) = 0.
    close(30)
    close(31)
    sim_ncl = NN
    npart = np
    if(sim_meshMe == MASTER_PE) then
        write(*,'("clumps and particles read in, sim_ncl =",i4,' &
            // '" npart =",i4)') sim_ncl,npart
    endif
end subroutine read_clumps
