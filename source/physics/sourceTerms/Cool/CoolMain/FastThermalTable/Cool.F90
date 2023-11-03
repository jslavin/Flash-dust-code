!!****if* source/physics/sourceTerms/Cool/CoolMain/APEC/Cool
!!
!! NAME
!!
!!  Cool
!!
!! SYNOPSIS
!!
!!  Cool(integer,intent(IN)  :: blockcount,
!!       integer,dimension(blockCount), intent(IN)  :: blocklist,
!!       real,intent(IN)  :: dt,
!!       real,intent(IN)  :: time)
!!
!! DESCRIPTION
!!     Implementation of the integration scheme by R. H. D. Townsend as a means
!!     of calculating the radiative cooling (ApJS 181, 391).
!!     Initial version by John Zuhone
!!
!! ARGUMENTS
!!
!!   blockcount : 
!!
!!   blocklist : 
!!
!!   dt : 
!!
!!   time : 
!!
!! AUTOGENROBODOC
!!
!!
!!***


subroutine Cool (blockCount, blockList, dt, time)

!==============================================================================
  use Cool_data, ONLY: cl_Nmax, cl_Zbar, &
       cl_Nt, cl_Abar, useCool, cl_rho, cl_gamma, cl_smalle, cl_meshMe, &
       cl_boltzmann, cl_AMU, cl_Xin  
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getBlkIndexLimits, &
       Grid_releaseBlkPtr, Grid_getDeltas
  use Eos_interface, ONLY : Eos_wrapped
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
       Multispecies_getSumFrac

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Flash_mpi.h"
  
  integer,intent(IN)  :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real,intent(IN) :: dt, time

  real,pointer,dimension(:,:,:,:) :: solnData
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: blockID
  real            :: dtmin, dtguess, t, abund, lcentral_cool, central_cool
  real            :: eid, eidold, ei, ek, dedt, rr, logr, de
  integer         :: i, j, k, nok, nbad, lb, ir, ierr
  real            :: cellVolume, cellCool
  real,allocatable, dimension(:) :: xCenter, yCenter, zCenter
  integer,dimension(MDIM)  :: dimSize
  logical :: gcell = .true.
  real,dimension(MDIM) :: del
  
  external cool_advance

!==============================================================================
  if(.not.useCool) return
  
  do lb = 1, blockCount
!!DEV: blockList(blockCount) in the next 2 lines looks really wrong - KW
     blockID=blockList(lb)
     call Grid_getBlkPtr(blockID,solnData)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              cl_rho = solnData(DENS_VAR,i,j,k)
              ! eid = eint * rho = 3/2*P
              eid    = cl_rho * solnData(EINT_VAR,i,j,k)
              eidold = eid
                 
              cl_gamma  = solnData(GAME_VAR,i,j,k) 
 
#if NSPECIES > 0
                  cl_Xin    = solnData(SPECIES_BEGIN:SPECIES_BEGIN+NSPECIES-1,i,j,k)
                  call Multispecies_getSumInv(A,cl_Abar,cl_Xin(:))
                  cl_Abar = 1.e0 / cl_Abar
                  call Multispecies_getSumFrac(Z,cl_Zbar,cl_Xin(:))
                  cl_Zbar = cl_Abar*cl_Zbar
#endif

              call cool_advance(eid, cl_rho, dt)

              solnData(COOL_VAR,i,j,k) = (eid - eidold)/dt

           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockID,solnData)
  end do
  !===============================================================================

  ! Now we update the internal and total energies

  do lb = 1, blockCount
     blockID=blockList(lb)
     call Grid_getBlkPtr(blockID,solnData)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              cl_rho = solnData(DENS_VAR,i,j,k)
              eid    = cl_rho * solnData(EINT_VAR,i,j,k)
              ek     = 0.5*(solnData(VELX_VAR,i,j,k)**2 + &
                   solnData(VELY_VAR,i,j,k)**2 + &
                   solnData(VELZ_VAR,i,j,k)**2)
              de = solnData(COOL_VAR,i,j,k)*dt
              eid = eid + de
              ei = max(eid/cl_rho, cl_smalle)   
              solnData(ENER_VAR,i,j,k) = ei + ek                                                                                               
              solnData(EINT_VAR,i,j,k) = ei                                                                                                    
           enddo
        enddo
     enddo
     ! Make all thermodynamic quantities consistent with each other.                                                                           
     call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
     call Grid_releaseBlkPtr(blockID,solnData)
  enddo

  return
end subroutine Cool

subroutine cool_advance(e, rho, dt)

!===============================================================================

  use Cool_data, ONLY : cl_T, cl_gamma, cl_boltzmann, cl_Abar, cl_chi, cl_XT, &
      cl_dL, cl_XTref, cl_lambda, cl_dene, cl_Y, cl_Nt, cl_Zbar, cl_AMU, &
      cl_smalle, cl_dLref, cl_alpha, cl_m, cl_nnh

  use ut_interpolationInterface, ONLY: ut_hunt
 
  implicit none
  
  real, intent(IN)       :: rho, dt
  real, intent(INOUT)    :: e
    
  integer            :: j, err
  real               :: dens, XT, Yold, Ynew, tcool, dlambda, chi

  real, external :: Yfunc, XTfunc
  logical :: lowT

  !------------------------------------------------------------------------------
    
  if (e < rho*cl_smalle) then
     return
  endif

  ! OK, find our spot in the table.
  ! e here is 1/(gamma - 1)*P = 3/2 * chi*k*n*T, XT = chi*T

  ! XT = cl_Abar*(cl_gamma-1.)*cl_AMU*(e/dens)/cl_boltzmann
  dens = rho/cl_m
  XT = (cl_gamma - 1.)*(e/dens)/cl_boltzmann
  ! If cl_Abar = m in amu then this assumes that chi is 1
  ! If cl_Abar = m/chi (m in amu) then we can assume a different chi - which
  ! is then combined with m to define Abar (i.e. eos_singleSpecisA)

  call ut_hunt(cl_XT, cl_Nt, XT, j)

  ! Compute the cooling rate.
  
     
  if (j > 0 .and. j < cl_Nt) then
     !chi = ((XT - cl_XT(j))*chi(j+1) + (cl_XT(j+1) - XT))/ &
     !    (cl_XT(j+1) - cl_XT(j))
     ! dlambda = (chi - 1)*Lambda
     dlambda = cl_dL(j)*(XT/cl_XT(j))**cl_alpha(j)
     tcool = e/(dens**2/cl_nnh*dlambda)
     Yold = Yfunc(XT, j)
     Ynew = Yold + (XT/cl_XTref)*(cl_dLref/dlambda)*dt/tcool
     XT = XTfunc(Ynew, j)
     e = XT*dens*cl_boltzmann/(cl_gamma - 1.)
  endif
  
  !=============================================================================
  
  return
end subroutine cool_advance

function Yfunc(XT, j)

  use Cool_data, ONLY: cl_XT, cl_dL, cl_alpha, cl_Nt, cl_Y, cl_lambda, &
       cl_dLref, cl_XTref

  implicit none

  real, intent(IN) :: XT
  integer, intent(IN) :: j

  real :: Yfunc

  Yfunc = cl_Y(j) + (1./(1. - cl_alpha(j))) * &
       (cl_dLref/cl_dL(j))*(cl_XT(j)/cl_XTref) * &
       (1. - (cl_XT(j)/XT)**(cl_alpha(j) - 1.))

  return

end function Yfunc

function XTfunc(Y, j)
	
  use Cool_data, ONLY: cl_XT, cl_dL, cl_alpha, cl_Y, cl_XTref, cl_dLref

  implicit none

  real, intent(IN) :: Y
  integer, intent(IN) :: j

  real :: XTfunc
	
  XTfunc = cl_XT(j)*(1. - (1. - cl_alpha(j))* &
      (cl_dL(j)/cl_dLref) * (cl_XTref/cl_XT(j)) * &
      (Y - cl_Y(j)))**(1./(1. - cl_alpha(j)))
		
  return
	
end function XTfunc
