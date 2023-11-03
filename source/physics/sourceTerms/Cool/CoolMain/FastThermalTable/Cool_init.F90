!!****if* source/physics/sourceTerms/Cool/CoolMain/APEC/Cool_init
!!
!! NAME
!!
!!  Cool_init
!!
!! SYNOPSIS
!!
!!  Cool_init()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   
!!
!! AUTOGENROBODOC
!!
!!
!!***


subroutine Cool_init()

  use Cool_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Simulation_interface, ONLY : Simulation_mapStrToInt

#include "constants.h"
#include "Flash.h"

  implicit none

  integer :: test_Nt, i

  real, external :: Yfunc

  !==============================================================================

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,cl_meshMe)

  call RuntimeParameters_get("cl_functionFile", cl_functionFile)
  call PhysicalConstants_get("proton mass", cl_AMU)
  call PhysicalConstants_get("Boltzmann", cl_boltzmann)
  call RuntimeParameters_get("smalle", cl_smalle)
  call RuntimeParameters_get("useCool",useCool)  
  call RuntimeParameters_get("useChi",useChi)  
  call RuntimeParameters_get("cl_m",cl_m)
  call RuntimeParameters_get("cl_nnh",cl_nnh)
  call RuntimeParameters_get("cl_Tcut",cl_Tcut)
#if NSPECIES == 0
  call RuntimeParameters_get("eos_singleSpeciesA",cl_Abar)
  call RuntimeParameters_get("eos_singleSpeciesZ",cl_Zbar)
#endif
  
  if (cl_meshMe .eq. MASTER_PE .and. .not. useCool) then
     write(6,*)'WARNING:  You have included the Cool unit but have set '
     write(6,*)'   the runtime parameter useCool to FALSE'
     write(6,*)'   No cooling will occur but Cool_init will continue.'
  end if
  
  ! Now allocate space for the arrays.
  
  if (useCool) then

     if (cl_meshMe .EQ. MASTER_PE) then
        print *, "Cool_init:  attempting to read cooling function file '", &   
        cl_functionFile(1:len_trim(cl_functionFile)), "'..."
     endif
     
     call readThermalTable(test_Nt)

     if (test_Nt < 1) then
        print *,cl_meshMe,":  Cool_init:  could not read '", cl_functionFile, "'!"  
        call Driver_abortFlash("Cool_init:  unable to read cooling function file")
     endif
     
     if (cl_meshMe .EQ. MASTER_PE) print *, "Cool_init:  successfully read file."

     do i = 1, cl_Nt-1
     	cl_alpha(i) = log(cl_dL(i+1)/cl_dL(i))/log(cl_XT(i+1)/cl_XT(i))
     enddo

     cl_dLref = cl_dL(cl_Nt)
     cl_XTref = cl_XT(cl_Nt)

     cl_Y(cl_Nt) = 0.0
     do i = cl_Nt-1, 1, -1
     	cl_Y(i) = cl_Y(i+1) - (1./(1. - cl_alpha(i))) * &
           (cl_dLref/cl_dL(i)*cl_XT(i)/cl_XTref) * &
           (1. - (cl_XT(i)/cl_XT(i+1))**(cl_alpha(i) - 1.))
     enddo
          
  end if

  return
end subroutine Cool_init
