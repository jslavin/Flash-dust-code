!!****if* source/physics/sourceTerms/Cool/CoolMain/CloudyJDS/Cool_init
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

  use Cool_data, ONLY :cl_fctn_file, cl_Nmax, cl_logT, cl_logL, cl_logG, &
       cl_dene, cl_m, cl_logXT, cl_dlogLdlogT, cl_dnedlogT, &
       cl_dlogGdlogT, cl_dlogTdlogXT, cl_N, cl_Xin, useCool, cl_smalle, &
       cl_Boltzmann, cl_meshMe, cl_dtmin, cl_smallt, cl_nnH, cl_useHeat, cl_Tcut
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

#include "constants.h"

  implicit none

  
  real    :: dlogTinv
  integer :: i

!==============================================================================

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,cl_meshMe)

  call RuntimeParameters_get("cl_functionFile", cl_fctn_file)
  call PhysicalConstants_get("Boltzmann", cl_Boltzmann)
  call RuntimeParameters_get("smalle", cl_smalle)
  call RuntimeParameters_get("smallt", cl_smallt)
  call RuntimeParameters_get("cl_m", cl_m)
  call RuntimeParameters_get("cl_nnH", cl_nnH)
  call RuntimeParameters_get("cl_Tcut", cl_Tcut)
  call RuntimeParameters_get("useCool",useCool)  
  call RuntimeParameters_get("cl_useHeat",cl_useHeat)  
  call RuntimeParameters_get("dtcoolmin",cl_dtmin)  
  if (.not. useCool) then
     write(6,*)'WARNING:  You have included the Cool unit but have set '
     write(6,*)'   the runtime parameter useCool to FALSE'
     write(6,*)'   No cooling will occur but Cool_init will continue.'
  end if

  
  ! Now allocate space for the arrays.
  
  if(useCool) then
     ! if MASTER_PE
     if (cl_meshMe .EQ. MASTER_PE) &
          print *, "Cool_init:  attempting to read cooling function file '", &
          cl_fctn_file(1:len_trim(cl_fctn_file)), "'..."

     call readCloudyTable(cl_fctn_file(1:len_trim(cl_fctn_file)), &
          cl_logT, cl_logL, cl_logG, cl_dene, cl_Nmax, CL_N)
     if (CL_N < 1) then
        print *,cl_meshMe,": Cool_init: could not read '", cl_fctn_file, "'!"  
        ! shows failure not on MASTER_PE
        call Driver_abortFlash("Cool_init: unable to read cooling function file")
     endif
     
    cl_logXT(1) = log10(1. + cl_dene(1)) + cl_logT(1)
     do i = 1, cl_N-1
        cl_logXT(i+1) = log10(1. + cl_dene(i+1)) + cl_logT(i+1)
        dlogTinv      = 1. / (cl_logT(i+1) - cl_logT(i))
        cl_dlogLdlogT(i) = (cl_logL(i+1) - cl_logL(i)) * dlogTinv
        cl_dlogGdlogT(i) = (cl_logG(i+1) - cl_logG(i)) * dlogTinv
        cl_dnedlogT(i)   = (cl_dene(i+1) - cl_dene(i)) * dlogTinv
        cl_dlogTdlogXT(i)   = (cl_logT(i+1) - cl_logT(i))/ &
            (cl_logXT(i+1) - cl_logXT(i))
     enddo
     
     if (cl_meshMe.EQ.MASTER_PE) print *, "Cool_init:  successfully read file."
     
  end if

  return
end subroutine Cool_init
