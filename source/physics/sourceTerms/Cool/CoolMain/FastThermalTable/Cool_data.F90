
module Cool_data

!==============================================================================

  implicit none
#include "Flash.h"

  integer, save :: cl_meshMe, cl_numProcs 

! Data tables containing cooling curve(s)

  character(len=80), save     :: cl_functionFile
  integer, parameter          :: cl_Nmax = 2000

  real, dimension(cl_Nmax), save :: cl_T, cl_lambda, cl_dene, cl_chi, &
      cl_XT, cl_dL, cl_Y, cl_alpha

  integer, save               :: cl_Nt
  real, save, dimension(NSPECIES) :: cl_Xin
  logical, save :: useCool, useChi

  real, save :: cl_rho, cl_gamma, cl_smalle, cl_boltzmann, cl_AMU, &
       cl_Abar, cl_Zbar, cl_m, cl_nnh, cl_Tcut

  real, save :: cl_XTref, cl_dLref
  
!==============================================================================

end module Cool_data
