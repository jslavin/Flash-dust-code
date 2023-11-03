
module Cool_data

!==============================================================================

  implicit none
#include "Flash.h"
! Data tables containing cooling curve(s)

  integer, save :: cl_meshMe, cl_numProcs 

  character(len=80), save :: cl_fctn_file
  integer, parameter :: cl_Nmax = 300
  real, dimension(cl_Nmax), save :: cl_logT, cl_logL, cl_logG, &
       cl_dene, cl_dlogLdlogT, cl_dlogGdlogT, cl_dnedlogT, cl_logXT, &
       cl_dlogTdlogXT
  integer, save :: cl_N
  real, save, dimension(NSPECIES) :: cl_Xin, cl_Abar
  logical, save :: useCool, cl_useHeat

  real, save :: cl_gamma, cl_smalle, cl_Boltzmann, cl_dtmin, cl_smallt, &
      cl_rho, cl_nnH, cl_m, cl_Tcut

!==============================================================================

end module Cool_data
