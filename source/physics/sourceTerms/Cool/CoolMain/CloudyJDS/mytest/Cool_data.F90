module Cool_data

!==============================================================================

  implicit none

! Data tables containing cooling curve(s)

  character(len=80), save :: cl_fctn_file
  integer, parameter :: cl_Nmax = 200
  real*8, dimension(cl_Nmax), save :: cl_logT, cl_logL, cl_logG, &
       cl_dene, cl_dlogLdlogT, cl_dlogGdlogT, cl_dnedlogT, cl_logXT, &
       cl_dlogTdlogXT
  integer, save :: cl_N

  real*8, save :: cl_gamma, cl_smalle, cl_Boltzmann, cl_dtmin, cl_smallt, &
      cl_rho, cl_nnH, cl_mubar

!==============================================================================

end module Cool_data
