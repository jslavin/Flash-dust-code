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
       cl_dene, cl_mubar, cl_logXT, cl_dlogLdlogT, cl_dnedlogT, &
       cl_dlogGdlogT, cl_dlogTdlogXT, cl_N, cl_smalle, &
       cl_Boltzmann, cl_nnH, cl_gamma

  implicit none
  
  real*8  :: dlogTinv
  integer :: i

!==============================================================================

    cl_fctn_file = '/export/slavin/FLASH/CLIC/SD_heatcoolcurve.txt'
    cl_Boltzmann = 1.38065E-16
    cl_smalle = 1.00E+9
    cl_mubar = 2.1716E-24
    cl_nnH = 1.099
    cl_gamma = 1.66666667

    call readCloudyTable(cl_fctn_file(1:len_trim(cl_fctn_file)), &
          cl_logT, cl_logL, cl_logG, cl_dene, cl_Nmax, CL_N)
     
    cl_logXT(1) = log10(1. + cl_dene(1)) + cl_logT(1)
    do i = 1, cl_N-1
        cl_logXT(i+1) = log10(1. + cl_dene(i+1)) + cl_logT(i+1)
        dlogTinv      = 1. / (cl_logT(i+1) - cl_logT(i))
        cl_dlogLdlogT(i) = (cl_logL(i+1) - cl_logL(i)) * dlogTinv
        cl_dlogGdlogT(i) = (cl_logG(i+1) - cl_logG(i)) * dlogTinv
        cl_dnedlogT(i)   = (cl_dene(i+1) - cl_dene(i)) * dlogTinv
        cl_dlogTdlogXT(i)   = (cl_logT(i+1) - cl_logT(i))/ &
            (cl_logXT(i+1) - cl_logXT(i))
        !write(*,*) 'i =',i,' cl_logT =',cl_logT(i),' cl_dene =',cl_dene(i)
        write(*,fmt='(i3,6(1PG12.5))') i,cl_logT(i),cl_dene(i), &
            cl_dlogLdlogT(i),cl_dlogGdlogT(i),cl_dnedlogT(i),cl_dlogTdlogXT(i)
    enddo
     
  return
end subroutine Cool_init
