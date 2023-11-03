!*******************************************************************************
!
! Routine:      readThermalTable
!
! Description:  Reads in a computed ascii table derived from the APEC model.
!
! Changed 1/9/2018 by JDS - now allow for comments in the table if the lines 
!     begin with # and are at the beginning of the file.
! Changed 3/15/2018 by JDS - now use table that includes ne/n - ratio of
!     electron density to nucleon density; scheme now modified to account for
!     variations in ionization

subroutine readThermalTable(Nt)

!===============================================================================

  use Cool_data, ONLY: cl_Nmax, cl_T, cl_lambda, cl_dene, cl_chi, cl_XT, &
      cl_dL, cl_Nt, cl_functionFile, cl_Tcut, useChi

  implicit none

  integer, intent(out) :: Nt
  
  real :: T, lambda, dene
  integer :: i, j, status
  character(len=72) :: line
  
  !=============================================================================
  
  ! First, determine whether the file can be opened and has the right format.
  
  open (1, file=cl_functionFile(1:len_trim(cl_functionFile)), iostat=status)
  if (status /= 0) then
     Nt = -1
     return
  endif
  
  Nt = 0
  do while (.true.)
      read(1,fmt='(A72)',iostat=status) line
      if (status == -1) exit
      if (status > 0) then
          Nt = -1
          close(1)
          return
      endif
      if (line(1:1) .ne. '#') then
          read(line, '(I4)') Nt
          exit
      endif
  enddo

  if ((Nt == 0) .or. (Nt > cl_Nmax)) then
     Nt = -1
     close(1)
     return
  endif

  j = 0
  do i = 1, Nt
     ! cooling files contain columns: log10(T), log10(L) + 23, ne/n
     if (useChi) then
         read(1,*) T, lambda, dene
     else
         read(1,*) T, lambda
         dene = 1.0909091
     endif
     if (T > cl_Tcut) then
         j = j + 1
         cl_T(j) = T
         cl_lambda(j) = lambda
         cl_dene(j) = dene
     else
         cycle
     endif
     cl_chi(j) = 1. + cl_dene(j)
     cl_T(j) = 10.**(cl_T(j))
     cl_lambda(j) = 10.**(cl_lambda(j) - 23.)
     cl_XT(j) = cl_chi(j)*cl_T(j)
     cl_dL(j) = cl_dene(j)*cl_lambda(j)
  enddo
  close(1)
  Nt = j
  cl_Nt = Nt
  
  !=====================================================================
  
  return
end subroutine readThermalTable
