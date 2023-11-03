!*******************************************************************************

! Routine:      readCloudyTable

! Description:  Reads in a table created by running Cloudy as a subroutine
!               to produce a plasma cooling and heating data.
!
! Parameters:   file              Name of file to read.
!               logT(:)           1D array returning log temperature in K.
!               logL(:)           Array returning cooling function L
!                                   in log(erg cm^3 s^-1)
!               logG(:)           Array returning heating function Gamma
!                                   in log(erg s^-1)
!               dene(:)           Array returning the ratio of electron number
!                                 density to nuclei density (dimensionless)
!               Nmax              Size of 1D arrays.
!               N                 Integer returning number of data rows found
!                                   in the file.  If > 0, this will be the
!                                   size of the arrays allocated for each of
!                                   the above pointers.  If < 0, there was a
!                                   problem and no data were read.
!

subroutine readCloudyTable (filename, logT, logL, logG, dene, Nmax, N)

!===============================================================================

  use Cool_data, ONLY : cl_useHeat, cl_meshMe
  use Driver_interface, ONLY : Driver_getMype

#include "constants.h"

  implicit none
  
  character(len=*), intent(IN) :: filename
  integer, intent(IN)   :: Nmax
  integer, intent(OUT)  :: N

  real, intent(INOUT), dimension(Nmax) :: logT, logL, logG, dene
  
  character (LEN=72) line
  integer :: i, nskip, status
  
!===============================================================================
  
  ! First, determine whether the file can be opened and has the right format.
  ! Count the data rows in the file.
  
  open(1, file=filename, iostat=status)
  if (status /= 0) then
      N = -1
      return
  endif
  
  N = 0
  nskip = 0
  do while (.true.)
      read(1,fmt='(A72)',iostat=status) line
      if (status == -1) exit
      if (status > 0) then
          N = -1
          close(1)
          return
      endif
      if (line(1:1).eq.'#') then
          nskip = nskip + 1
      else
          N = N + 1
      endif
  enddo
  
  close (1)
  
  if ((N == 0) .or. (N > Nmax)) then
      print *,'Problem with readCloudyTable: N =',N
      N = -1
      return
  endif
  ! These lines are to deal with files that have a line giving the no. of lines
  ! in the file on the first line after the header lines that start with #
  nskip = nskip + 1
  N = N - 1
  !! for cooling tables that don't have a line with the no. of lines after the
  !! header, should comment out those lines
  
  ! Next, free up any space that might have been allocated previously.
  
  !if (associated(logT))     deallocate (logT)
  !if (associated(logL))     deallocate (logL)
  !if (associated(logG))     deallocate (logG)
  !if (associated(dene))     deallocate (dene)
  
  ! Now allocate space for the arrays.
  
  !allocate (logT(N))
  !allocate (logL(N))
  !allocate (logG(N))
  !allocate (dene(N))
  
!===============================================================================
  
  ! Open the file again and read the data.
  
  open (1, file=filename, iostat=status)
  
  do i = 1, nskip    ! lines with # in first column contain header info
     read (1,*)
  enddo
  !call Driver_getMype(MESH_COMM,cl_meshMe)
  !if (cl_meshMe == MASTER_PE) write(*,*) 'N =',N
  do i = 1, N
     if(cl_useHeat) then
         read (1,*) logT(i), logL(i), logG(i), dene(i)
     else
         read (1,*) logT(i), logL(i), dene(i)
         !write(*,fmt='(i3,ES12.5,ES12.5,ES12.5)') i,logT(i), logL(i), dene(i)
         logG(i) = -10.
     endif

     ! values in table are multiplied by a factor of 1.E23
     logL(i) = logL(i) - 23.
     logG(i) = logG(i) - 23.
  enddo
  
  close (1)
  
!===============================================================================
  
  return
end subroutine readCloudyTable
