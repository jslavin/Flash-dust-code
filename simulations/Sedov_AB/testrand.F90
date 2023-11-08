program testrand
implicit none
intrinsic random_seed, random_number
integer, dimension(1) :: seed
integer :: n=10, i
real :: r
seed(1) = 314159265
!call random_seed(put=seed)
call random_seed()
do i=1,n
    call random_number(r)
    write(*,'("r =",F7.4)') r
enddo  
!integer size, seed(2), gseed(2), hiseed(2), zseed(2)
!real harvest(10)
!data seed /123456789, 987654321/
!data hiseed /-1, -1/
!data zseed /0, 0/
!call random_seed(SIZE=size)
!print *,"size ",size
!call random_seed(PUT=hiseed(1:size))
!call random_seed(GET=gseed(1:size))
!print *,"hiseed gseed", hiseed, gseed
!call random_seed(PUT=zseed(1:size))
!call random_seed(GET=gseed(1:size))
!print *,"zseed gseed ", zseed, gseed
!call random_seed(PUT=seed(1:size))
!call random_seed(GET=gseed(1:size))
!call random_number(HARVEST=harvest)
!print *, "seed gseed ", seed, gseed
!print *, "harvest"
!print *, harvest
!call random_seed(GET=gseed(1:size))
!print *,"gseed after harvest ", gseed
end program testrand
