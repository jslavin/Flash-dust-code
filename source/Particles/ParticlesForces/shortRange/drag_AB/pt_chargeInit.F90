!!****if* source/Particles/ParticlesForces/shortRange/drag/pt_chargeInit
!!
!! NAME
!!
!!  pt_chargeInit
!!
!!
!! SYNOPSIS
!!
!!  pt_yieldInit(charge_data_file)
!!
!! DESCRIPTION
!!
!!  Initializes the data for charging of dust grains
!!
!! ARGUMENTS
!!
!!     charge_data_file  name of text file containing tables of phi (potential
!!                       parameter) vs. Temperature, grain velocity and zeta
!!                       where zeta = G0*Qabs/n, G0 is UV field scaling 
!!                       relative to Draine's standard field, Qabs is the
!!                       absorption coefficient and n is the density
!!
!!***

subroutine pt_chargeInit(charge_data_file)

    use pt_chargeData

    implicit none
    character(len=*), intent(in) :: charge_data_file
    integer :: i,j,k
    real :: lgv

    nT = 51
    nv = 31
    nz = 6
    dTg = 0.1
    dvg = 0.1
    dzg = 0.5
    Tgmin = 3.0
    vgmin = 1.0
    zgmin = -0.5
    open(unit=99,file=charge_data_file,status='old')
    ! skip over header infomation
    do i = 1,3
        read(99,*)
    enddo
    ! Silicate tables
    do i=1,nz
        read(99,*)
        read(99,*)
        do j=1,nv
            read(99,*) lgv,(phiSi(k,j,i),k=1,nT)
            if(i.eq.1) vgrid(j) = lgv + 5.
        enddo
    enddo
    ! Carbonaceous tables
    do i=1,nz
        read(99,*)
        read(99,*)
        do j=1,nv
            read(99,*) lgv,(phiaC(k,j,i),k=1,nT)
        enddo
    enddo
    do i=1,nT
        Tgrid(i) = Tgmin + dTg*(i - 1.)
    enddo
    do i=1,nz
        zgrid(i) = zgmin + dzg*(i - 1.)
    enddo
    close(99)
end subroutine pt_chargeInit
