subroutine pt_sputYieldVel(ntype,vgrain,Tgas,ion,yieldvel)
!!     Subroutine to get the sputtering yield times grain velocity averaged 
!!     over the skewed maxwellian particle distribution for a given ion 
!!     (velocity, atomic mass and atomic number specified as input 
!!     parameters) as a function of target material (binding energy, atomic
!      mass, atomic number and sputtering parameter, K, given as input 
!!     parameters).
!!                 
!!     Parameters:
!!     ntype = index of grain type (aC=1, SiO4=2, ice=3, SiC=4, Fe=5, 
!!             diamond=6)
!!     vgrain = incident species velocity  (cm/s)
!!     ion = index of projectile ion species (H=1, He=2, O=3, C=4, N=5)
!!
!!     Returns:
!!     Yv  = sputtering yield Y ('atoms' per incident ion) times vgrain
!!         integrated over the particle velocity distribution
!!    
    use pt_yieldIntgrlData

    implicit none

    integer, intent(in) :: ntype,ion
    real, intent(in) :: vgrain,Tgas
    real, intent(out) :: yieldvel
    real :: yv11,yv12,yv21,yv22
    integer :: indexv1,indexv2,indexT1,indexT2
    real :: p,q

    !write(*,'("in pt_sputYieldVel: vtbl(1) =",ES12.5," dvtbl =",ES12.5)') &
    !    vtbl(1),dvtbl
    indexv1 = floor(log10(vgrain/vtbl(1))/dvtbl) + 1
    indexT1 = floor(log10(Tgas/Ttbl(1))/dTtbl) + 1
    indexv1 = max(min(indexv1,nvtbl),1)
    indexT1 = max(min(indexT1,nTtbl),1)
    if((indexT1 == nTtbl).or.(Tgas <= Ttbl(1))) then
        indexT2 = indexT1
        p = 1.
    else
        indexT2 = indexT1 + 1
        p = log10(Tgas/Ttbl(indexT1))/dTtbl
    endif
    if((indexv1 == nvtbl).or.(vgrain <= vtbl(1))) then
        indexv2 = indexv1
        q = 1.
    else
        indexv2 = indexv1 + 1
        q = log10(vgrain/vtbl(indexv1))/dvtbl
    endif
    if(ntype == 1) then
        if(ion == 1) then
            yv11 = Yv_aCH(indexT1,indexv1)
            yv12 = Yv_aCH(indexT1,indexv2)
            yv21 = Yv_aCH(indexT2,indexv1)
            yv22 = Yv_aCH(indexT2,indexv2)
        else if(ion == 2) then
            yv11 = Yv_aCHe(indexT1,indexv1)
            yv12 = Yv_aCHe(indexT1,indexv2)
            yv21 = Yv_aCHe(indexT2,indexv1)
            yv22 = Yv_aCHe(indexT2,indexv2)
        else if(ion == 3) then
            yv11 = Yv_aCO(indexT1,indexv1)
            yv12 = Yv_aCO(indexT1,indexv2)
            yv21 = Yv_aCO(indexT2,indexv1)
            yv22 = Yv_aCO(indexT2,indexv2)
        else if(ion == 4) then
            yv11 = Yv_aCC(indexT1,indexv1)
            yv12 = Yv_aCC(indexT1,indexv2)
            yv21 = Yv_aCC(indexT2,indexv1)
            yv22 = Yv_aCC(indexT2,indexv2)
        endif
    else if(ntype == 2) then
        if(ion == 1) then
            yv11 = Yv_SiH(indexT1,indexv1)
            yv12 = Yv_SiH(indexT1,indexv2)
            yv21 = Yv_SiH(indexT2,indexv1)
            yv22 = Yv_SiH(indexT2,indexv2)
        else if(ion == 2) then
            yv11 = Yv_SiHe(indexT1,indexv1)
            yv12 = Yv_SiHe(indexT1,indexv2)
            yv21 = Yv_SiHe(indexT2,indexv1)
            yv22 = Yv_SiHe(indexT2,indexv2)
        else if(ion == 3) then
            yv11 = Yv_SiO(indexT1,indexv1)
            yv12 = Yv_SiO(indexT1,indexv2)
            yv21 = Yv_SiO(indexT2,indexv1)
            yv22 = Yv_SiO(indexT2,indexv2)
        else if(ion == 4) then
            yv11 = Yv_SiC(indexT1,indexv1)
            yv12 = Yv_SiC(indexT1,indexv2)
            yv21 = Yv_SiC(indexT2,indexv1)
            yv22 = Yv_SiC(indexT2,indexv2)
        endif
    endif
    ! bilinear interpolation
    yieldvel = (yv11*(1. - p) + yv21*p)*(1. - q) + (yv12*(1. - p) + yv22*p)*q
    ! Yv_x values are log10s
    yieldvel = 10.**yieldvel
end subroutine pt_sputYieldVel
