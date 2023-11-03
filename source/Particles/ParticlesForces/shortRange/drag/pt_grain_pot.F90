subroutine pt_grain_pot(dens,temp,vel,ptype,phi)
    ! Interpolate grain potential parameter grid over zeta, temperature and
    ! velocity (phi = Z*e^2/a*k*T)
    !
    ! Input parameters:
    ! dens              density (H or electron) in cm^-3
    ! temp              gas temperature in K
    ! vel               relative gas-grain velocity in cm/s
    ! ptype             particle (grain) type, 1 => C, 2 => Si
    !
    ! Output:
    ! phi               grain potential parameter (dimensionless) = Ze^2/akT

    use pt_chargeData
    use Simulation_data, ONLY : sim_G0

    implicit none
    real, intent(in) :: dens,temp,vel
    integer, intent(in) :: ptype
    real, intent(out) :: phi
    integer :: zind,Tind,vind
    real :: pz, pt, pv, zeta
    real :: ph111, ph211, ph121, ph221, ph112, ph212, ph122, ph222

    zeta = sim_G0/dens
    Tind = floor((log10(temp) - Tgmin)/dTg) + 1
    vind = floor((log10(vel) - 5.0 - vgmin)/dvg) + 1
    zind = floor((log10(zeta) - zgmin)/dzg) + 1
    !write(*,'("temp =",ES11.3," vel =",ES11.3," zeta =",ES11.3)') temp,vel,zeta
    !write(*,'("Tind =",I3," vind =",I3," zind =",I3)') Tind,vind,zind
    !write(*,'("Tgrid(Tind) =",F5.2," vgrid(vind) =",F5.2," zgrid(zind) =",' &
    !    // 'F5.2)') Tgrid(Tind),vgrid(vind),zgrid(zind)
    if(Tind.le.0) then
        Tind = 1
        pt = 0.
    else if(Tind.ge.nT) then
        Tind = nT - 1
        pt = 1.
    else
        pt = (log10(temp) - Tgrid(Tind))/dTg
    endif

    if(vind.le.0) then
        vind = 1
        pv = 0.
    else if(vind.ge.nv) then
        vind = nv - 1
        pv = 1.
    else
        pv = (log10(vel) - vgrid(vind))/dvg
    endif

    if(zind.lt.1) then
        zind = 1
        pz = 0.
    else if(zind.ge.nz) then
        zind = nz - 1
        pz = 1.
    else
        pz = (log10(zeta) - zgrid(zind))/dzg
    endif

    if(ptype.eq.1) then
        ph111 = phiaC(Tind,vind,zind)
        ph211 = phiaC(Tind+1,vind,zind)
        ph121 = phiaC(Tind,vind+1,zind)
        ph221 = phiaC(Tind+1,vind+1,zind)
        ph112 = phiaC(Tind,vind,zind+1)
        ph212 = phiaC(Tind+1,vind,zind+1)
        ph122 = phiaC(Tind,vind+1,zind+1)
        ph222 = phiaC(Tind+1,vind+1,zind+1)
    else
        ph111 = phiSi(Tind,vind,zind)
        ph211 = phiSi(Tind+1,vind,zind)
        ph121 = phiSi(Tind,vind+1,zind)
        ph221 = phiSi(Tind+1,vind+1,zind)
        ph112 = phiSi(Tind,vind,zind+1)
        ph212 = phiSi(Tind+1,vind,zind+1)
        ph122 = phiSi(Tind,vind+1,zind+1)
        ph222 = phiSi(Tind+1,vind+1,zind+1)
    endif
    ! trilinear interpolation
    phi = ((ph111*(1. - pt) + ph211*pt)*(1. - pv) + &
        (ph121*(1. - pt) + ph221*pt)*pv)*(1. - pz) + &
        ((ph112*(1. - pt) + ph212*pt)*(1. - pv) + &
        (ph122*(1. - pt) + ph222*pt)*pv)*pz 
end subroutine pt_grain_pot
