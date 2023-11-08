!!****if* source/Particles/ParticlesForces/shortRange/drag/pt_yieldData
!!
!! NAME
!!
!!  pt_chargeData
!!
!! SYNOPSIS
!!
!!  use pt_chargeData 
!!
!!  DESCRIPTION
!!
!!  Stores the data for calculating the charge on a dust grain
!!  
!! PARAMETERS
!!
!!     nv               number of velocities in grid
!!     nT               number of temperatures in grid
!!     nz               number of zeta values in grid
!!     phiSi            phi (potential parameter) array for Silicate grains
!!     phiaC            phi (potential parameter) array for Carbonaceous grains
!!     Tgrid            Temperature grid in log T(K)
!!     vgrid            velocity grid in log v(cm/s)
!!     zgrid            zeta grid in log zeta (dimensionless)
!!     dTg              spacing for Temperature grid in log T(K)
!!     dvg              spacing for velocity grid in log v(cm/s)
!!     dzg              spacing for zeta grid in log zeta
!!     Tgmin            minimum Temperature in grid in log T(K)
!!     vgmin            minimum velocity in grid in log v(cm/s)
!!     zgmin            minimum zeta in grid in log zeta
!!
!! Note that zeta = G0*Qabs/n, G0 is scaling of FUV background relative to
!! Draine & Salpeter's standard value, Qabs is the absorption coefficient and
!! n is the density (actually H+ density or electron density, which we assume
!! are the same. We ignore Qabs (assume it's equal to 1, effectively), so we're
!! assuming that zeta = G0/n
!!
!!***
module pt_chargeData
   integer, save :: nv,nT,nz
   real, dimension(51,31,6), save :: phiSi,phiaC
   real, dimension(51), save :: Tgrid
   real, dimension(31), save :: vgrid
   real, dimension(6), save :: zgrid
   real, save :: dTg, dvg, dzg, Tgmin, vgmin, zgmin, G0
end module pt_chargeData
