!!****if* source/Particles/ParticlesForces/shortRange/drag/pt_yieldInit
!!
!! NAME
!!
!!  pt_yieldInit
!!
!!
!! SYNOPSIS
!!
!!  pt_yieldInit()
!!
!! DESCRIPTION
!!
!!  Initializes all the data needed for sputtering calculations.
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!     mproj            atomic mass (A) of gas species 'ion' (amu) 
!!     Zproj            atomic number (Z) of gas species 'ion'
!!     Cproj            charge of gas ion species 'ion'
!!     amu              amu in g
!!     nion             no. of different (projectile) ions included
!!     solidat          array including sputtering parameters
!!     sputparm(ntype,1) surface binding energy (U) of target (eV)
!!     sputparm(ntype,2) atomic mass of sputtered fragment (amu), mtarg
!!     sputparm(ntype,3) atomic number of sputtered fragment, Ztarg
!!     sputparm(ntype,4) free parameter "K" constant
!!
!!***

subroutine pt_yieldInit(dustdata,yielddata,yieldintdata)

    use pt_yieldIntgrlData
    use HDF5

    implicit none
    character(len=*), intent(IN) :: dustdata, yielddata, yieldintdata
    integer :: error,i,j
    real, dimension(5) :: depln, cosab

    integer(HID_T) :: file_id,group1_id,group2_id,group3_id,group4_id
    integer(HID_T) :: group31_id,group32_id,group33_id,group34_id
    integer(HID_T) :: group41_id,group42_id,group43_id,group44_id
    integer(HID_T) :: dset1_id,dset2_id
    integer(HID_T) :: dset31_id,dset32_id,dset33_id,dset34_id
    integer(HID_T) :: dset41_id,dset42_id,dset43_id,dset44_id
    character(LEN=13), parameter :: group1 = 'velocity grid'
    character(LEN=16), parameter :: group2 = 'temperature grid'
    character(LEN=9), parameter :: group3 = 'aC grains'
    character(LEN=11), parameter :: group4 = 'SiO4 grains'
    character(LEN=12), parameter :: group31 = 'H_sputtering'
    character(LEN=13), parameter :: group32 = 'He_sputtering'
    character(LEN=12), parameter :: group33 = 'O_sputtering'
    character(LEN=12), parameter :: group34 = 'C_sputtering'
    character(LEN=12), parameter :: group41 = 'H_sputtering'
    character(LEN=13), parameter :: group42 = 'He_sputtering'
    character(LEN=12), parameter :: group43 = 'O_sputtering'
    character(LEN=12), parameter :: group44 = 'C_sputtering'
    character(LEN=4), parameter :: dset1 = 'vtbl'
    character(LEN=4), parameter :: dset2 = 'Ttbl'
    character(LEN=6), parameter :: dset31 = 'Yv_aCH'
    character(LEN=7), parameter :: dset32 = 'Yv_aCHe'
    character(LEN=6), parameter :: dset33 = 'Yv_aCO'
    character(LEN=6), parameter :: dset34 = 'Yv_aCC'
    character(LEN=6), parameter :: dset41 = 'Yv_SiH'
    character(LEN=7), parameter :: dset42 = 'Yv_SiHe'
    character(LEN=6), parameter :: dset43 = 'Yv_SiO'
    character(LEN=6), parameter :: dset44 = 'Yv_SiC'

    character(LEN=5), parameter :: attrv1 = 'nvtbl'
    character(LEN=5), parameter :: attrv2 = 'dvtbl'
    character(LEN=5), parameter :: attrT1 = 'nTtbl'
    character(LEN=5), parameter :: attrT2 = 'dTtbl'
    integer(HID_T) :: attrv1_id, attrv2_id, attrT1_id, attrT2_id

    integer(HSIZE_T), dimension(1:1) :: dims
    integer(HSIZE_T), dimension(2) :: dims2

    nion = 6
    amu = 1.66053873E-24
    open(20,file=dustdata,status='old')
    ! Data for H+, He+, O+, C+, and N+
    read(20,'(5(F7.2))') (mproj(j),j=1,5)
    read(20,'(5(F7.2))') (Zproj(j),j=1,5)
    read(20,'(5(F7.2))') (Cproj(j),j=1,5)
    read(20,'(5(F7.2))') (depln(j),j=1,5)
    read(20,'(5(F7.2))') (cosab(j),j=1,5)
    close(20)
    ! For He++:
    mproj(6) = mproj(2)
    Zproj(6) = Zproj(2)
    Cproj(6) = 2.*Cproj(2)
    do i=1,5
        abund(i) = (1. - depln(i))*cosab(i)
    enddo
    ! get sputtering parameters U, mtarg, Ztarg, K
    !open(21,file='solidptrs',status='old')
    open(21,file=yielddata,status='old')
    ! Note: solidptrs also includes data on the solid particle mass density, 
    ! but we're inputting that via Simulation_init and pt_initPositions
    ! (also it includes vaporization data)
    do i=1,2
        read(21,'(47X,4(F7.2))') (sputparm(i,j),j=1,4)
    enddo
    close(21)

    ! Begin HDF5 section
    CALL h5open_f(error)
    CALL h5fopen_f(yieldintdata, H5F_ACC_RDONLY_F, file_id, error)

    ! velocity array data
    CALL h5gopen_f(file_id, group1, group1_id, error)
    CALL h5aopen_f(group1_id, attrv1, attrv1_id, error)
    CALL h5aread_f(attrv1_id, H5T_NATIVE_INTEGER, nvtbl, dims, error)
    CALL h5aclose_f(attrv1_id, error)

    CALL h5aopen_f(group1_id, attrv2, attrv2_id, error)
    CALL h5aread_f(attrv2_id, H5T_NATIVE_DOUBLE, dvtbl, dims, error)
    CALL h5aclose_f(attrv2_id, error)

    dims(1) = nvtbl
    !if(.not.allocated(vtbl)) allocate(vtbl(nvtbl))
    CALL h5dopen_f(group1_id, dset1, dset1_id, error)
    CALL h5dread_f(dset1_id, H5T_NATIVE_DOUBLE, vtbl, dims, error)
    CALL h5dclose_f(dset1_id, error)
    CALL h5gclose_f(group1_id, error)

    ! temperature array data
    CALL h5gopen_f(file_id, group2, group2_id, error)
    CALL h5aopen_f(group2_id, attrT1, attrT1_id, error)
    CALL h5aread_f(attrT1_id, H5T_NATIVE_INTEGER, nTtbl, dims, error)
    CALL h5aclose_f(attrT1_id, error)

    CALL h5aopen_f(group2_id, attrT2, attrT2_id, error)
    CALL h5aread_f(attrT2_id, H5T_NATIVE_DOUBLE, dTtbl, dims, error)
    CALL h5aclose_f(attrT2_id, error)

    !write(*,'("nvtbl =",i4," dvtbl =",F8.2)') nvtbl,dvtbl
    !write(*,'("nTtbl =",i4," dTtbl =",F8.2)') nTtbl,dTtbl

    dims(1) = nTtbl
    !if(.not.allocated(Ttbl)) allocate(Ttbl(nTtbl))
    CALL h5dopen_f(group2_id, dset2, dset2_id, error)
    CALL h5dread_f(dset2_id, H5T_NATIVE_DOUBLE, Ttbl, dims, error)
    CALL h5dclose_f(dset2_id, error)
    CALL h5gclose_f(group2_id, error)

    !! Carbonaceous grain data
    CALL h5gopen_f(file_id, group3, group3_id, error)
    dims2(1) = nTtbl
    dims2(2) = nvtbl
    !if(.not.allocated(Yv_aCH)) allocate(Yv_aCH(nTtbl,nvtbl))
    CALL h5gopen_f(group3_id, group31, group31_id, error)
    CALL h5dopen_f(group31_id, dset31, dset31_id, error)
    CALL h5dread_f(dset31_id, H5T_NATIVE_DOUBLE, Yv_aCH, dims2, error)
    CALL h5dclose_f(dset31_id, error)
    CALL h5gclose_f(group31_id, error)

    !if(.not.allocated(Yv_aCHe)) allocate(Yv_aCHe(nTtbl,nvtbl))
    CALL h5gopen_f(group3_id, group32, group32_id, error)
    CALL h5dopen_f(group32_id, dset32, dset32_id, error)
    CALL h5dread_f(dset32_id, H5T_NATIVE_DOUBLE, Yv_aCHe, dims2, error)
    CALL h5dclose_f(dset32_id, error)
    CALL h5gclose_f(group32_id, error)

    !if(.not.allocated(Yv_aCO)) allocate(Yv_aCO(nTtbl,nvtbl))
    CALL h5gopen_f(group3_id, group33, group33_id, error)
    CALL h5dopen_f(group33_id, dset33, dset33_id, error)
    CALL h5dread_f(dset33_id, H5T_NATIVE_DOUBLE, Yv_aCO, dims2, error)
    CALL h5dclose_f(dset33_id, error)
    CALL h5gclose_f(group33_id, error)

    !if(.not.allocated(Yv_aCC)) allocate(Yv_aCC(nTtbl,nvtbl))
    CALL h5gopen_f(group3_id, group34, group34_id, error)
    CALL h5dopen_f(group34_id, dset34, dset34_id, error)
    CALL h5dread_f(dset34_id, H5T_NATIVE_DOUBLE, Yv_aCC, dims2, error)
    CALL h5dclose_f(dset34_id, error)
    CALL h5gclose_f(group34_id, error)

    CALL h5gclose_f(group3_id, error)

    !! Silicate grain data
    CALL h5gopen_f(file_id, group4, group4_id, error)

    !if(.not.allocated(Yv_SiH)) allocate(Yv_SiH(nTtbl,nvtbl))
    CALL h5gopen_f(group4_id, group41, group41_id, error)
    CALL h5dopen_f(group41_id, dset41, dset41_id, error)
    CALL h5dread_f(dset41_id, H5T_NATIVE_DOUBLE, Yv_SiH, dims2, error)
    CALL h5dclose_f(dset41_id, error)
    CALL h5gclose_f(group41_id, error)

    !if(.not.allocated(Yv_SiHe)) allocate(Yv_SiHe(nTtbl,nvtbl))
    CALL h5gopen_f(group4_id, group42, group42_id, error)
    CALL h5dopen_f(group42_id, dset42, dset42_id, error)
    CALL h5dread_f(dset42_id, H5T_NATIVE_DOUBLE, Yv_SiHe, dims2, error)
    CALL h5dclose_f(dset42_id, error)
    CALL h5gclose_f(group42_id, error)

    !if(.not.allocated(Yv_SiO)) allocate(Yv_SiO(nTtbl,nvtbl))
    CALL h5gopen_f(group4_id, group43, group43_id, error)
    CALL h5dopen_f(group43_id, dset43, dset43_id, error)
    CALL h5dread_f(dset43_id, H5T_NATIVE_DOUBLE, Yv_SiO, dims2, error)
    CALL h5dclose_f(dset43_id, error)
    CALL h5gclose_f(group43_id, error)

    !if(.not.allocated(Yv_SiC)) allocate(Yv_SiC(nTtbl,nvtbl))
    CALL h5gopen_f(group4_id, group44, group44_id, error)
    CALL h5dopen_f(group44_id, dset44, dset44_id, error)
    CALL h5dread_f(dset44_id, H5T_NATIVE_DOUBLE, Yv_SiC, dims2, error)
    CALL h5dclose_f(dset44_id, error)
    CALL h5gclose_f(group44_id, error)

    CALL h5gclose_f(group4_id, error)

    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

end subroutine pt_yieldInit
