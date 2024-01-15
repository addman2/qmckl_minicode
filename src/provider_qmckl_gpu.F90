subroutine provider_qmckl_gpu(point&
                           &, ang_mom&
                           &, par&
                           &, values&
                           &, multiplicity&
                           &, gl)

    use qmckl_gpu_helpers

    implicit none

    integer(kind=4), intent(in) :: ang_mom
    integer(kind=4), intent(in) :: multiplicity
    integer(kind=4), intent(in) :: gl

    real(kind=8), target, dimension(3), intent(in) :: point
    real(kind=8), intent(in) :: par

    ! If gl is 0 only values are computed
    ! If gl is 1 values gradients and laplacians are computed
    ! this requires 5 times more space
    real(kind=8), dimension(multiplicity * (1 + gl * 4)), intent(out), target :: values

#ifdef _QMCKL_GPU

    integer(kind=4) :: rc

    ! Here we are defining local varibales to hold the data that will be passed

    real(kind=8), target, dimension(1) :: charges
    real(kind=8), target, dimension(3,1) :: positions
    real(kind=8), target, dimension(1) :: shell_factor
    real(kind=8), target, dimension(1) :: exponent
    real(kind=8), target, dimension(1) :: coefficient
    real(kind=8), target, dimension(1) :: prim_factor
    real(kind=8), target, dimension(:), allocatable :: ao_factor

    integer(kind=8), target, dimension(1) :: nucleus_shell_num
    integer(kind=8), target, dimension(1) :: nucleus_index
    integer(kind=4), target, dimension(1) :: shell_ang_mom
    integer(kind=8), target, dimension(1) :: shell_prim_num
    integer(kind=8), target, dimension(1) :: shell_prim_index

    integer(kind=4), dimension(1) :: array

    !##############################################################################
    !#                                                                            #
    !#  Set up the QMCkl context.                                                 #
    !#                                                                            #
    !##############################################################################

    qmckl_gpu_device_id = 0

    ! First test OpenMP offload

    array = 1
    !$omp target data map(tofrom:array)
    array = 2
    !$omp end target data
    if (array(1) /= 1) then
        print *, "OpenMP offload failed"
        stop 1
    else
        print *, "OpenMP offload succeeded"
    end if

    call qmckl_gpu_init()

    !$omp target data map(to:charges,positions,shell_factor&
    !$omp                  &,exponent,coefficient,prim_factor&
    !$omp                  &,nucleus_shell_num,nucleus_index&
    !$omp                  &,shell_ang_mom,shell_prim_num&
    !$omp                  &,shell_prim_index,point)

    rc = qmckl_set_nucleus_num_device(qmckl_gpu_ctx, 1_8)
    if (rc /= 0) then
        print *, "qmckl_set_nucleus_num failed"
        stop 1
    else
        print *, "qmckl_set_nucleus_num succeeded"
    end if

    charges(1) = 1.0_8
    !$omp target update to(charges)

    !$omp target data use_device_ptr(charges)
    rc = qmckl_set_nucleus_charge_device(qmckl_gpu_ctx, c_loc(charges), 1_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_nucleus_charge failed"
        stop 1
    else
        print *, "qmckl_set_nucleus_charge succeeded"
    end if

    ! Check if the nucleus charge is set correctly

    charges(1) = 0.0_8
    !$omp target update to(charges)

    !$omp target data use_device_ptr(charges)
    rc = qmckl_get_nucleus_charge_device(qmckl_gpu_ctx, c_loc(charges), 1_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_get_nucleus_charge failed"
        stop 1
    else
        print *, "qmckl_get_nucleus_charge succeeded"
    end if

    !$omp target update from(charges)
    if (charges(1) /= 1.0_8) then
        print *, "qmckl_set_nucleus_charge failed 2"
        stop 1
    else
        print *, "qmckl_set_nucleus_charge succeeded 2"
    end if

    positions(:,1) = (/ 0.0_8, 0.0_8, 0.0_8 /)
    !$omp target update to(positions)

    !$omp target data use_device_ptr(positions)
    rc = qmckl_set_nucleus_coord_device(qmckl_gpu_ctx, "N", c_loc(positions), 3_8 * size(positions))
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_nucleus_cood failed"
        stop 1
    else
        print *, "qmckl_set_nucleus_cood succeeded"
    end if

    !##############################################################################
    !#                                                                            #
    !#  Set up the basis set.                                                     #
    !#                                                                            #
    !##############################################################################

    rc = qmckl_set_ao_basis_type_device(qmckl_gpu_ctx, "G")
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_type failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_type succeeded"
    end if

    rc = qmckl_set_ao_basis_shell_num_device(qmckl_gpu_ctx, 1_8)
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_shell_num failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_shell_num succeeded"
    end if

    rc = qmckl_set_ao_basis_prim_num_device(qmckl_gpu_ctx, 1_8)
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_prim_num failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_prim_num succeeded"
    end if

    nucleus_shell_num(1) = 1
    !$omp target update to(nucleus_shell_num)

    !$omp target data use_device_ptr(nucleus_shell_num)
    rc = qmckl_set_ao_basis_nucleus_shell_num_device(qmckl_gpu_ctx, c_loc(nucleus_shell_num), 1_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_nucleus_shell_num failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_nucleus_shell_num succeeded"
    end if

    nucleus_index(1) = 0
    !$omp target update to(nucleus_index)

    !$omp target data use_device_ptr(nucleus_index)
    rc = qmckl_set_ao_basis_nucleus_index_device(qmckl_gpu_ctx, c_loc(nucleus_index), 1_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_nucleus_index failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_nucleus_index succeeded"
    end if

    shell_ang_mom(1) = ang_mom
    !$omp target update to(shell_ang_mom)

    !$omp target data use_device_ptr(shell_ang_mom)
    rc = qmckl_set_ao_basis_shell_ang_mom_device(qmckl_gpu_ctx, c_loc(shell_ang_mom), 1_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_shell_ang_mom failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_shell_ang_mom succeeded"
    end if

    shell_prim_num(1) = 1
    !$omp target update to(shell_prim_num)
    
    !$omp target data use_device_ptr(shell_prim_num)
    rc = qmckl_set_ao_basis_shell_prim_num_device(qmckl_gpu_ctx, c_loc(shell_prim_num), 1_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_shell_prim_num failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_shell_prim_num succeeded"
    end if

    shell_prim_index(1) = 0
    !$omp target update to(shell_prim_index)

    !$omp target data use_device_ptr(shell_prim_index)
    rc = qmckl_set_ao_basis_shell_prim_index_device(qmckl_gpu_ctx, c_loc(shell_prim_index), 1_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_shell_prim_index failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_shell_prim_index succeeded"
    end if

    shell_factor(1) = 1.0_8
    !$omp target update to(shell_factor)

    !$omp target data use_device_ptr(shell_factor)
    rc = qmckl_set_ao_basis_shell_factor_device(qmckl_gpu_ctx, c_loc(shell_factor), 1_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_shell_factor failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_shell_factor succeeded"
    end if

    exponent(1) = par
    !$omp target update to(exponent)

    !$omp target data use_device_ptr(exponent)
    rc = qmckl_set_ao_basis_exponent_device(qmckl_gpu_ctx, c_loc(exponent), 1_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_exponent failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_exponent succeeded"
    end if

    coefficient(1) = 1.0_8
    !$omp target update to(coefficient)

    !$omp target data use_device_ptr(coefficient)
    rc = qmckl_set_ao_basis_coefficient_device(qmckl_gpu_ctx, c_loc(coefficient), 1_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_coefficient failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_coefficient succeeded"
    end if

    prim_factor(1) = 1.0_8
    !$omp target update to(prim_factor)
    
    !$omp target data use_device_ptr(prim_factor)
    rc = qmckl_set_ao_basis_prim_factor_device(qmckl_gpu_ctx, c_loc(prim_factor), 1_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_prim_factor failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_prim_factor succeeded"
    end if

    allocate(ao_factor(multiplicity))
    ao_factor = 1.0_8
    !$omp target data map(to:ao_factor)

    rc = qmckl_set_ao_basis_ao_num_device(qmckl_gpu_ctx, 1_8 * multiplicity)
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_ao_num failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_ao_num succeeded"
    end if

    !$omp target data use_device_ptr(ao_factor)
    rc = qmckl_set_ao_basis_ao_factor_device(qmckl_gpu_ctx, c_loc(ao_factor), 1_8 * multiplicity)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_ao_factor failed"
        print *, "Error code", rc
        stop 1
    else
        print *, "qmckl_set_ao_basis_ao_factor succeeded"
    end if

    !$omp end target data

    !##############################################################################
    !#                                                                            #
    !#  Now calculate values of the AO basis                                      #
    !#                                                                            #
    !##############################################################################

    !$omp target data use_device_ptr(point)
    rc = qmckl_set_point_device(qmckl_gpu_ctx&
                             &, "N"&
                             &, 1_8&
                             &, c_loc(point)&
                             &, 3_8)
    !$omp end target data
    if (rc /= 0) then
        print *, "qmckl_set_point failed"
        print *, "Error code", rc
        stop 1
    else
        print *, "qmckl_set_point succeeded"
    end if

    if (gl.eq.0) then 
        !$omp target data map(from:values(1:multiplicity))

        !$omp target data use_device_ptr(values)
        rc = qmckl_get_ao_basis_ao_value_device(qmckl_gpu_ctx&
                                              &, c_loc(values)&
                                              &, 1_8 * multiplicity)
        !$omp end target data
        !$omp end target data
        if (rc /= 0) then
            print *, "qmckl_get_ao_basis_ao_value failed"
            print *, "Error code", rc
            stop 1
        else
            print *, "qmckl_get_ao_basis_ao_value succeeded"
        end if
    else
        !$omp target data map(from:values(1:multiplicity*5))

        !$omp target data use_device_ptr(values)
        rc = qmckl_get_ao_basis_ao_vgl_device(qmckl_gpu_ctx&
                                            &, c_loc(values)&
                                            &, 5_8 * multiplicity)
        !$omp end target data
        !$omp end target data
        if (rc /= 0) then
            print *, "qmckl_get_ao_basis_ao_vgl failed"
            print *, "Error code", rc
            stop 1
        else
            print *, "qmckl_get_ao_basis_ao_vgl succeeded"
        end if
    end if

    !$omp end target data

    call qmckl_gpu_finalize()

#endif

end subroutine provider_qmckl_gpu
