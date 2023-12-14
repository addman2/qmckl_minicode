program qmckl_example_scratch

    use qmckl_helpers

    real(kind=8), dimension(1) :: charges
    real(kind=8), dimension(3,1) :: positions
    real(kind=8), dimension(1) :: shell_factor
    real(kind=8), dimension(1) :: exponent
    real(kind=8), dimension(1) :: coefficient
    real(kind=8), dimension(1) :: prim_factor
    real(kind=8), dimension(1) :: ao_factor
    real(kind=8), dimension(:), allocatable :: values
    real(kind=8), dimension(1,3) :: point

    integer(kind=8), dimension(1) :: nucleus_shell_num
    integer(kind=8), dimension(1) :: nucleus_index
    integer(kind=4), dimension(1) :: shell_ang_mom
    integer(kind=8), dimension(1) :: shell_prim_num
    integer(kind=8), dimension(1) :: shell_prim_index

    logical*1 :: T = .true.

    logical :: rcb

    call qmckl_init()

    rc = qmckl_set_nucleus_num(qmckl_ctx, 1_8)
    if (rc /= 0) then
        print *, "qmckl_set_nucleus_num failed"
        stop 1
    else
        print *, "qmckl_set_nucleus_num succeeded"
    end if

    charges(1) = 1.0_8
    rc = qmckl_set_nucleus_charge(qmckl_ctx, charges, 1_8 * size(charges))
    if (rc /= 0) then
        print *, "qmckl_set_nucleus_charge failed"
        stop 1
    else
        print *, "qmckl_set_nucleus_charge succeeded"
    end if

    positions(:,1) = (/ 0.0_8, 0.0_8, 0.0_8 /)
    rc = qmckl_set_nucleus_coord(qmckl_ctx, "N", positions, 3_8 * size(positions))
    if (rc /= 0) then
        print *, "qmckl_set_nucleus_cood failed"
        stop 1
    else
        print *, "qmckl_set_nucleus_cood succeeded"
    end if

    !rcb = qmckl_nucleus_provided(qmckl_ctx)
    !print *, "QMCkl nucleus provided", qmckl_nucleus_provided(qmckl_ctx)

!##############################################################################
!#                                                                            #
!#  Set up the basis set.                                                     #
!#                                                                            #
!##############################################################################

    rc = qmckl_set_ao_basis_type(qmckl_ctx, "G")
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_type failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_type succeeded"
    end if

    rc = qmckl_set_ao_basis_shell_num(qmckl_ctx, 1_8)
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_shell_num failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_shell_num succeeded"
    end if

    rc = qmckl_set_ao_basis_prim_num(qmckl_ctx, 1_8)
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_prim_num failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_prim_num succeeded"
    end if

    nucleus_shell_num(1) = 1
    rc = qmckl_set_ao_basis_nucleus_shell_num(qmckl_ctx, nucleus_shell_num, 1_8 * size(nucleus_shell_num))
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_nucleus_shell_num failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_nucleus_shell_num succeeded"
    end if

    nucleus_index(1) = 0
    rc = qmckl_set_ao_basis_nucleus_index(qmckl_ctx, nucleus_index, 1_8 * size(nucleus_index))
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_nucleus_index failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_nucleus_index succeeded"
    end if

    shell_ang_mom(1) = 0
    rc = qmckl_set_ao_basis_shell_ang_mom(qmckl_ctx, shell_ang_mom, 1_8 * size(shell_ang_mom))
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_shell_ang_mom failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_shell_ang_mom succeeded"
    end if

    shell_prim_num(1) = 1
    rc = qmckl_set_ao_basis_shell_prim_num(qmckl_ctx, shell_prim_num, 1_8 * size(shell_prim_num))
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_shell_prim_num failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_shell_prim_num succeeded"
    end if

    shell_prim_index(1) = 1
    rc = qmckl_set_ao_basis_shell_prim_index(qmckl_ctx, shell_prim_index, 1_8 * size(shell_prim_index))
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_shell_prim_index failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_shell_prim_index succeeded"
    end if

    shell_factor(1) = 1.0_8
    rc = qmckl_set_ao_basis_shell_factor(qmckl_ctx, shell_factor, 1_8 * size(shell_factor))
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_shell_factor failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_shell_factor succeeded"
    end if

    exponent(1) = 1.0_8
    rc = qmckl_set_ao_basis_exponent(qmckl_ctx, exponent, 1_8 * size(exponent))
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_exponent failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_exponent succeeded"
    end if

    coefficient(1) = 1.0_8
    rc = qmckl_set_ao_basis_coefficient(qmckl_ctx, coefficient, 1_8 * size(coefficient))
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_coefficient failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_coefficient succeeded"
    end if

    prim_factor(1) = 1.0_8
    rc = qmckl_set_ao_basis_prim_factor(qmckl_ctx, prim_factor, 1_8 * size(prim_factor))
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_prim_factor failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_prim_factor succeeded"
    end if

    rc = qmckl_set_ao_basis_ao_num(qmckl_ctx, 1_8 * size(ao_factor))
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_ao_num failed"
        stop 1
    else
        print *, "qmckl_set_ao_basis_ao_num succeeded"
    end if

    ao_factor = 1.0_8
    rc = qmckl_set_ao_basis_ao_factor(qmckl_ctx, ao_factor, 1_8 * size(ao_factor))
    if (rc /= 0) then
        print *, "qmckl_set_ao_basis_ao_factor failed"
        print *, "Error code", rc
        stop 1
    else
        ! Print success in green
        print *, "qmckl_set_ao_basis_ao_factor succeeded"
    end if

    !rc = qmckl_set_ao_basis_cartesian(qmckl_ctx, T)
    !if (rc /= 0) then
    !    print *, "qmckl_set_ao_basis_cartesian failed"
    !    print *, "Error code", rc
    !    stop 1
    !else
    !    print *, "qmckl_set_ao_basis_cartesian succeeded"
    !end if

    !rc = qmckl_get_ao_basis_cartesian(qmckl_ctx, T)
    !if (rc /= 0) then
    !    print *, "qmckl_get_ao_basis_cartesian failed"
    !    print *, "Error code", rc
    !    stop 1
    !else
    !    print *, "qmckl_get_ao_basis_cartesian succeeded"
    !    print *, "Cartesian", T
    !end if

    !##############################################################################
    !#                                                                            #
    !#  Now calculate values of the AO basis                                      #
    !#                                                                            #
    !##############################################################################

    point(1,:) = (/ 0.1_8, 1.0_8, 1.0_8 /)

    print *, "Point", point

    rc = qmckl_set_point(qmckl_ctx&
                       &, "N"&
                       &, 3_8&
                       &, point&
                       &, 10_8 * size(point))

    if (rc /= 0) then
        print *, "qmckl_set_point failed"
        print *, "Error code", rc
        stop 1
    else
        print *, "qmckl_set_point succeeded"
    end if
                   
    allocate(values(5))

    rc = qmckl_get_ao_basis_ao_value_inplace(qmckl_ctx&
                                   &, values&
                                   &, 1_8 * size(values))
    if (rc /= 0) then
        print *, "qmckl_get_ao_basis_ao_value failed"
        print *, "Error code", rc
        stop 1
    else
        print *, "qmckl_get_ao_basis_ao_value succeeded"
    end if

    print *, "AO value", values

    deallocate(values)

    call qmckl_finalize()

end program qmckl_example_scratch
