program qmckl_reader

    use qmckl_helpers
    use, intrinsic :: iso_c_binding

    implicit none

    character(len=1024) :: args
    character(len=1) :: ao_basis_type

    integer :: i, count

    integer*8, parameter :: max_size = 128
    integer*8 :: ao_basis_ao_num
    integer*8 :: ao_basis_shell_num
    integer*8 :: ao_basis_prim_num
    integer*8 :: p3
    integer*8 :: ao_basis_nucleus_shell_num(max_size)
    integer*8 :: ao_basis_nucleus_index(max_size)
    integer*4 :: ao_basis_shell_ang_mom(max_size)
    integer*8 :: ao_basis_shell_prim_num(max_size)
    integer*8 :: ao_basis_shell_prim_index(max_size)
    integer*4 :: rc

    real*8 :: ao_basis_shell_factor(max_size)
    real*8 :: ao_basis_exponent(max_size)
    real*8 :: ao_basis_coefficient(max_size)
    real*8 :: ao_basis_prim_factor(max_size)
    real*8 :: point(1,3)
    real*8, dimension(:), allocatable :: values
    real(c_double), allocatable, dimension(:) :: ao_basis_ao_factor
    
    call qmckl_init()

    ao_basis_nucleus_shell_num = 0
    ao_basis_nucleus_index = 0
    ao_basis_shell_ang_mom = 0
    ao_basis_shell_prim_num = 0
    ao_basis_shell_prim_index = 0
    ao_basis_shell_factor = 0.0_8
    ao_basis_exponent = 0.0_8
    ao_basis_coefficient = 0.0_8
    ao_basis_prim_factor = 0.0_8


    call get_command_argument(1,args)

    rc = qmckl_trexio_read(qmckl_ctx&
                        &, trim(adjustl(args))&
                        &, 1_8 * len(trim(adjustl(args))))
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to load trexio file"
        stop 1
    else
        write (6,*) "Trexio file loaded"
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_type(qmckl_ctx, ao_basis_type)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis type"
        stop 1
    else
        write (6,*) "AO basis type: ", ao_basis_type
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_shell_num(qmckl_ctx, ao_basis_shell_num)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis shell num"
        stop 1
    else
        write (6,*) "AO basis shell num: ", ao_basis_shell_num
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_prim_num(qmckl_ctx, ao_basis_prim_num)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis prim num"
        stop 1
    else
        write (6,*) "AO basis prim num: ", ao_basis_prim_num
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_nucleus_shell_num(qmckl_ctx&
                                            &, ao_basis_nucleus_shell_num&
                                            &, max_size)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis nucleus shell num"
        stop 1
    else
        write (6,*) "AO basis nucleus shell num: "
        do i = 1, max_size, 8
            write (6, '((i3): (8i5))') i, ao_basis_nucleus_shell_num(i:i+7)
        end do
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_nucleus_index(qmckl_ctx&
                                        &, ao_basis_nucleus_index&
                                        &, max_size)

    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis nucleus index"
        stop 1
    else
        write (6,*) "AO basis nucleus index: "
        do i = 1, max_size, 8
            write (6, '((i3),":", (8i5))') i, ao_basis_nucleus_index(i:i+7)
        end do
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_shell_ang_mom(qmckl_ctx&
                                        &, ao_basis_shell_ang_mom&
                                        &, max_size)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis shell ang mom"
        stop 1
    else
        write (6,*) "AO basis shell ang mom: "
        do i = 1, max_size, 8
            write (6, '((i3),":", (8i5))') i, ao_basis_shell_ang_mom(i:i+7)
        end do
    end if

    print *, ""
    print *, "###############################################"
    print *, ""
    
    rc = qmckl_get_ao_basis_shell_prim_num(qmckl_ctx&
                                        &, ao_basis_shell_prim_num&
                                        &, max_size)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis shell prim num"
        stop 1
    else
        write (6,*) "AO basis shell prim num: "
        do i = 1, max_size, 8
            write (6, '((i3),":", (8i5))') i, ao_basis_shell_prim_num(i:i+7)
        end do
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_shell_prim_index(qmckl_ctx&
                                            &, ao_basis_shell_prim_index&
                                            &, max_size)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis shell prim index"
        stop 1
    else
        write (6,*) "AO basis shell prim index: "
        do i = 1, max_size, 8
            write (6, '((i3),":", (8i5))') i, ao_basis_shell_prim_index(i:i+7)
        end do
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_shell_factor(qmckl_ctx&
                                    &, ao_basis_shell_factor&
                                    &, max_size)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis shell factor"
        stop 1
    else
        write (6,*) "AO basis shell factor: "
        do i = 1, max_size, 8
            write (6, '((i3),":", (8f8.3))') i, ao_basis_shell_factor(i:i+7)
        end do
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_exponent(qmckl_ctx&
                                   &, ao_basis_exponent&
                                   &, max_size)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis exponent"
        stop 1
    else
        write (6,*) "AO basis exponents: "
        do i = 1, max_size, 8
            write (6, '((i3),":", (8f8.3))') i, ao_basis_exponent(i:i+7)
        end do
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_coefficient(qmckl_ctx&
                                       &, ao_basis_coefficient&
                                       &, max_size)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis coefficient"
        stop 1
    else
        write (6,*) "AO basis coefficient: "
        do i = 1, max_size, 8
            write (6, '((i3),":", (8f8.3))') i, ao_basis_coefficient(i:i+7)
        end do
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_prim_factor(qmckl_ctx&
                                       &, ao_basis_prim_factor&
                                       &, max_size)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis prim factor"
        stop 1
    else
        write (6,*) "AO basis prim factor: "
        do i = 1, max_size, 8
            write (6, '((i3),":", (8f8.3))') i, ao_basis_prim_factor(i:i+7)
        end do
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    rc = qmckl_get_ao_basis_ao_num(qmckl_ctx, ao_basis_ao_num)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis ao num"
        stop 1
    else
        write (6,*) "AO basis ao num: ", ao_basis_ao_num
    end if

    print *, ""
    print *, "###############################################"
    print *, ""

    allocate(ao_basis_ao_factor(ao_basis_ao_num))

    rc = qmckl_get_ao_basis_ao_factor(qmckl_ctx&
                                    &, ao_basis_ao_factor&
                                    &, ao_basis_ao_num)
    if (rc.ne.QMCKL_SUCCESS) then
        write (0, *) "Unable to get ao basis ao factor"
        write (0, *) "Parameter 3: ", p3
        write (0, *) "Error code: ", rc
        stop 1
    else
        write (6,*) "AO basis ao factor: "
        do i = 1, ao_basis_ao_num, 8
            write (6, '((i3),":", (8f8.3))') i, ao_basis_ao_factor(i:i+7)
        end do
    end if

    ! Count non zero values in ao_basis_ao_factor
    count = 0
    do i = 1, ao_basis_ao_num
        if (abs(ao_basis_ao_factor(i)) > 0.00001_8) then
            count = count + 1
        end if
    end do

    print *, ""
    print *, " Number of non zero values in ao_basis_ao_factor: ", count
    print *, ""

    point(1,:) = (/ 0.1_8, 1.0_8, 1.0_8 /)

    print *, "Point", point

    rc = qmckl_set_point(qmckl_ctx&
                       &, "N"&
                       &, ao_basis_ao_num&
                       &, point&
                       &, 3_8 * size(point))

    if (rc /= 0) then
        print *, "qmckl_set_point failed"
        print *, "Error code", rc
        stop 1
    else
        print *, "qmckl_set_point succeeded"
    end if
                   
    allocate(values(ao_basis_ao_num*10))
    values = 0.0_8

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

    print *, "AO values:"
    do i = 1, ao_basis_ao_num, 8
        write (6, '((i3),":", (8f8.3))') i, values(i:i+7)
    end do

    deallocate(values)

    deallocate(ao_basis_ao_factor)

    call qmckl_finalize()

end program qmckl_reader
