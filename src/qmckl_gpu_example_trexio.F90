#ifdef _OLD_QMCKL_GPU_INTERFACE
#define PTR_C(argument) c_loc(argument)
#else
#define PTR_C(argument) argument
#endif

program qmckl_gpu_example_trexio

    use qmckl_gpu_helpers
    use qmckl_helpers

    implicit none

    integer(kind=4) :: rc, ii, num_points, p_
    integer(kind=8), target :: num_ao
    character(len=256), target :: trexiofile

    real(kind=8), target, dimension(:,:), allocatable :: points&
                                                      &, values_cpu&
                                                      &, values_gpu

#ifdef _QMCKL_GPU

    num_points = 3
    p_ = 1

!#####################################################
!#                                                   #
!#  Initialize QMCkl and QMCkl GPU context           #
!#                                                   #
!#####################################################

    call qmckl_init()

    if (qmckl_ctx .eq. 0) then
        print *, "Error creating QMCkl context"
        stop 1
    else
        print *, "QMCkl context created"
    end if

    call qmckl_gpu_init()

    if (qmckl_gpu_ctx .eq. 0) then
        print *, "Error creating QMCkl GPU context"
        stop 1
    else
        print *, "QMCkl GPU context created"
    end if

!#####################################################
!#                                                   #
!#  Read TREXIO file                                 #
!#                                                   #
!#####################################################

    trexiofile = "h.hdf5"

    rc = qmckl_trexio_read(&
                         &  qmckl_ctx&
                         &, trexiofile&
                         &, 1_8*len(trim(trexiofile)))

    if (rc .ne. 0) then
        print *, "Error reading TREXIO file"
        stop 1
    else
        print *, "TREXIO file read"
    end if

    rc = qmckl_trexio_read_device(&
                               &  qmckl_gpu_ctx&
                               &, c_loc(trexiofile)&
                               &, 1_8*len(trim(trexiofile)))

    if (rc .ne. 0) then
        print *, "Error reading TREXIO file"
        stop 1
    else
        print *, "TREXIO file read"
    end if

!#####################################################
!#                                                   #
!#  Set points                                       #
!#                                                   #
!#####################################################

    allocate(points(3, num_points))

    !call random_number(points)
    points(1, :) = 0.0d0
    points(2, :) = 0.0d0
    points(3, :) = 0.7d0

    print *, "Points:"
    do ii = 1, num_points
        print *, "point ", ii
        print *, points(:,ii)
    end do

    rc = qmckl_set_point(&
                     &  qmckl_ctx&
                     &, "N"&
                     &, 1_8 * p_&
                     &, points&
                     &, 3_8*num_points)

    if (rc .ne. 0) then
        print *, "qmckl_set_point failed"
        print *, "Error code", rc
        stop 1
    else
        print *, "qmckl_set_point succeeded"
    end if

    call setup_points(num_points, "nat", points)

    deallocate(points)

!#####################################################
!#                                                   #
!#  Get values                                       #
!#                                                   #
!#####################################################

    rc = qmckl_get_ao_basis_ao_num_device(&
                                       &  qmckl_gpu_ctx&
                                       &, PTR_C(num_ao))

    if (rc .ne. 0) then
        print *, "qmckl_get_ao_basis_ao_num failed"
        print *, "Error code", rc
        stop 1
    else
        print *, "qmckl_get_ao_basis_ao_num succeeded"
    end if

    print *, "Number of AO", num_ao

    allocate(values_gpu(num_ao,num_points))
    allocate(values_cpu(num_ao,num_points))

    values_gpu = 0.0d0
    values_cpu = 0.0d0

    rc = qmckl_get_ao_basis_ao_value_inplace(&
                                          &  qmckl_ctx&
                                          &, values_cpu&
                                          &, 1_8 * num_ao * num_points)

    if (rc .ne. 0) then
        print *, "qmckl_get_ao_basis_ao_value failed"
        print *, "Error code", rc
        stop 1
    else
        print *, "qmckl_get_ao_basis_ao_value succeeded"
    end if

    !$omp target data map(tofrom:values_gpu)

    !$omp target
    values_gpu = 0.0d0
    !$omp end target

    !$omp target data use_device_ptr(values_gpu)
    rc = qmckl_get_ao_basis_ao_value_device(&
                                         &  qmckl_gpu_ctx&
                                         &, PTR_C(values_gpu)&
                                         &, 1_8 * num_ao * num_points)
    !$omp end target data
    if (rc .ne. 0) then
        print *, "qmckl_get_ao_basis_ao_value failed"
        print *, "Error code", rc
        stop 1
    else
        print *, "qmckl_get_ao_basis_ao_value succeeded"
    end if

    !$omp end target data

    print *, "AO value:"
    do ii = 1, num_points
        print *, "point ", ii
        write(*, '("H1 cpu: ", 1F6.3, " | ", 3F6.3)') values_cpu(1:4,ii)
        write(*, '("H1 gpu: ", 1F6.3, " | ", 3F6.3)') values_gpu(1:4,ii)
        write(*, '("H2 cpu: ", 1F6.3, " | ", 3F6.3)') values_cpu(5:8,ii)
        write(*, '("H2 gpu: ", 1F6.3, " | ", 3F6.3)') values_gpu(5:8,ii)
    end do

    deallocate(values_cpu)
    deallocate(values_gpu)

    call qmckl_finalize()
    call qmckl_gpu_finalize()

#endif

contains


    subroutine setup_points(num_points, typ, points)

        integer(kind=4), intent(in) :: num_points
        character(len=*), intent(in) :: typ
        real(kind=8), dimension(:,:), intent(in) :: points(3, num_points)

        type(c_ptr) :: points_device

        ! I forgot what type is context therefore
        ! I am passing ti from main program

        if (trim(typ).eq."omp") then
            !$omp target data map(to:points)
            !$omp target data use_device_ptr(points)
            rc = qmckl_set_point_device(&
                                     &  qmckl_gpu_ctx&
                                     &, "N"&
                                     &, 1_8 * p_&
                                     &, PTR_C(points)&
                                     &, 3_8*num_points)
            !$omp end target data
            !$omp end target data

        else if (trim(typ).eq."host") then
            rc =  qmckl_set_point_device_from_host(&
                                                &  qmckl_gpu_ctx&
                                                &, "N"&
                                                &, 1_8 * p_&
                                                &, points&
                                                &, 3_8*num_points)
        else if (trim(typ).eq."nat") then
            points_device = qmckl_malloc_device(qmckl_gpu_ctx, 8_8*3*num_points)

            if (points_device .eq. c_null_ptr) then
                print *, "Error allocating device memory"
                stop 1
            else
                print *, "Device memory allocated"
            end if

            rc = qmckl_memcpy_H2D_double(qmckl_gpu_ctx, points_device, points, 8_8*3*num_points)

            if (rc .ne. 0) then
                print *, "Error copying data to device"
                stop 1
            else
                print *, "Data copied to device"
            end if

            rc = qmckl_set_point_device(&
                                     &  qmckl_gpu_ctx&
                                     &, "N"&
                                     &, 1_8 * p_&
                                     &, points_device&
                                     &, 3_8*num_points)

            if (rc .ne. 0) then
                print *, "qmckl_set_point failed"
                print *, "Error code", rc
                stop 1
            else
                print *, "qmckl_set_point succeeded"
            end if

            rc = qmckl_free_device(qmckl_gpu_ctx, points_device)

            if (rc .ne. 0) then
                print *, "Error freeing device memory"
                stop 1
            else
                print *, "Device memory freed"
            end if

        else
            write(*,*) "Unknown type"
            stop 1
        end if

        if (rc .ne. 0) then
            print *, "qmckl_set_point failed"
            print *, "Error code", rc
            stop 1
        else
            print *, "qmckl_set_point succeeded"
        end if

    end subroutine setup_points

end program qmckl_gpu_example_trexio
