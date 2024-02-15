#ifdef _OLD_QMCKL_GPU_INTERFACE
#define PTR_C(argument) c_loc(argument)
#else
#define PTR_C(argument) argument
#endif

program qmckl_gpu_example_trexio

    use qmckl_gpu_helpers

    implicit none

    integer(kind=4) :: rc, ii, num_points
    integer(kind=8), target :: num_ao
    character(len=256), target :: trexiofile

    real(kind=8), target, dimension(:,:), allocatable :: points, values

#ifdef _QMCKL_GPU

    num_points = 3

    call qmckl_gpu_init()

    if (qmckl_gpu_ctx .eq. 0) then
        print *, "Error creating QMCkl GPU context"
        stop 1
    else
        print *, "QMCkl GPU context created"
    end if

    trexiofile = "h.hdf5"

    rc = qmckl_trexio_read_device(qmckl_gpu_ctx&
                               &, c_loc(trexiofile)&
                               &, 1_8*len(trim(trexiofile)))

    if (rc .ne. 0) then
        print *, "Error reading TREXIO file"
        stop 1
    else
        print *, "TREXIO file read"
    end if

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

    !$omp target data map(to:points)
    !$omp target data use_device_ptr(points)
    rc = qmckl_set_point_device(qmckl_gpu_ctx&
                             &, "N"&
                             &, 1_8&
                             &, PTR_C(points)&
                             &, 3_8*2)
    !$omp end target data

    if (rc .ne. 0) then
        print *, "qmckl_set_point failed"
        print *, "Error code", rc
        stop 1
    else
        print *, "qmckl_set_point succeeded"
    end if

    !$omp end target data

    deallocate(points)

    rc = qmckl_get_ao_basis_ao_num_device(qmckl_gpu_ctx&
                                       &, PTR_C(num_ao))

    if (rc .ne. 0) then
        print *, "qmckl_get_ao_basis_ao_num failed"
        print *, "Error code", rc
        stop 1
    else
        print *, "qmckl_get_ao_basis_ao_num succeeded"
    end if

    print *, "Number of AO", num_ao

    allocate(values(num_ao,num_points))

    values = 0.0d0

    !$omp target data map(tofrom:values)
    !$omp target data use_device_ptr(values)
    rc = qmckl_get_ao_basis_ao_value_device(qmckl_gpu_ctx&
                                         &, PTR_C(values)&
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
        print *, values(:,ii)
    end do

    deallocate(values)

    call qmckl_gpu_finalize()

#endif
end program qmckl_gpu_example_trexio
