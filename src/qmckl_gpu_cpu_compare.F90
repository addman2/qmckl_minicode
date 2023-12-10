program qmckl_gpu_cpu_compare

    use qmckl_helpers
    use qmckl_gpu_helpers

    implicit none

    integer*4 :: rc
    character(len=1024), target :: args
    type(c_ptr) :: helper_ptr
    character(len=1), pointer :: helper_ptr

    qmckl_gpu_device_id = 1

    call qmckl_init()
    call qmckl_gpu_init()

    call get_command_argument(1,args)

    rc = qmckl_trexio_read(qmckl_ctx&
                        &, trim(args)&
                        &, 1_8 * len_trim(args))

    if (rc.ne.QMCKL_SUCCESS) then
        print *, "Error reading file"
        stop
    else
        print *, "File read successfully (CPU)"
    end if

    ! Cast the pointer to a c_ptr

    helper_ptr = c_loc(args(1:1))
    !helper_ptr => helper_ptr(1:1)

    rc = qmckl_trexio_read_device(qmckl_gpu_ctx&
                        &, trim(args)&
                        &, 1_8 * len_trim(args))

    if (rc.ne.QMCKL_SUCCESS) then
        print *, "Error reading file"
        print *, "rc = ", rc
        stop
    else
        print *, "File read successfully (GPU)"
    end if

    call qmckl_gpu_finalize()
    call qmckl_finalize()

end program qmckl_gpu_cpu_compare
