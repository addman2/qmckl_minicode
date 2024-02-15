module qmckl_gpu_helpers

    use qmckl_gpu_f
    use qmckl
    use iso_c_binding

    implicit none
#ifdef _QMCKL_GPU
    integer(kind=qmckl_context_device) :: qmckl_gpu_ctx
    integer(c_int32_t) :: qmckl_gpu_device_id

contains

    subroutine qmckl_gpu_init()
        use qmckl_gpu_f

        implicit none

        integer(kind=4) :: rc

        qmckl_gpu_ctx = qmckl_context_create_device(qmckl_gpu_device_id)
        if (qmckl_gpu_ctx.eq.0_8) then
            write (0,*) "QMCKL GPU context is a null pointer, but it should never happen"
            stop 1
        else
            write (6,*) "QMCKL GPU context created"
        end if

    end subroutine qmckl_gpu_init

    subroutine qmckl_gpu_finalize()

        use qmckl_gpu_f

        implicit none

        integer(kind=4) :: rc

        rc = qmckl_context_destroy_device(qmckl_gpu_ctx)
        if (rc.ne.QMCKL_SUCCESS) then
            write (0, *) "Unable to destroy QMCkl GPU context"
        else
            write (6, *) "QMCKL GPU context destroyed"
        end if

    end subroutine qmckl_gpu_finalize
#endif

end module qmckl_gpu_helpers

