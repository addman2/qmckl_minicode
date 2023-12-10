module qmckl_helpers

    use qmckl

    implicit none

    integer(kind=8) :: qmckl_ctx

contains

    subroutine qmckl_init()
        use qmckl

        implicit none

        integer(kind=qmckl_exit_code) :: rc

        qmckl_ctx = qmckl_context_create()
        if (qmckl_ctx.eq.0_8) then
            write (0,*) "QMCKL context is a null pointer, but it should never happen"
            stop 1
        else
            write (6,*) "QMCKL context created"
        end if
        
    end subroutine qmckl_init

    subroutine qmckl_finalize()
        use qmckl

        implicit none

        integer(kind=qmckl_exit_code) :: rc

        rc = qmckl_context_destroy(qmckl_ctx)
        if (rc.ne.QMCKL_SUCCESS) then
            write (0, *) "Unable to destroy QMCkl context"
        end if

    end subroutine qmckl_finalize
end module qmckl_helpers

