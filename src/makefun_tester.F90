program makefun_tester

    implicit none

    integer*4 :: iorb, indt, indpar, indorb, i0, iflagnorm_unused&
            &, indtm, typec, nelskip
    real*8, dimension(:), allocatable :: dd, zeta, r
    real*8, dimension(:,:), allocatable :: z, rmu, distp

    iorb = 16_4

    allocate(dd(1))
    allocate(zeta(1))
    allocate(r(1))

    if (allocated(dd)) deallocate(dd)
    if (allocated(zeta)) deallocate(zeta)
    if (allocated(r)) deallocate(r)

contains


    subroutine test

    end subroutine test

end program makefun_tester
