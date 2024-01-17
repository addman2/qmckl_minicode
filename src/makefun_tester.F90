program makefun_tester

    implicit none

    integer*4 :: iorb, indt, indpar, indorb, i0, iflagnorm_unused&
            &, indtm, typec, nelskip, num_lines
    real*8, dimension(:), allocatable :: dd, zeta, r
    real*8, dimension(:,:), allocatable :: z, rmu, distp

    integer*4, dimension(:), allocatable :: iorbs


    ! Load parameters form file parameters.csv
    open(unit=10, file='parameters.csv', status='old', action='read')
    ! Skip header
    read(10,*)
    ! Count number of lines
    num_lines = 0
    do
        read(10,*,iostat=iflagnorm_unused)
        if (iflagnorm_unused /= 0) exit
        num_lines = num_lines + 1
    end do
    ! Allocate arrays
    allocate(iorbs(num_lines))
    ! Rewind file
    rewind(10)
    ! Skip header
    read(10,*)
    ! Read data
    do iorb = 1, num_lines
        read(10,*) iorbs(iorb)
    end do

    print *, iorbs

    allocate(dd(1))
    allocate(zeta(1))
    allocate(r(1))

    if (allocated(iorbs)) deallocate(iorbs)

    if (allocated(dd)) deallocate(dd)
    if (allocated(zeta)) deallocate(zeta)
    if (allocated(r)) deallocate(r)

contains


    subroutine test

    end subroutine test

end program makefun_tester
