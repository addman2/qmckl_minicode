program makefun_tester

    implicit none

    integer*4                      :: iorb, indt, indpar, indorb, i0, iflagnorm_unused        &
           &, indtm, typec, nelskip, num_lines, ii
    real*8, dimension(:), allocatable :: dd, zeta, r
    real*8, dimension(:,:), allocatable :: z, rmu, distp

    integer*4, dimension(:), allocatable :: iorbs, multiplicities, npars
    character*80                   :: dummy_1, dummy_2, dummy_3, dummy_4

    indtm = 5
    indt = 10

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
    allocate(multiplicities(num_lines))
    allocate(npars(num_lines))
    ! Rewind file
    rewind(10)
    ! Skip header
    read(10,*)
    ! Read data
    do ii = 1, num_lines
        read(10,*) iorbs(ii), dummy_1, dummy_2, dummy_3, dummy_4, multiplicities(ii), npars(ii)
    end do
    
    do ii = 1, num_lines
        write(*,*) iorbs(ii), multiplicities(ii), npars(ii)
    
        if (allocated(dd)) deallocate(dd)
        allocate(dd(npars(ii)))
    
        if (allocated(rmu)) deallocate(rmu)
        allocate(rmu(3,0:indtm))
    
        if (allocated(r)) deallocate(r)
        allocate(r(0:indtm))
    
        if (allocated(z)) deallocate(z)
        allocate(z(multiplicities(ii),0:indt+4))
    
        ! Generate random points
        call random_number(rmu)
        rmu = rmu * 2.0d0 - 1.0d0
    
        r = sqrt(sum(rmu**2, dim=1))
    
        ! Randomize parameters
        call random_number(dd)

        i0 = 0
        indtmin = 0
        typec = 0
        indpar = 0
        indorb = 0
        indshell = 0
    
        call makefun(iorbs(ii),indt,i0,indtmin,indtm,typec,indpar      &
               &,indorb,indshell,nelskip,z,dd,zeta,r,rmu,distp    &
               &,iflagnorm_unused,cr)
    
        deallocate(z)
        deallocate(r)
        deallocate(dd)
        deallocate(rmu)
    end do
    
    allocate(zeta(1))
    allocate(r(1))
    
    if (allocated(iorbs)) deallocate(iorbs)
    if (allocated(multiplicities)) deallocate(multiplicities)
    if (allocated(npars)) deallocate(npars)
    
    if (allocated(zeta)) deallocate(zeta)
    if (allocated(r)) deallocate(r)

contains


subroutine test

end subroutine test

end program makefun_tester
