#define MAX_FILES 65536
#define MAX_FILE_LENGTH 30

program makefun_tester

    implicit none

    interface
        subroutine list_directory(path, path_length, files, num_files, prefix, prefix_length) bind(C, name="list_directory")
            use, intrinsic :: iso_c_binding
            character(c_char), intent(in) :: path(*)
            character(c_char), intent(in) :: prefix(*)
            character(c_char), intent(out) :: files(MAX_FILE_LENGTH, MAX_FILES)
            integer(c_int), intent(in), value :: prefix_length
            integer(c_int), intent(in), value :: path_length
            integer(c_int), intent(out) :: num_files
        end subroutine list_directory
    end interface

    integer*4 :: iorb, indt, indpar, indorb, indshell, i0, iflagnorm_unused&
              &, indtm, indtmin, typec, nelskip, num_lines, ii
    real*8 :: cr 
    real*8, dimension(:), allocatable :: dd, zeta, r
    real*8, dimension(:,:), allocatable :: z, rmu, distp

    integer*4, dimension(:), allocatable :: iorbs, multiplicities, npars, failed_test
    character*80 :: dummy_1, dummy_2, dummy_3, dummy_4

    character*12000 :: error_message = "", tmp

    logical :: failed

    call initialize()

    cr = 0.0d0
    iflagnorm_unused = 0
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
    allocate(failed_test(num_lines))

    ! -1 means that test not failed
    failed_test = -1

    ! Rewind file
    rewind(10)
    ! Skip header
    read(10,*)
    ! Read data
    do ii = 1, num_lines
        read(10,*) iorbs(ii), dummy_1, dummy_2, dummy_3, dummy_4, multiplicities(ii), npars(ii)
    end do
    close(10)

    multiplicities(2) = multiplicities(2) + 1
    multiplicities(3) = multiplicities(3) + 1
    multiplicities(4) = multiplicities(4) + 1
    
    do ii = 1, num_lines
        failed = .false.
        write(*,'("Checking orbital index ", I3)') iorbs(ii)
        write(*, *)
        write(*, *) "random_string(10) = ", random_string(10)
    
        if (allocated(dd)) deallocate(dd)
        allocate(dd(npars(ii)))
    
        if (allocated(rmu)) deallocate(rmu)
        allocate(rmu(3,0:indtm))
    
        if (allocated(r)) deallocate(r)
        allocate(r(0:indtm))
    
        if (allocated(z)) deallocate(z)
        allocate(z(multiplicities(ii),0:indt+4))

        if (allocated(distp)) deallocate(distp)
        allocate(distp(0:indt+4,20))
    
        !call check_index_movement()
        !call create_single_value()
        call test_single_value()

        if (failed) then
            failed_test(ii) = 1
            tmp = trim(error_message)
            write(error_message, '(A,A,I3.3)') trim(tmp), ",", iorbs(ii)
        end if
    
        deallocate(distp)
        deallocate(z)
        deallocate(r)
        deallocate(dd)
        deallocate(rmu)

        write(*, *)
        write(*, '("Result = ", L1)') .not. failed
        write(*, *)
        write(*, '("######################################")')
    end do
    
    allocate(zeta(1))
    allocate(r(1))
    
    if (allocated(iorbs)) deallocate(iorbs)
    if (allocated(multiplicities)) deallocate(multiplicities)
    if (allocated(npars)) deallocate(npars)
    if (allocated(failed_test)) deallocate(failed_test)
    
    if (allocated(zeta)) deallocate(zeta)
    if (allocated(r)) deallocate(r)

    print *, "Failed tests: ", trim(error_message)
    if (len_trim(error_message) > 0) then
        stop 1
    else
        stop 0
    end if

contains

function random_string(length) result(str)

    implicit none

    integer, intent(in) :: length
    character(len=length) :: str
    integer :: i
    real :: r
    character(len=29) :: letters = 'abcdefghijklmnopqrstuvwxyz148'

    do i = 1, length
        ! Generate random number between 1 and 29
        call random_number(r)
        r = r * 29 + 1
        str(i:i) = letters(int(r):int(r))
    end do

end function random_string

subroutine initialize

    implicit none
    logical :: data_dir_exists

    ! Initialize random number generator
    call random_seed()

    ! Check if data directory exists, and create it if not

    inquire(file='data', exist=data_dir_exists)
    if (.not. data_dir_exists) then
        ! Create directory TODO: make it cross-platform
        call system('mkdir data')
    end if

end subroutine initialize

subroutine run_random

    implicit none

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
    z = 0.0d0
    
    call makefun(iorbs(ii),indt,i0,indtmin,indtm,typec,indpar      &
               &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
               &,iflagnorm_unused,cr)

    print *, "indorb = ", indorb
    print *, "indshell = ", indshell
    print *, "indpar = ", indpar

end subroutine run_random

subroutine check_index_movement

    implicit none

    failed = .false.

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
    z = 0.0d0
    
    call makefun(iorbs(ii),indt,i0,indtmin,indtm,typec,indpar      &
               &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
               &,iflagnorm_unused,cr)

    write(*,'("indorb  = ",I3)', advance="no") indorb
    if (indorb /= multiplicities(ii)) then
        write(*,'(" FAILED")')
        failed = .true.
    else
        write(*,'(" OK")')
    end if

    write(*,'("indshell = ",I3)', advance="no") indshell
    if (indshell /= multiplicities(ii)) then
        write(*,'(" FAILED")')
        failed = .true.
    else
        write(*,'(" OK")')
    end if

    write(*,'("indpar  = ",I3)', advance="no") indpar
    if (indpar /= npars(ii)) then
        write(*,'(" FAILED")')
        failed = .true.
    else
        write(*,'(" OK")')
    end if

end subroutine check_index_movement

subroutine create_single_value

    implicit none
    character*80 :: filename

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
    typec = 1
    z = 0.0d0

    call makefun(iorbs(ii),indt,i0,indtmin,0,typec,indpar      &
               &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
               &,iflagnorm_unused,cr)

    ! In a binary file data/sv.dat we store the following data:
    ! 1. iorbs(ii)
    ! 2. dd array
    ! 3. rmu array

    write(filename, '(A,I3.3,A,A,A,A)') 'data/sv_', iorbs(ii), "_", random_string(15), '.dat'
    open(unit=20, file=filename, form='unformatted', access='sequential', status='new', action='write')
    write(20) dd
    write(20) rmu(:,0)
    write(20) z(0:multiplicities(ii), 0)
    close(20)

end subroutine create_single_value

subroutine test_single_value

    implicit none

    character :: files(MAX_FILE_LENGTH, MAX_FILES)
    character*MAX_FILE_LENGTH :: filename
    integer :: num_files, i, j, iorb_test, ios
    real :: z_test(0:multiplicities(ii), 0)

    call list_directory('data', 4, files, num_files, 'sv_', 3)
    print *, "num_files = ", num_files
    
    do i = 1, num_files
        do j = 1, MAX_FILE_LENGTH
            filename(j:j) = files(j,i)
        end do
        write(*, '(I3," ",A)') i, adjustl(trim(filename))
        
        read(filename(4:6), '(I3)', iostat=ios ) iorb_test
        if (ios /= 0) then
            continue
        end if
        if (iorb_test /= iorbs(ii)) then
            continue
        end if
        open(unit=20, file='data/'//trim(filename), form='unformatted', access='sequential', status='old', action='read')
        do j = 1, 1
            read(20) dd(j)
        end do
        read(20) rmu(1,0)
        read(20) rmu(2,0)
        read(20) rmu(3,0)
        !read(20) z_test(0:multiplicities(ii), 0)
        close(20)

        !write(*, '(A)') "dd = ", dd

        !i0 = 0
        !indtmin = 0
        !typec = 0
        !indpar = 0
        !indorb = 0
        !indshell = 0
        !z = 0.0d0

        !call makefun(iorbs(ii),indt,i0,indtmin,0,typec,indpar      &
        !           &,indorb,indshell,multiplicities(ii),z,dd,zeta,r,rmu,distp    &
        !           &,iflagnorm_unused,cr)



        !write(*, '(A)') "z_test = ", z_test
        !write(*, '(A)') "z = ", z(0:multiplicities(ii), 0)

    end do

end subroutine test_single_value

end program makefun_tester
