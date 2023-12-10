program test_makefun_qmckl
    use qmckl

    implicit none

#ifdef _QMCKL_

    integer(kind=4), parameter :: OPENMP_TEST = 1
    integer(kind=4), parameter :: QMCKL_TEST = 2
    integer(kind=4), parameter :: QMCKL_GPU_TEST = 3
    integer(kind=4), parameter :: MAKEFUN_TEST = 4

    integer(kind=4) :: ii, operation = 0
    integer(kind=8) :: qmckl_ctx
    integer(kind=qmckl_exit_code) :: rc
    character(len=30), dimension(:), allocatable :: args

    call init()

    select case (operation)
    case (OPENMP_TEST)
        call run_openmp_test()
    case (QMCKL_TEST)
        call run_qmckl
    case (MAKEFUN_TEST)
        call run_makefun
    case default
        write (0,*) "Please provide the operation to test"
    end select

    call finish()

contains

    subroutine init()
            
        implicit none

        qmckl_ctx = qmckl_context_create()
        if (qmckl_ctx.eq.0_8) then
            write (0,*) "QMCKL context is a null pointer, but it should never happen"
            stop 1
        else
            write (6,*) "QMCKL context created"
        end if

        ! Assumes one argument: name of the trexio file

        ! If the argument is not provided, the program will crash
        if (command_argument_count().lt.0) then
            write (0,*) "Please provide arguments"
            stop 1
        else 
            allocate(args(command_argument_count()))
        end if

        ! Load arguments
        do ii = 1, command_argument_count()
            call get_command_argument(ii,args(ii))
        end do

        if (command_argument_count().gt.1) then
          if (args(2).eq.'qmckl') then
              operation = QMCKL_TEST
          end if
          if (args(2).eq.'openmp') then
              operation = OPENMP_TEST
          end if
          if (args(2).eq.'makefun') then
              operation = MAKEFUN_TEST
          end if
        end if

        rc = qmckl_trexio_read(qmckl_ctx&
                            &, trim(adjustl(args(1)))&
                            &, 1_8 * len(trim(adjustl(args(1)))))
        if (rc.ne.QMCKL_SUCCESS) then
            write (0, *) "Unable to load trexio file"
            stop 1
        else
            write (6,*) "Trexio file loaded"
        end if

    end subroutine init

    subroutine finish()

        implicit none

        rc = qmckl_context_destroy(qmckl_ctx)
        if (rc.ne.QMCKL_SUCCESS) then
            write (0, *) "Unable to destroy QMCkl context"
        end if

        if (allocated(args)) deallocate(args)

        write (6,*) "Makefun test completed successfully"

    end subroutine finish

    subroutine run_makefun()

        use makefun_pars

        implicit none

        call setup_parameters(6)
        call randomize_electrons(1)

        !call makefun(16, indt, i0, indtmin, indtm, typec, indpar, indorb, indshell,nelskip, z, parameters, zeta, r, rmu, distp, iflagnorm_unused, cr) 

        call finish_parameters()

    end subroutine run_makefun

    subroutine run_qmckl()

        implicit none
        
        integer(kind=4) :: ii
        integer(kind=8) :: ao_num, mo_num
        real(kind=8), allocatable, dimension(:) :: z
        real(kind=8), allocatable, dimension(:,:) :: electron

        allocate(electron(3,1))

        ! Randomize the electron position
        call random_number(electron)

        rc = qmckl_get_ao_basis_ao_num(qmckl_ctx, ao_num)
        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting ao_num', rc, qmckl_ctx, ao_num
            call abort()
        end if

        write (6,*) "ao_num = ", ao_num

        if (ao_num.eq.0) then
            rc = qmckl_get_mo_basis_mo_num(qmckl_ctx, mo_num)
            if (rc /= QMCKL_SUCCESS) then
                print *, 'Error getting ao_num', rc, qmckl_ctx, mo_num
                call abort()
            end if
        end if

        write (6,*) "mo_num = ", mo_num

        allocate(z(ao_num))

        do ii = 1, ao_num
            z(ii) = 0.0_8
        end do

        rc = qmckl_set_point(qmckl_ctx, 'N', 1_8, electron, 3_8)

        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error setting electron coords'
            call abort()
        end if

        rc = qmckl_get_ao_basis_ao_value_inplace(                        &
               &qmckl_ctx,                                               &
               &z,                                                       &
               &ao_num)

        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting AOs from QMCkl'
            call abort()
        end if

        print *, "z = ", z

        if (allocated(electron)) deallocate(electron) 
        if (allocated(z)) deallocate(z) 

    end subroutine run_qmckl

    subroutine test_makefun_orbital(orb)

        implicit none

        integer(kind=4), intent(in) :: orb

    end subroutine test_makefun_orbital

    subroutine test_qmckl_gpu(orb)

        implicit none

        integer(kind=4), intent(in) :: orb

#ifdef _QMCKL_GPU

        integer(kind=8) :: ao_num, mo_num

        real(kind=8), allocatable, dimension(:,:) :: electron
        real(kind=8), allocatable, dimension(:) :: ao_values, mo_values

        allocate(electron(3,1))
        ! Randomize the electron position
        call random_number(electron)

        !#######################################################################
        !#                                                                     #
        !#  Test QMCkl GPU                                                     #
        !#                                                                     #
        !#  Part regarding loading the trexio file is missing here. It should  #
        !#  be part if the init/finish subroutines.                            #
        !#                                                                     #
        !#######################################################################

        write (6,*) "Testing QMCkl GPU"

        !#######################################################################
        !#                                                                     #
        !# First get ao_num and mo_num                                         #
        !# these can be obtained from qmckl or qmckl_gpu                       #
        !#                                                                     #
        !#######################################################################
        
        rc = qmckl_gpu_get_ao_basis_ao_num(qmckl_ctx, ao_num)
        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting ao_num', rc, qmckl_ctx, ao_num
            call abort()
        end if

        write (6,*) "ao_num = ", ao_num

        rc = qmckl_gpu_get_mo_basis_mo_num(qmckl_ctx, mo_num)
        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting mo_num', rc, qmckl_ctx, mo_num
            call abort()
        end if

        write (6,*) "mo_num = ", mo_num

        !#######################################################################
        !#                                                                     #
        !# Is it possible to get them from qmckl_gpu?                          #
        !#                                                                     #
        !# rc = qmckl_gpu_get_ao_basis_ao_num(qmckl_ctx, ao_num)               #
        !# if (rc /= QMCKL_SUCCESS) then                                       #
        !#     print *, 'Error getting ao_num', rc, qmckl_ctx, ao_num          #
        !#     call abort()                                                    #
        !# end if                                                              #
        !#                                                                     #
        !# rc = qmckl_gpu_get_mo_basis_mo_num(qmckl_ctx, mo_num)               #
        !# if (rc /= QMCKL_SUCCESS) then                                       #
        !#     print *, 'Error getting mo_num', rc, qmckl_ctx, mo_num          #
        !#     call abort()                                                    #
        !# end if                                                              #
        !#                                                                     #
        !# Or is it even necessary? Does it make sense to have two contexts?   #
        !#                                                                     #
        !#######################################################################

        !#######################################################################
        !#                                                                     #
        !# Now allocate the arrays                                             #
        !#                                                                     #
        !#######################################################################

        allocate(ao_values(ao_num))
        allocate(mo_values(mo_num))

        !#######################################################################
        !#                                                                     #
        !# Now offload the data to the GPU                                     #
        !#                                                                     #
        !#######################################################################

!$omp target enter data map(to:electron) map(from:ao_values,mo_values)

!$omp target data use_device_ptr(electron)
        rc = qmckl_gpu_set_point(qmckl_gpu_ctx, 'N', 1_8, electron, 3_8)

        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error setting electron coords'
            call abort()
        end if

!$omp target data use_device_ptr(ao_values)
        rc = qmckl_gpu_get_ao_basis_ao_value_inplace(                    &
               &qmckl_gpu_ctx,                                           &
               &ao_values,                                               &
               &ao_num)

        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting AOs from QMCkl'
            call abort()
        end if

!$omp target data use_device_ptr(mo_values)
        rc = qmckl_gpu_get_mo_basis_mo_value_inplace(                    &
               &qmckl_gpu_ctx,                                           &
               &mo_values,                                               &
               &mo_num)

        if (rc /= QMCKL_SUCCESS) then
            print *, 'Error getting MOs from QMCkl'
            call abort()
        end if
!$omp end data

        print *, "ao_values = ", ao_values
        print *, "mo_values = ", mo_values

        if (allocated(electron)) deallocate(electron) 
        if (allocated(ao_values)) deallocate(ao_values) 
        if (allocated(mo_values)) deallocate(mo_values)
        
        !#######################################################################
        !#                                                                     #
        !#  By looking at the documentation I understand the arrays are        #
        !#  assumedto be in the device memory. While numbers like ao_num and   #
        !#  mo_num are scalars and are passed by value on CPU.                 #
        !#                                                                     #
        !#  End of test QMCkl GPU                                              #
        !#                                                                     #
        !#######################################################################

#endif

    end subroutine test_qmckl_gpu

    subroutine run_openmp_test()

        implicit none

        integer(kind=4) :: array(1)

        array(1) = 0
!$omp target data map(tofrom:array)
        array(1) = 1
!$omp end target data
        if (array(1).ne.0) then
            write (0,*) "OpenMP test failed"
            stop 1
        else
            write (6,*) "OpenMP test passed"
        end if

    end subroutine run_openmp_test
    
#endif

end program test_makefun_qmckl
