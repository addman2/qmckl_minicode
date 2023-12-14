program comparator

    use qmckl_gpu_helpers
    use qmckl_helpers

    implicit none

    real(kind=8), dimension(3) :: point
    real(kind=8) :: par
    real(kind=8), dimension(:), allocatable :: values_1, values_2, values_3

    integer(kind=4) :: ang_mom, multiplicity

    ang_mom = 3
    multiplicity = (ang_mom + 1) * (ang_mom + 2) / 2

    allocate(values_1(multiplicity))
    allocate(values_2(multiplicity))
    allocate(values_3(multiplicity))

    ! Setup point
    point = (/ 0.1d0, 0.2d0, 0.3d0 /)

    ! Setup parameters
    par = 1.0d0

    call provider_qmckl(point&
                     &, ang_mom&
                     &, par&
                     &, values_1&
                     &, multiplicity)

    call provider_qmckl_gpu(point&
                     &, ang_mom&
                     &, par&
                     &, values_2&
                     &, multiplicity)

    call provider_makefun(point&
                     &, ang_mom&
                     &, par&
                     &, values_3&
                     &, multiplicity)

    ! Print point with 3 decimal spaces
    write(*, *)  
    write(*, *) "#############################################"
    write(*, *)  
    write(*, '("Point = ", 3F7.3)') point
    write(*, '("Ang mom = ", I3)') ang_mom
    write(*, '("Par = ", F7.3)') par
    write(*, '("Values qmckl = ", 16F7.3)') values_1
    write(*, '("Values qmckl gpu = ", 16F7.3)') values_2
    write(*, '("Values makefun = ", 16F7.3)') values_3
    write(*, *)  
    write(*, *) "#############################################"
    write(*, *)  

    if (allocated(values_1)) deallocate(values_1)
    if (allocated(values_2)) deallocate(values_2)
    if (allocated(values_3)) deallocate(values_3)

end program comparator
