program comparator

    use qmckl_gpu_helpers
    use qmckl_helpers

    implicit none

    real(kind=8), dimension(3) :: point
    real(kind=8) :: par
    real(kind=8), dimension(:), allocatable :: values

    integer(kind=4) :: ang_mom, multiplicity

    ang_mom = 1
    multiplicity = (ang_mom + 1) * (ang_mom + 2) / 2

    allocate(values(multiplicity))

    ! Setup point
    point = (/ 0.1d0, 0.2d0, 0.3d0 /)

    ! Setup parameters
    par = 1.0d0

    call provider_qmckl(point&
                     &, ang_mom&
                     &, par&
                     &, values&
                     &, multiplicity)

    ! Print point with 3 decimal spaces
    write(*, '("Point = ", 3F7.3)') point
    write(*, '("Ang mom = ", I7)') ang_mom
    write(*, '("Par = ", F7.3)') par
    write(*, '("Values = ", 16F7.3)') values

    if (allocated(values)) then
        deallocate(values)
    end if

end program comparator
