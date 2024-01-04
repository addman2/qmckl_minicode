program comparator

    use qmckl_gpu_helpers
    use qmckl_helpers

    implicit none

    real(kind=8), dimension(3) :: point
    real(kind=8) :: par
    real(kind=8), dimension(:), allocatable :: values_1, values_2, values_3

    integer(kind=4) :: ang_mom, multiplicity, gl

    gl = 1_4
    ang_mom = 0
    multiplicity = (ang_mom + 1) * (ang_mom + 2) / 2

    allocate(values_1((gl * 4 + 1) * multiplicity))
    allocate(values_2((gl * 4 + 1) * multiplicity))
    allocate(values_3((gl * 4 + 1) * multiplicity))

    ! Clear values
    values_1 = 0.0d0
    values_2 = 0.0d0
    values_3 = 0.0d0

    ! Setup point
    point = (/ 0.1d0, 0.2d0, 0.3d0 /)

    ! Setup parameters
    par = 2.0d0

    call provider_qmckl(point&
                     &, ang_mom&
                     &, par&
                     &, values_1&
                     &, multiplicity&
                     &, gl)

    call provider_qmckl_gpu(point&
                     &, ang_mom&
                     &, par&
                     &, values_2&
                     &, multiplicity&
                     &, gl)

    call provider_makefun(point&
                     &, ang_mom&
                     &, par&
                     &, values_3&
                     &, multiplicity&
                     &, gl)

    write(*, *)  
    write(*, *) "#############################################"
    write(*, *)  
    write(*, '("Point   = ", 3F7.3)') point
    write(*, '("Ang mom = ", I3)') ang_mom
    write(*, '("Par     = ", F7.3)') par
    write(*, *)
    !write(*, '("Values qmckl     = ", 16F7.3)') values_1
    !write(*, '("Values qmckl gpu = ", 16F7.3)') values_2
    !write(*, '("Values makefun   = ", 16F7.3)') values_3
    
    call show_values(gl, multiplicity, values_1, "qmckl")
    call show_values(gl, multiplicity, values_2, "gpu qmckl")
    call show_values(gl, multiplicity, values_3, "makefun")

    write(*, *)  
    write(*, *) "#############################################"
    write(*, *)  

    if (allocated(values_1)) deallocate(values_1)
    if (allocated(values_2)) deallocate(values_2)
    if (allocated(values_3)) deallocate(values_3)

end program comparator

subroutine show_values(gl, multiplicity, values, name)

    implicit none

    integer(kind=4), intent(in) :: multiplicity, gl
    character(len=*), intent(in) :: name
    real(kind=8), dimension( ((4 * gl) + 1) * multiplicity ), intent(in) :: values

    integer(kind=4) :: ii

    if (gl.eq.0_4) then
        write(*, '("Values ", A10, " = ", 16F7.3)') adjustl(name), values
    else
        write(*, '("Values ", A10)') adjustl(name)
        do ii = 1, multiplicity
            write(*, '(I3, ": ", F7.3, " | ", 3F7.3, " | ", F7.3)') ii&
                                                                 &, values(ii)&
                                                                 &, values(1 * multiplicity + ii)&
                                                                 &, values(2 * multiplicity + ii)&
                                                                 &, values(3 * multiplicity + ii)&
                                                                 &, values(4 * multiplicity + ii)
        end do
    end if
        
end subroutine show_values
