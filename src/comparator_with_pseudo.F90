program comparator

    use qmckl_gpu_helpers
    use qmckl_helpers

    implicit none

    real(kind=8), dimension(:,:), allocatable :: point
    real(kind=8) :: par
    real(kind=8), dimension(:), allocatable :: values_1, values_2, values_3

    integer(kind=4) :: ang_mom, multiplicity, gl, indt, indtm

    gl = 1_4
    indt = 8_4
    indtm = 6_4
    ang_mom = 3
    multiplicity = (ang_mom + 1) * (ang_mom + 2) / 2

    allocate(values_1((gl * 4 + indt) * multiplicity))
    allocate(values_2((gl * 4 + indt) * multiplicity))
    allocate(values_3((gl * 4 + indt) * multiplicity))
    allocate(point(3, 0:indtm))

    ! Clear values
    values_1 = 0.0d0
    values_2 = 0.0d0
    values_3 = 0.0d0

    ! Setup seed 1
    !call random_init(.true., .true.)

    ! Setup point
    call random_number(point)
    point = point * 2.0d0 - 1.0d0

    ! Setup parameters
    par = 2.0d0

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
    
    call provider_pseudo_makefun(indt, indtm, point, ang_mom, par, gl, multiplicity, values_3)
    call provider_pseudo_qmckl(indt, indtm, point, ang_mom, par, gl, multiplicity, values_1)

    if (ang_mom.lt.2) then
        call provider_pseudo_makefun_old(indt, indtm, point, ang_mom, par, gl, multiplicity, values_2)
        call show_values(gl, indt, multiplicity, values_2, "makefun old")
    end if
    call show_values(gl, indt, multiplicity, values_3, "makefun")
    call show_values(gl, indt, multiplicity, values_1, "qmckl")

    write(*, *)  
    write(*, *) "#############################################"
    write(*, *)  

    if (allocated(values_1)) deallocate(values_1)
    if (allocated(values_2)) deallocate(values_2)
    if (allocated(values_3)) deallocate(values_3)

end program comparator

subroutine show_values(gl, indt, multiplicity, values, name)

    implicit none

    integer(kind=4), intent(in) :: multiplicity, gl, indt
    character(len=*), intent(in) :: name
    real(kind=8), dimension( ((4 * gl) + indt) * multiplicity ), intent(in) :: values

    integer(kind=4) :: ii, jj
    character(len=100) ::format_string

    format_string = " "
    write(unit=format_string, fmt='(A, I1, A)') '(I3, ": ", ',indt, 'F7.3, " | ", 3F7.3, " | ", F7.3)'

    if (gl.eq.0_4) then
        write(*, '("Values ", A10, " = ", 16F7.3)') adjustl(name), values
    else
        write(*, '("Values ", A10)') adjustl(name)
        do ii = 1, multiplicity
            write(*, format_string) ii&
                                 &, (values(ii + jj*multiplicity), jj = 0, indt - 1)&
                                 &, values((indt + 0) * multiplicity + ii)&
                                 &, values((indt + 1) * multiplicity + ii)&
                                 &, values((indt + 2) * multiplicity + ii)&
                                 &, values((indt + 3) * multiplicity + ii)
        end do
    end if
        
end subroutine show_values
