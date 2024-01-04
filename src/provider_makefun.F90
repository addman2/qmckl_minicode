subroutine plot(ii, jj, kk)

    ! Small routine to show the powers

    implicit none

    integer(kind=4), intent(in) :: ii, jj, kk
    integer(kind=4) :: ii2, jj2, kk2

    do ii2 = 1, ii
        write(* ,fmt="(a)", advance="no") "X"
    end do
    do jj2 = 1, jj
        write(* ,fmt="(a)", advance="no") "Y"
    end do
    do kk2 = 1, kk
        write(* ,fmt="(a)", advance="no") "Z"
    end do
    write(*,*) " "

end subroutine plot

subroutine provider_makefun(point&
                       &, ang_mom&
                       &, par&
                       &, values&
                       &, multiplicity&
                       &, gl)

    implicit none

    integer(kind=4), intent(in) :: ang_mom
    integer(kind=4), intent(in) :: multiplicity
    integer(kind=4), intent(in) :: gl

    real(kind=8), dimension(3), intent(in) :: point
    real(kind=8), intent(in) :: par
    ! If gl is 0 only values are computed
    ! If gl is 1 values gradients and laplacians are computed
    ! this requires 5 times more space
    real(kind=8), dimension(multiplicity * (1 + gl * 4)), intent(out) :: values

    ! Local variables

    integer(kind=4), parameter :: max_power = 10

    integer(kind=4) :: ii, jj, kk, count

    real(kind=8), dimension(3,0:max_power) :: powers
    real(kind=8) :: radial_v
    real(kind=8), dimension(3) :: radial_g
    real(kind=8) :: radial_l

    ! Initialize helper values

    powers(:,0) = 1.0d0

    do ii = 1, max_power
        powers(1, ii) = powers(1, ii-1) * point(1)
        powers(2, ii) = powers(2, ii-1) * point(2)
        powers(3, ii) = powers(3, ii-1) * point(3)
    end do

    radial_v = exp(-par * (point(1) * point(1) + point(2) * point(2) + point(3) * point(3)))

    do ii = 1, size(values)
        values(ii) = 0.712705470354990_8
        values(ii) = values(ii) * par ** 0.75
        values(ii) = values(ii) * (8.0_8 * par) ** (ang_mom)
    end do

    ! Initialize gradients and laplacians

    if (gl == 1) then
        radial_g(1) = -2.0d0 * par * point(1) * radial_v
        radial_g(2) = -2.0d0 * par * point(2) * radial_v
        radial_g(3) = -2.0d0 * par * point(3) * radial_v
        radial_l = par * (4.0d0 * par * (point(1) * point(1) + point(2) * point(2) + point(3) * point(3)) - 6.0d0) * radial_v
    end if

    count = 1
    do ii = ang_mom, 0, -1
        do jj = ang_mom - ii, 0, -1
            kk = ang_mom - ii - jj
            !call plot(ii, jj, kk)
            values(count) = values(count) * powers(1, ii)
            values(count) = values(count) * powers(2, jj)
            values(count) = values(count) * powers(3, kk)
            count = count + 1
        end do
    end do

    if (gl == 1) then
        ! Solve ang_mom = 0, 1 separately
        if (ang_mom == 0) then
            values(2) = radial_g(1)
            values(3) = radial_g(2)
            values(4) = radial_g(3)
            values(5) = radial_l
        else if (ang_mom == 1) then
            values(4) = radial_g(1) * point(1) + radial_v
            values(7) = radial_g(2) * point(1)
            values(10) = radial_g(3) * point(1)

            values(5) = radial_g(1) * point(2)
            values(8) = radial_g(2) * point(2) + radial_v
            values(11) = radial_g(3) * point(2)

            values(6) = radial_g(1) * point(3)
            values(9) = radial_g(2) * point(3)
            values(12) = radial_g(3) * point(3) + radial_v

            values(13) = radial_l * point(1) + 2.0d0 * radial_g(1)
            values(14) = radial_l * point(2) + 2.0d0 * radial_g(2)
            values(15) = radial_l * point(3) + 2.0d0 * radial_g(3)
        else
            count = 1
            do ii = ang_mom, 0, -1
                do jj = ang_mom - ii, 0, -1
                    kk = ang_mom - ii - jj
                    !call plot(ii, jj, kk)
                    ! First store polynomial part into respective places
                    ! Then solve full laplacian using using lower derivatives
                    ! Then do the same thing for gradients
                    ! Then finally the values
                    values(1 * multiplicity + count) = values(1 * multiplicity + count) * powers(1, ii-1)
                    values(1 * multiplicity + count) = values(1 * multiplicity + count) * powers(2, jj)
                    values(1 * multiplicity + count) = values(1 * multiplicity + count) * powers(3, kk)
                    values(1 * multiplicity + count) = values(1 * multiplicity + count) * ii

                    values(2 * multiplicity + count) = values(2 * multiplicity + count) * powers(1, ii)
                    values(2 * multiplicity + count) = values(2 * multiplicity + count) * powers(2, jj-1)
                    values(2 * multiplicity + count) = values(2 * multiplicity + count) * powers(3, kk)
                    values(2 * multiplicity + count) = values(2 * multiplicity + count) * jj

                    values(3 * multiplicity + count) = values(3 * multiplicity + count) * powers(1, ii)
                    values(3 * multiplicity + count) = values(3 * multiplicity + count) * powers(2, jj)
                    values(3 * multiplicity + count) = values(3 * multiplicity + count) * powers(3, kk-1)
                    values(3 * multiplicity + count) = values(3 * multiplicity + count) * kk

                    values(4 * multiplicity + count) = powers(1, ii-2) * powers(2, jj) * powers(3, kk) * ii * (ii-1)&
                                                    &+ powers(1, ii) * powers(2, jj-2) * powers(3, kk) * jj * (jj-1)&
                                                    &+ powers(1, ii) * powers(2, jj) * powers(3, kk-2) * kk * (kk-1)

                    
                    ! All polynomial parts are now stored
                    ! Now solve laplacian
                    values(4 * multiplicity + count) = values(4 * multiplicity + count) * radial_v &
                                                   & + 2.0_8 * values(1 * multiplicity + count) * radial_g(1) &
                                                   & + 2.0_8 * values(2 * multiplicity + count) * radial_g(2) &
                                                   & + 2.0_8 * values(3 * multiplicity + count) * radial_g(3) &
                                                   & + values(count) * radial_l

                    ! Now solve gradients
                    values(1 * multiplicity + count) = values(1 * multiplicity + count) * radial_v &
                                                   & + values(count) * radial_g(1)
                    values(2 * multiplicity + count) = values(2 * multiplicity + count) * radial_v &
                                                   & + values(count) * radial_g(2)
                    values(3 * multiplicity + count) = values(3 * multiplicity + count) * radial_v &
                                                   & + values(count) * radial_g(3)
                    count = count + 1
                end do
            end do
        end if

        ! Multiply by radial part for values
        do ii = 1, multiplicity
            values(ii) = values(ii) * radial_v
        end do

    end if

end subroutine provider_makefun
