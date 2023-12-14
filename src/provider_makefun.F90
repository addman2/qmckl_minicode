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
                       &, multiplicity)

    implicit none

    integer(kind=4), intent(in) :: ang_mom
    integer(kind=4), intent(in) :: multiplicity

    real(kind=8), dimension(3), intent(in) :: point
    real(kind=8), intent(in) :: par
    real(kind=8), dimension(multiplicity), intent(out) :: values

    ! Local variables

    integer(kind=4), parameter :: max_power = 10

    integer(kind=4) :: ii, jj, kk, count

    real(kind=8), dimension(3,0:max_power) :: powers

    ! Initialize powers

    powers(:,0) = 1.0d0

    do ii = 1, max_power
        powers(1, ii) = powers(1, ii-1) * point(1)
        powers(2, ii) = powers(2, ii-1) * point(2)
        powers(3, ii) = powers(3, ii-1) * point(3)
    end do

    do ii = 1, multiplicity
        values(ii) = exp(-par * (point(1) * point(1) + point(2) * point(2) + point(3) * point(3)))
    end do

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

end subroutine provider_makefun
