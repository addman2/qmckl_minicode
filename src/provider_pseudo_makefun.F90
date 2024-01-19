subroutine provider_pseudo_makefun(indt&
                                &, indtm&
                                &, points&
                                &, ang_mom&
                                &, par&
                                &, gl&
                                &, multiplicity&
                                &, values)

    implicit none

    integer*4, intent(in) :: indt, indtm, ang_mom, multiplicity, gl
    real*8, intent(in) :: par
    real*8, intent(in) :: points(3,indtm)

    real*8, intent(inout) :: values((indt + 4 * gl) * multiplicity)
    
    real*8 :: zeta(1), dd(1)&
           &, r(0:indtm), rmu(3, 0:indtm), distp(0:indtm, 20)

    real*8, dimension(:, :), allocatable :: z

    integer*4 :: i0, iopt, indtmin, ii, jj, kk
    integer*4 :: indpar, indorb, indshell
    integer*4 :: iflagnorm_unused, cr, nelskip

    print *, "indt", indt
    print *, "indtm", indtm
    print *, "ang_mom", ang_mom
    print *, "par", par
    print *, "gl", gl
    print *, "multiplicity", multiplicity

    iopt = ang_mom + 90
    i0 = 0_4
    indtmin = 0_4
    cr = 0_4
    iflagnorm_unused = 0_4

    indpar = 0_4
    indorb = 0_4
    indshell = 0_4
    nelskip = multiplicity
    allocate(z(nelskip, 0:indt+4))
    z = 0

    dd(1) = par

    rmu = points
    r = sqrt(sum(rmu**2, dim=1))
    
    call makefun(iopt&
              &, indt&
              &, i0&
              &, indtmin&
              &, indtm&
              &, 0_4&
              &, indpar&
              &, indorb&
              &, indshell&
              &, nelskip&
              &, z&
              &, dd&
              &, zeta&
              &, r&
              &, rmu&
              &, distp&
              &, iflagnorm_unused&
              &, cr)

    do ii = 1, indtm
        do jj = 1, multiplicity
            values(jj + multiplicity*(ii-1)) = z(jj, ii-1)
        end do
    end do

    do jj = 1, multiplicity
        values(jj + multiplicity*(indt + 0)) = z(jj, indt + 1)
        values(jj + multiplicity*(indt + 1)) = z(jj, indt + 2)
        values(jj + multiplicity*(indt + 2)) = z(jj, indt + 3)
        values(jj + multiplicity*(indt + 3)) = z(jj, indt + 4)
    end do

    deallocate(z)

end subroutine provider_pseudo_makefun
