
    indshellp=indshell+1
    indorbp=indorb+1
    dd1=dd(indpar+1)

    multiplicity = (iopt - 90 + 2) * (iopt - 90 + 1) / 2

    powers(:,0) = 1.0d0

    do i = 1, max_power
        do k = i0, indtm
            powers(1, i, k) = powers(1, i-1, k) * rmu(1, k)
            powers(2, i, k) = powers(2, i-1, k) * rmu(2, k)
            powers(3, i, k) = powers(3, i-1, k) * rmu(3, k)
        end do
    end do

    do k = indtmin, indtm
        distp(k,1) = dexp(-1.0_8 * dd1 * r(k) * r(k))
    end do

    c = 0.712705470354990_8 * dd1 ** 0.75_8
    do k = i0, indtm
        do i = i, multiplicity
            z(indorbp + i - 1, k) = c
        end do
    end do

    ! Initialize gradients and laplacians (radial part)

    if (typec .ne. 0) then
        distp(indtm + 1, 1) = -2.0d0 * dd1 * rmu(1, k) * distp(1, 1)
        distp(indtm + 2, 1) = -2.0d0 * dd1 * rmu(2, k) * distp(1, 1)
        distp(indtm + 3, 1) = -2.0d0 * dd1 * rmu(3, k) * distp(1, 1)
        distp(indtm + 4, 1) = dd1 * (4.0d0 * dd1 * (r(1) * r(1)) - 6.0d0) * distp(indtm + 1, 1)
    end if

    count = 0
    do i = (iopt - 90), 0, -1
        do j = (iopt - 90) - i, 0, -1
            k = (iopt - 90) - i - j
            z(indorbp + count, k) = z(indorbp + count, k) * powers(1, i, k)
            z(indorbp + count, k) = z(indorbp + count, k) * powers(2, j, k)
            z(indorbp + count, k) = z(indorbp + count, k) * powers(3, k, k)
            count = count + 1
        end do
    end do

    if (typec .ne. 0) then
        ! Solve ang_mom = 0, 1 separately
        if (iopt - 90 .eq. 0) then
            z(indorbp, indtm + 1) = distp(indtm + 1, 1)
            z(indorbp, indtm + 2) = distp(indtm + 2, 1)
            z(indorbp, indtm + 3) = distp(indtm + 3, 1)
            z(indorbp, indtm + 4) = distp(indtm + 4, 1)
        else if (iopt - 90 .eq. 1) then
            z(indorbp    , indtm + 1) = distp(indtm + 1, 1) * rmu(1, i0) + distp(1, 1)
            z(indorbp    , indtm + 2) = distp(indtm + 2, 1) * rmu(1, i0)
            z(indorbp    , indtm + 3) = distp(indtm + 3, 1) * rmu(1, i0)

            z(indorbp + 1, indtm + 1) = distp(indtm + 1, 1) * rmu(2, i0)
            z(indorbp + 1, indtm + 2) = distp(indtm + 2, 1) * rmu(2, i0) + distp(1, 1)
            z(indorbp + 1, indtm + 3) = distp(indtm + 3, 1) * rmu(2, i0)

            z(indorbp + 2, indtm + 1) = distp(indtm + 1, 1) * rmu(3, i0)
            z(indorbp + 2, indtm + 2) = distp(indtm + 2, 1) * rmu(3, i0)
            z(indorbp + 2, indtm + 3) = distp(indtm + 3, 1) * rmu(3, i0) + distp(1, 1)

            z(indorbp    , indtm + 4) = distp(indtm + 4, 1) * rmu(1, i0) + 2.0d0 * distp(indtm + 1, 1)
            z(indorbp + 1, indtm + 4) = distp(indtm + 4, 1) * rmu(2, i0) + 2.0d0 * distp(indtm + 2, 1)
            z(indorbp + 2, indtm + 4) = distp(indtm + 4, 1) * rmu(3, i0) + 2.0d0 * distp(indtm + 3, 1)
        else
            count = 1
            do i = (iopt - 90), 0, -1
                do j = (iopt - 90) - i, 0, -1
                    k = (iopt - 90) - i - j

                    ! First store polynomial part into respective places
                    ! Then solve full laplacian using using lower derivatives
                    ! Then do the same thing for gradients
                    ! Then finally the values

                    z(indorbp + count - 1, indtm + 1) = z(indorbp + count - 1, indtm + 1) * powers(1, i-1)
                    z(indorbp + count - 1, indtm + 1) = z(indorbp + count - 1, indtm + 1) * powers(2, j)
                    z(indorbp + count - 1, indtm + 1) = z(indorbp + count - 1, indtm + 1) * powers(3, k)
                    z(indorbp + count - 1, indtm + 1) = z(indorbp + count - 1, indtm + 1) * i

                    z(indorbp + count - 1, indtm + 2) = z(indorbp + count - 1, indtm + 2) * powers(1, i)
                    z(indorbp + count - 1, indtm + 2) = z(indorbp + count - 1, indtm + 2) * powers(2, j-1)
                    z(indorbp + count - 1, indtm + 2) = z(indorbp + count - 1, indtm + 2) * powers(3, k)
                    z(indorbp + count - 1, indtm + 2) = z(indorbp + count - 1, indtm + 2) * j

                    z(indorbp + count - 1, indtm + 3) = z(indorbp + count - 1, indtm + 3) * powers(1, i)
                    z(indorbp + count - 1, indtm + 3) = z(indorbp + count - 1, indtm + 3) * powers(2, j)
                    z(indorbp + count - 1, indtm + 3) = z(indorbp + count - 1, indtm + 3) * powers(3, k-1)
                    z(indorbp + count - 1, indtm + 3) = z(indorbp + count - 1, indtm + 3) * k

                    z(indorbp + count - 1, indtm + 4) = powers(1, i-2) * powers(2, j) * powers(3, k) * i * (i-1)&
                                                    & + powers(1, i) * powers(2, j-2) * powers(3, k) * j * (j-1)&
                                                    & + powers(1, i) * powers(2, j) * powers(3, k-2) * k * (k-1)

                    
                    ! All polynomial parts are now stored
                    ! Now solve laplacian
                    z(indorbp + count - 1, indtm + 4) = z(indorbp + count - 1, indtm + 4) * distp(1, 0) &
                                                    & + 2.0_8 * z(indorbp + count - 1, indtm + 1) * distp(1, indtm + 1) &
                                                    & + 2.0_8 * z(indorbp + count - 1, indtm + 2) * distp(1, indtm + 2) &
                                                    & + 2.0_8 * z(indorbp + count - 1, indtm + 3) * distp(1, indtm + 3) &
                                                    & + z(indorbp + count - 1, indtmin) * distp(1, indtm + 4)

                    ! Now solve gradients
                    z(indorbp + count - 1, indtm + 1) = z(indorbp + count - 1, indtm + 1) * distp(count, 0) &
                                                    & + z(indorbp, i0) * distp(1, indtm + 1)
                    z(indorbp + count - 1, indtm + 2) = z(indorbp + count - 1, indtm + 2) * distp(count, 0) &
                                                    & + z(indorbp, i0) * distp(1, indtm + 2)
                    z(indorbp + count - 1, indtm + 3) = z(indorbp + count - 1, indtm + 3) * distp(count, 0) &
                                                    & + z(indorbp, i0) * distp(1, indtm + 3)
                    count = count + 1
                end do
            end do
        end if

        ! Multiply by radial part for values
        do i = 1, multiplicity
            do k = i0, indtm
                z(indorbp + i - 1, k) = z(indorbp + i - 1, k) * distp(k, i)
            end do
        end do

    end if

    indpar=indpar + 1
    indshell=indshell + multiplicity
    indorb=indorbp
