
    indshellp=indshell+1
    indorbp=indorb+1
    dd1=dd(indpar+1)

    multiplicity = (iopt - 90 + 2) * (iopt - 90 + 1) / 2

    powers(:,0,:) = 1.0d0

    do i = 1, max_power
        do k = indtmin, indtm
            powers(1, i, k) = powers(1, i-1, k) * rmu(1, k)
            powers(2, i, k) = powers(2, i-1, k) * rmu(2, k)
            powers(3, i, k) = powers(3, i-1, k) * rmu(3, k)
        end do
    end do

    do k = i0, indtm
        distp(k,1) = dexp(-1.0_8 * dd1 * r(k) * r(k))
    end do

    c = 0.712705470354990_8 * dd1 ** 0.75_8! * 2.829
    if (iopt - 90 .ne. 0) then
        c = c * (8_4 * dd1) ** ((iopt - 90)/2.0_8)
    end if
    do k = i0, indtm
        count = 0
        do ii = (iopt - 90), 0, -1
            do jj = (iopt - 90) - ii, 0, -1
                kk = (iopt - 90) - ii - jj
                z(indorbp + count, k) = c
                rp1 = 1.0_8
                do i = ii + 1, 2 * ii
                    rp1 = rp1 * i
                end do
                z(indorbp + count, k) = z(indorbp + count, k) / dsqrt(rp1)
                rp1 = 1.0_8
                do i = jj + 1, 2 * jj
                    rp1 = rp1 * i
                end do
                z(indorbp + count, k) = z(indorbp + count, k) / dsqrt(rp1)
                rp1 = 1.0_8
                do i = kk + 1, 2 * kk
                    rp1 = rp1 * i
                end do
                z(indorbp + count, k) = z(indorbp + count, k) / dsqrt(rp1)
                count = count + 1
            end do
        end do
    end do

    ! We need to calculate it again for derivatives, it could not be done in previous loop because of case if i0 /= indtmin
    if (typec .ne. 1) then
        count = 0
        do ii = (iopt - 90), 0, -1
            do jj = (iopt - 90) - ii, 0, -1
                kk = (iopt - 90) - ii - jj
                z(indorbp + count, indt + 1) = c
                z(indorbp + count, indt + 2) = c
                z(indorbp + count, indt + 3) = c
                z(indorbp + count, indt + 4) = c
                rp1 = 1.0_8
                do i = ii + 1, 2 * ii
                    rp1 = rp1 * i
                end do
                z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) / dsqrt(rp1)
                z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) / dsqrt(rp1)
                z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) / dsqrt(rp1)
                z(indorbp + count, indt + 4) = z(indorbp + count, indt + 4) / dsqrt(rp1)
                rp1 = 1.0_8
                do i = jj + 1, 2 * jj
                    rp1 = rp1 * i
                end do
                z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) / dsqrt(rp1)
                z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) / dsqrt(rp1)
                z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) / dsqrt(rp1)
                z(indorbp + count, indt + 4) = z(indorbp + count, indt + 4) / dsqrt(rp1)
                rp1 = 1.0_8
                do i = kk + 1, 2 * kk
                    rp1 = rp1 * i
                end do
                z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) / dsqrt(rp1)
                z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) / dsqrt(rp1)
                z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) / dsqrt(rp1)
                z(indorbp + count, indt + 4) = z(indorbp + count, indt + 4) / dsqrt(rp1)
                count = count + 1
            end do
        end do
    end if

    ! Initialize gradients and laplacians (radial part)

    if (typec .ne. 1) then
        distp(indt + 1, 1) = -2.0d0 * dd1 * rmu(1, 0) * distp(0, 1) * c
        distp(indt + 2, 1) = -2.0d0 * dd1 * rmu(2, 0) * distp(0, 1) * c
        distp(indt + 3, 1) = -2.0d0 * dd1 * rmu(3, 0) * distp(0, 1) * c
        distp(indt + 4, 1) = dd1 * (4.0d0 * dd1 * (r(0) * r(0)) - 6.0d0) * distp(0, 1) * c
    end if

    do k = i0, indtm
        count = 0
        do ii = (iopt - 90), 0, -1
            do jj = (iopt - 90) - ii, 0, -1
                kk = (iopt - 90) - ii - jj
                z(indorbp + count, k) = z(indorbp + count, k) * powers(1, ii, k)
                z(indorbp + count, k) = z(indorbp + count, k) * powers(2, jj, k)
                z(indorbp + count, k) = z(indorbp + count, k) * powers(3, kk, k)
                count = count + 1
            end do
        end do
    end do

    if (typec .ne. 1) then
        ! Solve ang_mom = 0, 1 separately
        if (iopt - 90 .eq. 0) then
            z(indorbp, indt + 1) = distp(indt + 1, 1) 
            z(indorbp, indt + 2) = distp(indt + 2, 1)
            z(indorbp, indt + 3) = distp(indt + 3, 1)
            z(indorbp, indt + 4) = distp(indt + 4, 1)
        else if (iopt - 90 .eq. 1) then
            rp1 = dsqrt(2.0_8)
            z(indorbp    , indt + 1) = (distp(indt + 1, 1) * rmu(1, indtmin) + c * distp(0, 1)) / rp1
            z(indorbp    , indt + 2) = (distp(indt + 2, 1) * rmu(1, indtmin)) / rp1
            z(indorbp    , indt + 3) = (distp(indt + 3, 1) * rmu(1, indtmin)) / rp1

            z(indorbp + 1, indt + 1) = (distp(indt + 1, 1) * rmu(2, indtmin)) / rp1
            z(indorbp + 1, indt + 2) = (distp(indt + 2, 1) * rmu(2, indtmin) + c * distp(0, 1)) / rp1
            z(indorbp + 1, indt + 3) = (distp(indt + 3, 1) * rmu(2, indtmin)) / rp1

            z(indorbp + 2, indt + 1) = (distp(indt + 1, 1) * rmu(3, indtmin)) / rp1
            z(indorbp + 2, indt + 2) = (distp(indt + 2, 1) * rmu(3, indtmin)) / rp1
            z(indorbp + 2, indt + 3) = (distp(indt + 3, 1) * rmu(3, indtmin) + c * distp(0, 1)) / rp1

            z(indorbp    , indt + 4) = (distp(indt + 4, 1) * rmu(1, indtmin) + 2.0d0 * distp(indt + 1, 1)) / rp1
            z(indorbp + 1, indt + 4) = (distp(indt + 4, 1) * rmu(2, indtmin) + 2.0d0 * distp(indt + 2, 1)) / rp1
            z(indorbp + 2, indt + 4) = (distp(indt + 4, 1) * rmu(3, indtmin) + 2.0d0 * distp(indt + 3, 1)) / rp1
        else
            count = 0
            do ii = (iopt - 90), 0, -1
                do jj = (iopt - 90) - ii, 0, -1
                    kk = (iopt - 90) - ii - jj

                    ! First store polynomial part into respective places
                    ! Then solve full laplacian using using lower derivatives
                    ! Then do the same thing for gradients
                    ! Then finally the values

                   z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) * powers(1, ii-1, 0)
                   z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) * powers(2, jj, 0)
                   z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) * powers(3, kk, 0)
                   z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) * ii

                   z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) * powers(1, ii, 0)
                   z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) * powers(2, jj-1, 0)
                   z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) * powers(3, kk, 0)
                   z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) * jj

                   z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) * powers(1, ii, 0)
                   z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) * powers(2, jj, 0)
                   z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) * powers(3, kk-1, 0)
                   z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) * kk

                   z(indorbp + count, indt + 4) = powers(1, ii-2, 0) * powers(2, jj, 0) * powers(3, kk, 0) * ii * (ii-1)&
                                              & + powers(1, ii, 0) * powers(2, jj-2, 0) * powers(3, kk, 0) * jj * (jj-1)&
                                              & + powers(1, ii, 0) * powers(2, jj, 0) * powers(3, kk-2, 0) * kk * (kk-1)

                    
                    ! All polynomial parts are now stored
                    ! Now solve laplacian
                    z(indorbp + count, indt + 4) = z(indorbp + count, indtm + 4) * distp(1, 0) &
                                               & + 2.0_8 * z(indorbp + count, indt + 1) * distp(1, indt + 1) &
                                               & + 2.0_8 * z(indorbp + count, indt + 2) * distp(1, indt + 2) &
                                               & + 2.0_8 * z(indorbp + count, indt + 3) * distp(1, indt + 3) &
                                               & + z(indorbp + count, indtmin) * distp(1, indt + 4)

                    ! Now solve gradients
                    z(indorbp + count, indt + 1) = z(indorbp + count, indt + 1) * distp(count + 1, 0) &
                                               & + z(indorbp, i0) * distp(1, indt + 1)
                    z(indorbp + count, indt + 2) = z(indorbp + count, indt + 2) * distp(count + 1, 0) &
                                               & + z(indorbp, i0) * distp(1, indt + 2)
                    z(indorbp + count, indt + 3) = z(indorbp + count, indt + 3) * distp(count + 1, 0) &
                                               & + z(indorbp, i0) * distp(1, indt + 3)
                    count = count + 1
                end do
            end do
        end if

    end if

    ! Multiply by radial part for values
    do ii = 1, multiplicity
        do kk = i0, indtm
            z(indorbp + ii - 1, kk) = z(indorbp + ii - 1, kk) * distp(kk, 1)
        end do
    end do

    indpar=indpar + 1
    indshell=indshell + multiplicity
    indorb=indorbp
