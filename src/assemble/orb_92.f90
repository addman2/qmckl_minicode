  ! d orbitals
  ! R(r)= exp(-alpha r^2)
  ! each gaussian term is normalized

  indparp=indpar+1
  dd1=dd(indparp)

  c=dd1**1.75d0*1.64592278064948967213d0
  
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)*r(k))
  end do

  do i=indtmin,indtm
    distp(i,2)=rmu(1,i)*rmu(1,i)
    distp(i,3)=rmu(1,i)*rmu(2,i)
    distp(i,4)=rmu(1,i)*rmu(3,i)
    distp(i,5)=rmu(2,i)*rmu(2,i)
    distp(i,6)=rmu(2,i)*rmu(3,i)
    distp(i,7)=rmu(3,i)*rmu(3,i)
  end do

  do ic=1,6
    indorbp=indorb+ic
    do k=i0,indtm
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)
    end do
  end do

  if(typec.ne.1) then
  !if(.false.) then

    fun0=distp(0,1)
    fun=-2.0_8*dd1
    fun2=dd1*(4.0_8 * r(0) * r(0) * dd1 - 6.0_8)* distp(0, 1)

    do ic=1,6
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=distp(0,1+ic)*fun
        if (i.eq.1) then
            z(indorbp, indt+1) = 0.0 ! (z(indorbp, indt+1) +&
                                !&(2*rmu(1,indtmin) + rmu(2,indtmin) + rmu(3,indtmin)))* fun0
        end if
        if (i.eq.2) then
            z(indorbp, indt+2) = 0.0 ! (z(indorbp, indt+2) +&
                                !&(rmu(1,indtmin) + 2*rmu(2,indtmin) + rmu(3,indtmin)))* fun0
        end if
        if (i.eq.3) then
            z(indorbp, indt+3) = 0.0 ! (z(indorbp, indt+3) +&
                                !&(rmu(1,indtmin) + rmu(2,indtmin) + 2*rmu(3,indtmin)))* fun0
        end if
      end do
      z(indorbp,indt+4)=0.0 !distp(0,1+ic)*fun2
    end do

  end if

  indpar=indpar+1
  indshell=indshell+6
  indorb=indorbp

