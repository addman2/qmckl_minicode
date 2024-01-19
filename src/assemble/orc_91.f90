  !
  ! Cartesian p orbital

  dd1=dd(indpar+1)

  c=dd1**1.25d0*1.42541094070998d0

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)*r(k))
  end do

  indorbp=indorb+1
  do i=i0,indtm
    z(indorbp,i)=rmu(3,i)*distp(i,1)
  end do

  indorbp=indorb+2
  do i=i0,indtm
    z(indorbp,i)=rmu(1,i)*distp(i,1)
  end do

  indorbp=indorb+3
  do i=i0,indtm
    z(indorbp,i)=rmu(2,i)*distp(i,1)
  end do

  if(typec.ne.1) then
    fun0=distp(0,1)
    fun=-2.d0*dd1*distp(0,1)
    fun2=fun*(1.d0-2.d0*dd1*r(0)*r(0))

    indorbp=indorb+1
    do i=1,3
      z(indorbp,indt+i)=rmu(3,0)*rmu(i,0)*fun
    end do
    z(indorbp,indt+3)=z(indorbp,indt+3)+fun0
    z(indorbp,indt+4)=rmu(3,0)*(4.d0*fun+fun2)

    indorbp=indorb+2
    do i=1,3
      z(indorbp,indt+i)=rmu(1,0)*rmu(i,0)*fun
    end do
    z(indorbp,indt+1)=z(indorbp,indt+1)+fun0
    z(indorbp,indt+4)=rmu(1,0)*(4.d0*fun+fun2)

    indorbp=indorb+3
    do i=1,3
      z(indorbp,indt+i)=rmu(2,0)*rmu(i,0)*fun
    end do
    z(indorbp,indt+2)=z(indorbp,indt+2)+fun0
    z(indorbp,indt+4)=rmu(2,0)*(4.d0*fun+fun2)
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp

