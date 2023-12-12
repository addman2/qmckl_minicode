


  dd1=dd(indpar+1)

  !\print *, "i0, indtmin, indtm"
  !\print *,i0, indtmin, indtm


  !        if(iflagnorm.gt.2) then
  !        c=dsqrt(2.d0)*pi**(-0.75d0)*(2.d0*dd1)**1.25d0
  c=dd1**1.25d0*1.42541094070998d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,1)
    end do
    !           endif
  end do

  if(typec.ne.1) then
    fun0=distp(0,1)
    fun=-2.d0*dd1*distp(0,1)
    fun2=fun*(1.d0-2.d0*dd1*r(0)**2)

    !              indorbp=indorb

    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun
      end do
      z(indorbp,indt+ic)=z(indorbp,indt+ic)+fun0
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
      !                 endif
    end do
    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp

