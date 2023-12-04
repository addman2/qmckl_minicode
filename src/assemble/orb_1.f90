


  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  dd1=dd(indpar+1)
  !           if(iflagnorm.gt.2) then
  !           c=dd1*dsqrt(dd1)/dsqrt(pi)
  c=dd1*dsqrt(dd1)*0.56418958354775628695d0
  !           endif

  indorbp=indorb+1
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)
  end do

  if(typec.ne.1) then
    fun=-dd1*distp(0,1)


    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do

    z(indorbp,indt+4)=(-2.d0*dd1/r(0)+dd1**2)                        &
        *distp(0,1)


  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp


  ! 1s double Z with cusp cond
