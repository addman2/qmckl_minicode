  indparp=indpar+1

  dd1=dd(indparp)

  c=dd1**2.25d0*1.47215808929909374563d0
  
  !multiplicity=(iopt+1)*(iopt+2)/2

  ! distp(:,1) contains the radial part of the wavefunction
  ! distp(:,2-10) contains the angular part of the wavefunction

  do k=indtmin, indtm
    distp(k,1)=c*dexp(-dd1*r(k)*r(k))
    i = 2
    do ii=1, 3
      do jj=ii, 3
        do kk=jj, 3
          distp(k,i)=rmu(ii,k)*rmu(jj,k)*rmu(kk,k)
          i=i+1
        end do
      end do
    end do  
  end do

  do ic=1, 10
    indorbp=indorb+ic
    do k=i0,indtm
      z(indorbp,k)=distp(k,1)*distp(k,1+ic)
    end do
  end do

  if (typec.ne.1) then
  !if(.false.) then
    
    fun0=distp(0,1)
    fun=-2.0_8*dd1
    fun2=2.0_8*dd1*(2.0_8*dd1*r(0)*r(0)-3.0_8)

    ! Add derivative of the radial part multiplied by the angular part
    indorbp=indorb

    ! xxx
    indorbp=indorbp+1
    z(indorbp,indtm+1)=3*rmu(1,1)*rmu(1,1)
    z(indorbp,indtm+2)=0
    z(indorbp,indtm+3)=0
    ! Calculate the second derivative of the radial part
    z(indorbp,indtm+4)=6*rmu(1,1)
    ! Add mixed term

    ! xxy
    indorbp=indorbp+1
    z(indorbp,indtm+1)=2*rmu(1,1)*rmu(2,1)
    z(indorbp,indtm+2)=rmu(1,1)*rmu(1,1)
    z(indorbp,indtm+3)=0
    z(indorbp,indtm+4)=2*rmu(2,1)

    ! xxz
    indorbp=indorbp+1
    z(indorbp,indtm+1)=2*rmu(1,1)*rmu(2,1)
    z(indorbp,indtm+2)=0
    z(indorbp,indtm+3)=rmu(1,1)*rmu(1,1)
    z(indorbp,indtm+4)=2*rmu(3,1)

    ! xyy
    indorbp=indorbp+1
    z(indorbp,indtm+1)=rmu(2,1)*rmu(2,1)
    z(indorbp,indtm+2)=2*rmu(1,1)*rmu(2,1)
    z(indorbp,indtm+3)=0
    z(indorbp,indtm+4)=2*rmu(1,1)

    ! xyz
    indorbp=indorbp+1
    z(indorbp,indtm+1)=rmu(2,1)*rmu(3,1)
    z(indorbp,indtm+2)=rmu(1,1)*rmu(3,1)
    z(indorbp,indtm+3)=rmu(1,1)*rmu(2,1)
    z(indorbp,indtm+4)=0.0_8

    ! xzz
    indorbp=indorbp+1
    z(indorbp,indtm+1)=rmu(3,1)*rmu(3,1)
    z(indorbp,indtm+2)=0
    z(indorbp,indtm+3)=2*rmu(1,1)*rmu(3,1)
    z(indorbp,indtm+4)=2*rmu(1,1)

    ! yyy
    indorbp=indorbp+1
    z(indorbp,indtm+1)=0
    z(indorbp,indtm+2)=3*rmu(2,1)*rmu(2,1)
    z(indorbp,indtm+3)=0
    z(indorbp,indtm+4)=6*rmu(2,1)

    ! yyz
    indorbp=indorbp+1
    z(indorbp,indtm+1)=0
    z(indorbp,indtm+2)=2*rmu(2,1)*rmu(3,1)
    z(indorbp,indtm+3)=rmu(2,1)*rmu(2,1)
    z(indorbp,indtm+4)=6*rmu(3,1)

    ! yzz
    indorbp=indorbp+1
    z(indorbp,indtm+1)=0
    z(indorbp,indtm+2)=rmu(3,1)*rmu(3,1)
    z(indorbp,indtm+3)=2*rmu(2,1)*rmu(3,1)
    z(indorbp,indtm+4)=2*rmu(2,1)

    ! zzz
    indorbp=indorbp+1
    z(indorbp,indtm+1)=0
    z(indorbp,indtm+2)=0
    z(indorbp,indtm+3)=3*rmu(3,1)*rmu(3,1)
    z(indorbp,indtm+4)=6*rmu(3,1)

    ! P''(r)*R(r) + 2*P'(r)*R'(r) + P(r)*R''(r) =
    ! R(r)*[P''(r) + 2*P'(r)*G(r) + P(r)*L(r)] =
    !
    ! Now in indtm+4 the P''(r) is stored

    do ic=1,10
      indorbp=indorb+ic
      ! Add mixed term 2*R'(r)*P'(r). First derivatives are stored pristine in z(indorbp,indtm+1:indtm+3)
      do ii = 1, 3
        z(indorbp,indtm+4)=z(indorbp,indtm+4) + dd1*rmu(ii,1)*z(indorbp,indtm+ii)
      end do
      ! Add mixed term P(r)*L(r)
      ! P is stored in distp(:,1)
      ! L is stared in fun2
      z(indorbp,indtm+4)=z(indorbp,indtm+4) + distp(ic,i0)*fun2
      ! now multiply by R(r)
      z(indorbp,indtm+4)=z(indorbp,indtm+4)*fun0
      do ii = 1, 3
        z(indorbp,indtm+ii)=(fun*distp(ic, i0) + z(indorbp,indtm+ii))*fun0
      end do
    end do

  end if

  indpar=indpar+1
  indshell=indshell+10
  indorb=indorbp

