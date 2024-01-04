
  indparp=indpar+1
  dd1=dd(indparp)
  c=dd1**2.75d0*1.11284691281640568826d0
  
  do k=indtmin, indtm
    distp(k,1)=c*dexp(-dd1*r(k)*r(k))
    pows(1,0,k) = 1.0_8
    pows(2,0,k) = 1.0_8
    pows(3,0,k) = 1.0_8
    do i=1, 4
      pows(1,i,k) = pows(1,i-1,k) * rmu(1,k)
      pows(2,i,k) = pows(2,i-1,k) * rmu(2,k)
      pows(3,i,k) = pows(3,i-1,k) * rmu(3,k)
    end do
  end do

  indorbp=indorb
  do k=i0,indtm
    do ii=4, 0, -1
      do jj=4-ii, 0, -1
        kk = 4 - (ii + jj)
        z(indorbp,k)=distp(0,1)*pows(1,ii,k)*pows(2,jj,k)*pows(3,kk,k)
        indorbp = indorbp + 1
      end do
    end do
  end do

  if (typec.eq.1) then
  !if(.false.) then
    
    fun0=distp(0,1)
    fun=-2.0_8*dd1
    fun2=2.0_8*dd1*(2.0_8*dd1*r(0)*r(0)-3.0_8)

    ! Add derivative of the radial part multiplied by the angular part
    indorbp=indorb

    do ic=1,15
      indorbp=indorb+ic
      do ii = 1, 4
        z(indorbp,indtm+ii)=0.0_8
      end do
    end do

  end if

  indpar=indpar+1
  indshell=indshell+15
  indorb=indorb+15

