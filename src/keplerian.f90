!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Module containing various public constants      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module keplerian
use geometry
use constants
  contains
subroutine solve_kep_eq(e,amean,fan)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Solves the Kepler's equation and provides the    !
  !  true anomaly                                    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(kind=16), intent(in) :: e, amean
  real(kind=16), intent(out) :: fan
  real(kind=16) :: eec
  integer :: iii
  eec=0.d0
  do iii=1,100
  eec=amean+e*sin(eec)
  end do
  fan=2.d0*atan(sqrt((1.d0+e)/(1.d0-e))*tan(eec/2.d0))
end subroutine solve_kep_eq

subroutine indeb(x, a, e, bom, om, inc, fan, mu)
  real(kind=16), intent(in) :: a, e, bom, om, inc, fan, mu
  real(kind=16), intent(out), dimension(6) :: x
  real(kind=16) :: r, rdot, rfdot, yy, zz, yydot, zzdot
  r=a*(1.d0-e*e)/(1.d0+e*cos(fan))
  rdot=(sqrt(mu)/sqrt(a*(1.d0-e*e)))*e*sin(fan)
  rfdot=(sqrt(mu)/sqrt(a*(1.d0-e*e)))*(1.d0+e*cos(fan))
  x(1)=r*(cos(bom)*cos(om+fan)-sin(bom)*sin(om+fan)*cos(inc))
  yy=r*(sin(bom)*cos(om+fan)+cos(bom)*sin(om+fan)*cos(inc))
  zz=r*sin(om+fan)*sin(inc)
  x(2)=yy
  x(3)=zz
  x(4)=rdot*(cos(bom)*cos(om+fan)-sin(bom)*sin(om+fan)*cos(inc)) &
  &-rfdot*(cos(bom)*sin(om+fan)+sin(bom)*cos(om+fan)*cos(inc))
  yydot=rdot*(sin(bom)*cos(om+fan)+cos(bom)*sin(om+fan)*cos(inc)) &
  &-rfdot*(sin(bom)*sin(om+fan)-cos(bom)*cos(om+fan)*cos(inc))
  zzdot=rdot*sin(om+fan)*sin(inc)+rfdot*cos(om+fan)*sin(inc)
  x(5)=yydot
  x(6)=zzdot
  return
end subroutine indeb

subroutine orbital(x, orb,amuE)
  real(kind=16), intent(in) :: amuE
  real(kind=16), intent(in), dimension(6) :: x
  real(kind=16), intent(out), dimension(7) :: orb
  real(kind=16), dimension(3) :: Vr, Vv, Vh, Ve, Vn, Vi, Vj, Vk
  real(kind=16), dimension(3) :: VhovSh, VnovSn, auxvec
  real(kind=16) :: Sr, Sv, Sh, Se, Sa, Si, Sn, SbOm, aux, Som, Sf
  real(kind=16) :: SEE, SM, SMM, St, tiny
  tiny=1.d-8
  Vi(1)=1.d0
  Vi(2)=0.d0
  Vi(3)=0.d0
  Vj(1)=0.d0
  Vj(2)=1.d0
  Vj(3)=0.d0
  Vk(1)=0.d0
  Vk(2)=0.d0
  Vk(3)=1.d0
  Vr(1)=x(1)
  Vr(2)=x(2)
  Vr(3)=x(3)
  Sr=sqrt(x(1)**2+x(2)**2+x(3)**2)
  Vv(1)=x(4)
  Vv(2)=x(5)
  Vv(3)=x(6)
  Sv=sqrt(x(4)**2+x(5)**2+x(6)**2)
  ! Vector aligned along the angular momentum
  call vector_product(Vr, Vv, Vh)
  Sh=sqrt(Vh(1)**2+Vh(2)**2+Vh(3)**2)
  call vector_product(Vv, Vh, Ve)
  Ve(1)=Ve(1)/amuE-Vr(1)/Sr
  Ve(2)=Ve(2)/amuE-Vr(2)/Sr
  Ve(3)=Ve(3)/amuE-Vr(3)/Sr
  Se=sqrt(Ve(1)**2+Ve(2)**2+Ve(3)**2)
  if (Se .eq. 0.d0) then
    Se=Se+tiny
  end if
  Sa=dot_product(Vh,Vh)/(amuE*(1-Se**2))
  VhovSh(1)=Vh(1)/Sh
  VhovSh(2)=Vh(2)/Sh
  VhovSh(3)=Vh(3)/Sh
  Si=acos(dot_product(Vk,VhovSh))
  if (Si .eq. 0.d0) then
    Si=Si+tiny
  end if
  call vector_product(Vk, Vh, Vn)
  Sn=sqrt(dot_product(Vn, Vn))
  if (Sn .eq. 0.d0) then
    Sn=Sn+tiny
  end if


  VnovSn(1)=Vn(1)/Sn
  VnovSn(2)=Vn(2)/Sn
  VnovSn(3)=Vn(3)/Sn

  SbOm=acos(dot_product(Vi,VnovSn))
  SbOm=mod(SbOm, pi2)

  if (dot_product(Vn,Vj) .lt. 0.d0) then
    SbOm=pi2 -SbOm
    SbOm=mod(SbOm, pi2)
  end if
  ! SbOm=atan2(dot_product(Vj,VnovSn),dot_product(Vi,VnovSn))
  ! TBC
  aux=dot_product(Vn,Ve)/(Sn*Se)
  call vector_product(VhovSh,VnovSn,auxvec)
  Som=atan2(dot_product(Ve,auxvec/sqrt(dot_product(auxvec,auxvec))),aux)

  if (aux .gt.1.d0) then
    aux=1.d0
  end if
  if (aux .lt. -1.d0) then
    aux=-1.d0
  end if
  Som=acos(aux)
  Som=mod(Som, pi2)
  if (dot_product(Ve,Vk) .lt. 0.d0) then
    Som=pi2 -Som
    Som=mod(Som, pi2)
  end if
  aux=dot_product(Ve,Vr)/(Se*Sr)
  if (aux .gt.1.d0) then
     aux=1.d0
  end if
  if (aux .lt. -1.d0) then
    aux=-1.d0
  end if
  Sf=acos(aux)
  if (dot_product(Vr,Vv) .lt. 0.d0) then
    Sf=pi2 -Sf
  end if

  SEE=acos((Se+cos(Sf))/(1.d0+Se*cos(Sf)))
  if ((pi2/2.d0 .lt. Sf) .and. (Sf .lt. pi2)) then
    SEE=pi2 -SEE
  end if

  SM=SEE-Se*sin(SEE)
  SMM=sqrt(amuE/(Sa**3))
  St=SM/SMM

  orb(1)=Sa
  orb(2)=Se
  orb(3)=Si/degree
  orb(4)=Som/degree
  orb(5)=SbOm/degree
  orb(6)=SM/degree
  orb(7)=St


  return
end subroutine orbital

end module keplerian
