module integration
  use constants

  implicit none
 
    contains
  subroutine adams(fun,xx, ff, tsid, h, iter, eps, conv)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  The Adams method                                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    external fun
    real(kind=16), dimension(:,:), intent(inout) :: xx, ff !(10,dim)
    real(kind=16), intent(in) :: tsid, h, eps
    integer, intent(in) :: iter
    real(kind=16),dimension(1), intent(out) :: conv
    integer :: jj,kk, dim
    integer, dimension(2) :: sha
    
    sha=shape(xx)
    dim=sha(2)

    call adba(fun,xx, ff, tsid, h)
    call admo(fun,xx, ff, tsid, h, iter, eps, conv)
  
    do jj=1,9
    do kk=1,dim
    xx(jj,kk)=xx(jj+1,kk)
    ff(jj,kk)=ff(jj+1,kk)
    end do
    end do
  
    return
  end subroutine adams

  subroutine adba(fun,xx, ff, tsid, h)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Adams-Bashfort                                  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We do not want to solve the variational equation for now, so we
    ! instead focus only on a simple propagation
    external fun
    real(kind=16), intent(in) :: tsid, h
    integer :: ii,jj, dim
    real(kind=16), dimension (:,:), intent(inout) :: xx,ff
    real(kind=16), allocatable, dimension(:) :: ff10, xx10 !dim
    real(kind=16), dimension(9) :: a
    real(kind=16) :: xtemp
    integer, dimension(2) :: sha
  
    a=(/14097247.d0,-43125206.d0, 95476786.d0,-139855262.d0,137968480.d0,  &
    & -91172642.d0,38833486.d0,-9664106.d0,1070017.d0/)
    sha=shape(xx)
    dim=sha(2)
    allocate(ff10(dim),xx10(dim))
    do ii=1,dim
    xtemp=0.0d0
    do jj=1,9
    xtemp=xtemp+(h/3628800.d0)*a(jj)*ff(10-jj,ii)
    end do
    xx(10,ii)=xtemp+xx(9,ii)
    end do
  
    do ii=1,dim
    xx10(ii)=xx(10, ii)
    end do
  
    call fun(tsid, xx10, ff10)
  
    do ii=1,dim
    ff(10, ii)=ff10(ii)
    end do
  
    return
  end subroutine adba
  
  subroutine admo(fun,xx, ff, tsid, h, iter, eps, conv)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Adams-Moulton                                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    external fun
    real(kind=16), dimension(:,:), intent(inout) :: xx, ff !10,dim
    real(kind=16), intent(in) :: tsid, h, eps
    integer, intent(in) :: iter
    real(kind=16), allocatable, dimension(:) :: ff10, xx10
    real(kind=16), allocatable, dimension(:) :: xold, xtemp
    real(kind=16), intent(inout), dimension(1) :: conv
    integer :: jj, ii, it, dim
    real(kind=16), dimension(9) :: b
    integer, dimension(2) :: sha
    b=(/1070017.d0,4467094.d0, -4604594.d0,5595358.d0,-5033120.d0,  &
    & 3146338.d0,-1291214.d0,312874.d0,-33953.d0/)
    sha=shape(xx)
    dim=sha(2)
    allocate(ff10(dim),xx10(dim),xold(dim),xtemp(dim))
    do ii=1,dim  
      xtemp(ii)=0.0d0
      xold(ii)=0.0d0
      xx10(ii)=0.0d0
    end do
    do ii=1,dim    
      xtemp(ii) = xtemp(ii) + xx(9,ii) 
      do jj=1,8
        xtemp(ii)=xtemp(ii)+(h/3628800.d0)*b(jj+1)*ff(10-jj,ii)
      end do
    end do
  
    do it=1,iter
      do ii=1,dim
        xold(ii)=xx(10,ii)
        xx(10,ii)=xtemp(ii)+(h/3628800.d0)*b(1)*ff(10,ii)
      end do
      do ii=1,dim
        xx10(ii)=xx(10, ii)
      end do
      call fun(tsid+h, xx10, ff10)
      do ii=1,dim
        ff(10, ii)=ff10(ii)
      end do
      if (hconv(xx10,xold, eps))  then
        conv(1)=1.0
        !	write(*,*) 'Convergence at iteration ',it,'/', iter
        go to 25
      else
        conv(1)=0.d0
      	write(*,*) ' Adams-Moulton not converge at time = ', tsid
      end if
    end do  ! end iter
    25   continue
    return
  end subroutine admo

  function hconv(xx10,xold, eps)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !            hconv                                 !
    !This function returns 1.d0 (true) if the relative !
    ! difference between the current and previous esti-!
    ! mate of xx10 is less than eps when |xx10| and    !
    ! |xold| are greater than 1 or if the absolute     !
    ! difference of xx10 and xold is less than eps when!
    ! |xx10| and |xold| are less than 1. If neither    !
    ! condition is true, the function returns 0.d0     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=16), dimension(:) :: xx10, xold 
    real(kind=16) :: xxnorm, xonorm, xn
    real(kind=16) :: eps, eps2, maxxn!, hconv
    integer :: i,dim
    logical hconv
    dim=size(xx10)
    eps2=eps
    xxnorm=0.0d0
    xonorm=0.0d0
    xn=0.0d0
    maxxn=0.0d0
    do i=1,dim
      xxnorm=xxnorm+xx10(i)**2
      xonorm=xonorm+xold(i)**2
      xn=xn + (xx10(i)-xold(i))**2
      if (maxxn .lt. abs(xx10(i)-xold(i))) then
        maxxn=abs(xx10(i)-xold(i))
      end if
    end do
    xxnorm=sqrt(xxnorm)
    xonorm=sqrt(xonorm)
    xn=sqrt(xn)
    if((xxnorm .gt. 1.d0) .and. (xonorm .gt. 1.d0)) then
      eps2=eps*xxnorm
    end if
    if(xn .lt. eps2) then
    !	if(maxxn .lt. eps2) then
      ! hconv=1.d0
      hconv=.true.
    else
      ! hconv=0.d0
      hconv=.false.
    end if
    return
  end function hconv

  subroutine ainit(fun,x, xx, ff, h, time0, prel_steps)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Initialize using Runge-Kutta                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    external fun
    real(kind=16), dimension(:), intent(inout) :: x
    integer, intent(in) :: prel_steps
    real(kind=16), intent(in) :: h, time0
    real(kind=16), dimension(:,:), intent(out) :: xx, ff
    real(kind=16), allocatable,dimension(:) :: xtemp, ftemp
    real(kind=16) :: hr, tsid, ait
    integer :: i, ii, jj, it, itMax ,dim

    dim=size(x)
    allocate(xtemp(dim),ftemp(dim))
    ! allocate(xx(10,dim),ff(10,dim))

    do ii=1,9
      do jj=1,dim
        xx(ii,jj)=0.0d0
        ff(ii,jj)=0.0d0
      end do
    end do
    hr=-h/dble(prel_steps) !10000.d0
    tsid=time0
    do ii=1,dim
      xtemp(ii)=x(ii)
    end do
    itMax=9
    do it=1,itMax
      call progress_bar(it,itMax)
      ait=dble(it-1)
      call fun(tsid, xtemp, ftemp)
      do ii=1,dim
        xx(10-it, ii)=xtemp(ii)
        ff(10-it, ii)=ftemp(ii)
      end do! end ii
      do i=1,prel_steps
        call RK6(fun,xtemp, tsid, hr)
        tsid=time0-ait*h+hr*dble(i)
      end do ! end i
    end do ! end it
    return
  end subroutine ainit
  
  ! subroutine ABM(fun,x,tsid,h)
  !   external :: fun
  !   real(kind=16), dimension(:), intent(inout) :: x
  !   real(kind=16), intent(in) :: tsid, h
  !   ! real(kind=16) :: 

  !   call ainit(fun,x,xx,ff,h,tsi)
  ! endsubroutine

  subroutine RK6(fun,x, tsid, h)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Runge-Kutta 6th Order, 7 steps                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    external :: fun
    real(kind=16), dimension(:), intent(inout) :: x
    real(kind=16), intent(in) :: tsid, h
    real(kind=16), allocatable, dimension(:) :: z
    real(kind=16), allocatable, dimension(:) :: bone, btwo, bthree, bfour, bfive
    real(kind=16), allocatable, dimension(:) :: bsix, bseven
    real(kind=16), allocatable, dimension(:) :: b2temp, b3temp, b4temp, b5temp
    real(kind=16), allocatable, dimension(:) :: b6temp, b7temp
    integer :: ii, dim
    
    dim=size(x)
    allocate(z(dim))
    allocate(bone(dim), btwo(dim), bthree(dim), bfour(dim), bfive(dim),bsix(dim),bseven(dim))
    allocate(b2temp(dim), b3temp(dim), b4temp(dim), b5temp(dim),b6temp(dim),b7temp(dim))

    do ii=1,dim
      z(ii)=x(ii)
    end do
    call fun(tsid, z, bone)
    do ii=1,dim
      b2temp(ii)=z(ii)+(h/3.d0)*bone(ii)
    end do
    call fun(tsid+h/3.d0, b2temp, btwo)
    do ii=1,dim
      b3temp(ii)=z(ii)+(2.d0*h/3.d0)*btwo(ii)
    end do
    call fun(tsid+2.d0*h/3.d0, b3temp, bthree)
    do ii=1,dim
      b4temp(ii)=z(ii)+h*(bone(ii)+4.d0*btwo(ii)-bthree(ii))/12.d0
    end do
    call fun(tsid+h/3.d0, b4temp, bfour)
    do ii=1,dim
      b5temp(ii)=z(ii)+h*(-bone(ii)+18.d0*btwo(ii) &
    &     -3.d0*bthree(ii)-6.d0*bfour(ii))/16.d0
    end do
    call fun(tsid+h/2.d0, b5temp, bfive)
    do ii=1,dim
      b6temp(ii)=z(ii)+h*(9.d0*btwo(ii) &
    &     -3.d0*bthree(ii)-6.d0*bfour(ii)+4.d0*bfive(ii))/8.d0
    end do
    call fun(tsid+h/2.d0, b6temp, bsix)
    do ii=1,dim
      b7temp(ii)=z(ii)+h*(9.d0*bone(ii)-36.d0*btwo(ii) &
    &     +63.d0*bthree(ii)+72.d0*bfour(ii)-64.d0*bfive(ii))/44.d0
    end do
    call fun(tsid+h, b7temp, bseven)
    do ii=1,dim
      x(ii)=z(ii)+(h*(11.d0*bone(ii)+81.d0*bthree(ii)+81.d0*bfour(ii)&
      &    -32.d0*bfive(ii)-32.d0*bsix(ii)+11.d0*bseven(ii)))/120.d0
    end do
    return
  end subroutine RK6
  
  subroutine var_RK6(fun,var_fun,x,v, tsid, h)
    external :: fun, var_fun
    real(kind=16), dimension(:), intent(inout) :: x,v
    real(kind=16) :: tsid, h
    real(kind=16), allocatable, dimension(:) :: z
    real(kind=16), allocatable, dimension(:) :: bone, btwo, bthree, bfour, bfive
    real(kind=16), allocatable, dimension(:) :: bsix, bseven
    real(kind=16), allocatable, dimension(:) :: b2temp, b3temp, b4temp, b5temp
    real(kind=16), allocatable, dimension(:) :: b6temp, b7temp
    real(kind=16), allocatable, dimension(:) :: cone, ctwo, cthree, cfour, cfive
    real(kind=16), allocatable, dimension(:) :: csix, cseven
    real(kind=16), allocatable, dimension(:) :: ctemp
    integer :: ii, dim
    dim=size(x)
    allocate(z(dim))
    allocate(bone(dim), btwo(dim), bthree(dim), bfour(dim), bfive(dim),bsix(dim),bseven(dim))
    allocate(b2temp(dim), b3temp(dim), b4temp(dim), b5temp(dim),b6temp(dim),b7temp(dim))
    allocate(cone(dim), ctwo(dim), cthree(dim), cfour(dim), cfive(dim),csix(dim),cseven(dim),ctemp(dim))
    ! First compute the solution of the ODE
    do ii=1,2
      z(ii)=x(ii)
    end do
    call fun(tsid, z, bone)
    do ii=1,2
      b2temp(ii)=z(ii)+(h/3.d0)*bone(ii)
    end do
    call fun(tsid+h/3.d0, b2temp, btwo)
    do ii=1,2
      b3temp(ii)=z(ii)+(2.d0*h/3.d0)*btwo(ii)
    end do
    call fun(tsid+2.d0*h/3.d0, b3temp, bthree)
    do ii=1,2
      b4temp(ii)=z(ii)+h*(bone(ii)+4.d0*btwo(ii)-bthree(ii))/12.d0
    end do
    call fun(tsid+h/3.d0, b4temp, bfour)
    do ii=1,2
      b5temp(ii)=z(ii)+h*(-bone(ii)+18.d0*btwo(ii) &
    &     -3.d0*bthree(ii)-6.d0*bfour(ii))/16.d0
    end do
    call fun(tsid+h/2.d0, b5temp, bfive)
    do ii=1,2
      b6temp(ii)=z(ii)+h*(9.d0*btwo(ii) &
    &     -3.d0*bthree(ii)-6.d0*bfour(ii)+4.d0*bfive(ii))/8.d0
    end do
    call fun(tsid+h/2.d0, b6temp, bsix)
    do ii=1,2
      b7temp(ii)=z(ii)+h*(9.d0*bone(ii)-36.d0*btwo(ii) &
    &     +63.d0*bthree(ii)+72.d0*bfour(ii)-64.d0*bfive(ii))/44.d0
    end do
    call fun(tsid+h, b7temp, bseven)
    do ii=1,2
      x(ii)=z(ii)+(h*(11.d0*bone(ii)+81.d0*bthree(ii)+81.d0*bfour(ii)&
      &    -32.d0*bfive(ii)-32.d0*bsix(ii)+11.d0*bseven(ii)))/120.d0
    end do
    ! Then compute the solution of the variational equation using
    ! the auxiliary solutions computed in the previous step
    ! var_fun(tsid,x,v,dv)
    call var_fun(tsid, z, v, cone)
    do ii=1,2
      ctemp(ii)=v(ii)+(h/3.d0)*cone(ii)
    end do
    call var_fun(tsid+h/3.d0, b2temp, ctemp, ctwo)
    do ii=1,2
      ctemp(ii)=v(ii)+(2.d0*h/3.d0)*ctwo(ii)
    end do
    call var_fun(tsid+2.d0*h/3.d0, b3temp,ctemp, cthree)
    do ii=1,2
      ctemp(ii)=v(ii)+h*(cone(ii)+4.d0*ctwo(ii)-cthree(ii))/12.d0
    end do
    call var_fun(tsid+h/3.d0, b4temp, ctemp,cfour)
    do ii=1,2
      ctemp(ii)=v(ii)+h*(-cone(ii)+18.d0*ctwo(ii) &
    &     -3.d0*cthree(ii)-6.d0*cfour(ii))/16.d0
    end do
    call var_fun(tsid+h/2.d0, b5temp,ctemp, cfive)
    do ii=1,2
      ctemp(ii)=v(ii)+h*(9.d0*ctwo(ii) &
    &     -3.d0*cthree(ii)-6.d0*cfour(ii)+4.d0*cfive(ii))/8.d0
    end do
    call var_fun(tsid+h/2.d0, b6temp,ctemp, csix)
    do ii=1,2
      ctemp(ii)=v(ii)+h*(9.d0*cone(ii)-36.d0*ctwo(ii) &
    &     +63.d0*cthree(ii)+72.d0*cfour(ii)-64.d0*cfive(ii))/44.d0
    end do
    call var_fun(tsid+h, b7temp,ctemp, cseven)
    do ii=1,2
      v(ii)=v(ii)+(h*(11.d0*cone(ii)+81.d0*cthree(ii)+81.d0*cfour(ii)&
      &    -32.d0*cfive(ii)-32.d0*csix(ii)+11.d0*cseven(ii)))/120.d0
    end do
    return
  end subroutine var_RK6

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  end module integration
  