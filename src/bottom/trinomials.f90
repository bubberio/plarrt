!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Module for trinomials handling                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module trinomials
implicit none 
public trin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definitions of custom data type                  !
! for handling trinomials                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! CHECKKK
type trin
integer :: deg
real(kind=16), allocatable, dimension(:) :: coeff
contains 
procedure :: trinloc
procedure :: TriMult, TriSum, TriDiff
procedure :: TriScal, TriNormSq
end type
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function for indexing auxiliary                  !
! vector of coefficients                           !
! input: i,j,k - integers associated to a trinomial!
!            n - degree of the trinomial           !
! output: ind - index corresponding to (i,j,k) in  !
!               auxiliary vector.                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function indexing(i,j,k,n) result (ind)
    integer i,j,k,n,ind
    ind = (j+k)*(j+k+1)/2 + k+1
end function indexing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to allocate auxiliary vector          !
! input: d - degree of the trinomial to initialize !
! output: t - initialized trinomial	               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine trinloc(t, d)
    integer d
    class(trin) :: t
    t%deg=d
    !	write(*,*) "Before Allocation"
    allocate (t%coeff((d+1)*(d+2)/2))
    !	write(*,*) "After Allocation"
    !	t%coeff=0.0d0
end subroutine trinloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Implement multiplication between trinomials      !
! input:  a,b - trinomials                         !
! output: c - trinomial                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function TriMult(a,b) result (c)
    class(trin), intent(in) :: a,b
    type(trin) :: c
    integer i,j,k,r,s,t,n1,n2,n3,inda,indb,indc ! Auxiliary indices
    n1= a%deg
    n2= b%deg
    n3= n1+n2
    i=n1
    call trinloc(c, n3)
    c%coeff=0.0
    do while (i .ge. 0)
    k=0
    do while (k .le. n1-i)
    j = n1 - i - k
    inda =indexing(i,j,k,n1)
    r=n2
    do while (r .ge. 0)
    t=0
    do while (t .le. n2-r)
    s = n2 - r - t
    indb =indexing(r,s,t,n2)
    indc=indexing(i+r,j+s,k+t,n3)
    c%coeff(indc)=c%coeff(indc) +a%coeff(inda)*b%coeff(indb)
    t = t+1
    end do
    r=r-1
    end do
    k = k+1
    end do
    i=i-1
    end do
end function TriMult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function to compute the multiplication between a !
! trinomial and a scalar.                          !
! input:  a - trinomial                            !
!         b - scalar                               !
! output: c - trinomial	                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function TriScal(a,b) result (c)
    class(trin), intent(in)  :: a
    real(kind=16), intent(in) :: b
    type(trin) :: c
    type(trin) :: b_trin

    ! Auxiliary trinomial
    call trinloc(b_trin, 0)
    b_trin%coeff(1)=b

    c=TriMult(a,b_trin)
    end function TriScal

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Function to compute the sum between trinomials.  !
    ! input:  a,b - trinomials                         !
    ! output: c - trinomial                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function TriSum(a,b) result (c)
    class(trin), intent(in)  :: a,b
    type(trin) :: c

    ! Auxiliary indices

    if (a%deg .ne. b%deg) then
    write(*,*) "Error(+): Trinomials have different degree"
    write(*,*) "deg a = ", a%deg ,"deg b = ", b%deg
    return
    end if
    call trinloc(c, a%deg)
    c%coeff= a%coeff+b%coeff
    return 
end function TriSum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function to compute the difference               !
! between trinomials.                              !
! input:  a,b - trinomials                         !
! output: c - trinomial                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function TriDiff(a,b) result (c)
    class(trin), intent(in)  :: a,b
    type(trin) :: c
    ! Auxiliary indices
    if (a%deg .ne. b%deg) then
    write(*,*) "Error(-): Trinomials have different degree"
    write(*,*) "deg a = ", a%deg ,"deg b = ", b%deg 
    return
    end if
    call trinloc(c, a%deg)

    c%coeff= a%coeff-b%coeff
end function TriDiff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function to compute the norm of a given trinomial! 
! input:  a,b,c - trinomials                       !
! output: d - trinomial                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function TriNormSq(a,b,c) result (d)
    class(trin), intent(in)  :: a,b,c
    type(trin) :: d
    if (a%deg .ne. b%deg) then
    write(*,*) "Error(||): Trinomials have different degree"
    write(*,*) "deg a = ", a%deg ,"deg b = ", b%deg 
    return
    end if

    call trinloc(d, 2*(a%deg))
    d=TriSum(TriSum(TriMult(a,a),TriMult(b,b)),TriMult(c,c))
end function TriNormSq

end module trinomials
