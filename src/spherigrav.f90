!     Definitions of the required module for handling polyhedron dynamics
!**************************************************************************
    module spherigrav
    use constants
    use geometry
    implicit none
    
    contains ! compute_V_and_W, spheriacc_rot, spheripot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_V_and_W(pos, a, nMax, V, W)
    real(kind=16), dimension(3), intent(in) :: pos
    real(kind=16), intent(in) :: a
    integer, intent(in) :: nMax
    real(kind=16), allocatable, dimension(:,:), intent(out) :: V, W
    real(kind=16) :: r, x, y, z, coef
    integer :: i, j

    x=pos(1)
    y=pos(2)
    z=pos(3)
    r=sqrt(x**2+y**2+z**2)
    coef= a/r**2
    allocate(V(nMax+1,nMax+1))
    allocate(W(nMax+1,nMax+1))

    !    V=0.0d0*V
    !    W=0.0d0*W
    do i=1,nMax+1
    do j=1,nMax+1
    V(i,j)=0.0d0
    W(i,j)=0.0d0
    end do
    end do

    ! NOTE i=n+1, j=m+1, where n,m appear in the formulas 3.29,30,31
    ! Formula 3.31
    
    ! n=m=0
    V(1,1)=a/r
    W(1,1)=0.0d0
    ! n=1, m=0
    ! Formula 3.30
    V(2,1)=z*coef*V(1,1)
    W(2,1)=0.0d0    

    ! m=0, n>=2, i.e. i>=3
    do i=3,nMax+1

    ! Formula 3.30
    V(i,1)=(real(2*i-3)/real(i-1))*z*coef*V(i-1,1) -(real(i-2)/real(i-1))*coef*a*V(i-2,1)
    W(i,1)=0.0d0    ! As noted in Montenbruck

    end do


    do j=2,nMax+1
    ! Formula 3.29
    V(j,j)=real(2*j-3)*(x*V(j-1,j-1)-y*W(j-1,j-1))*coef
    W(j,j)=real(2*j-3)*(x*W(j-1,j-1)+y*V(j-1,j-1))*coef
    
    do i=j+1,nMax+1
        ! Formula 3.30 (Recall that n=i-1, m=j-1)
    !    V(i,j)= ((2*i-3)/(i-j))*z*coef*V(i-1,j) - ((i+j-3)/(i-j))*coef*a*V(i-2,j)
    !    W(i,j)= ((2*i-3)/(i-j))*z*coef*W(i-1,j) - ((i+j-3)/(i-j))*coef*a*W(i-2,j)
        V(i,j)= (real(2*i-3)/real(i-j)*z*V(i-1,j) - real(i+j-3)/real(i-j)*a*V(i-2,j))*coef
        W(i,j)= (real(2*i-3)/real(i-j)*z*W(i-1,j) - real(i+j-3)/real(i-j)*a*W(i-2,j))*coef
    end do

    end do

end subroutine compute_V_and_W
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spheriacc_rot(pos,a,M,C,S,rate,time,acc)
    ! acc is the resulting acceleration
    real(kind=16), dimension(3), intent(in) :: pos
    real(kind=16), intent(in) :: a, M, rate, time
    real(kind=16), allocatable, dimension(:,:), intent(in) :: C, S
    real(kind=16), dimension(3), intent(out) :: acc
    real(kind=16), allocatable, dimension(:,:) :: V, W
    integer :: i, j, nMax, nMax2
    real(kind=16), dimension(3) :: pos_aux, acc_aux
    real(kind=16) :: arg, arg2, acc_x, acc_y, acc_z, coef

    do i=1,3
    acc(i)=0.0d0
    acc_aux(i)=0.0d0*G*M/(sqrt(pos(1)**2+pos(2)**2+pos(3)**2))**3*pos(i)
    pos_aux(i)=0.d0
    end do
    
    nMax=size(C,dim=1)-1
    nMax2=nMax+1
    arg=rate*time
    arg2=-1.0d0*arg
    coef=G*M/a**2
    ! Step 1: Rotate pos of an angle -arg around the z axis
    call rotate_around_z(pos,arg2,pos_aux)
    ! Step 2: Compute the auxiliary quantities V and W at position pos_aux
    call compute_V_and_W(pos_aux,a,nMax2,V,W)
    ! Step 3: Compute the Cartesian acceleration acc using formulas 3.32 and 3.33
    do i=1,nMax+1
    do j=1,i
    if (j .eq. 1) then ! i.e. m=j-1=0
        acc_x=-coef*C(i,1)*V(i+1,2)
        acc_y=-coef*C(i,1)*W(i+1,2)
    else if (j .gt. 1) then
        acc_x=coef*0.5d0*((-C(i,j)*V(i+1,j+1)-S(i,j)*W(i+1,j+1))&
        & + real((i-j+2)*(i-j+1))*(C(i,j)*V(i+1,j-1)+S(i,j)*W(i+1,j-1)))
    !    & + (fact(i-j+2)/fact(i-j))*(C(i,j)*V(i+1,j-1)+S(i,j)*W(i+1,j-1)))


        acc_y=coef*0.5d0*((-C(i,j)*W(i+1,j+1)+S(i,j)*V(i+1,j+1))&
        & + real((i-j+2)*(i-j+1))*(-C(i,j)*W(i+1,j-1)+S(i,j)*V(i+1,j-1)))
    !    & + (fact(i-j+2)/fact(i-j))*(-C(i,j)*W(i+1,j-1)+S(i,j)*V(i+1,j-1)))
    end if
    acc_z=coef*(i-j+1)*(-C(i,j)*V(i+1,j)-S(i,j)*W(i+1,j))

    acc_aux(1)=acc_aux(1)+acc_x
    acc_aux(2)=acc_aux(2)+acc_y
    acc_aux(3)=acc_aux(3)+acc_z

    end do
    end do
    ! Step 3: Rotate acc of angle +arg around the z axis
    call rotate_around_z(acc_aux,arg,acc)
    return 
end subroutine spheriacc_rot

subroutine spheriggm_rot(pos,a,M,C,S,rate,time,ggm)
    ! acc is the resulting acceleration
    real(kind=16), dimension(3), intent(in) :: pos
    real(kind=16), intent(in) :: a, M, rate, time
    real(kind=16), dimension(3,3),intent(out) :: ggm
    real(kind=16), allocatable, dimension(:,:), intent(in) :: C, S
    real(kind=16), allocatable, dimension(:,:) :: V, W
    integer :: i, j, nMax, nMax2, p, q
    real(kind=16), dimension(3) :: pos_aux, acc_aux
    real(kind=16), dimension(3,3) :: rot_aux, ggm_aux
    real(kind=16) :: arg, arg2, acc_x, acc_y, acc_z, coef

    do p=1,3
        do q=1,3
            ggm(p,q)=0.d0
        end do
    end do
    
    nMax=size(C,dim=1)-1
    nMax2=nMax+1
    arg=rate*time !rate*(time-t0)
    arg2=-1.0d0*arg
    coef=G*M/a**3
    ! Step 1: Rotate pos of an angle -arg around the z axis
    call rotate_around_z(pos,arg2,pos_aux)
    ! pos_aux(3)=pos(3)
    ! pos_aux(1)=cos(arg2)*pos(1)-sin(arg2)*pos(2)
    ! pos_aux(2)=cos(arg2)*pos(2)+sin(arg2)*pos(1)
    ! Step 2: Compute the auxiliary quantities V and W at position pos_aux
    call compute_V_and_W(pos_aux,a,nMax2,V,W)
    ! Step 3: Compute the Cartesian acceleration acc using formulas 3.32 and 3.33

    do i=1,nMax+1
        do j=1,i
            if (j .eq. 1) then ! i.e. m=j-1=0
                ggm(1,1)=ggm(1,1)+coef*C(i,1)*(V(i+2,3)-(i+2)*(i+1)*V(i+2,1))/2.d0
                ggm(2,2)=ggm(2,2)+coef*C(i,1)*(-V(i+2,3)-(i+2)*(i+1)*V(i+2,1))/2.d0
                ggm(3,3)=ggm(3,3)+coef*C(i,1)*(i+2)*(i+1)*V(i+2,1)
                ggm(1,2)=ggm(1,2)+coef*C(i,1)*W(i+2,3)/2.d0
                ggm(1,3)=ggm(1,3)+coef*C(i,1)*(i+1)*V(i+2,2)
                ggm(2,3)=ggm(2,3)+coef*C(i,1)*(i+1)*W(i+2,2)
            else if (j .eq. 2) then ! i.e. m=j-2=1
                ggm(1,1)=ggm(1,1)+(coef/4)*(C(i,2)*(V(i+2,4)-6*(i+1)*i*V(i+2,2))+S(i,2)*(W(i+2,4)+2*(i+1)*i*W(i+2,2)))
                ggm(2,2)=ggm(2,2)+(coef/4)*(C(i,2)*(-V(i+2,4)+2*(i+1)*i*V(i+2,2))+S(i,2)*(-W(i+2,4)-6*(i+1)*i*W(i+2,2)))
                ggm(3,3)=ggm(3,3)+coef*(i+1)*i*(C(i,2)*V(i+2,2)+S(i,2)*W(i+2,2))
                ggm(1,2)=ggm(1,2)+(coef/4)*(C(i,2)*(W(i+2,4)-(i+1)*i*W(i+2,2))+S(i,2)*(-V(i+2,4)-(i+1)*i*V(i+2,2)))
                ggm(1,3)=ggm(1,3)+(coef/2)*(C(i,2)*(i*V(i+2,3)-(i+2)*(i+1)*i*V(i+2,1))+S(i,2)*W(i+2,3))
                ggm(2,3)=ggm(2,3)+(coef/2)*(C(i,2)*W(i+2,3)+S(i,2)*(i*V(i+2,3)-(i+2)*(i+1)*i*V(i+2,1)))
            else 
                ggm(1,1)=ggm(1,1)+(coef/4)*(C(i,j)*(V(i+2,j+2)-2*(i-j+2)*(i-j+1)*V(i+2,j)+(i-j+2)*(i-j+1)*V(i+2,j-2)) &
                + S(i,j)*(W(i+2,j+2)-2*(i-j+2)*(i-j+1)*W(i+2,j)+(i-j+2)*(i-j+1)*W(i+2,j-2)))
                ggm(2,2)=(ggm(2,2)+coef/4)*(C(i,j)*(-V(i+2,j+2)-2*(i-j+2)*(i-j+1)*V(i+2,j)-(i-j+2)*(i-j+1)*V(i+2,j-2)) &
                + S(i,j)*(-W(i+2,j+2)-2*(i-j+2)*(i-j+1)*W(i+2,j)-(i-j+2)*(i-j+1)*W(i+2,j-2)))
                ggm(3,3)=ggm(3,3)+coef*(C(i,j)*V(i+2,j) + S(i,j)*W(i+2,j))
                ggm(1,2)=ggm(1,2)+(coef/4)*(C(i,j)*(W(i+2,j+2)-(i-j+4)*(i-j+3)*(i-j+2)*(i-j+1)*W(i+2,j-2))&
                 + S(i,j)*(-V(i+2,j+2)+(i-j+4)*(i-j+3)*(i-j+2)*(i-j+1)*W(i+2,j-2)*V(i+2,j-2)))
                ggm(1,3)=ggm(1,3)+(coef/2)*(C(i,j)*((i-j+1)*V(i+2,j+2)-(i-j+3)*(i-j+2)*(i-j+1)*V(i+2,j-1))&
                 + S(i,j)*((i-j+1)*W(i+2,j+1)-(i-j+3)*(i-j+2)*(i-j+1)*W(i+2,j-1)))
                ggm(2,3)=ggm(2,3)+(coef/2)*(C(i,j)*((i-j+1)*W(i+2,j+2)+(i-j+3)*(i-j+2)*(i-j+1)*W(i+2,j-1))&
                 + S(i,j)*(-(i-j+1)*V(i+2,j+1)+(i-j+3)*(i-j+2)*(i-j+1)*V(i+2,j-1)))
            end if
    
        ggm(2,1)=ggm(1,2)
        ggm(3,1)=ggm(1,3)
        ggm(3,2)=ggm(2,3)
    
        end do
        end do
    
    ! Step 3: Rotate acc of angle +arg around the z axis
    ! call rotate_around_z(acc_aux,arg,acc)
    rot_aux(1,1)=cos(arg)
    rot_aux(2,1)=sin(arg)
    rot_aux(3,3)=1.d0
    rot_aux(2,2)=rot_aux(1,1)
    rot_aux(1,2)=-rot_aux(2,1)
    rot_aux(1,3)=0.d0
    rot_aux(2,3)=0.d0
    rot_aux(3,2)=0.d0
    rot_aux(3,1)=0.d0
    ggm_aux=matmul(rot_aux,ggm)

    ggm=ggm_aux

    rot_aux(2,1)=sin(-arg)
    rot_aux(1,2)=-rot_aux(2,1)
    ggm_aux=matmul(ggm,rot_aux)
    ggm=ggm_aux

    return 
end subroutine spheriggm_rot

subroutine spheriacc_rot_new(time,state,acc)
    ! acc is the resulting acceleration
    real(kind=16), dimension(6), intent(in) :: state
    real(kind=16), intent(in) ::  time
    real(kind=16), dimension(6), intent(out) :: acc
    real(kind=16), allocatable, dimension(:,:) :: V, W
    integer :: i, j, nMax, nMax2
    real(kind=16), dimension(3) :: pos, vel, pos_aux, acc_aux, acc_rot
    real(kind=16) :: arg, arg2, acc_x, acc_y, acc_z, coef


    do i=1,3
        pos(i)=state(i)
        vel(i)=state(i+3)
    end do

    ! do i=1,3
    ! acc(i)=0.0d0
    ! acc_aux(i)=0.0d0*G*M/(sqrt(pos(1)**2+pos(2)**2+pos(3)**2))**3*pos(i)
    ! pos_aux(i)=0.d0
    ! end do
    
    nMax=size(C,dim=1)-1
    nMax2=nMax+1
    time0=0.d0
    arg=rate*(time-time0)
    arg2=-1.0d0*arg
    ! write(*,*) arg
    ! write(*,*) arg2
    ! write(*,*) rr
    coef=G*Ma/rr**2
    ! write(*,*) coef
    ! coef=mu/rr**2
    ! write(*,*) coef
    ! stop
    ! Step 1: Rotate pos of an angle -arg around the z axis
    call rotate_around_z(pos,arg2,pos_aux)
    ! Step 2: Compute the auxiliary quantities V and W at position pos_aux
    call compute_V_and_W(pos_aux,rr,nMax2,V,W)
    ! Step 3: Compute the Cartesian acceleration acc using formulas 3.32 and 3.33
    do i=1,nMax+1
        do j=1,i
            if (j .eq. 1) then ! i.e. m=j-1=0
                acc_x=-coef*C(i,1)*V(i+1,2)
                acc_y=-coef*C(i,1)*W(i+1,2)
            else if (j .gt. 1) then
                acc_x=coef*0.5d0*((-C(i,j)*V(i+1,j+1)-S(i,j)*W(i+1,j+1))&
                & + dble((i-j+2)*(i-j+1))*(C(i,j)*V(i+1,j-1)+S(i,j)*W(i+1,j-1)))
            !    & + (fact(i-j+2)/fact(i-j))*(C(i,j)*V(i+1,j-1)+S(i,j)*W(i+1,j-1)))

                acc_y=coef*0.5d0*((-C(i,j)*W(i+1,j+1)+S(i,j)*V(i+1,j+1))&
                & + dble((i-j+2)*(i-j+1))*(-C(i,j)*W(i+1,j-1)+S(i,j)*V(i+1,j-1)))
            !    & + (fact(i-j+2)/fact(i-j))*(-C(i,j)*W(i+1,j-1)+S(i,j)*V(i+1,j-1)))
            end if
            acc_z=coef*dble(i-j+1)*(-C(i,j)*V(i+1,j)-S(i,j)*W(i+1,j))

            acc_aux(1)=acc_aux(1)+acc_x
            acc_aux(2)=acc_aux(2)+acc_y
            acc_aux(3)=acc_aux(3)+acc_z
        end do
    end do
    ! Step 3: Rotate acc of angle +arg around the z axis
    call rotate_around_z(acc_aux,arg,acc_rot)

    do i=1,3
        acc(i)=vel(i)
        acc(i+3)=acc_rot(i)
    end do
    return 
end subroutine spheriacc_rot_new
end module spherigrav

