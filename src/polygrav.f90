!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Module containing routines necessary to        !
!   implement polyhedron dynamics                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module polygrav
    use trinomials
    use geometry
    use constants
    implicit none
    
    contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_Le(fp,v1,v2, L)
    real(kind=16), dimension(3), intent(in) :: fp, v1, v2
    real(kind=16), intent(out) :: L

    real(kind=16) :: e,a,b
    
    e=sqrt(dot_product(v1-v2,v1-v2))
    a=sqrt(dot_product(fp-v1,fp-v1))
    b=sqrt(dot_product(fp-v2,fp-v2))
    L=log((a+b+e)/(a+b-e))
    return
end subroutine get_Le
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_omega(fp,v1,v2,v3,omega)
    real(kind=16), dimension(3), intent(in) :: fp, v1, v2, v3
    real(kind=16), intent(out) :: omega
    real(kind=16), dimension(3) :: auxvec, rv1, rv2, rv3
    real(kind=16) :: num, den, r1, r2, r3
    
    rv1= v1-fp
    rv2= v2-fp
    rv3= v3-fp
    r1=sqrt(dot_product(rv1,rv1))
    r2=sqrt(dot_product(rv2,rv2))
    r3=sqrt(dot_product(rv3,rv3))
    
    call vector_product(rv2,rv3,auxvec)
    num=dot_product(rv1,auxvec)
    den= r1*r2*r3 + r1*dot_product(rv2,rv3) + r2*dot_product(rv3,rv1) &
    &  + r3*dot_product(rv1,rv2)
    omega= 2*atan2(num,den) 
    return
end subroutine get_omega
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine edge_dyads(i,ee,vv,nnff,dyad_E)
    integer, intent(in) :: i
    integer, allocatable, dimension(:,:),intent(in) :: ee
    real(kind=16), allocatable, dimension(:,:),intent(in) :: vv,nnff
    real(kind=16), dimension(3,3), intent(out) :: dyad_E
    integer :: iv1, iv2, iA, iB
    real(kind=16), dimension(3) :: v1, v2, nA, nB, nA12, nB21

    iv1=ee(i,1)
    iv2=ee(i,2)
    iA=ee(i,3)
    iB=ee(i,4)
    nA=nnff(iA,:)
    nB=nnff(iB,:)
    v1=vv(iv1,:)
    v2=vv(iv2,:)
    call vector_product(v2-v1, nA,nA12)
    nA12=nA12/sqrt(dot_product(nA12,nA12))
    call vector_product(v1-v2, nB,nB21)
    nB21=nB21/sqrt(dot_product(nB21,nB21))
    dyad_E=matmul(reshape(nA,(/3,1/)),reshape(nA12,(/1,3/))) &
        &+ matmul(reshape(nB,(/3,1/)),reshape(nB21,(/1,3/)))
    return
end subroutine edge_dyads	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine face_dyads(i,nnff,dyad_F)
    integer, intent(in) :: i
    real(kind=16), allocatable, dimension(:,:),intent(in) :: nnff
    real(kind=16), dimension(3,3), intent(out) :: dyad_F

    real(kind=16), dimension(3) :: nf
    nf=nnff(i,:)
    dyad_F=matmul(reshape(nf, (/3,1/)),reshape(nf, (/1,3/)))
    return
end subroutine face_dyads
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine poly_all(fp,vv,ff,ee,sigma,rate,time,pot,acc,ggmat,lap)
    real(kind=16), intent(in) :: sigma, time, rate
    real(kind=16),dimension(3), intent(in) :: fp
    real(kind=16), allocatable, dimension(:,:), intent(in) :: vv
    integer, allocatable, dimension(:,:), intent(in) :: ff
    integer, allocatable, dimension(:,:), intent(in) :: ee
    real(kind=16), intent(out) :: pot, lap
    real(kind=16),dimension(3), intent(out) :: acc
    real(kind=16),dimension(3,3), intent(out) :: ggmat
    integer :: auxf,auxe, i, j,k
    real(kind=16), allocatable, dimension(:,:) :: nnff
    real(kind=16), allocatable, dimension(:,:) :: vv_rot
    real(kind=16), dimension(3) :: v1, v2, v3, rf, re
    real(kind=16), dimension(3,3) :: dyad_F, dyad_E
    real(kind=16), dimension(3) :: acc_aux
    real(kind=16), dimension(1,1) :: pot_aux
    real(kind=16) :: omega, L, arg
    
    arg=rate*(time-t0)
    auxf= size(ff,dim=1)
    auxe= size(ee,dim=1)
    call rotate_vertices(vv,arg,vv_rot)
    call face_normals(vv_rot,ff,nnff) 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialization of output vectors 				   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pot=0.d0
    lap=0.d0
    do j=1,3
    acc(j)= 0.0d0
    end do

    do j=1,3
        do k=1,3
        ggmat(j,k)=0.d0
        end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Iterate over the faces
    do i=1, auxf
        v1=vv_rot(ff(i,1),:)
        v2=vv_rot(ff(i,2),:)
        v3=vv_rot(ff(i,3),:)
        rf=(v1+v2)/2 - fp !	rf=vv_rot(ff(i,1),:)-fp	
        call face_dyads(i,nnff,dyad_F)
        call get_omega(fp,v1,v2,v3,omega)
        pot_aux= matmul(reshape(rf, (/1,3/)),matmul(dyad_F,reshape(rf,(/3,1/))))
        pot= pot-0.5d0*G*sigma*pot_aux(1,1)*omega
        lap=lap-G*sigma*omega 
        acc_aux=reshape( matmul(dyad_F,reshape(rf,(/3,1/))),(/3/))
        do j=1,3
            acc(j)= acc(j) + G*sigma*acc_aux(j)*omega
        end do
        do j=1,3
            do k=1,3
            ggmat(j,k)=ggmat(j,k) - G*sigma*omega*dyad_F(j,k)
            end do
        end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Iterate over the edges
    do i=1, auxe
        v1=vv_rot(ee(i,1),:)
        v2=vv_rot(ee(i,2),:)
        re=(v1+v2)/2 -fp
        call get_Le(fp,v1,v2, L)
        call edge_dyads(i,ee,vv_rot,nnff,dyad_E)
        pot_aux = matmul(reshape(re, (/1,3/)),matmul(dyad_E,reshape(re,(/3,1/))))
        pot= pot+0.5d0*G*sigma*pot_aux(1,1)*L
        acc_aux=reshape( matmul(dyad_E,reshape(re,(/3,1/))),(/3/))
        do j=1,3
        acc(j)= acc(j) - G*sigma*acc_aux(j)*L
        end do

        do j=1,3
            do k=1,3
            ggmat(j,k)=ggmat(j,k) + G*sigma*L*dyad_E(j,k)
            end do
        end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    return 
end subroutine poly_all
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine polyacc_rot(fp,vv,ff,ee,sigma,rate,time,acc_out)
    ! acc is the resulting acceleration
    real(kind=16), intent(in) :: sigma, time, rate
    real(kind=16),dimension(3), intent(in) :: fp
    real(kind=16), allocatable, dimension(:,:), intent(in) :: vv
    integer, allocatable, dimension(:,:), intent(in) :: ff
    integer, allocatable, dimension(:,:), intent(in) :: ee
    real(kind=16),dimension(3), intent(out) :: acc_out
    real(kind=16),dimension(3) :: acc
    integer :: auxf,auxe, i, j
    real(kind=16), allocatable, dimension(:,:) :: nnff
    real(kind=16), allocatable, dimension(:,:) :: vv_rot
    real(kind=16), dimension(3) :: v1, v2, v3, rf, re
    real(kind=16), dimension(3,3) :: dyad_F, dyad_E
    real(kind=16), dimension(3) :: acc_aux
    real(kind=16) :: omega, L, arg
    
    arg=rate*(time-time0)
    auxf= size(ff,dim=1)
    call rotate_vertices(vv,arg,vv_rot)
    call face_normals(vv_rot,ff,nnff)

    do j=1,3
    acc(j)= 0.0d0
    end do

    do i=1, auxf
    v1=vv_rot(ff(i,1),:)
    v2=vv_rot(ff(i,2),:)
    v3=vv_rot(ff(i,3),:)
    rf=(v1+v2)/2 - fp !	rf=vv_rot(ff(i,1),:)-fp	
    call face_dyads(i,nnff,dyad_F)
    call get_omega(fp,v1,v2,v3,omega)
    acc_aux=reshape( matmul(dyad_F,reshape(rf,(/3,1/))),(/3/))
    do j=1,3
    acc(j)= acc(j) + G*sigma*acc_aux(j)*omega
    end do
    end do
    auxe= size(ee,dim=1)
    do i=1, auxe
    v1=vv_rot(ee(i,1),:)
    v2=vv_rot(ee(i,2),:)
    re=(v1+v2)/2 -fp
    call get_Le(fp,v1,v2, L)
    call edge_dyads(i,ee,vv_rot,nnff,dyad_E)
    acc_aux=reshape( matmul(dyad_E,reshape(re,(/3,1/))),(/3/))
    do j=1,3
    acc(j)= acc(j) - G*sigma*acc_aux(j)*L
    end do
    end do
    do j=1,3
    acc_out(j)=acc(j)
    end do
    return 
end subroutine polyacc_rot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ggm(fp,vv,ff,ee,sigma,rate,time,ggmat)
    ! acc is the resulting acceleration
    real(kind=16), intent(in) :: sigma, time, rate
    real(kind=16),dimension(3), intent(in) :: fp
    real(kind=16), allocatable, dimension(:,:), intent(in) :: vv
    integer, allocatable, dimension(:,:), intent(in) :: ff
    integer, allocatable, dimension(:,:), intent(in) :: ee
    real(kind=16),dimension(3,3), intent(out) :: ggmat
    integer :: auxf,auxe, i, j,k
    real(kind=16), allocatable, dimension(:,:) :: nnff
    real(kind=16), allocatable, dimension(:,:) :: vv_rot
    real(kind=16), dimension(3) :: v1, v2, v3
    real(kind=16), dimension(3,3) :: dyad_F, dyad_E
    real(kind=16) :: omega, L, arg
    
    arg=rate*(time-t0)
    auxf= size(ff,dim=1)
    call rotate_vertices(vv,arg,vv_rot)
    call face_normals(vv_rot,ff,nnff)

    do j=1,3
        do k=1,3
        ggmat(j,k)=0.d0
        end do
    end do

    do i=1, auxf
        v1=vv_rot(ff(i,1),:)
        v2=vv_rot(ff(i,2),:)
        v3=vv_rot(ff(i,3),:)
        call face_dyads(i,nnff,dyad_F)
        call get_omega(fp,v1,v2,v3,omega)
        do j=1,3
            do k=1,3
            ggmat(j,k)=ggmat(j,k) - G*sigma*omega*dyad_F(j,k)
            end do
        end do

    end do
    auxe= size(ee,dim=1)
    do i=1, auxe
        v1=vv_rot(ee(i,1),:)
        v2=vv_rot(ee(i,2),:)
        call get_Le(fp,v1,v2, L)
        call edge_dyads(i,ee,vv_rot,nnff,dyad_E)
        do j=1,3
            do k=1,3
            ggmat(j,k)=ggmat(j,k) + G*sigma*L*dyad_E(j,k)
            end do
        end do
    end do

    return 
end subroutine ggm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine polypot(sigma,fp,vv,ff,pot)
    real(kind=16), intent(in) :: sigma
    real(kind=16),dimension(3), intent(in) :: fp
    real(kind=16), allocatable, dimension(:,:), intent(in) :: vv
    integer, allocatable, dimension(:,:), intent(in) :: ff
    real(kind=16),intent(out) :: pot
    real(kind=16), allocatable, dimension(:,:) :: nnff
    integer, allocatable, dimension(:,:) :: ee
    integer :: aux, i
    real(kind=16),dimension(3) :: nf, v1, v2, v3, rf, re
    real(kind=16), dimension(3,3) :: dyad_F, dyad_E
    real(kind=16), dimension(1,1) :: pot_aux	
    real(kind=16) :: omega,L	

    aux= size(ff,dim=1)
    call face_normals(vv,ff,nnff)
    call extract_edges(ff,ee)
    
    pot=0.0d0
    do i=1, aux
        rf=vv(ff(i,1),:)-fp	
        v1=vv(ff(i,1),:)
        v2=vv(ff(i,2),:)
        v3=vv(ff(i,3),:)
        call face_dyads(i,nnff,dyad_F)
        call get_omega(fp,v1,v2,v3,omega)
        pot_aux= matmul(reshape(rf, (/1,3/)),matmul(dyad_F,reshape(rf,(/3,1/))))
        pot= pot-0.5d0*G*sigma*pot_aux(1,1)*omega
    end do

    aux= size(ee,dim=1)
    do i=1, aux
        re=vv(ee(i,1),:)-fp	
        v1=vv(ee(i,1),:)
        v2=vv(ee(i,2),:)
        call get_Le(fp,v1,v2, L)
        call edge_dyads(i,ee,vv,nnff,dyad_E) 
        pot_aux = matmul(reshape(re, (/1,3/)),matmul(dyad_E,reshape(re,(/3,1/))))
        pot= pot+0.5d0*G*sigma*pot_aux(1,1)*L
    end do
end subroutine polypot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine polylap(sigma,fp,vv,ff,lap)
    ! lap is the resulting laplacian
    real(kind=16), intent(in) :: sigma
    real(kind=16), intent(out) :: lap
    real(kind=16),dimension(3), intent(in) :: fp
    real(kind=16), allocatable, dimension(:,:), intent(in) :: vv
    integer, allocatable, dimension(:,:), intent(in) :: ff
    real(kind=16),dimension(3) ::  v1, v2, v3
    real(kind=16) :: omega
    integer :: aux,i
    lap=0.0d0
    aux= size(ff,dim=1)
    do i=1, aux
    v1=vv(ff(i,1),:)
    v2=vv(ff(i,2),:)
    v3=vv(ff(i,3),:)
    call get_omega(fp,v1,v2,v3,omega)
    lap=lap-G*sigma*omega 
    end do
end subroutine polylap
end module polygrav
    
    