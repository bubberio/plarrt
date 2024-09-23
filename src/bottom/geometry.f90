!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Module containing geometric tools       	     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module geometry
implicit none
contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to extract vertices and faces data    !
! from an obj file                         	       !
! input: address - address of the file             !
!         nVert - number of vertices       	       !
!         nFaces - number of faces        	       !
! output: vv - matrix of vertices coordinates      !
!         ff - matrix of faces             	       !
!         a - radius of Brillouin Sphere  	       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine extract_ver_and_fac(address,nVert,nFaces,sca,vv,ff,a)
    ! Also computes the reference radius!
    implicit none 
    character(len=100), intent(in) :: address
    character(len=1) :: ident
    integer, intent(in) :: nVert, nFaces
    real(kind=16), intent(in) :: sca
    real(kind=16), allocatable, dimension(:,:), intent(out) :: vv
    integer, allocatable, dimension(:,:), intent(out) :: ff
    real(kind=16), intent(out) :: a
    integer :: i_v, i_f, v1, v2, v3
    real(kind=16) :: a_aux, x,y,z
    allocate (vv(nVert,3))
    allocate (ff(nFaces,3))

    open (10, file=address, form='formatted')
    ! Initialize indexes for do while loops
    i_v = 1
    i_f = 1
    ! Store the data from the .obj file into matrices
    ! Each row of vv contains the coordinates of a single vertex
    read(10,*) ident
    if (ident .eq. '#') then
    read(10,*) ident
    else
    close(10)
    open (10, file=address, form='formatted')
    end if
    do while (i_v .le. nVert)
        read(10,*) ident, x, y, z
        if (i_v .eq. 1) then
            a=sqrt(x**2 + y**2 + z**2 )
            a=a*sca
        end if
        a_aux=sqrt(x**2 + y**2 + z**2 )
        a_aux=a_aux*sca
        vv(i_v,1) = x*sca
        vv(i_v,2) = y*sca
        vv(i_v,3) = z*sca
        i_v = i_v+1
        if (a_aux .ge. a) then
            a=a_aux
        end if
    end do

    do while (i_f .le. nFaces)
        read(10,*) ident, v1, v2, v3
        ff(i_f,1) = v1
        ff(i_f,2) = v2
        ff(i_f,3) = v3
        i_f = i_f+1
    end do
    close(10)
end subroutine extract_ver_and_fac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to extract the edges in a (nEdges,4)  !
! matrix, including the faces              	       !
! input:  ff - matrix of faces             	       !
! output: ee - matrix of edges            	       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine extract_edges(ff,ee)
    implicit none
    integer, allocatable, dimension(:,:), intent(in) :: ff
    integer, allocatable, dimension(:,:), intent(out) :: ee
    integer, allocatable, dimension(:,:) :: ee_aux
    integer :: aux, i, j, k, l,last, last_i, n, m

    aux=size(ff,dim=1)
    allocate(ee_aux(3*aux,3))
    allocate(ee(3*aux/2,4))

    do i=1, aux
        j=ff(i,1)
        k=ff(i,2)
        l=ff(i,3)

        ee_aux(3*(i-1)+1,1)= j
        ee_aux(3*(i-1)+1,2)= k
        ee_aux(3*(i-1)+1,3)= i

        ee_aux(3*(i-1)+2,1)= k
        ee_aux(3*(i-1)+2,2)= l
        ee_aux(3*(i-1)+2,3)= i

        ee_aux(3*(i-1)+3,1)= l
        ee_aux(3*(i-1)+3,2)= j
        ee_aux(3*(i-1)+3,3)= i
    end do

    ee(1,1)= ee_aux(1,1)
    ee(1,2)= ee_aux(1,2)
    ee(1,3)= ee_aux(1,3)

    ee(2,1)= ee_aux(2,1)
    ee(2,2)= ee_aux(2,2)
    ee(2,3)= ee_aux(2,3)

    ee(3,1)= ee_aux(3,1)
    ee(3,2)= ee_aux(3,2)
    ee(3,3)= ee_aux(3,3)

    last_i=3
    do i=4, 3*aux
        last=last_i
        j=ee_aux(i,1)
        k=ee_aux(i,2)
        do n=1 , last
            if (k .eq. ee(n,1) .and. (j .eq. ee(n,2))) then
                ee(n,4)= int(ceiling(i/3.0))
                exit
            end if
            if (n .eq. last) then
                ee(last_i+1,1)= j
                ee(last_i+1,2)= k
                ee(last_i+1,3)= int(ceiling(i/3.0))
                last_i=last_i+1
            end if
        end do
    end do
end subroutine extract_edges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to extract the face normals           !
! input: vv - matrix of vertices coordinates       !
!         ff - matrix of faces             	       !
! output: ee - matrix of edges            	       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine face_normals(vv,ff,nnff)
    implicit none 
    real(kind=16), allocatable, dimension(:,:),intent(in) :: vv
    integer, allocatable, dimension(:,:), intent(in) :: ff
    real(kind=16), allocatable, dimension(:,:),intent(out) :: nnff
    real(kind=16), dimension(3) :: v1, v2, v3, nf
    integer :: i

    allocate(nnff(size(ff,dim=1),3))

    do i=1, size(ff,dim=1)
    v1=vv(ff(i,1),:)
    v2=vv(ff(i,2),:)
    v3=vv(ff(i,3),:)
    call vector_product(v2-v1,v3-v2,nf)
    nf=nf/sqrt(dot_product(nf,nf))
    nnff(i,:)=nf
    end do

end subroutine face_normals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to shift the .obj file in the         !
! computed center of mass                          !
! input: new_address - where to save shifted file  !
!        com - vector including coordinates of com !
!        sca - scaling parameter                   !
! output: vv - matrix of vertices coordinates      !
!         ff - matrix of faces             	       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine shiftCoM(new_address,com,vv,ff,sca)
    implicit none
    real(kind=16), dimension(3), intent(in) :: com
    CHARACTER(LEN=100) :: new_address
    integer :: nVert, nFaces, i_v,i_f,v1,v2,v3
    ! Data array to be evaluated
    real(kind=16), allocatable, dimension(:,:),intent(inout) :: vv 
    integer, allocatable, dimension(:,:),intent(inout) :: ff 
    real(kind=16) sca ! Should be an input?

    v1=1
    v2=1
    v3=1
    nVert=size(vv,dim=1)
    nFaces=size(ff,dim=1)
    open (3, file=new_address, status='unknown', form='formatted')

    write(3,*) '#', nVert  ! Write the number of vertices in the new file
    write(3,*) '#', nFaces ! Write the number of faces in the new file

    ! Initialize indexes for do while loops
    i_v = 1
    i_f = 1
    ! Store the data from the .obj file into matrices
    ! Each row of vv contains the coordinates of a single vertex
    do while (i_v .le. nVert)

    vv(i_v,1) = vv(i_v,1) - com(1)
    vv(i_v,2) = vv(i_v,2) - com(2)
    vv(i_v,3) = vv(i_v,3) - com(3)
    write(3,*) 'v', vv(i_v,1)/sca, vv(i_v,2)/sca, vv(i_v,3)/sca
    i_v = i_v+1
    end do
    ! Each row of ff contains the labels for 
    ! the vertices which define a simplex
    do while (i_f .le. nFaces)
    v1 = ff(i_f,1)
    v2 = ff(i_f,2)
    v3 = ff(i_f,3)
    write(3,*) 'f', v1, v2, v3
    i_f = i_f+1
    end do

    close(3)

end subroutine shiftCoM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to shift the .obj file in the         !
! computed center of mass                          !
! input:  vv - matrix of vertices coordinates      !
!         ff - matrix of faces             	       !
! output: ee - matrix of edges            	       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine count_v_and_f(address,nVert,nFaces)
    !!!!!!!!!!!!! SOME MISTAKES HERE (Shouldn't)
    implicit none
    CHARACTER(LEN=100) :: address
    CHARACTER(LEN=1) :: ident 
    integer :: i_v,i_f,v1,v2,v3,IOstatus
    integer :: nVert,nFaces
    real(kind=8) :: x,y,z
    IOstatus=0
    open (3, file=address, form='formatted')
    read(3,*) ident, x,y,z
    i_v=1
    do while (ident .eq. 'v')
    read(3,*,IOSTAT=IOstatus) ident, x,y,z
    if (ident .eq. 'v') then
    i_v=i_v+1
    end if
    end do
    nVert=i_v
    i_f=1

    do while (ident .eq. 'f')
    if (IOstatus .ge. 0) then
    read(3,*,IOSTAT=IOstatus) ident, v1,v2,v3
    i_f=i_f+1
    else
    exit
    end if

    end do
    nFaces=i_f-1
    close(3)
end subroutine count_v_and_f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to compute the vector product         !
! between two three-dimensional vectors            !
! input: a, b - 3d vectors                         !
! output: c - output vector            	           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vector_product(a,b,c)
    implicit none 
    real(kind=16), dimension(3), intent(in) :: a,b
    real(kind=16), dimension(3), intent(out) :: c

    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
    return
end subroutine vector_product

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to rotate a given vector around the   !
! z-axis by a given angle                          !
! input: vec - 3d vector                           !
!        arg - rotation angle in radians           !
! output: rotvec - output vector                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rotate_around_z(vec,arg,rotvec)
    real(kind=16), dimension(3), intent(in) :: vec
    real(kind=16), intent(in) :: arg
    real(kind=16), dimension(3), intent(out) :: rotvec

    rotvec(1)=cos(arg)*vec(1)-sin(arg)*vec(2)
    rotvec(2)=cos(arg)*vec(2)+sin(arg)*vec(1)
    rotvec(3)=vec(3)
    return
end subroutine rotate_around_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to rotate all the vectors in an .obj  !
! around the z-axis by a given angle               !
! input: vv - matrix of vertices coordinates       !
!        arg - rotation angle in radians           !
! output: vv_out - output matrix                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rotate_vertices(vv,arg,vv_out)
    implicit none 
    real(kind=16), intent(in) :: arg
    real(kind=16), allocatable, dimension(:,:), intent(in) :: vv
    real(kind=16), allocatable, dimension(:,:), intent(out) :: vv_out
    integer :: i
    real(kind=16) :: vx, vy
    allocate(vv_out(size(vv,dim=1),3))
    do i=1,size(vv,dim=1)
    vx=vv(i,1)
    vy=vv(i,2)
    vv_out(i,1)=cos(arg)*vx - sin(arg)*vy
    vv_out(i,2)=cos(arg)*vy + sin(arg)*vx
    vv_out(i,3)=vv(i,3)
    end do
    return
end subroutine

end module geometry
