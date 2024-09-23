module SHcomputer
    use trinomials
    use constants
    use geometry
    use OMP_LIB
    implicit none 
        
        contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine denormalize_C_and_S(C,S)
    real(kind=16), allocatable, dimension (:,:), intent(inout) :: C, S
    integer :: nMax, i, j, n, m
    real(kind=16) :: Nnm
    nMax=size(C,dim=1)-1	
    do i=1,nMax+1
        do j=1,i
            n= i-1
            m= j-1
            if (m .eq. 0) then
            Nnm=sqrt(real(2*n+1)*fact(n-m)/fact(n+m))
            else
            Nnm=sqrt(2*real(2*n+1)*fact(n-m)/fact(n+m))
            end if
            C(i,j)=Nnm*C(i,j)
            S(i,j)=Nnm*S(i,j)
            ! write(*,*) C(i,j), S(i,j)
        end do
    end do

end subroutine denormalize_C_and_S
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_C_and_S(C,S,output_address,flag_denorm,flag_no_sigma)
    real(kind=16), allocatable, dimension (:,:), intent(in) :: C, S
    character(len=100), intent(in) :: output_address
    logical, intent(in) :: flag_denorm, flag_no_sigma
    character(len=30) :: fmt
    integer :: n, m, nMax
    nMax=size(C,dim=1)-1
    fmt='(I3,2X,I3,2(2X,ES18.7))'
    open(1,file=output_address,status='unknown', form='formatted')
    if (flag_no_sigma) then
        write(*,*) "No available info on mass or density."
    end if
    if (flag_denorm) then
        write(1,*) "Not Normalized"
    else
        write(1,*) "Normalized"
    end if
    write(1,*) "n:		m:		C:		S:"
    write(1,*) "***********************************************"
    do n=1,nMax+1
        do m=1,n
        write(1,fmt) n-1, m-1, C(n,m), S(n,m)
        end do
    end do
    close(1)
end subroutine write_C_and_S
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_C_and_S_and_rr(C,S,rr,Ma,output_address,flag_denorm,flag_no_sigma)
    real(kind=16), allocatable, dimension (:,:), intent(in) :: C, S
    real(kind=16), intent(in) :: rr, Ma ! reference radius
    character(len=100), intent(in) :: output_address
    logical, intent(in) :: flag_denorm, flag_no_sigma
    character(len=30) :: fmt
    integer :: n, m, nMax
    nMax=size(C,dim=1)-1
    fmt='(I3,2X,I3,2(2X,ES18.7))'
    open(1,file=output_address,status='unknown', form='formatted')
    if (flag_no_sigma) then
        write(*,*) "No available info on mass or density."
    end if
    if (flag_denorm) then
        write(1,*) "Not Normalized"
    else
        write(1,*) "Normalized"
    end if
    write(1,"(A1,2X,F16.6)") "r",rr
    write(1,"(A5,2X,F32.6)") "Mass:", Ma
    write(1,*) "n:		m:		C:		S:"
    do n=1,nMax+1
        do m=1,n
        write(1,fmt) n-1, m-1, C(n,m), S(n,m)
        end do
    end do
    close(1)
end subroutine write_C_and_S_and_rr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_C_and_S(input_address,nMax,C,S)
    character(len=100), intent(in) :: input_address
    integer, intent(in) :: nMax
    real(kind=16), allocatable, dimension (:,:), intent(out) :: C, S
    character(len=30) :: fmt
    integer :: n, m, IOstatus
    real(kind=16) :: Caux, Saux
    IOstatus=0
    allocate(C(nMax+1,nMax+1))
    allocate(S(nMax+1,nMax+1))
    open(10, file=input_address, form='formatted')
    read(10,*,IOSTAT=IOstatus)
    read(10,*,IOSTAT=IOstatus)
    read(10,*,IOSTAT=IOstatus) n, m, Caux, Saux
    C(n+1,m+1)=Caux
    S(n+1,m+1)=Saux
    do while  (n .le. nMax)
    read(10,*,IOSTAT=IOstatus) n, m, Caux, Saux
    if (IOstatus .ge. 0) then
    C(n+1,m+1)=Caux
    S(n+1,m+1)=Saux
    else
    go to 55
    end if
    end do
    55 	continue
    close(10)
end subroutine read_C_and_S
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_C_and_S_and_rr(input_address,nMax,C,S,rr,Ma)
    character(len=100), intent(in) :: input_address
    integer, intent(in) :: nMax
    real(kind=16), allocatable, dimension (:,:), intent(out) :: C, S
    real(kind=16),  intent(out) :: rr, Ma
    character(len=30) :: fmt
    character(len=4) :: chma
    character(len=1) :: hash
    integer :: n, m, IOstatus
    real(kind=16) :: Caux, Saux
    IOstatus=0
    allocate(C(nMax+1,nMax+1))
    allocate(S(nMax+1,nMax+1))
    open(10, file=input_address, form='formatted')
    read(10,*,IOSTAT=IOstatus)
    read(10,*,IOSTAT=IOstatus) hash, rr
    read(10,*,IOSTAT=IOstatus) chma, Ma
    read(10,*,IOSTAT=IOstatus)
    read(10,*,IOSTAT=IOstatus) n, m, Caux, Saux
    C(n+1,m+1)=Caux
    S(n+1,m+1)=Saux
    do while  (n .le. nMax)
    read(10,*,IOSTAT=IOstatus) n, m, Caux, Saux
    if (IOstatus .ge. 0) then
    C(n+1,m+1)=Caux
    S(n+1,m+1)=Saux
    else
    go to 55
    end if
    end do
    55 	continue
    close(10)

end subroutine read_C_and_S_and_rr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_C_and_S(vv,ff, sigma ,a,Ma,nMax,C,S)
    real(kind=16), allocatable, dimension (:,:), intent(in) :: vv
    integer, allocatable, dimension (:,:), intent(in) :: ff
    real(kind=16), intent(in) :: a, sigma, Ma
    integer, intent(in) :: nMax
    real(kind=16), allocatable, dimension (:,:), intent(out) :: C, S
    real(kind=16), allocatable, dimension(:) :: mix
    integer :: n, m, i, j, k, inda, i_f, nFaces, v1, v2, v3
    real(kind=16) :: C00, Cnm, Snm, detJ
        type(trin) :: xx, yy, zz
    type(trin),dimension(2) :: arrCS
    
    nFaces=size(ff,dim=1)
    allocate(C(nMax+1,nMax+1))
    allocate(S(nMax+1,nMax+1))
    C=0.0d0
    S=0.0d0
    n=0
    do while (n .le. nMax)
        if (.not.(nMax .eq. 0)) then
            call progress_bar(n,nMax)
        end if
    ! Create mixing vector (actually table)
        allocate(mix((n+1)*(n+2)/2))

        i=n
        do while (i .ge. 0)
            k=0
            do while (k .le. n-i)
                j = n - i - k
                inda =indexing(i,j,k,n)
                mix(inda)=fact(i)*fact(j)*fact(k)/fact(n+3)
                k = k+1
            end do
            i=i-1
        end do
        ! Initialize m to 0
        m=0
    ! Loop over m
        do while (m .le. n)
        Cnm=0.0d0
        Snm=0.0d0
    ! allocation of arrCS and the basic deg 1 trinomials xx,yy,zz

        call trinloc(xx,1)
        call trinloc(yy,1)
        call trinloc(zz,1)
        call trinloc(arrCS(1),n)
        call trinloc(arrCS(2),n)
        do i_f=1, nFaces

            v1=ff(i_f,1)
            v2=ff(i_f,2)
            v3=ff(i_f,3)

            xx%coeff(1)=vv(v1,1)
            xx%coeff(2)=vv(v2,1)
            xx%coeff(3)=vv(v3,1)

            yy%coeff(1)=vv(v1,2)
            yy%coeff(2)=vv(v2,2)
            yy%coeff(3)=vv(v3,2)

            zz%coeff(1)=vv(v1,3)
            zz%coeff(2)=vv(v2,3)
            zz%coeff(3)=vv(v3,3)
            call deter_sub(xx%coeff,yy%coeff,zz%coeff,detJ)
            arrCS=compSmCandS(n, m,Ma, a,xx,yy,zz)
            Cnm=Cnm+detJ*dot_product(mix,arrCS(1)%coeff)
            Snm=Snm+detJ*dot_product(mix,arrCS(2)%coeff)
        end do
        deallocate(arrCS(1)%coeff,arrCS(2)%coeff)
        deallocate(xx%coeff,yy%coeff,zz%coeff)
        C(n+1,m+1)= C(n+1,m+1)+sigma*Cnm
        S(n+1,m+1)= S(n+1,m+1)+sigma*Snm

        m=m+1
        end do ! over m
    m=0
    n=n+1
    deallocate(mix)
    end do !over n
    if (nMax .ge. 1) then ! In order for this subroutine to work with compute_volume
        C00=C(1,1)
    else
        C00=1.d0
    end if

    do n=1,nMax+1
    do m=1,n
    C(n,m)=C(n,m)/C00
    S(n,m)=S(n,m)/C00

    ! This block is introduced to remove extremely small numbers from the coefficients
    if (abs(C(n,m)) .lt. 1.d-16) then
    C(n,m)=0.0d0
    end if
    if (abs(S(n,m)) .lt. 1.d-16) then
    S(n,m)=0.0d0
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end do
    end do
    write(*,*) 
    return
end subroutine compute_C_and_S
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive function compSmCandS(n, m,Ma, a, x, y,z) result (arrCS)
    implicit none
    real(kind=16) a, Ma,cost_aux,cost_aux2
    integer,intent(in) :: n,m
    class(trin), intent(in)  :: x, y, z 
    type(trin) :: rSquared,c,s, c_aux, s_aux, c_aux2, s_aux2
    type(trin),  dimension(2) ::auxArrCS,auxArrCS2	!allocatable,

    type(trin),dimension(2) :: arrCS
    
    !	write(*,*) "Allocating rSquared"	
    call trinloc(rSquared, 2)!*(x%deg))
    rSquared = TriNormSq(x,y,z)
    !	write(*,*) "Allocating small c and s"	
    call trinloc(c, n)
    call trinloc(s, n)
    call trinloc(arrCS(1),n)
    call trinloc(arrCS(2),n)

    ! Sectorial
    if (n.eq.m) then 
        if (n.eq.0) then
    !		write(*,*) "Enter n=m=0"
            c%coeff(1) = 1/Ma
            s%coeff(1) = 0.0
        else if (n.eq.1) then
            c = TriScal(x,1/(sqrt(3.0)*Ma*a))
            s = TriScal(y,1/(sqrt(3.0)*Ma*a))
        else if (n.gt.1) then
        call trinloc(c_aux,n-1)
        call trinloc(s_aux,n-1)
            auxArrCS = compSmCandS(n-1, m-1, Ma, a, x, y, z)
            c_aux=auxArrCS(1)
            s_aux=auxArrCS(2)
    cost_aux=((2*real(n) - 1) / sqrt(2*real(n)*(2*real(n) + 1))) / a
    c = TriScal(TriDiff(TriMult(x,c_aux),TriMult(y,s_aux)),cost_aux) 
    s = TriScal(TriSum(TriMult(y,c_aux),TriMult(x,s_aux)),cost_aux)

        end if
    ! Subdiagonal	
    else if (m.eq.n-1) then
    call trinloc(c_aux,n-1)
    call trinloc(s_aux,n-1)
        auxArrCS = compSmCandS(n-1, n-1, Ma, a, x, y, z)
        c_aux=auxArrCS(1)
        s_aux=auxArrCS(2)
        cost_aux=((2*real(n) - 1) / (sqrt(2*real(n) + 1)))/a
            c = TriScal(TriMult(z, c_aux),cost_aux)
            s = TriScal(TriMult(z, s_aux),cost_aux)
    else
    call trinloc(c_aux,n-1)
    call trinloc(s_aux,n-1)
    call trinloc(c_aux2,n-2)
    call trinloc(s_aux2,n-2)
        auxArrCS = compSmCandS(n-1, m, Ma, a, x, y, z)
        auxArrCS2 = compSmCandS(n-2, m, Ma, a, x, y, z)
        c_aux=auxArrCS(1)
        s_aux=auxArrCS(2)
        c_aux2=auxArrCS2(1)
        s_aux2=auxArrCS2(2)

    cost_aux=(2*real(n)-1)*sqrt((2*real(n) - 1)/&
        &((2*real(n) + 1) * (real(n) + real(m)) * (real(n) -real(m))))/a

    cost_aux2=sqrt(((2*real(n) - 3)*(real(n) + real(m) - 1)*&
        &(real(n) - real(m) - 1))/&
    &((2*real(n)+1) * (real(n) + real(m))*(real(n) - real(m))))/a**2

    c = TriDiff(TriScal(TriMult(z, c_aux),cost_aux),&
    &TriScal(TriMult(rSquared, c_aux2),cost_aux2))
    s = TriDiff(TriScal(TriMult(z, s_aux),cost_aux),&
    &TriScal(TriMult(rSquared, s_aux2),cost_aux2))
    end if

    arrCS(1)=c
    arrCS(2)=s	

end function compSmCandS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_volume(vv,ff,rr,vol) 
    ! Computation of the volume by means of the SH coefficients
    ! Actually one does not need the density
    real(kind=16), allocatable, dimension(:,:),intent(in) :: vv
    integer, allocatable, dimension(:,:),intent(in) :: ff
    real(kind=16), intent(in) :: rr
    real(kind=16), allocatable, dimension (:,:) :: C, S
    real(kind=16), intent(out) :: vol
    real(kind=16) :: ma=1.0d0, sig=1.0d0 ! Ma
    integer ::nmax=0
    ! type(trin),dimension(2) :: arrCS
    
    call compute_C_and_S(vv,ff, sig ,rr,ma,nmax,C,S)
    vol=C(1,1)
    return
end subroutine compute_volume

end module SHcomputer