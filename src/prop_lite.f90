!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	Program to propagate various type of dynamics                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module prop_lite
    use constants
    use SHcomputer
    use dynamics
    use geometry
    use trinomials
    use polygrav
    use keplerian
    use integration
    use SHcomputer

contains
subroutine prop_sd_ham(orb_in, totyears,stepdays)
    real(kind=16), dimension(7), intent(in) :: orb_in
    real(kind=16), intent(in) :: totyears, stepdays
    character(len=100) :: orb_fmt,orb_fmt2
    real(kind=16), dimension(4) :: x_ham ! secular state vector
    real(kind=16) :: tsid, mjd   ! Time variables
    real(kind=16) :: h  !sca= scale, h=stepsize, eps=small number
    real(kind=16) ::  a, a_ham, e, inc, om, bom, amean, fan, lam, Mstar
    ! real(kind=16) start_time, end_time	! start and end times
    real(kind=16), dimension(7) :: orb
    integer :: jsrp0, ksrp0
    real(kind=16), dimension(10,4) :: xx, ff
    ! real(kind=16), dimension(1) :: conv

    integer ::  totsteps, nmax, l   ! integers for iterations, counting the verts 
    a=orb_in(1)
    e=orb_in(2)
    inc=orb_in(3)
    om=orb_in(4)
    bom=orb_in(5)
    amean=orb_in(6)
    fan=orb_in(7)

    if (dyn_vec(2)) then
        write(*,*) "Testing SRP toy model"
        write(*,*) "j = ", jsrp
        write(*,*) "k = ", ksrp
        call init_srp_con()
        a2m=(are2ms*1d-6)/(ageo**2)
        ttime=(365.242196d0*time0)/(36525.6363d0*366.242196d0*pi2)
        Mstar=(357.5256+35999.049*ttime)*degree
        write(*,*) "resang = ", resang0
        om= resang0 - (jsrp*bom + ksrp*(Mstar/degree+omstardeg))
        om=mod(om,360.d0)
        if (om .lt. 0) then
            om=om+360.d0
        end if
    end if
    jsrp0=jsrp
    ksrp0=ksrp
    inc=inc*degree
    om=om*degree
    bom=bom*degree
    amean=amean*degree


    ! and faces of poly, steps of integration and maximum degree for SH, N and l ADD DESCRIPTION 
    open(1,file='test_out/ham_states_sd.plt',  status='unknown',  form='formatted')
    open(2,file='test_out/ham_orb_states_sd.plt',  status='unknown',  form='formatted')

    a_ham=a*1.d-3/ageo ! 
    ! CAREFUL!!! Here we are dealing with scaled variables, mu=1
    ! and the unit of time is different
    rL=sqrt(a_ham) ! mu=1 with scaled variables
    x_ham(1)=rL*sqrt(1-e**2)  ! it's the rG Delauney element
    x_ham(2)=rL*sqrt(1-e**2)*cos(inc) ! if inc is given in radians
    x_ham(3)=om
    x_ham(4)=bom


    ! Hamiltonian states (1) and Orbital elements (2) (only 4, not included: a and M)
    orb_fmt= '(F15.5,4(2X,F15.8))'
    orb_fmt2= '(F15.5,2X,F18.8,4(2X,F15.10))'
    tsid= time0!0.d0
    tdays=(tsid-time0)/pi2
    ! Hailtonian states: rG, rH, om, bom, with parameter rL
    write(1,orb_fmt) tdays, x_ham(1), x_ham(2), x_ham(3), x_ham(4)
    if (dyn_vec(3)) then
        write(2,orb_fmt2) tdays, a/1.d3, e, inc/degree, om/degree, resang0
    else
        write(2,orb_fmt2) tdays, a/1.d3, e, inc/degree, om/degree, bom/degree!a/1.d3, amean/degree, lam/degree
    end if
        ! write(*,*) 'Step of integration: ', h
    h=pi2*stepdays
    totsteps=int(totyears*366.2419/stepdays)
    ! write(*,*) totyears*366.2419d0,totsteps, h
    write(*,*) 'Total number of steps: ', totsteps
    dateju=mjd+2400000.5d0
    ! We set the initial time to 0 because there are no third bodies whose position is relevant
    tsid=time0
    ! Start integration
    if (flag_ABM) then
        write(*,*) "Initializing ABM method..."
        call ainit(sd_ham_dyn,x_ham,xx,ff,h,time0,prel_steps)
    end if
    write(*,*) 
    write(*,*) "Starting integration..."
    do i=1, totsteps
        call progress_bar(i, totsteps)
        ! call RK6(sd_ham_dyn,x_ham,tsid,h)

        if (.not. flag_ABM) then
            call RK6(sd_ham_dyn,x_ham,tsid,h)
        else
            call adams(sd_ham_dyn,xx,ff,tsid,h,iter,eps,conv)
            do j=1,4
                x_ham(j)=xx(8,j)
            end do
        end if

        ! Writing the output to file
        e=sqrt(1-(x_ham(1)/rL)**2)
        inc= acos(x_ham(2)/x_ham(1))
        om=x_ham(3)
        bom=x_ham(4)
        om=mod(om,pi2)
        bom=mod(bom,pi2)
        ! om and bom are considered modulo 2pi, and we chose them to belong
        ! to the interval [0,2pi)
        do while ((om .lt. 0.0d0) .or. (om .gt. pi2))
            if (om .lt. 0.d0) then
                om= om+pi2
            else if (om .gt. pi2) then
                om=om-pi2
            end if
        end do
        do while ((bom .lt. 0) .or. (bom .gt. pi2))
            if (bom .lt. 0.d0) then
                bom= bom+pi2
            else if (bom .gt. pi2) then
                bom=bom-pi2
            end if
        end do
        x_ham(3)=om
        x_ham(4)=bom
        tsid=tsid+h
        tdays=(tsid-time0)/pi2
        write(1,orb_fmt) tdays, x_ham(1), x_ham(2), x_ham(3)/degree, x_ham(4)/degree
        if (dyn_vec(2)) then
            ttime=(365.242196d0*tsid)/(36525.6363d0*366.242196d0*pi2)
            Mstar=(357.5256+35999.049*ttime)*degree
            lam = om+ jsrp0*bom + ksrp0*(Mstar+omstar)
            lam=mod(lam,pi2)
            do while ((lam .lt. 0) .or. (lam .gt. pi2))
                if (lam .lt. 0.d0) then
                    lam= lam+pi2
                else if (lam .gt. pi2) then
                    lam=lam-pi2
                end if
            end do
            write(2,orb_fmt2) tdays, a/1.d3, e, inc/degree, om/degree, lam/degree
        else
            write(2,orb_fmt2) tdays, a/1.d3, e, inc/degree, om/degree, bom/degree!a/1.d3, amean/degree, lam/degree
        end if
    end do
    close(1)
    close(2)
end subroutine prop_sd_ham

subroutine prop_sd_cart(orb_in, totdays, stepminutes,nmax)
    real(kind=16), dimension(7), intent(in) :: orb_in
    real(kind=16), intent(in) :: totdays, stepminutes
    integer, intent(in) :: nmax
    character(len=100) :: orb_fmt,orb_fmt2, sh_address
    real(kind=16), dimension(6) :: x    ! Cartesian state vector
    real(kind=16) ::  tsid   ! Time variables
    real(kind=16) :: h
    real(kind=16) :: a, e, inc, om, bom, amean, fan, lam
    ! real(kind=16) start_time, end_time	! start and end times
    real(kind=16), dimension(7) :: orb
    real(kind=16), dimension(10,6) :: xx, ff
    ! real(kind=16), dimension(1) :: conv
    ! real(kind=16) ttime0,Mstar0 !Effects for third body
    integer ::  totsteps, l   ! integers for iterations, counting the verts 
    a=orb_in(1)
    e=orb_in(2)
    inc=orb_in(3)
    om=orb_in(4)
    bom=orb_in(5)
    amean=orb_in(6)

    mu=G*maE


    if (dyn_vec(2)) then
        write(*,*) "Testing SRP toy model"
        write(*,*) "j = ", jsrp
        write(*,*) "k = ", ksrp
        ! ttime=(365.242196d0*time0)/(36525.6363d0*366.242196d0*pi2)
        ttime=time0/day/36525.d0
        Mstar=357.5256+35999.049*ttime
        write(*,*) "resang = ", resang0
        om= resang0 - (jsrp*bom + ksrp*(Mstar+omstardeg))
        om=mod(om,360.d0)
        if (om .lt. 0) then
            om=om+360.d0
        end if
    end if

    inc=inc*degree
    om=om*degree
    bom=bom*degree
    amean=amean*degree

    open(1,file='test_out/car_states_sd.plt',  status='unknown',  form='formatted')
    if (flag_full_sh) then ! Earth with full SH coefficients
        open(2,file='test_out/SH_orb_states_sd.plt',  status='unknown',  form='formatted')
    else 
        open(2,file='test_out/orb_states_sd.plt',  status='unknown',  form='formatted')
    end if


    !   Solve Kepler's equation					 
    !   eec is the eccentric anomaly and fan is the true anomaly	 !
    write(*,*) 'Solving Kepler Equation'
    call solve_kep_eq(e,amean, fan)
    write(*,*) 'Done.'
    ! Initial position and velocity of the Satellite (or Debris) 		 !
    ! (in cartesian coordinates)							 !
    write(*,*) 'Converting to Cartesian coordinates'
    call indeb(x, a, e, bom, om, inc, fan,mu) ! Careful, indeb outputs angles in INSERT
    write(*,*) 'Done.'
    ! write(*,*) sqrt(x(1)**2+x(2)**2+x(3)**2)
    ! lam= fan - rate*time0 + om + bom !stroboscopic mean node
    ! Cartesian states (1) and Orbital elements (2)
    orb_fmt= '(F15.6,6(2X,F15.8))'
    orb_fmt2= '(F15.6,2X,F15.10,5(2X,F15.10))'
    tsid=0.d0 !time0
    tdays=(tsid-time0)/day
    write(1,orb_fmt) tdays, x(1), x(2), x(3), x(4), x(5), x(6)  
    if (dyn_vec(2)) then
        write(2,orb_fmt2) tdays, a/1.d3, e, inc/degree, om/degree, bom/degree, resang0
        else
            write(2,orb_fmt2) tdays, a/1.d3, e, inc/degree, om/degree, bom/degree, amean/degree!, lam/degree
        end if
    h=stepminutes*60.0d0
    totsteps=int(totdays*day/h)
    write(*,*) 'Total number of steps: ', totsteps
    ! write(*,*) "Earth's GM = ", mu
    ! time0= 0.d0 ! We set the initial time to 0 because there are no third bodies whose position is relevant
    tsid=time0
    ! Load SH coefficients if Full SH dynamics is engaged
    if (flag_full_sh) then ! Earth with full SH coefficients
        sh_address= 'spher_harm/SH_Earth_N_denorm.txt'
        write(*,*) 'Reading SH coefficients from file, up to degree', nmax
        call read_C_and_S_and_rr(sh_address,nmax,C,S,rr,Ma)
        Ma=maE
        write(*,*) 'Done reading.'
        write(*,*) "Ref. Rad. = ", rr, "[m]"
        write(*,*) "Mass = ", Ma, maE
        write(*,*) "GM = ", G*maE
        write(*,*) "mu = ", mu        
        write(*,*) "rate = ", rate
    end if
    
    if (flag_ABM) then
        write(*,*) "Initializing ABM method..."
        call ainit(sd_car_dyn,x,xx,ff,h,time0,prel_steps)
    end if
    ! Start integration
    write(*,*) 
    write(*,*) "Starting integration..."
    do i=1, totsteps
        call progress_bar(i, totsteps)
        ! call RK6(sd_ham_dyn,x_ham,tsid,h)

        if (.not. flag_ABM) then
            call RK6(sd_car_dyn,x,tsid,h)
        else
            call adams(sd_car_dyn,xx,ff,tsid,h,iter,eps,conv)
            do j=1,6
                x(j)=xx(8,j)
            end do
        end if

    ! write(*,*) "Starting integration..."
    ! do i=1, totsteps
    !     call progress_bar(i, totsteps)
    !     call RK6(sd_car_dyn,x,tsid,h)
        ! Conversion from Cartesian Coordinates to orbital elements
        call orbital(x, orb,mu)
        tdays=(tsid-time0)/day
        if (mod(i,50).eq.0) then
            write(1,orb_fmt) tdays, x(1), x(2), x(3), x(4), x(5), x(6)
            if (dyn_vec(2)) then
                ttime=tsid/day/36525.d0
                Mstar=357.5256d0+35999.049d0*ttime
                lam = orb(4)+ jsrp*orb(5) + ksrp*(Mstar+omstardeg)
                lam=mod(lam,360.d0)
                if (lam .lt. 0) then
                    lam=lam+360.d0
                end if
                write(2,orb_fmt2) tdays, orb(1)/1.d3, orb(2), orb(3), orb(4), orb(5), lam
            else
                write(2,orb_fmt2) tdays, orb(1)/1.d3, orb(2), orb(3), orb(4), orb(5), orb(6)!, lam
            end if
        end if
        if (tdays .ge. totdays) then
            ! write(*,*) "Exceeding totdays..."
            go to 342
        end if
        tsid=tsid+h
    end do
        342     continue
    write(*,*)

    close(1)
    close(2)

end subroutine prop_sd_cart

subroutine prop_gen_poly(orb_in, tot_time, h, filename,period,sca)
    real(kind=16), dimension(7), intent(in) :: orb_in
    real(kind=16), intent(in) :: tot_time, h, period, sca
    character(len=40), intent(in) :: filename
    character(len=1) ::  hash
    character(len=100) :: orb_fmt,orb_fmt2
    character(len=100) :: address, sh_address, shifted_address
    ! real(kind=16), dimension(3) :: fp, com 	! point at which the poly acc or pot is computed, com = coordinate of com
    real(kind=16), dimension(6) :: x    ! Cartesian state vector
    real(kind=16), dimension(10,6) :: xx, ff
    real(kind=16) :: tsid, tout   ! Time variables
    real(kind=16) :: lap ! sca= scale, h=stepsize, eps=small number, mu=Gm product, lap=laplacian
    real(kind=16) ::  a, e, inc, om, bom, amean, fan, lam, mu
    real(kind=16) :: vol, sig!, rr		! period and volume of the polyhedron
    ! real(kind=16) start_time, end_time	! start and end times
    real(kind=16), dimension(7) :: orb
    real(kind=16), dimension(3) :: pos_aux
	! real(kind=16), dimension(1) :: conv
    integer :: totsteps, i, j
    integer(kind=4) nVert, nFaces


    a=orb_in(1)
    e=orb_in(2)
    inc=orb_in(3)
    om=orb_in(4)
    bom=orb_in(5)
    amean=orb_in(6)

    inc=inc*degree
    om=om*degree
    bom=bom*degree
    amean=amean*degree

    ! Open the output files
    open(1,file='test_out/poly_states_ast.plt',  status='unknown',  form='formatted')
    open(2,file='test_out/poly_orb_states_ast.plt',  status='unknown',  form='formatted') 
    address= "shape_models/"//trim(filename)//'.obj'
    shifted_address = address
    ! shifted_address = 'shift_poly_mod/'//trim(filename)//'_shift.obj'
    ! other polyhedron data
    rate=pi2/period 
    Ma=1.0d0        ! This value is set to one for convenience
    !      EXTRACTION OF VERTICES AND FACES
    ! open the poly .obj located in "address"
    open (3, file=address, form='formatted')
    read(3,*) hash  ! Read the number of vertices
    close(3)
    ! Note that we request a special form for the .obj, in the first line we want an hash (#)
    ! followed by the number of vertices first and faces after
    ! Future versions should include also different types of .obj files
    if (hash .ne. 'v') then
        open (3, file=address, form='formatted')
        read(3,*) hash, nVert
        read(3,*) hash, nFaces 
    else
        ! if we don't know how many vertices and faces we can count them using
        ! the subroutine count_v_and_f from polygrav
        call count_v_and_f(address,nVert,nFaces)
        open (3, file=address, form='formatted')
    end if
    close(3)
    !Allocate matrices of vertices coords and vertices composing a simplex !
    vol=0.0d0
    ! write(*,*) sca
    write(*,*) 'Extraction of vertices, faces, edges and normals in progress.'
    call extract_ver_and_fac(address,nVert,nFaces,sca,ver,fac,rr)   ! Also computes the reference radius!
    call extract_edges(fac,edg)
    call face_normals(ver,fac,ffnn)
    write(*,*) 'Done.'
    write(*,*) 'Computation of volume in progress.'
    call compute_volume(ver,fac,rr,vol)
    write(*,*) 'Done.'
    write(*,*) 'Vol = ', vol, ' [m**3], or', vol/sca**3, ' [km**3]'
    Ma=vol*sigma
    write(*,*) 'Mass = ', Ma, ' [kg]'
    write(*,*) 'RefRad = ', rr, ' [m]'
    mu=G*Ma
    !   Solve Kepler's equation	 
    !   eec is the eccentric anomaly and fan is the true anomaly	 !
    write(*,*) 'Solving Kepler Equation'
    call solve_kep_eq(e,amean, fan)
    write(*,*) 'Done.'
    ! Initial position and velocity of the Satellite (or Debris)!
    ! (in cartesian coordinates)
    write(*,*) 'Converting to Cartesian coordinates'
    call indeb(x, a, e, bom, om, inc, fan,mu) ! Careful, indeb outputs angles in INSERT
    write(*,*) 'Done.'
	! lam= fan - rate*time0 + om + bom !stroboscopic mean node
    ! Cartesian states (1) and Orbital elements (2)
    orb_fmt= '(F18.6,6(2X,F18.8))'
    orb_fmt2= '(F15.6,2X,F15.10,5(2X,F15.10))'
    tsid=0.d0 !time0
    tdays=(tsid-time0)/day
    write(1,orb_fmt) tdays, x(1), x(2), x(3), x(4), x(5), x(6)
    write(2,orb_fmt2) tdays, a/1.d3, e, inc/degree, om/degree, bom/degree, fan/degree !, lam/degree
    ! h=stepminutes*60.d0
    totsteps=int(tot_time/h)
    write(*,*) 'Total number of steps: ', totsteps
    write(*,*) "Asteroid's GM = ", G*Ma
    time0= 0.d0 ! We set the initial time to 0 because there are no third bodies whose position is relevant
    tsid=time0

    if (flag_ABM) then
        write(*,*) "Initializing ABM method..."
        call ainit(poly_dyn,x,xx,ff,h,time0,prel_steps)
    end if

    ! Start integration
    write(*,*) "Starting integration..."
    do i=1, totsteps
        call progress_bar(i, totsteps)
        ! Check if the object is inside the asteroid using the Laplacian
        ! pos_aux(1)=x(1)
        ! pos_aux(2)=x(2)
        ! pos_aux(3)=x(3)
        ! call polylap(sigma,pos_aux,ver,fac,lap)
        ! if (lap/(G*sigma) .lt. -pi2) then
        !     write(*,*) lap/(G*sigma)
        !     write(*,*) "INSIDE THE ASTEROID!"
        !     goto 342
        ! end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (.not. flag_ABM) then
            call RK6(poly_dyn,x,tsid,h)
        else
            call adams(poly_dyn,xx,ff,tsid,h,iter,eps,conv)
            do j=1,6
                x(j)=xx(8,j)
            end do
        end if
        ! Conversion from Cartesian Coordinates to orbital elements
        call orbital(x, orb,mu)
        tdays=(tsid-time0)/day ! Remove abs in the future
        ! write(*,*) tdays, orb(1)/1.d3, orb(2), orb(3), orb(4), orb(5), orb(6)
        ! stop
        if (mod(i,50) .eq. 0) then
            write(1,orb_fmt) tdays, x(1)/1.d3, x(2)/1.d3, x(3)/1.d3, x(4)/1.d3, x(5)/1.d3, x(6)/1.d3
            write(2,orb_fmt2) tdays, orb(1)/1.d3, orb(2), orb(3), orb(4), orb(5), orb(6)
        end if
        if (tdays*day .ge. tot_time) then
            go to 420
        end if
        tsid=tsid+h
    end do
     
    420     continue
    close(1)
    close(2)
end subroutine prop_gen_poly

subroutine prop_gen_SH(orb_in, tot_time, h, filename, period, nmax)
    real(kind=16), dimension(7), intent(in) :: orb_in
    real(kind=16), intent(in) :: tot_time, h, period
    character(len=40), intent(in) :: filename
    integer, intent(in) :: nmax
    character(len=100) :: sh_address
    character(len=1) ::  hash
    character(len=10) :: nmax_char
    character(len=100) :: orb_fmt,orb_fmt2
    real(kind=16), dimension(6) :: x    ! Cartesian state vector
    real(kind=16), dimension(10,6) :: xx, ff
    real(kind=16) :: tsid, tout   ! Time variables
    real(kind=16) :: vol! volume (refRad è shared_data)
    real(kind=16) :: nre, a, e, inc, om, bom, amean, fan, lam, mu, theta
    ! real(kind=16) start_time, end_time	! start and end times
    real(kind=16), dimension(7) :: orb
	! real(kind=16), dimension(1) :: conv
    real(kind=16) relerr, abserr
    ! real(kind=16) ttime0,Mstar0 !Effects for third body
    integer ::  totsteps, l   ! integers for iterations, counting the verts 

    a=orb_in(1)
    e=orb_in(2)
    inc=orb_in(3)
    om=orb_in(4)
    bom=orb_in(5)
    amean=orb_in(6)

    inc=inc*degree
    om=om*degree
    bom=bom*degree
    amean=amean*degree
    rate=pi2/period 

    ! Open the output files
    write(nmax_char,'(I3.3)') nmax_aux
    open(1,file='test_out/'//trim(nmax_char)//'_SH_states_ast.plt',  status='unknown',  form='formatted')
    open(2,file='test_out/'//trim(nmax_char)//'_SH_orb_states_ast.plt',  status='unknown',  form='formatted') 

    ! dyn_type="SH"
    ! SH Dynamics auxiliary data
    sh_address= 'spher_harm/SH_'//trim(filename)//'.txt'
    write(*,*) 'Reading SH coefficients from file, up to degree', nmax_aux
    call read_C_and_S_and_rr(sh_address,nmax_aux,C,S,rr,Ma)
    write(*,*) 'Done reading.'
    write(*,*) "Ref. Rad. = ", rr, "[m]"
    write(*,*) "Mass = ", Ma

    !   Solve Kepler's equation					 
    !   eec is the eccentric anomaly and fan is the true anomaly	 !
    write(*,*) 'Solving Kepler Equation'
    call solve_kep_eq(e,amean, fan)
    write(*,*) 'Done.'
    ! Initial position and velocity of the Satellite (or Debris) 		 !
    ! (in cartesian coordinates)							 !
    write(*,*) 'Converting to Cartesian coordinates'
    call indeb(x, a, e, bom, om, inc, fan,G*Ma) ! Careful, indeb outputs angles in INSERT
    write(*,*) 'Done.'

    ! Cartesian states (1) and Orbital elements (2)
    orb_fmt= '(F15.6,6(2X,F15.8))'
    orb_fmt2= '(F15.6,2X,F15.10,5(2X,F15.10))'
    tsid=time0
    ! time0=0.d0
    tdays=(tsid-time0)/day
    write(1,orb_fmt) tdays, x(1), x(2), x(3), x(4), x(5), x(6)
    write(2,orb_fmt2) tdays, a/1.d3, e, inc/degree, om/degree, bom/degree, amean/degree!, lam/degree
    ! h=stepminutes*60.d0
    totsteps=int(tot_time/h)
    write(*,*) 'Total number of steps: ', totsteps
    write(*,*) "Asteroid's GM = ", G*Ma

    if (flag_ABM) then
        write(*,*) "Initializing ABM method..."
        call ainit(SH_dyn,x,xx,ff,h,time0,prel_steps)
    end if

    ! Start integration
    write(*,*) "Starting integration..."
    do i=1, totsteps
        call progress_bar(i, totsteps)

        if (.not. flag_ABM) then
            call RK6(SH_dyn,x,tsid,h)
        else
            call adams(SH_dyn,xx,ff,tsid,h,iter,eps,conv)
            do j=1,6
                x(j)=xx(8,j)
            end do
        end if

        ! Conversion from Cartesian Coordinates to orbital elements
        call orbital(x, orb, G*Ma)
        tdays=(tsid-time0)/day 
        theta=rate*tsid/degree
        lam=(orb(6)-2.d0*theta+2.d0*orb(4)+2.d0*orb(5))/2.d0
        lam=mod(lam,360.d0)
        ! do while ((lam .lt. 0.0d0) .or. (lam .gt. 360.0d0))
        ! if (lam .lt. 0.d0) then
        ! lam= lam+360.0d0
        ! else if (lam .gt. 360.0d0) then
        ! lam=lam-360.0d0
        ! end if
        ! end do
        if (mod(i,50).eq.0) then
            write(1,orb_fmt) tdays, x(1), x(2), x(3), x(4), x(5), x(6)
            write(2,orb_fmt2) tdays, orb(1)/1.d3, orb(2), orb(3), orb(4), orb(5), orb(6)!, lam
        end if
        if (tdays*day .ge. tot_time) then
            go to 342
        end if
        tsid=tsid+h
    end do

    342     continue
    close(1)
    close(2)
end subroutine prop_gen_SH
end module

 