!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	Program to propagate various type of dynamics                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module var_lite
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

subroutine var_gen_SH(orb_in, tot_time, h, filename, period, nmax,delta,v0)
    real(kind=16), dimension(7), intent(in) :: orb_in
    real(kind=16), dimension(6), intent(in) :: v0
    real(kind=16), intent(in) :: tot_time, h, period,delta
    character(len=40), intent(in) :: filename
    integer, intent(in) :: nmax
    character(len=100) :: sh_address
    character(len=1) ::  hash
    character(len=10) :: nmax_char
    character(len=100) :: orb_fmt,orb_fmt2, fli_fmt
    real(kind=16), dimension(6) :: x,v  ! Cartesian state vector
    real(kind=16), dimension(10,6) :: xx, ff
    real(kind=16) :: tsid, tout   ! Time variables
    real(kind=16) :: vol! volume (refRad è shared_data)
    real(kind=16) :: nre, a, e, inc, om, bom, amean, fan, lam, theta
    real(kind=16) :: e0, inc0, bom0, amean0, fan0, ha, hsig
    integer :: iX,iY
    ! real(kind=16) start_time, end_time	! start and end times
    real(kind=16), dimension(7) :: orb
	! real(kind=16), dimension(1) :: conv
    real(kind=16) suppnorm
    ! real(kind=16) ttime0,Mstar0 !Effects for third body
    integer ::  totsteps, l   ! integers for iterations, counting the verts 

    ! Open the output files
    ! write(nmax_char,'(I3.3)') nmax_aux
    open(2,file='test_out/fli/fli_sigma_a.plt',  status='unknown',  form='formatted')
    fli_fmt= '(F15.8,2(2X,F15.8))'
    ! SH Dynamics auxiliary data
    totsteps=int(tot_time/h) ! Remember that tot_time and time_steps are given in seconds
    
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Save the initial conditions, so to re use them at!
    ! every iteration in the grid                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    e0=e
    bom0=bom
    inc0=inc
    om0=om
    amean0=amean
    
    !   Solve Kepler's equation					 
    !   eec is the eccentric anomaly and fan is the true anomaly	 !
    write(*,*) 'Solving Kepler Equation'
    call solve_kep_eq(e0,amean0, fan0)
    write(*,*) 'Done.'
    v=v0
    ! Initialize tangent vector
    pnorm=0.d0
    do ii=1,6
        pnorm=pnorm+v(ii)**2
    end do    
    
    pnorm=sqrt(pnorm)
    do ii=1,6
        v(ii)=v(ii)/pnorm
    end do    
    ! v_0=v
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Preparation                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (reso_type .eq. '11') then
        write(*,*) 'Testing Resonance 1:1'
        a=(mu/(rate)**2)**(1.d0/3.d0)
        write(*,*) 'Semi-major axis for 1:1 resonance of ', filename
        write(*,*) 'Located at: [km]', a/1.d3
    else  if (reso_type .eq. '21') then
        write(*,*) 'Testing Resonance 2:1'
        a=(mu/(2.d0*rate)**2)**(1.d0/3.d0)
        write(*,*) 'Semi-major axis for 2:1 resonance of ', filename
        write(*,*) 'Located at: [km]', a/1.d3
    end if

    write(*,*) "Other initial conditions:"
    write(*,*) "eccentricity=", e0
    write(*,*) "inclination=", inc0/degree
    write(*,*) "om=", om0/degree
    write(*,*) "Omega=", bom0/degree
    !    call compute_amplitude_mmr(a,delta)
    amin= a-1.d0*delta/2
    amax=a+1.d0*delta/2
    ha= (amax-amin)/grid
    hsig= 360.d0/grid
    write(*,*) 'Amplitude =' , delta/1.d3, " km"
    write(*,*) 'amin=', amin, ' amax=', amax
    write(*,*) 'grid=',grid,', hsig=', hsig, ', ha=', ha
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Computation of the FLIs                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call cpu_time(start_time)
    do iY= 1, grid+1 ! varies a
    !   call progress_bar(iY,grid+1)
      do iX= 1, grid+1 ! varies lam
        ! call progress_bar(iX,grid+1)
        tsid=time0
        a= amin+(iY-1)*ha
        lam = (iX-1)*hsig*degree
        if (reso_type .eq. '11') then
            amean= lam - bom0 - om0 + rate*time0
            ! amean=mod(amean,pi2)
            ! amean0=amean
        else if (reso_type .eq. '21') then
            amean= 2*lam - om0 -2.d0*bom0 + 2.d0*rate*time0
            amean=mod(amean,pi2)
            amean0=amean
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Initial position and velocity of the Satellite    !
        !(or Debris) in cartesian coordinates              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   write(*,*) 'Converting to Cartesian coordinates'
        call indeb(x, a, e0, bom0, om, inc0, amean,mu)
        !   write(*,*) 'Done.'
        v=v0
        suppnorm= -10.d0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Propagation starts here
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !tsid=time0
        tsid=0.d0
        do i=1, totsteps
          ! call progress_bar(i,totsteps)
          taux=tsid    
          aux_car_state=x
          call RK6(SH_dyn,x,tsid,h)        
          tsid=taux
          call RK6(var_SH,v,tsid,h)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Computation of the FLIs                          !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          pnorm=0.d0
          do ii=1,6
            pnorm=pnorm+v(ii)**2
          end do
          pnorm=0.5d0*log10(pnorm)
          suppnorm = max(abs(pnorm),suppnorm)
          if (i .lt. 50) then
            suppnorm=0.d0
          end if
          if ( suppnorm .gt. maxFLI) then
            ! write(*,*) 'FLI TOO BIG'
            ! stop
            goto  200
          end if
          if ( i .gt. max_steps) then
           ! write(*,*) 'Iteration over 100000'
            ! stop
            goto  200
          end if
          tsid=tout
        end do  ! over i=1,N
  
        200     continue
    
        if (suppnorm .gt. maxFLI) then
          suppnorm=maxFLI
        end if

        write(2,fli_fmt) lam/degree, a/1.d3, suppnorm
        suppnorm=-10.d0
      end do !iX
    !   write(2,fli_fmt) 
    end do !iY
    call cpu_time(end_time)
    write(*,*)
    write(*,*) 'Computation of FLIs done in ', dble(end_time-start_time)/60.d0, ' [min]'
    close(2)
end subroutine var_gen_SH
end module

 