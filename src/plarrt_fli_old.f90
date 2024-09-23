!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Program to propagate around a given          !
!   (possibly rotating) polyhedron                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program plarrt_fli
  use constants
  use shared_data
  use geometry
  use trinomials
  use keplerian
  use integration
  use dynamics
  use spherigrav
  use srp_constants
  use SHcomputer

  ! Declaration of characters
  character(len=30) :: filename
  character(len=1) :: bull
  character(len=100) :: fli_fmt
  character(len=100) :: address, output_address, shifted_address, orb_fmt, sh_address
  ! Declaration of reals (including vectors)
  real(kind=16), dimension(6) :: x,v,x_0, v_0
  real(kind=16), dimension(3) :: fp
  real(kind=16), dimension(10,6) :: xx,ff
  real(kind=16) ::  time0, taux, tsid, tout, vol
  real(kind=16) :: totdays, stepminutes, totdays_fli
  real(kind=16) :: sca, h, eps, lap, maxFLI, eps_var
  real(kind=16) :: ares, delta
  real(kind=16), dimension(7) :: orb
  real(kind=16) :: a, e, inc, om, bom, lam, fan, nre
  real(kind=16) :: e0, inc0, bom0, amean0, fan0
  real(kind=16), dimension(1) :: conv
  real(kind=16) relerr, abserr, relerr_var, abserr_var, ha, hsig
  real(kind=16) start_time, end_time
  real(kind=16), allocatable, dimension(:) :: work, work_var
  ! Declaration of integers
  integer(kind=4), dimension(5) :: iwork, iwork_var
  integer(kind=4) iflag, iflag_var
  integer :: iter, prel_steps, totsteps, grid, i, ii
  integer :: fu, rc, N, l, counter,iX,iY, tot_count
  integer :: nVert, nFaces, nmax, max_steps
  ! Declaration of logicals
  logical :: flag_e, flag_auto, shift_flag, flag_res11, flag_res21

  namelist /INIT/ filename, period,sigma, sca, nmax

  namelist /ORBEL/ nre, a, e, inc, om, bom, amean

  namelist /INTEG/ dyn_type, totdays, stepminutes, time0, eps,v    
  
  namelist /FLI/  grid, maxFLI, delta, flag_res11, flag_res21,totdays_fli,eps_var, max_steps

  namelist /FLAGS/ flag_res, shift_flag

  open (action='read', file='namelist.nml', iostat=rc, newunit=fu)
  read (nml=INIT, iostat=rc, unit=fu)
  read (nml=ORBEL, iostat=rc, unit=fu)
  read (nml=INTEG, iostat=rc, unit=fu)
  read (nml=FLI, iostat=rc, unit=fu)
  read (nml=FLAGS, iostat=rc, unit=fu)
  close(fu)

  allocate(work(100+21*neqn))
  allocate(work_var(100+21*neqn))
  sigma=sigma*1000 !Originally sigma is provided in g/cm**3
  rate=pi2/period
  Ma=1.0d0
  address= "poly_mod/" // trim(filename)//'.obj'
  output_address= 'SH_text/SH_'//trim(filename)//'.txt'
  shifted_address = 'shift_poly_mod/'//trim(filename)//'_shift.obj'
    
  !    if (flag_e) then
  !    open(2,file='test_out/fli/fli_sigma_e.plt',  status='unknown',  form='formatted')
  !    else
  open(2,file='test_out/fli/fli_sigma_a.plt',  status='unknown',  form='formatted')
  !    end if
  fli_fmt= '(F15.8,2(2X,F15.8))'
  if (dyn_type .eq. 'poly') then
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      EXTRACTION OF VERTICES AND FACES                  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (shift_flag) then
      continue
      else 
      address=shifted_address
    end if
    open (3, file=address, form='formatted')
    read(3,*) bull  ! Read the number of vertices
    close(3)
    if (bull .ne. 'v') then
      open (3, file=address, form='formatted')
      read(3,*) bull, nVert
      read(3,*) bull, nFaces 
    else
      call count_v_and_f(address,nVert,nFaces)
      ! write(*,*) nVert, nFaces
      open (3, file=address, form='formatted')
    end if
    ! Read the number of faces
    close(3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Allocate matrices of vertices coords           !
    !   and vertices composing a simplex               ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'Extraction of vertices, faces, edges and normals in progress.'
    call extract_ver_and_fac(address,nVert,nFaces,sca,ver,fac,rr)
    call extract_edges(fac,edg)
    call face_normals(ver,fac,ffnn)
    write(*,*) 'Done.'
    write(*,*) 'Computation of volume in progress.'
    call compute_volume(ver,fac,rr,vol)
    write(*,*) 'Done.'
    write(*,*) 'Vol = ', vol, ' [m**3]'
    write(*,*) 'Mass = ', vol*sigma, ' [kg]'
    Ma=vol*sigma
    mu= G*vol*sigma
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else if (dyn_type .eq. 'SH') then
    sh_address= 'spher_harm/SH_'//trim(filename)//'.txt'
    write(*,*) 'Reading SH coefficients from file, up to degree', nmax
    call read_C_and_S_and_rr(sh_address,nmax,C,S,rr,Ma)
    write(*,*) 'Done reading.'
    write(*,*) "Ref. Rad. = ", rr
    write(*,*) "Mass = ", Ma
  end if
  mu=G*Ma
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      TIME                                        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  N=int(totdays_fli)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The step of integration. "stepminutes" represents!
  ! the number of minutes:                           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  h=stepminutes*60.0d0
  !    write(*,*) 'Stepminutes: ', stepminutes        
  write(*,*) 'Total time of integration per point on the grid: ', totdays_fli
  write(*,*) 'Step of integration: ', h
  time0=0.0d0
  totsteps=int(totdays*day/h)
  write(*,*) 'Total number of steps per dots on the grid: ', totsteps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Initial position of debris                  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  bom=bom*degree
  inc=inc*degree
  om=om*degree
  amean=amean*degree
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Save the initial conditions, so to re use them at!
  ! every iteration in the grid                      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  e0=e
  bom0=bom
  inc0=inc
  om0=om

  delta=delta*sca 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Solve the Kepler's equation    eec is the       !
  !  eccentric anomaly and fan is the true anomaly   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call solve_kep_eq(e0,amean0, fan0)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Tolerances                                  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  relerr=eps !eps*1.d5
  abserr=eps*1d-1
  
  relerr_var=eps_var
  abserr_var=eps_var*1d-1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Normalization of the tangent vector         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pnorm=0.d0
  do ii=1,6
      pnorm=pnorm+v(ii)**2
  end do    
  
  pnorm=sqrt(pnorm)
  do ii=1,6
      v(ii)=v(ii)/pnorm
  end do    

  v_0=v
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Preparation                                      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if (flag_res11) then
    write(*,*) 'Testing Resonance 1:1'
    a=(mu/(rate)**2)**(1.d0/3.d0)
    write(*,*) 'Semi-major axis for 1:1 resonance of ', filename
    write(*,*) 'Located at: [km]', a/sca
  else if (flag_res21) then
    write(*,*) 'Testing Resonance 2:1'
    a=(mu/(2.d0*rate)**2)**(1.d0/3.d0)
    write(*,*) 'Semi-major axis for 1:1 resonance of ', filename
    write(*,*) 'Located at: [km]', a/sca
  end if

  write(*,*) "Other initial conditions:"
  write(*,*) "inclination=", inc/degree
  write(*,*) "om=", om0/degree
  write(*,*) "Omega=", bom0/degree
  !    call compute_amplitude_mmr(a,delta)
  amin= a-1.d0*delta/2
  amax=a+1.d0*delta/2
  ha= (amax-amin)/grid
  hsig= 180.d0/grid
  write(*,*) 'Amplitude =' , delta/sca, " km"

  !    orb_fmt= '(F15.8,6(2X,F15.8))'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Computation of the FLIs                          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call cpu_time(start_time)
  do iY= 1, grid+1 ! varies a
    call progress_bar(iY,grid+1)
    do iX= 1, grid+1 ! varies lam
      ! call progress_bar(iX,grid+1)
      iflag=1
      iflag_var=1
      tsid=time0
      a= amin+(iY-1)*ha
      lam = (iX-1)*hsig*degree
      amean= lam - bom0 - om0 + rate*time0
      amean=mod(amean,pi2)
      amean0=amean
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Initial position and velocity of the Satellite    !
      !(or Debris) in cartesian coordinates              !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !   write(*,*) 'Converting to Cartesian coordinates'
      call indeb(x, a, e0, bom0, om, inc0, amean0,mu)
      !   write(*,*) 'Done.'
      v=v_0
      suppnorm= -10.d0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Propagation starts here
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tsid=time0
      do i=1, totsteps
        ! call progress_bar(i,totsteps)
        taux=tsid    
        tout=tsid+h
        aux_car_state=x
        call ode(vel_and_acc, neqn, x, tsid, tout, relerr, abserr, iflag, work, iwork)        
        ! if (lap_check .lt. -pi2) then
        !   write(*,*) 'Collision!'
        !   go to 200
        ! end if
        tsid=taux
        if (dyn_type .eq. 'poly') then
          call ode(var_poly_old, neqn, v, tsid, tout, relerr_var, abserr_var, iflag_var, work_var, iwork_var)
        else if (dyn_type .eq. 'SH') then 
          call ode(var_SH, neqn, v, tsid, tout, relerr_var, abserr_var, iflag_var, work_var, iwork_var)
        end if
        if ((iflag .ne. 2) .or. (iflag_var .ne. 2)) then
          write(*,*) 
          write(*,*) iflag
          write(*,*) iflag_var
          write(*,*) aux_car_state
          write(*,*) 
          write(*,*) i, iX, iY
          suppnorm=maxFLI
          goto  200
        end if
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
      write(2,fli_fmt) lam/degree, a/sca, suppnorm
      suppnorm=-10.d0
    end do !iX
    write(2,fli_fmt) 
  end do !iY
  call cpu_time(end_time)
  write(*,*)
  write(*,*) 'Computation of FLIs done in ', dble(end_time-start_time)/60.d0, ' [min]'
  close(2)
end
    