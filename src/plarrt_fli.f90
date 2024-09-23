!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Program to propagate around a given          !
!   (possibly rotating) polyhedron                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program plarrt_fli
  use constants
  ! use shared_data
  ! use geometry
  use trinomials
  use keplerian
  ! use integration
  ! use dynamics
  ! use spherigrav
  ! use srp_constants
  use SHcomputer
  ! use prop_lite
  use var_lite

  ! Declaration of characters
  character(len=40) :: filename, orb_fmt
  real(kind=16) :: totdays, stepminutes 
  real(kind=16) :: totyears, stepdays, period
  real(kind=16) :: tot_time, time_step, h,tdays !h=stepsize
  real(kind=16), dimension(7) :: orb
  real(kind=16) :: nre, a,a_ham, e, inc, om, bom, amean
  character(len=5), dimension(2) :: type_time
  real(kind=16) :: aux_time
  character(len=5) :: aux_char
  character(len=100) :: fli_fmt, hash
  character(len=100) :: address, output_address, shifted_address, sh_address
  ! Declaration of reals (including vectors)
  real(kind=16), dimension(6) :: x,v,x_0, v_0
  real(kind=16), dimension(3) :: fp
  real(kind=16), dimension(10,6) :: xx,ff
  real(kind=16) ::  taux, tsid, tout, vol
  real(kind=16) :: sca, lap,  eps_var
  real(kind=16) :: ares, delta
  real(kind=16) :: e0, inc0, bom0, amean0, fan0
  real(kind=16) relerr, abserr, relerr_var, abserr_var, ha, hsig
  real(kind=16) start_time, end_time
  real(kind=16), allocatable, dimension(:) :: work, work_var
  ! Declaration of integers
  integer(kind=4), dimension(5) :: iwork, iwork_var
  integer(kind=4) iflag, iflag_var
  integer ::   totsteps, i, ii
  integer :: fu, rc, N, l, counter,iX,iY, tot_count
  integer :: nVert, nFaces, nmax!, max_steps
  ! Declaration of logicals
  logical :: flag_e, flag_auto, shift_flag

  namelist /OBJECT/ filename, sca, sigma, period   ! read filename, part of the address as "poly_mod/filename.obj", the period of the rotation, the density, the scale and the maximum degree of SH
  !nre
  namelist /SH/ nmax

  namelist /ORBEL/ a, e, inc, om, bom, amean  ! read the initial orbital elements for a propagation

  namelist /TIME/ time0, tot_time, time_step, type_time

  namelist /DYN/ dyn_vec, flag_full_sh, flag_alt_j2, flag_ABM
  
  namelist /FLI/  grid, delta, reso_type, max_steps,maxFLI,dyn_type,v_0

  ! namelist /FLAGS/ flag_res, shift_flag

  open (action='read', file='namelist.nml', iostat=rc, newunit=fu)
  read (nml=OBJECT, iostat=rc, unit=fu)
  read (nml=SH, iostat=rc, unit=fu)
  read (nml=ORBEL, iostat=rc, unit=fu)
  read (nml=TIME, iostat=rc, unit=fu)
  read (nml=DYN, iostat=rc, unit=fu)
  read (nml=FLI, iostat=rc, unit=fu)
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
  
  orb(1)=a
  orb(2)=e
  orb(3)=inc
  orb(4)=om
  orb(5)=bom
  orb(6)=amean
  orb(7)=fan

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      TIME                                        !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) 'Total time of integration:'
  write(*,*) tot_time, type_time(1)
  write(*,*) time_step, type_time(2)
  aux_char=type_time(1)
  call time_converter(tot_time,aux_char,aux_time)
  tot_time=aux_time
  aux_char=type_time(2)
  call time_converter(time_step,aux_char,aux_time)
  time_step=aux_time
  write(*,*) 'Total time of integration per point on the grid:'
  write(*,*) tot_time, "s"
  write(*,*) 'Time step:'
  write(*,*) time_step, "s"
  time0=0.0d0
  totsteps=int(tot_time/time_step)
  write(*,*) 'Total number of steps per dots on the grid: ', totsteps

  if (dyn_type .eq. 'poly') then
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !      EXTRACTION OF VERTICES AND FACES                  !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! if (shift_flag) then
    !   continue
    !   else 
    !   address=shifted_address
    ! end if
    ! open (3, file=address, form='formatted')
    ! read(3,*) hash  ! Read the number of vertices
    ! close(3)
    ! if (hash .ne. 'v') then
    !   open (3, file=address, form='formatted')
    !   read(3,*) bull, nVert
    !   read(3,*) bull, nFaces 
    ! else
    !   call count_v_and_f(address,nVert,nFaces)
    !   ! write(*,*) nVert, nFaces
    !   open (3, file=address, form='formatted')
    ! end if
    ! ! Read the number of faces
    ! close(3)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !   Allocate matrices of vertices coords           !
    ! !   and vertices composing a simplex               ! 
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write(*,*) 'Extraction of vertices, faces, edges and normals in progress.'
    ! call extract_ver_and_fac(address,nVert,nFaces,sca,ver,fac,rr)
    ! call extract_edges(fac,edg)
    ! call face_normals(ver,fac,ffnn)
    ! write(*,*) 'Done.'
    ! write(*,*) 'Computation of volume in progress.'
    ! call compute_volume(ver,fac,rr,vol)
    ! write(*,*) 'Done.'
    ! write(*,*) 'Vol = ', vol, ' [m**3]'
    ! write(*,*) 'Mass = ', vol*sigma, ' [kg]'
    ! Ma=vol*sigma
    ! mu= G*vol*sigma
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else if (dyn_type .eq. 'SH') then
    sh_address= 'spher_harm/SH_'//trim(filename)//'.txt'
    write(*,*) 'Reading SH coefficients from file, up to degree', nmax
    call read_C_and_S_and_rr(sh_address,nmax,C,S,rr,Ma)
    write(*,*) 'Done reading.'
    write(*,*) "Ref. Rad. = ", rr/1.d3,' [km]'
    write(*,*) "Mass = ", Ma
  end if
  mu=G*Ma
  
  if (dyn_type .eq. 'SH') then
    call var_gen_SH(orb, tot_time,time_step,filename,period,nmax,delta,v_0)
  end if


end
    