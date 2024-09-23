!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Module containing various public constants	   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module constants
    ! generic useful constants
    real(kind=16) :: pi2=8.d0*datan(1.d0)
    real(kind=16) :: degree=8.d0*datan(1.d0)/360.0d0
    real(kind=16) :: G=6.674d0*10.0**(-11) !m3⋅kg−1⋅s−2 gravitational constant
    real(kind=16) :: day=86400.d0 ! length of day in seconds
    real(kind=16) :: t0=0.0d0 !initial time in seconds
    real(kind=16) :: eJ2=1.08262602d-3      ! Earth_J2
    real(kind=16) :: ageo=42164.1696d0      ! semi-major axis of geostationary orbits in km
    real(kind=16) :: rE=6378.1370d0         ! Earth radius in km
    real(kind=16) :: istardeg=23.4392795  
    real(kind=16) :: omstardeg=282.94d0     ! argument of perigee of the Sun as viewed from Earth
    real(kind=16) :: crpr=4.56d-3           ![kg/(km*s^2)] Constant for SRP propagation
    real(kind=16) :: crpr_car=4.56d-6           ![kg/(m*s^2)] Constant for SRP propagation
    real(kind=16) :: muE=1.d0               ! value of mu_earth with modified units
    real(kind=16) :: maE=5.972d24, maS=1.989d30, maM=7.34767309d22
    real(kind=16) :: astrounit=149597870700.d0  !astronomical unit in meters
    integer :: neqn=6
    integer :: neqnsec=4
    integer :: neqpen=2
    integer :: prel_steps=10000
    integer :: nmax_aux,max_steps
    character(len=30) :: reso_type

    ! Declaration of shared variables
    real(kind=16) :: mu
    real(kind=16) :: ut                     ! value of the unit of time in seconds
    real(kind=16) :: rL, n_s, maxFLI,lon
    real(kind=16), dimension(4) :: aux_state 
    real(kind=16), dimension(6) :: aux_car_state
    real(kind=16) :: are2ms, a2m !to be initialized in main program with (are2ms*1d-6)/(ageo**2) where are2ms is given in m**2/kg are2ms should be read from the namelist
    integer :: jsrp, ksrp, jsrp0, ksrp0 ! SRP toy model labels
    real(kind=16) :: omstar
    real(kind=16) :: istar, time0, resang0
    real(kind=16), dimension(2) :: aux_state_pend
    real(kind=16), allocatable, dimension(:,:) :: ver, ffnn
    integer, allocatable, dimension(:,:) :: fac
    integer, allocatable, dimension(:,:) :: edg
    real(kind=16) :: rate, sigma, rr, Ma, lap_check
    character(len=30) :: dyn_type
    real(kind=16), allocatable, dimension(:,:) :: C, S
    real(kind=16), dimension(6,6) :: var_mat
    logical, dimension(4) :: dyn_vec ! (J2,SRP,Sun,Moon)
    logical :: flag_full_sh, flag_all_tm, flag_alt_j2, flag_ABM
    integer :: iter=5000
    integer :: grid
    real(kind=16) :: eps =1.d-9
    real(kind=16),dimension(1) :: conv
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Routine which prints a progressbar on screen	   !
    !  Input: iter - current iteration number          !
    !         max_iter - maximum number of iterations  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine progress_bar(itera, max_iter)
        use, intrinsic :: iso_fortran_env, only: output_unit
        integer, intent(in) :: itera, max_iter
        integer          :: uu
        character(len=1) :: bar, back, dot
        back = char(8)
        bar  = '='
        dot  = ' '
        write(output_unit,'(256a1)', advance='no') (back, uu =1,30+10)
        flush(output_unit)
        write(output_unit,'(1x, 1i3,1a1,2x,1a1,256a1,1a1,256a1,1a1)', advance='no') 100*itera/max_iter,'%','[', &
        (bar, uu =1,30*itera/max_iter), '>', (dot, uu=1,(30-30*itera/max_iter)), ']'
        flush(output_unit)
    end subroutine progress_bar
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Function which computes the factorial of an     !
    !  integer number                                  !
    !  input: n - natural number                       !
    !  output: f - factorial                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive function fact(n) result (f)
        integer, intent(in) :: n
        real(kind=16) :: f
    
        if (n .gt. 0) then
        f= dble(n)*fact(n-1)
        else
        f=1.d0
        end if
    
    end function fact
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Routine which computes the determinant of the   !
    !  matrix whose column are the vectors x,y,z       !
    !  input: x,y,z - columns of the matrix            !
    !  output: detJ - output                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine deter_sub(x,y,z,detJ)
        real(kind=16), dimension(3),intent(in):: x,y,z
        real(kind=16),intent(out) :: detJ
    
        detJ= real(x(1)*y(2)*z(3)+x(2)*y(3)*z(1)+x(3)*y(1)*z(2)-(&
        &x(1)*y(3)*z(2) + x(2)*y(1)*z(3)+x(3)*y(2)*z(1)),16)
    end subroutine deter_sub

    subroutine time_converter(time,type,time_out)
        real(kind=16), intent(in) :: time
        character(len=5), intent(in) :: type
        real(kind=16), intent(out) :: time_out
        if (type .eq. "s") then
            time_out=time
        else if (type .eq. "m") then
            time_out=time*60.d0
        else if (type .eq. "h") then
            time_out=time*3600.d0
        else if (type .eq. "d") then
            time_out=time*day
        else if (type .eq. "y") then
            time_out=time*365.25d0*day
        end if

        return
    end subroutine time_converter

    subroutine time_converter_inverse(time,type,time_out)
        real(kind=16), intent(in) :: time
        real(kind=16), intent(out) :: time_out
        character(len=5), intent(in) :: type

        if (type .eq. "s") then
            time_out=time
        else if (type .eq. "m") then
            time_out=time/60.d0
        else if (type .eq. "h") then
            time_out=time/3600.d0
        else if (type .eq. "d") then
            time_out=time/day
        else if (type .eq. "y") then
            time_out=time/(365.25d0*day)
        end if
        
        return
    end subroutine time_converter_inverse

    subroutine init_srp_con()
        istar=istardeg*degree
        omstar=omstardeg*degree
        ut=(365.24219d0/366.24219d0)*day/pi2
        crpr=4.56d-3*ageo*ut**2
        n_s=35999.049*365.242196d0/(36525.6363d0*366.242196d0)/360.d0 
        ! rE=rE/ageo
        return
    end subroutine init_srp_con

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sun(time, xs)
    real(kind=16), intent(out), dimension(4) :: xs
    real(kind=16), intent(in) :: time
    real(kind=16) :: eps, time_aux,ttime, anmns, slam, rs

    time_aux=time
    eps=23.4392911d0*degree
    ! ttime=(365.242196d0*time_aux)/(36525.d0*366.242196d0*pi2)
    ttime=(time/day)/36525.d0

    anmns=357.5256d0*degree+35999.049d0*degree*ttime

    slam=282.94d0*degree+anmns+(6892.d0/3600.d0)*degree*sin(anmns) &
  &  + (72.d0/3600.d0)*degree*sin(2.d0*anmns)

    rs= (149.619d0-2.499d0*cos(anmns) &
  &  -0.021d0*cos(2.d0*anmns))!*23.71681950544094d0
    rs=rs*1.d9!*ageo*1000.d0

    xs(1)=rs*cos(slam)
    xs(2)=rs*sin(slam)*cos(eps)
    xs(3)=rs*sin(slam)*sin(eps)
    xs(4)=rs
    return
end subroutine sun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine moon(time,xm)
    real(kind=16), intent(out), dimension(4) :: xm
    real(kind=16), intent(in) :: time
    real(kind=16) :: eps,ttime,ool0, ool, oolpr, ff, d
    real(kind=16) :: oolam, betmo, rm

     eps=23.4392911d0*degree

     ttime=(time/day)/36525.d0


    ! L_0 in Montenbruck
        ool0=(218.31617d0+481267.88088d0*ttime &
    &  -4.06d0/3600.d0*ttime**2)*degree
    ! l in M.
        ool=(134.96292d0+477198.86753d0*ttime)*degree
    ! l' in M.
        oolpr=(357.52543d0+35999.04944d0*ttime)*degree
    ! F in M.
        ff=(93.27283d0+483202.01873d0*ttime)*degree
    ! D in M.
        d=(297.85027d0+445267.11135d0*ttime)*degree

    oolam=ool0+(1.d0/3600.d0)*degree*(22640.d0*sin(ool)    &
  &  +769.d0*sin(2.d0*ool)-4586.d0*sin(ool-2.d0*d)   &
  &  +2370.d0*sin(2.d0*d)-668.d0*sin(oolpr)-412.d0*sin(2.d0*ff)  &
  &  -212.d0*sin(2.d0*ool-2.d0*d)-206.d0*sin(ool+oolpr-2.d0*d) &
  &  +192.d0*sin(ool+2.d0*d)-165.d0*sin(oolpr-2.d0*d) &
  &  +148.d0*sin(ool-oolpr)-125.d0*sin(d) &
  &  -110.d0*sin(ool+oolpr)-55.d0*sin(2.d0*ff-2.d0*d)) 



     betmo=  (1.d0/3600.d0)*degree*(18520.d0*sin(ff+oolam-ool0 &
  &  +(1.d0/3600.d0)*degree*412.d0*sin(2.d0*ff) &
  & + (1.d0/3600.d0)*degree*541.d0*sin(oolpr)) &
  &  - 526.d0*sin(ff-2.d0*d) &
  &  +44.d0*sin(ool+ff-2.d0*d)-31.d0*sin(-ool+ff-2.d0*d) &
  &  -25.d0*sin(-2.d0*ool+ff)-23.d0*sin(oolpr+ff-2.d0*d) &
  &  +21.d0*sin(-ool+ff)+11.d0*sin(-oolpr+ff-2.d0*d))


    rm=(385000.d0-20905.d0*cos(ool)-3699.d0*cos(2.d0*d-ool) &
  &  -2956.d0*cos(2.d0*d)-570*cos(2.d0*ool)+246.d0*cos(2.d0*ool &
  &  -2.d0*d)-205.d0*cos(oolpr-2.d0*d)-171.d0*cos(ool+2.d0*d) &
  &  -152.d0*cos(ool+oolpr-2.d0*d))

    rm=rm*1.d3
    xm(1)=rm*cos(oolam)*cos(betmo)
    xm(2)=rm*sin(oolam)*cos(betmo)*cos(eps) &
  &  - rm*sin(betmo)*sin(eps)
    xm(3)=rm*sin(oolam)*cos(betmo)*sin(eps) &
  &  + rm*sin(betmo)*cos(eps)
    xm(4)=rm
    return

end subroutine moon


end module constants
