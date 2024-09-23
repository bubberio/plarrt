!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the RHS of some space dynamics problems   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module dynamics
    use geometry
    use constants
    use polygrav
    use spherigrav
    implicit none
    contains

subroutine sd_car_dyn(time,state,acc_car)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Secular Space Debris dynamics                    !
    !  input: state - Ham. state of the satellite  	   !
    !         time - time at which to compute accel.   !
    ! output: acc_ham - secular acceleration due to J2 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
    real(kind=16), intent(in), dimension(6) :: state
    real(kind=16), intent(in) :: time
    real(kind=16), intent(out), dimension(6) :: acc_car
    real(kind=16), dimension(6) :: acc_aux, accJ2, acc22
    real(kind=16), dimension(6) :: accSRP, accSun, accMoon, accSH
    integer :: i
    do i=1,6
        acc_aux(i)=0.d0
        acc_car(i)=0.d0
    end do

    if (flag_full_sh) then ! Earth with full SH coefficients
        call SH_dyn(time,state,accSH)
        do i=1,6
            acc_aux(i)=acc_aux(i)+accSH(i)
        end do
    else 
        if (dyn_vec(1)) then ! Dynamics with at least J2=-C20
            if (flag_alt_j2) then
                call car_Earth_22(time,state,accJ2)
            else
                call car_Earth_J2(time, state, accJ2)
            end if
            do i=1,6
                acc_aux(i)=acc_aux(i)+accJ2(i)
            end do
        end if
    end if

    if (dyn_vec(2)) then
        call car_SRP(time,state,accSRP)
        do i=4,6
            acc_aux(i)=acc_aux(i)+accSRP(i)
        end do
    end if
    if (dyn_vec(3)) then
        call car_Sun(time,state,accSun)
        do i=1,6
            acc_aux(i)=acc_aux(i)+accSun(i)
        end do
    end if
    if (dyn_vec(4)) then
        call car_Moon(time,state,accMoon)
        do i=1,6
            acc_aux(i)=acc_aux(i)+accMoon(i)
        end do
    end if
    do i=1,6
        acc_car(i)=acc_aux(i)
    end do

    return
end subroutine sd_car_dyn
! We first collect Cartesian dynamics (and their variational equation if necessary)
subroutine car_Earth_J2(time,state,  accJ2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Routine which compute the Cartesian acceleration !
    ! due to the Earth's J2	                           !
    !  input: state - Cartesian state of the satellite !
    !         time - time at which to compute accel.   !
    ! output: accJ2 - Cartesian acceleration due to J2 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=16), intent(in), dimension(6) :: state
    real(kind=16), intent(in) :: time
    real(kind=16), intent(out), dimension(6) :: accJ2
    real(kind=16) ::C20, auxC20, aux
    real(kind=16) ::x,y,z,r,reN
    integer :: i
    reN=rE*1.d3 ! rE is stored in constants.f90
    C20=-eJ2    ! eJ2 is stored in constants.f90
    x=state(1)
    y=state(2)
    z=state(3)
    r=sqrt(x**2+y**2+z**2)
    aux=-mu/r**3! mu is declared in constants.f90
                ! it should be initialized before calling  
                ! this subroutine
    auxC20=-aux*C20*((reN/r)**2)*(1.5d0 - 7.5d0*(z/r)**2)

    do i=1,6
        accJ2(i)=0.0d0
    end do
    
    accJ2(1)=state(4)
    accJ2(2)=state(5)
    accJ2(3)=state(6)
    accJ2(4)= (aux+auxC20)*x
    accJ2(5)= (aux+auxC20)*y
    accJ2(6)= (aux+auxC20-3*aux*C20*((reN/r)**2))*z
    return
end subroutine car_Earth_J2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine car_Earth_22( time,state, accJ22)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Routine which compute the Cartesian acceleration !
    ! due to the Earth's C22 and S22 SH coefficients.  !
    !  input: state - Cartesian state of the satellite !
    !         time - time at which to compute accel.   !
    ! output: accJ2 - Cartesian acceleration due to J2 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(kind=16), intent(in), dimension(6) :: state
	real(kind=16), intent(in) :: time
	real(kind=16), intent(out), dimension(6) :: accJ22
	real(kind=16) :: C20,C22, S22, CSp, CSm, aux, auxC20
	real(kind=16) :: x, y, z,r, theta,reN
	integer :: i
    reN=rE*1.d3
	x=state(1)
	y=state(2)
	z=state(3)
	r=sqrt(x**2+y**2+z**2)
    theta=rate*time
    C20=-eJ2 
    C22=1.5746153d-6
    S22=-9.0387279d-7
	CSm=C22*cos(2*theta)-S22*sin(2*theta)
	CSp=C22*sin(2*theta)+S22*cos(2*theta)
	aux=-mu/r**3
    auxC20=-aux*((reN/r)**2)*(C20*(1.5d0-7.5d0*(z/r)**2) &
    & + (15.d0/(r**2))*(CSm*(y**2-x**2)-2*x*y*CSp))
	! auxC20=-aux*C20*((rE/r)**2)*(1.5d0-7.5d0*(z/r)**2 &
    ! & + (15/(r**2))*(CSm*(y**2-x**2)-2*x*y*CSp))
    ! old version - Check if this is wrong
	do i=1,6
		accJ22(i)=0.0d0
	end do
	accJ22(1)=state(4)
	accJ22(2)=state(5)
	accJ22(3)=state(6)
	accJ22(4)= (aux+auxC20)*x -aux*((reN/r)**2)*6*(CSm*x+CSp*y)
	accJ22(5)= (aux+auxC20)*y -aux*((reN/r)**2)*6*(CSp*x-CSm*y)
	accJ22(6)= (aux+auxC20)*z -aux*((reN/r)**2)*C20*3.d0*z
	return
    ! old version
    ! accJ22(4)= (aux+auxC20)*x -aux*C20*((rE/r)**2)*6*(CSm*x+CSp*y)
	! accJ22(5)= (aux+auxC20)*y -aux*C20*((rE/r)**2)*6*(CSm*x-CSp*y)
endsubroutine car_Earth_22
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine car_SRP(time,state,accSRP)
    real(kind=16), intent(in), dimension(6) :: state
    real(kind=16), intent(in) :: time
    real(kind=16), intent(out), dimension(6) :: accSRP
    real(kind=16), dimension(4) :: pos_sun
    real(kind=16), dimension(3) ::  pos, diff
    real(kind=16) :: diffnorm, coef, time_aux
    integer :: i
    
    time_aux=time
    do i=1,3
        accSRP(i)=0.d0
    end do
    call sun(time_aux,pos_sun)

    do i=1,3
        diff(i)=state(i)-pos_sun(i)
    end do
    diffnorm=sqrt(diff(1)**2+diff(2)**2+diff(3)**2)
    coef=crpr_car*astrounit**2*are2ms/diffnorm**3
    do i=1,3
        accSRP(i)=0.d0
        accSRP(i+3)=coef*diff(i)
    end do
    return
end subroutine car_SRP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine car_Sun(time,state,accSun)
    real(kind=16), intent(in), dimension(6) :: state
    real(kind=16), intent(in) :: time
    real(kind=16), intent(out), dimension(6) :: accSun
    real(kind=16), dimension(4) :: pos_sun
    real(kind=16), dimension(3) :: pos, diff
    real(kind=16) :: diffnorm, coef, time_aux, norm_pos_sun
    integer :: i
    
    time_aux=time
    do i=1,3
        accSun(i)=0.d0
    end do
    call sun(time_aux,pos_sun)
    norm_pos_sun=pos_sun(4)
    do i=1,3
        diff(i)=state(i)-pos_sun(i)
    end do
    diffnorm=sqrt(diff(1)**2+diff(2)**2+diff(3)**2)
    coef=-G*maS/diffnorm**3
    do i=1,3
        accSun(i)=0.d0
        accSun(i+3)=coef*diff(i)-G*maS*pos_sun(i)/norm_pos_sun**3
    end do
    return
end subroutine car_Sun
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine car_Moon(time,state,accMoon)
    real(kind=16), intent(in), dimension(6) :: state
    real(kind=16), intent(in) :: time
    real(kind=16), intent(out), dimension(6) :: accMoon
    real(kind=16), dimension(4) :: pos_moon
    real(kind=16), dimension(3) :: pos, diff
    real(kind=16) :: diffnorm, coef, time_aux, norm_pos_moon
    integer :: i
    
    time_aux=time
    do i=1,3
        accMoon(i)=0.d0
    end do
    call moon(time_aux,pos_moon)
    norm_pos_moon=pos_moon(4)
    do i=1,3
        diff(i)=state(i)-pos_moon(i)
    end do
    diffnorm=sqrt(diff(1)**2+diff(2)**2+diff(3)**2)
    coef=-G*maM/diffnorm**3
    do i=1,3
        accMoon(i)=0.d0
        accMoon(i+3)=coef*diff(i)-G*maM*pos_moon(i)/norm_pos_moon**3
    end do
    return
end subroutine car_Moon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sd_ham_dyn(time,state,acc_ham)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Secular Space Debris dynamics                    !
    !  input: state - Ham. state of the satellite  	   !
    !         time - time at which to compute accel.   !
    ! output: acc_ham - secular acceleration due to J2 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
    real(kind=16), intent(in), dimension(4) :: state
    real(kind=16), intent(in) :: time
    real(kind=16), intent(out), dimension(4) :: acc_ham
    real(kind=16), dimension(4) :: acc_aux, accJ2, accSRP, accSun, accMoon
    integer :: i
    do i=1,4
        acc_aux(i)=0.d0
        acc_ham(i)=0.d0
    end do
    if (dyn_vec(1)) then
        call secular_Earth_J2(time,state,accJ2)
        do i=1,4
            acc_aux(i)=acc_aux(i)+accJ2(i)
        end do   
    end if
    if (dyn_vec(2)) then
        if (flag_all_tm) then
        call secular_SRP_all(time,state, accSRP)
        else
            call secular_SRP_jk(time,state, accSRP)
        end if
        do i=1,4
            acc_aux(i)=acc_aux(i)+accSRP(i)
        end do  
    end if

    if (dyn_vec(3)) then
        ! DA_Sun
    end if
    if (dyn_vec(4)) then
        ! DA_Moon
    end if
    do i=1,4
        acc_ham(i)=acc_aux(i)
    end do
    return
end subroutine sd_ham_dyn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following dynamics are secular and Hamiltonian
subroutine secular_Earth_J2( time, x, accJ2)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Secular J2 dynamics                              !
	!  input: x - partial Ham. state of the satellite  !
	!         time - time at which to compute accel.   !
	!         rL - first Delauney momentum (const.)    !
	! output: accJ2 - secular acceleration due to J2   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	real(kind=16), intent(in), dimension(4) :: x
	real(kind=16), intent(in) :: time
	real(kind=16), intent(out), dimension(4) :: accJ2
	real(kind=16) :: rH, rG, J2,reN
	integer :: i
    ! rL is declared in constants.f90 and should be 
    ! initialized upon calling this subroutine
	reN=rE/ageo
	J2=eJ2
	rG=x(1)
	rH=x(2)
	do i=1,4
		accJ2(i)=0.0d0
	end do
	accJ2(1)=0.d0
	accJ2(2)=0.d0
	accJ2(3)= (-3*J2*reN**2*(rG**2 - 5*rH**2))/(4.*rG**6*rL**3)
	accJ2(4)= (-3*J2*reN**2*rH)/(2.*rG**5*rL**3)
	return
endsubroutine secular_Earth_J2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine secular_SRP_jk( time,x, accSRP)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Secular SRP dynamics, depending on j and k       !
    ! (global constants)                               !
    !  input: x - partial Ham. state of the satellite  !
    !         time - time at which to compute accel.   !
    !         rL - first Delauney momentum (const.)    !
    ! output: varJ2 - RHS of the variational equation  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=16), intent(in), dimension(4) :: x
      real(kind=16), intent(in) :: time
      real(kind=16), intent(out), dimension(4) :: accSRP
      real(kind=16) :: ttime, Mstar, rH, rG, resang, om, bom
      real(kind=16) :: coef, f, auxg, auxh, auxc, auxjh
      integer :: i
      ! rL is declared in constants.f90 and should be 
      ! initialized upon calling this subroutine
      coef= (3.d0/2.d0)*crpr*a2m*rL**2/muE 
      ! write(*,*) are2ms, a2m
      ! stop
      ! crpr should be initialized using the 
      !muE should be 1 in the scaled variables.
      call fjk(jsrp,ksrp,f)
  
      ttime=(365.242196d0*time)/(36525.6363d0*366.242196d0*pi2)
        Mstar=(357.5256+35999.049*ttime)*degree !pi2/360.d0  
   !	Mstar=(35999.049*ttime)*degree !*pi2/360.d0
      rG=x(1)
      rH=x(2)
      om=x(3)
      bom=x(4)
      resang=om + jsrp*bom + ksrp*(Mstar+omstar)
  
      auxh=sqrt(1.d0 - (rH/rG)**2)
      auxjh=(1.d0 + jsrp*rH/rG)
      auxg=sqrt(1.d0 - (rG/rL)**2)
      auxc=coef*f
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initializing the accSRP vector to zero.          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,4
          accSRP(i)=0.0d0
      end do
  
      if (jsrp .eq. 0) then
          accSRP(1)=-auxc*auxh*auxg*sin(resang)
          accSRP(2)=0.d0
          accSRP(3)=auxc*((rH**2/rG**3)*auxg/auxh - (rG/rL**2)*auxh/auxg)*cos(resang)
          accSRP(4)=-auxc*(rH/rG**2)*(auxg/auxh)*cos(resang)
      else 
          accSRP(1)=-auxc*auxjh*auxg*sin(resang)
          accSRP(2)=-auxc*auxjh*auxg*real(jsrp)*sin(resang)
          accSRP(3)=-auxc*(rH*jsrp*auxg/rG**2 + (rG/rL**2)*auxjh/auxg)*cos(resang)
          accSRP(4)=auxc*jsrp*(auxg/rG)*cos(resang)
      end if
      return
endsubroutine secular_SRP_jk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine secular_SRP_all( time,x, accSRP)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Secular SRP dynamics, depending on j and k       !
    ! (global constants)                               !
    !  input: x - partial Ham. state of the satellite  !
    !         time - time at which to compute accel.   !
    !         rL - first Delauney momentum (const.)    !
    ! output: varJ2 - RHS of the variational equation  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=16), intent(in), dimension(4) :: x
      real(kind=16), intent(in) :: time
      real(kind=16), intent(out), dimension(4) :: accSRP
      real(kind=16),  dimension(4) :: acc_aux
      integer :: i
      do i=1,4
        accSRP(i)=0.d0
      end do 
      jsrp=1
      ksrp=-1
      call secular_SRP_jk( time,x, acc_aux)
      do i=1,4
        accSRP(i)=accSRP(i)+acc_aux(i)
      end do  
      ksrp=1
      call secular_SRP_jk( time,x, acc_aux)
      do i=1,4
        accSRP(i)=accSRP(i)+acc_aux(i)
      end do 
      jsrp=-1
      ksrp=-1
      call secular_SRP_jk( time,x, acc_aux)
      do i=1,4
        accSRP(i)=accSRP(i)+acc_aux(i)
      end do 
      ksrp=1
      call secular_SRP_jk( time,x, acc_aux)
      do i=1,4
        accSRP(i)=accSRP(i)+acc_aux(i)
      end do 
      jsrp=0
      ksrp=-1
      call secular_SRP_jk( time,x, acc_aux)
      do i=1,4
        accSRP(i)=accSRP(i)+acc_aux(i)
      end do 
      ksrp=1
      call secular_SRP_jk( time,x, acc_aux)
      do i=1,4
        accSRP(i)=accSRP(i)+acc_aux(i)
      end do 
      return
endsubroutine secular_SRP_all
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine DA_Sun( time,x, accSun)
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ! Secular SRP dynamics, depending on j and k       !
!     ! (global constants)                               !
!     !  input: x - partial Ham. state of the satellite  !
!     !         time - time at which to compute accel.   !
!     !         rL - first Delauney momentum (const.)    !
!     ! output: varJ2 - RHS of the variational equation  !
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       real(kind=16), intent(in), dimension(4) :: x
!       real(kind=16), intent(in) :: time
!       real(kind=16), intent(out), dimension(4) :: accSRP
!       real(kind=16),  dimension(4) :: acc_aux
!       integer :: i
!       do i=1,4
!         accSRP(i)=0.d0
!       end do 

!       return
! endsubroutine DA_Sun
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine DA_Moon( time,x, accMoon)
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ! Secular SRP dynamics, depending on j and k       !
!     ! (global constants)                               !
!     !  input: x - partial Ham. state of the satellite  !
!     !         time - time at which to compute accel.   !
!     !         rL - first Delauney momentum (const.)    !
!     ! output: varJ2 - RHS of the variational equation  !
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       real(kind=16), intent(in), dimension(4) :: x
!       real(kind=16), intent(in) :: time
!       real(kind=16), intent(out), dimension(4) :: accSRP
!       real(kind=16),  dimension(4) :: acc_aux
!       integer :: i
!       do i=1,4
!         accSRP(i)=0.d0
!       end do 
!       jsrp=1
!       ksrp=-1
!       call secular_SRP_jk( time,x, acc_aux)
!       do i=1,4
!         accSRP(i)=accSRP(i)+acc_aux(i)
!       end do  
!       ksrp=1
!       call secular_SRP_jk( time,x, acc_aux)
!       do i=1,4
!         accSRP(i)=accSRP(i)+acc_aux(i)
!       end do 
!       jsrp=-1
!       ksrp=-1
!       call secular_SRP_jk( time,x, acc_aux)
!       do i=1,4
!         accSRP(i)=accSRP(i)+acc_aux(i)
!       end do 
!       ksrp=1
!       call secular_SRP_jk( time,x, acc_aux)
!       do i=1,4
!         accSRP(i)=accSRP(i)+acc_aux(i)
!       end do 
!       jsrp=0
!       ksrp=-1
!       call secular_SRP_jk( time,x, acc_aux)
!       do i=1,4
!         accSRP(i)=accSRP(i)+acc_aux(i)
!       end do 
!       ksrp=1
!       call secular_SRP_jk( time,x, acc_aux)
!       do i=1,4
!         accSRP(i)=accSRP(i)+acc_aux(i)
!       end do 
!       return
! endsubroutine DA_Moon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine poly_dyn(time,x,accpoly)
    real(kind=16), intent(in), dimension(6) :: x
	real(kind=16), intent(in) :: time
	real(kind=16), intent(out), dimension(6) :: accpoly
	! real(kind=16) :: pot
    real(kind=16), dimension(3) :: pos,vel,acc
    ! real(kind=16),dimension(3,3) :: ggmat
	integer :: j
    do j=1,3
        pos(j)=x(j)
        vel(j)=x(j+3)
        acc(j)=0.0d0
      end do	
    ! call poly_all(pos,ver,fac, edg, sigma,rate,time,pot,acc,ggmat,lap_check)
    call polyacc_rot(pos, ver, fac, edg, sigma, rate,time,acc)
    
    do j=1,3
        accpoly(j)=vel(j)
        accpoly(j+3)=acc(j)
      end do
end subroutine poly_dyn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SH_dyn(time,x,accSH)
    real(kind=16), intent(in), dimension(6) :: x
	real(kind=16), intent(in) :: time
	real(kind=16), intent(out), dimension(6) :: accSH
    real(kind=16), dimension(3) :: pos,vel,acc
	integer :: j
    do j=1,3
        pos(j)=x(j)
        vel(j)=x(j+3)
        acc(j)=0.0d0
      end do	
      call spheriacc_rot(pos,rr,Ma,C,S,rate,time,acc)
      do j=1,3
        accSH(j)=vel(j)
        accSH(j+3)=acc(j)
      end do
end subroutine SH_dyn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The next routine should be used together with the one
! concerning the SRP secular effect
subroutine var_secular_Earth_J2(time,x, v, varJ2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RHS of the variational equation for the 		   !
    ! Secular J2 dynamics                              !
    !  input: x - partial Ham. state of the satellite  !
    !         v - tangent vector                       !
    !         time - time at which to compute accel.   !
    !         rL - first Delauney momentum (const.)    !
    ! output: varJ2 - RHS of the variational equation  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=16), intent(in), dimension(4) :: x, v
      real(kind=16), intent(in) :: time
      real(kind=16), intent(out), dimension(4) :: varJ2
      real(kind=16) :: rH, rG, J2, reN
      real(kind=16) :: df1dG,df1dH,df1dom,df1dbo
      real(kind=16) :: df2dG,df2dH,df2dom,df2dbo
      real(kind=16) :: df3dG,df3dH,df3dom,df3dbo
      real(kind=16) :: df4dG,df4dH,df4dom,df4dbo
      integer :: i
      ! rL is declared in constants.f90 and should be 
      ! initialized upon calling this subroutine
      rG=x(1)
      rH=x(2)
      J2=eJ2
      reN=rE/ageo
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initializing the varJ2 vector to zero.           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,4
          varJ2(i)=0.0d0
      end do
  
      df1dG= 0.d0
      df1dH= 0.d0
      df1dom= 0.d0
      df1dbo= 0.d0
  
      df2dG= 0.d0
      df2dH= 0.d0
      df2dom= df1dbo
      df2dbo= 0.d0
  
    df3dG= (3*J2*(2*rG**2-15*rH**2)*reN**2)/(2.d0*rG**7*rL**3)
    !(-3*J2*reN**2)/(2.d0*rG**5*rL**3) + (9*J2*reN**2*(G**2 - 5*rH**2))/ &
    !     &  (2.d0*rG**7*rL**3)
      df3dH= (15*J2*reN**2*rH)/(2.d0*rG**6*rL**3)
      df3dom= -df1dG
      df3dbo= -df2dG
  
      df4dG=df3dH
      df4dH= (-3*J2*reN**2)/(2.d0*rG**5*rL**3)
      df4dom=-df1dH
      df4dbo=-df2dH
      
      varJ2(1)=(df1dG*v(1)+df1dH*v(2)+df1dom*v(3)+df1dbo*v(4))
      varJ2(2)=(df2dG*v(1)+df2dH*v(2)+df2dom*v(3)+df2dbo*v(4))
      varJ2(3)=(df3dG*v(1)+df3dH*v(2)+df3dom*v(3)+df3dbo*v(4))
      varJ2(4)=(df4dG*v(1)+df4dH*v(2)+df4dom*v(3)+df4dbo*v(4))
  
      return	
end subroutine var_secular_Earth_J2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine var_secular_SRP_jk(time,x, v, varSRP)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RHS of the variational equation for the          !
    ! Secular SRP dynamics, depending on j and k       !
    ! (global constants)                               !
    !  input: x - partial Ham. state of the satellite  !
    !         v - tangent vector                       !
    !         time - time at which to compute accel.   !
    !         rL - first Delauney momentum (const.)    !
    ! output: varJ2 - RHS of the variational equation  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=16), intent(in), dimension(4) :: x, v
    real(kind=16), intent(in) :: time
    real(kind=16), intent(out), dimension(4) :: varSRP
    real(kind=16) :: ttime, Mstar, rH, rG, om, bom, resang
    real(kind=16) :: auxh, auxjh, auxg, auxc, coef, f
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=16) :: df1dG,df1dH,df1dom,df1dbo
    real(kind=16) :: df2dG,df2dH,df2dom,df2dbo
    real(kind=16) :: df3dG,df3dH,df3dom,df3dbo
    real(kind=16) :: df4dG,df4dH,df4dom,df4dbo
    integer :: i
    coef= -(3.d0/2.d0)*crpr*a2m*rL**2!/muE !muE should be 1.
    ttime=(365.242196d0*time)/(36525.6363d0*366.242196d0*pi2)
    !	Mstar=(357.5256+35999.049*ttime)*pi2/360.d0  
    Mstar=(35999.049*ttime)*pi2/360.d0
    rG=x(1)
    rH=x(2)
    om=x(3)
    bom=x(4)
    resang=om + real(jsrp)*bom + real(ksrp)*(Mstar+omstar)
    call fjk(jsrp,ksrp,f)
    auxh=sqrt(1.d0 - (rH/rG)**2)
    auxjh=(1.d0 + jsrp*rH/rG)
    auxg=sqrt(1.d0 - (rG/rL)**2)
    auxc=coef*f
    do i=1,4
        varSRP(i)=0.0d0
    end do
    if (jsrp .eq. 0) then
        df1dG= -auxc*sin(resang)*((rH**2/rG**3)*auxg/auxh - (rG/rL**2)*auxh/auxg)
        df1dH= auxc*rH*sin(resang)*(auxg/auxh)/rG**2
        df1dom= -auxc*auxh*auxg*cos(resang)
        df1dbo= 0.d0

        df2dG= 0.d0
        df2dH= 0.d0
        df2dom= df1dbo
        df2dbo= 0.d0

        !df3dG= auxc*cos(resang)*(-2*(3*rG**4+(rH*rL)**2)*(rH*rL)**2 +(rH**2+rL**2)*&
        !	& (rG**4 + 3*(rH*rL)**2))/((rG**4)*(rG**2-rH**2)*(rG**2-rL**2)*auxh*auxg*rL**2)

        df3dG= auxc*cos(resang)*((-(auxh**4*rG**8) - auxg**2*auxh**4*rG**6*rL**2 - &
        &    2*auxg**2*auxh**2*rG**4*rH**2*rL**2 - &
        &    3*auxg**4*auxh**2*rG**2*rH**2*rL**4 - auxg**4*rH**4*rL**4)/ &
        &  (auxg**3*auxh**3*rG**6*rL**4))

        !df3dH= -auxc*cos(resang)*(rH/rG**3)*(rG**4-2*(rG*rL)**2+(rH*rL)**2)/&
        !	& ((rG**2-rH**2)*auxh*auxg*rL**2)

        df3dH= auxc*cos(resang)*((auxh**2*rG**4*rH + 2*auxg**2*auxh**2*rG**2*rH*rL**2 + &
        &    auxg**2*rH**3*rL**2)/(auxg*auxh**3*rG**5*rL**2))
        df3dom= -df1dG
        df3dbo= -df2dG

        df4dG= df3dH
        df4dH= -auxc*cos(resang)*auxg*(auxh**2*rG**2 + rH**2)/(auxh**3*rG**4)
        df4dom= -df1dH
        df4dbo= -df2dH
    else

        df1dG= auxc*sin(resang)*(rG + jsrp*rH*(1.d0 + (rL*auxg/rG)**2))/(auxg*rL**2)
        df1dH= -auxc*auxg*jsrp*sin(resang)/rG
        df1dom=-auxc*auxjh*auxg*cos(resang)
        df1dbo=-auxc*auxjh*auxg*jsrp*cos(resang)

        df2dG= auxc*jsrp*sin(resang)*(rG + jsrp*rH*(1.d0 + (rL*auxg/rG)**2))/(auxg*rL**2)
        df2dH= auxc*auxg*jsrp*jsrp*sin(resang)/rG
        df2dom= df1dbo
        df2dbo= -auxc*auxjh*auxg*cos(resang)*jsrp**2

        df3dG= -auxc*cos(resang)*(rG**2+ (auxg*rL)**2)*(auxjh*rG**3 &
    & -2.d0*rH*jsrp*(auxg*rL)**2)/((rL**4)*(rG*auxg)**3)

        df3dH= -auxc*cos(resang)*jsrp*(1+ (auxg*rL/rG)**2)/(auxg*rL**2)
        df3dom= -df1dG
        df3dbo= -df2dG

        df4dG= df3dH
        df4dH=0.d0
        df4dom= -df1dH
        df4dbo= -df2dH
    end if
    varSRP(1)=(df1dG*v(1)+df1dH*v(2)+df1dom*v(3)+df1dbo*v(4))
    varSRP(2)=(df2dG*v(1)+df2dH*v(2)+df2dom*v(3)+df2dbo*v(4))
    varSRP(3)=(df3dG*v(1)+df3dH*v(2)+df3dom*v(3)+df3dbo*v(4))
    varSRP(4)=(df4dG*v(1)+df4dH*v(2)+df4dom*v(3)+df4dbo*v(4))
    return

endsubroutine var_secular_SRP_jk
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fjk(j,k,f)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Routine which produces the inclination functions !
    ! which appear in the SRP toy models               !
    !  input: j,k - indices  which label the toy models!
    ! output: f - value of the inclination function j,k!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer,intent(in) :: j,k 
      !No declaration because I am trying with global variables
      real(kind=16), intent(out) :: f
      istar=istardeg*degree
      f=0.d0
  
      if (jsrp .eq. 0) then
          f=-real(k)*cos(istar/2.d0)*sin(istar/2.d0)
      else 
          if (j*k .gt. 0) then
              f= 0.5d0*(sin(istar/2.d0))**2
          else if (j*k .lt. 0) then
              f= 0.5d0*(cos(istar/2.d0))**2
          end if
      end if
      return
end subroutine fjk
! Dynamics of the pendulum for the examples of Chapter 2.
subroutine pendulum(time,x, dx)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Simple Pendulum dynamics                         !
    !  input: x - (theta,dtheta) state of the pendulum !
    !         time - time at which to compute accel.   !
    !         l - length of the pendulum (const.)      !
    ! output: dx - vector field for the simple pendulum!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
    real(kind=16), intent(in), dimension(2) :: x
    real(kind=16), intent(in) :: time
    real(kind=16), intent(out), dimension(2) :: dx
    real(kind=16) :: theta, dtheta, grav
    integer :: i
    grav=9.81d0
    theta=x(1)
    dtheta=x(2)
    do i=1,2
        dx(i)=0.0d0
    end do
    dx(1)=dtheta
    dx(2)=-(grav/lon)*sin(theta) ! lon is declared in constants
    ! it should be initialized upon call of this function
    return
endsubroutine pendulum 

subroutine var_pendulum(time,x,v, dvar)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Simple Pendulum variational equation             !
	!  input: x - (theta,dtheta) state of the pendulum !
	!         v - tangent vector					   !
	!         time - time at which to compute accel.   !
	!         l - length of the pendulum (const.)      !
	! output: dvar - output of one iter. of var. eq.   !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	real(kind=16), intent(in), dimension(2) :: x,v
	real(kind=16), intent(in) :: time
	real(kind=16), intent(out), dimension(2) :: dvar
	real(kind=16) :: theta, grav
	integer :: i
	grav=9.81d0
	theta=x(1)
	! write(*,*) theta
	do i=1,2
		dvar(i)=0.0d0
	end do
	! write(*,*) 'Inside VAR subroutine'
	dvar(1)=v(2)
	dvar(2)=-v(1)*(grav/lon)*cos(theta)
	return
endsubroutine var_pendulum

subroutine var_SH(time,v,dv)
    real(kind=16), intent(in), dimension(6) :: v
	real(kind=16), intent(in) :: time
	real(kind=16), intent(out), dimension(6) :: dv
    real(kind=16), dimension(3) :: aux_v, aux_dv, pos
    real(kind=16), dimension(3,3) :: ggm
	integer :: j
    do j=1,3
        aux_v(j)=v(j)
        pos(j)=aux_car_state(j)
    end do
      call spheriggm_rot(pos,rr,Ma,C,S,rate,time,ggm)
      do j=1,3
        aux_dv(j)=v(1)*ggm(j,1)+v(2)*ggm(j,2)+v(3)*ggm(j,3)
      end do
    !   aux_dv=matmul(ggm,aux_v)
      do j=1,3
        dv(j)=v(j+3)
        dv(j+3)=aux_dv(j)
      end do
end subroutine var_SH
end module dynamics