!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	General propagator                                                 !
!   To Be Set in namelist.nml:                                         !
!       . /OBJECT/ filename, sca, sigma, period                        !
!       . /SH/ nmax                                                    !
!       . /ORBEL/ a, e, inc, om, bom, amean                            ! 
!       . /TIME/ time0, tot_time, time_step, type_time                 !
!       . /DYN/ dyn_vec, flag_full_sh, flag_alt_j2, flag_ABM           !
!       . /SRP/  are2ms, jsrp, ksrp, resang0, flag_all_tm              !
!       . /FLAGS_PROP/ flag_sd, flag_ham, flag_poly                    !
!       Description:                                                   !
!       - filename: name of the asteroid .obj file                     !
!       - sca: scale at which the .obj model is given (1=m,1000=km,ecc.)
!       - sigma: mean density of the body, given in [g/cm^3]           !
!       - period: rotation period of the body, given in seconds        !
!                                                                      !
!       - nmax: maximum degree of spherical harmonics to be used       !
!           (requires  flag_full_sh=.true.)                            !
!                                                                      !
!       - a,e,inc,om,bom,amean: initial orbital elements               !
!                                                                      !
!       - time0: initial time                                          !
!       - tot_time: total time of integration given in (see type_time) !
!       - time_step: time step given in (see type_time)                !
!       - type_time: provides unit of measure of tot_time and time_step!
!           (for example "y","d" implies a tot_time given in years and !
!           a time_step given in days)                                 !
!                                                                      !
!       - dyn_vec: 4-dim boolean vector to switch on/off given dynamics!
!           (J2,SRP,Sun,Moon)                                          !
!       - flag_full_sh: only for Earth dynamics, include all available !
!           SH coefficients to model Earth's gravity field             !
!       - flag_alt_j2: include C22 and S22 to the Earth J2 dynamics    !
!           (requires flag_sd=.true. )                                 !
!       - flag_ABM: activate the Adams-Bashfort-Moulton prop. scheme   !
!                                                                      !
!       - are2ms: area-to-mass ratio, given in [m^2/kg]                !
!       - jsrp, ksrp: SRP toy-models references                        !
!       - resang0: initial resonant angle, defined by jsrp and ksrp    !
!       - flag_all_tm: include all the toy models in the dynamics      !
!                                                                      !
!       - flag_sd: select Earth space debris dynamics                  !
!       - flag_ham: select Hamiltonian formulation of the dynamics     !
!           (available only for the space debris case flag_sd=.true.)  !
!       - flag_poly: select polyhedron dynamics                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program plarrt_prop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 	call to use modules
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! use OMP_LIB
    use constants
    ! use geometry
    ! use trinomials
    ! use polygrav ! to compute the RHS of polyhedron dynamics
    use keplerian
    ! use integration
    use prop_lite
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 	Variables declaration
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! strings
    character(len=40) :: filename
    ! numbers
    ! real(kind=16) :: time0, tsid, tout  ! Time variables
    real(kind=16) :: totdays, stepminutes 
    real(kind=16) :: totyears, stepdays
    real(kind=16) :: tot_time, time_step, h !h=stepsize
    real(kind=16) :: sca !, eps	! sca= scale, eps=small number
    real(kind=16) :: nre, a,a_ham, e, inc, om, bom, amean, fan
    real(kind=16) :: period, vol        ! period and volume of the polyhedron
    real(kind=16), dimension(7) :: orb
    real(kind=16), dimension(6) :: x
    real(kind=16), dimension(10,6) :: xx, ff
    integer :: fu, rc   ! integers to be used as auxiliary quantities
    integer :: totsteps, nmax, N, l ! integers for iterations, counting the verts 
    ! and faces of poly, steps of integration and maximum degree for SH, N and l ADD DESCRIPTION 
    ! flags
    logical :: flag_sd,flag_ham, flag_J2, flag_SRP, flag_poly!flags for the different option
    ! The available flags should be: flag_J2, flag_SRP, flag_3rd, flag_poly
    character(len=5), dimension(2) :: type_time
    real(kind=16) :: aux_time
    character(len=5) :: aux_char

     ! Definitions of where to find specific data in the namelist
        namelist /OBJECT/ filename, sca, sigma, period   ! read filename, part of the address as "poly_mod/filename.obj", the period of the rotation, the density, the scale and the maximum degree of SH
        !nre
        namelist /SH/ nmax

        namelist /ORBEL/ a, e, inc, om, bom, amean  ! read the initial orbital elements for a propagation
    
        namelist /TIME/ time0, tot_time, time_step, type_time

        namelist /DYN/ dyn_vec, flag_full_sh, flag_alt_j2, flag_ABM
        
        namelist /SRP/  are2ms, jsrp, ksrp, resang0, flag_all_tm

        namelist /FLAGS_PROP/ flag_sd, flag_ham, flag_poly! read what the program should do
    ! Opening and reading from the namelist
        open (action='read', file='namelist.nml', iostat=rc, newunit=fu)
        read (nml=OBJECT, iostat=rc, unit=fu)
        read (nml=SH, iostat=rc, unit=fu)           
        read (nml=ORBEL, iostat=rc, unit=fu)
        read (nml=TIME, iostat=rc, unit=fu) 
        read (nml=DYN, iostat=rc, unit=fu)
        read (nml=SRP, iostat=rc, unit=fu)
        read (nml=FLAGS_PROP, iostat=rc, unit=fu)
        close(fu)

        orb(1)=a
        orb(2)=e
        orb(3)=inc
        orb(4)=om
        orb(5)=bom
        orb(6)=amean
        orb(7)=fan

        nmax_aux=nmax

        write(*,*) 'Total time of integration:'
        write(*,*) tot_time, type_time(1)
        write(*,*) time_step, type_time(2)
        aux_char=type_time(1)
        call time_converter(tot_time,aux_char,aux_time)
        tot_time=aux_time
        aux_char=type_time(2)
        call time_converter(time_step,aux_char,aux_time)
        time_step=aux_time
        write(*,*) 'Total time of integration in seconds:'
        write(*,*) tot_time, "s"
        write(*,*) time_step, "s"
if (flag_sd) then
    rate=pi2/86400.d0
    if (flag_ham) then
        write(*,*) "Executing Hamiltonian Space Debris Propagation"
        aux_char="y"
        call time_converter_inverse(tot_time,aux_char,totyears)
        aux_char="d"
        call time_converter_inverse(time_step,aux_char,stepdays)
        ! write(*,*) 'Hamiltonian times:'
        ! write(*,*) totyears, "y"
        ! write(*,*) stepdays, "d"
        call prop_sd_ham(orb, totyears,stepdays)
    else
        write(*,*) "Executing Cartesian Space Debris Propagation"
        aux_char="d"
        call time_converter_inverse(tot_time,aux_char,totdays)
        aux_char="m"
        call time_converter_inverse(time_step,aux_char,stepminutes)
        ! write(*,*) 'Cartesian times:'
        ! write(*,*) totdays, "d"
        ! write(*,*) stepminutes, "m"
        call prop_sd_cart(orb, totdays,stepminutes,nmax)
    end if
else
    sigma=sigma*1000.d0
    if (flag_poly) then
        write(*,*) "Executing Polyhedron Asteroid Propagation"
        write(*,*) "Polyhedron model: ", trim(filename),".obj"
        call prop_gen_poly(orb, tot_time,time_step, filename, period, sca)
    else
        write(*,*) "Executing SH Asteroid Propagation"
        call prop_gen_SH(orb, tot_time, time_step, filename, period, nmax) 
    end if
end if 
end program
    