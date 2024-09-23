!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	Program to make various types of comparisons                       !
!   To Be Set in namelist.nml:                                         !
!       . /OBJECT/ filename, sca, sigma, period                        !
!       . /SH/ nmax                                                    !
!       . /ORBEL/ a, e, inc, om, bom, amean                            ! 
!       . /TIME/ time0, tot_time, time_step, type_time                 !
!       . /DYN/ dyn_vec, flag_full_sh, flag_alt_j2, flag_ABM           !
!       . /SRP/  are2ms, jsrp, ksrp, resang0, flag_all_tm              !
!       . /FLAGS_PROP/ flag_sd, flag_ham, flag_poly                    !
!       . /COMPARE/ what                                               !
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
!                                                                      !
!       - what: the type of comparison that we are considering         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program plarrt_compare
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 	call to use modules
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! use OMP_LIB
    use constants
    ! use geometry
    ! use trinomials
    ! use polygrav ! to compute the RHS of polyhedron dynamics
    use keplerian
    ! use integration
    use prop_lite
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 	Variables declaration
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! strings
    character(len=40) :: filename, what,orb_fmt
    ! numbers
    ! real(kind=16) :: time0, tsid, tout  ! Time variables
    real(kind=16) :: totdays, stepminutes 
    real(kind=16) :: totyears, stepdays
    real(kind=16) :: tot_time, time_step, h,tdays !h=stepsize
    real(kind=16) :: sca !, eps	! sca= scale, eps=small number
    real(kind=16) :: nre, a,a_ham, e, inc, om, bom, amean
    real(kind=16) :: a1, e1, inc1, om1, bom1, amean1,a2, e2, inc2, om2, bom2, amean2
    real(kind=16) :: period, vol        ! period and volume of the polyhedron
    character(len=10) :: nmax_char
    real(kind=16), dimension(7) :: orb
    integer :: fu, rc   ! integers to be used as auxiliary quantities
    integer :: totsteps, nmax, N, l ! integers for iterations, counting the verts 
    ! and faces of poly, steps of integration and maximum degree for SH, N and l ADD DESCRIPTION 
    ! flags
    logical :: flag_sd,flag_ham, flag_J2, flag_SRP, flag_poly!flags for the different option
    logical :: abort_poly
    ! The available flags should be: flag_J2, flag_SRP, flag_3rd, flag_poly
    character(len=5), dimension(2) :: type_time
    real(kind=16) :: aux_time
    character(len=5) :: aux_char
    integer :: iostat1, iostat2
    real :: start, finish !aux fort time
     ! Definitions of where to find specific data in the namelist
        namelist /OBJECT/ filename, sca, sigma, period   ! read filename, part of the address as "poly_mod/filename.obj", the period of the rotation, the density, the scale and the maximum degree of SH
        !nre
        namelist /SH/ nmax

        namelist /ORBEL/ a, e, inc, om, bom, amean  ! read the initial orbital elements for a propagation

        namelist /TIME/ time0, tot_time, time_step, type_time

        namelist /DYN/ dyn_vec, flag_full_sh, flag_alt_j2, flag_ABM
        
        namelist /SRP/  are2ms, jsrp, ksrp, resang0, flag_all_tm

        namelist /FLAGS_PROP/ flag_sd, flag_ham, flag_poly! read what the program should do

        namelist /COMPARE/ what, abort_poly
    ! Opening and reading from the namelist
        open (action='read', file='namelist.nml', iostat=rc, newunit=fu)
        read (nml=OBJECT, iostat=rc, unit=fu)
        read (nml=SH, iostat=rc, unit=fu)           
        read (nml=ORBEL, iostat=rc, unit=fu)
        read (nml=TIME, iostat=rc, unit=fu) 
        read (nml=DYN, iostat=rc, unit=fu)
        read (nml=SRP, iostat=rc, unit=fu)
        read (nml=FLAGS_PROP, iostat=rc, unit=fu)
        read (nml=COMPARE, iostat=rc, unit=fu)
        close(fu)

        orb(1)=a
        orb(2)=e
        orb(3)=inc
        orb(4)=om
        orb(5)=bom
        orb(6)=amean
        orb(7)=fan

        jsrp0=jsrp
        ksrp0=ksrp
        nmax_aux=nmax
        write(*,*) 'Total time of integration:'
        write(*,*) tot_time, type_time(1)
        write(*,*) time_step, type_time(2)
        aux_char=type_time(1)
        call time_converter(tot_time,aux_char,aux_time)
        tot_time=aux_time ! Now tot_time is in seconds
        aux_char=type_time(2)
        call time_converter(time_step,aux_char,aux_time)
        time_step=aux_time ! Now time_step is in seconds
        write(*,*) 'Total time of integration in seconds:'
        write(*,*) tot_time, "s"
        write(*,*) 'Time step in seconds:'
        write(*,*) time_step, "s"
        if (what .eq. "sec") then
            ! propagate using Hamiltonian
            write(*,*) "Executing Hamiltonian Space Debris Propagation"
            aux_char="y"
            call time_converter_inverse(tot_time,aux_char,totyears)
            aux_char="d"
            call time_converter_inverse(time_step,aux_char,stepdays)
            call prop_sd_ham(orb, totyears,stepdays)
            jsrp=jsrp0
            ksrp=ksrp0
            ! propagate using Cartesian
            write(*,*)
            write(*,*) "Executing Cartesian Space Debris Propagation"
            aux_char="d"
            call time_converter_inverse(tot_time,aux_char,totdays)
            stepminutes=stepdays!*3.d0
            call prop_sd_cart(orb, totdays,stepminutes,nmax)
        end if

        if (what .eq. "polySH") then
            ! propagate using associated SH coefficients
            call cpu_time(start)
            write(*,*) "Executing SH Asteroid Propagation"
            write(*,*) 
            call prop_gen_SH(orb, tot_time, time_step, filename, period, nmax)
            call cpu_time(finish)
            write(*,*) "Elapsed: ", finish-start, " seconds"
            ! propagate using polyhedron dynamics
            sigma=sigma*1000.d0
            write(*,*) 
            if (.not. abort_poly) then
                call cpu_time(start)
                write(*,*) "Executing Polyhedron Asteroid Propagation"
                write(*,*) "Polyhedron model: ", trim(filename),".obj"
                call prop_gen_poly(orb, tot_time,time_step, filename, period, sca)
                call cpu_time(finish)
                write(*,*) "Elapsed: ", finish-start, " seconds"
            else 
                write(*,*) "Skipping polyhedron propagation"
            end if
            write(*,*)
            write(*,*) "Creating comparison file"
            write(nmax_char,'(I3.3)') nmax
            open(1,file='test_out/'//trim(nmax_char)//'_SH_orb_states_ast.plt', form='formatted')
            open(2,file='test_out/poly_orb_states_ast.plt', form='formatted')
            open(3,file='test_out/polyvSH.plt',  status='unknown',  form='formatted')
            orb_fmt= '((F15.6),6(2X,F15.10))'
            iostat1=0
            iostat2=0
            read(1,*,IOSTAT=iostat1) tdays, a1, e1, inc1, om1, bom1, amean1
            read(2,*,IOSTAT=iostat2) tdays, a2, e2, inc2, om2, bom2, amean2
            do while (iostat1 .eq. 0)
            write(3,orb_fmt) tdays, (a1-a2)*1.d3,e1-e2,inc1-inc2,om1-om2,bom1-bom2,amean1-amean2
            read(1,*,IOSTAT=iostat1) tdays, a1, e1, inc1, om1, bom1, amean1
            read(2,*,IOSTAT=iostat2) tdays, a2, e2, inc2, om2, bom2, amean2
            end do
            close(1)
            close(2)
            close(3)
        end if 
        if (what .eq. "multiSH") then
            nmax_aux=2
            do while (nmax_aux .le. nmax)
                call cpu_time(start)
                write(*,*) "Executing SH Asteroid Propagation - nmax=", nmax_aux
                write(*,*) 
                call prop_gen_SH(orb, tot_time, time_step, filename, period, nmax_aux)
                call cpu_time(finish)
                write(*,*) "Elapsed: ", finish-start, " seconds"
                nmax_aux=nmax_aux+1
            end do
        end if
end program
    