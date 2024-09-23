!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	Program to handle polyhedra, shift to CoM and compute SH coeff.    !
!   To Be Set in namelist.nml:                                         !
!       . /OBJECT/ filename, sca, sigma                                !
!       . /SH/ nmax                                                    !
!       . /FLAGS_PH/ flag_already_shifted, flag_gen, flag_shift_com,   !
!           flag_denorm, flag_no_sigma                                 !
!       Description:                                                   !
!       - filename: name of the asteroid .obj file                     !
!       - sca: scale at which the .obj model is given (1=m,1000=km,ecc.)
!       - sigma: mean density of the body, given in [g/cm^3]           !
!                                                                      !
!       - nmax: maximum degree of spherical harmonics to be used       !
!           (requires  flag_full_sh=.true.)                            !
!                                                                      !
!       - flag_already_shifted: select Earth space debris dynamics     !
!       - flag_gen: select Hamiltonian formulation of the dynamics     !
!           (available only for the space debris case flag_sd=.true.)  !
!       - flag_shift_com: select polyhedron dynamics                   !
!       - flag_denorm: select polyhedron dynamics                      !
!       - flag_no_sigma: select polyhedron dynamics                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program plarrt_poly_handler
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   call to use modules                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use constants
    use SHcomputer
    use geometry
    use trinomials
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	Variables declaration                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! strings
    character(len=30) :: filename
    character(len=1) ::  hash
    character(len=100) :: address, output_address, shifted_address
    real(kind=16), dimension(3) ::  com 	! com = coordinate of com
    real(kind=16) :: sca, h	! sca= scale, h=stepsize, eps=small number, mu=Gm product
    real(kind=16) :: vol		! volume of the polyhedron
    real(kind=16) start_time, end_time	! start and end times
    integer :: nVert, nFaces, totsteps, nmax, nmax1, N, l	! integers for iterations, counting the verts 
    ! and faces of poly, steps of integration and maximum degree for SH, N and l ADD DESCRIPTION 
    logical :: flag_already_shifted, flag_gen, flag_shift_com, flag_denorm, flag_par,flag_no_sigma!flags for the different option
    ! The available flags should be: flag_SH, flag_shift (also flag_denorm)
         integer :: fu, rc	! integers to be used as auxiliary quantities
    ! Definitions of where to find specific data in the namelist
    ! read filename, part of the address as "poly_mod/filename.obj",
    ! the density, the scale and the maximum degree of SH
        namelist /OBJECT/ filename, sca,sigma
        namelist /SH/ nmax
    ! read fags that tell what the program should do
        namelist /FLAGS_PH/ flag_already_shifted, flag_gen, flag_shift_com, flag_denorm, flag_no_sigma
    ! Opening and reading from the namelist
        open (action='read', file='namelist.nml', iostat=rc, newunit=fu)
        read (nml=OBJECT, iostat=rc, unit=fu)
        read (nml=SH, iostat=rc, unit=fu)
        read (nml=FLAGS_PH, iostat=rc, unit=fu)
        close(fu)
        Ma=1.0d0		! This value is set to one for convenience
    ! allocation
        address= "shape_models/"//trim(filename)//'.obj'
        output_address= 'spher_harm/SH_'//trim(filename)//'.txt'
        shifted_address = 'shape_models/'//trim(filename)//'_shift.obj'
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !      EXTRACTION OF VERTICES AND FACES			   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! if the object has already been shifted in a 
        ! previous run change address to the shifted one
        if (flag_already_shifted) then
            address=shifted_address
        else 
            continue
        end if
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
            !	write(*,*) nVert, nFaces
            open (3, file=address, form='formatted')
        end if
        close(3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Allocate matrices of vertices coords and vertices composing a simplex !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        vol=0.0d0
        write(*,*) 'Scale = ', sca, "[m]"
        write(*,*) 'Total number of vertices: ', nVert
        write(*,*) 'Total number of faces: ', nFaces
        write(*,*) 'Reading from: ', address
        write(*,*) 'Extraction of vertices, faces, edges and normals in progress.'
        call extract_ver_and_fac(address,nVert,nFaces,sca,ver,fac,rr) 	! Also computes the reference radius!
        call extract_edges(fac,edg)
        call face_normals(ver,fac,ffnn)
        write(*,*) 'Done.'
        write(*,*) 'Computation of volume in progress.'
        call compute_volume(ver,fac,rr,vol)
        if (flag_no_sigma) then
            sigma=1.d0
            else
           sigma=sigma*1000.d0 ! Originally the density sigma is provided in g/cm**3
            end if  
        write(*,*) 'Done.'
        write(*,*) 'Vol = ', vol, ' [m**3]'
        Ma=vol*sigma
        write(*,*) 'Mass = ', Ma, ' [kg]'
        write(*,*) 'RefRad = ', rr/sca, ' [km]'
    
        if ((flag_shift_com) .and. .not.(flag_already_shifted)) then
            nmax1=1
            write(*,*) 'Computing Centre of Mass.'
            call cpu_time(start_time)
            call compute_C_and_S(ver,fac,sigma,rr,Ma,nmax1,C,S)
            call cpu_time(end_time)
            write(*,*) 'Done in ', end_time-start_time, ' [s].'
            ! The centre of mass is defined in terms of the spherical harmonics as
            ! See Milani & Gronchi
            com(1)=C(2,2)*rr
            com(2)=S(2,2)*rr
            com(3)=C(2,1)*rr
            write(*,*) 'Coordinates of the c.o.m. :', com
            call shiftCoM(shifted_address,com,ver,fac,sca)
            write(*,*) 'Done. Shifted .obj file can be found in ',trim(shifted_address)
            if (flag_gen) then
                write(*,*) 'NEW Extraction of vertices, faces, edges and normals in progress.'
                call extract_ver_and_fac(shifted_address,nVert,nFaces,sca,ver,fac,rr) 	! Also computes the reference radius!
                call extract_edges(fac,edg)
                call face_normals(ver,fac,ffnn)
                write(*,*) 'Done.'
                write(*,*) 'RefRad = ', rr/sca, ' [km]'
            else
                stop
            end if
        else
            continue
        end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !	Normalized Harmonic Coefficients computation   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    100		continue
        if (flag_gen) then
        write(*,*) 'Computing SH coefficients.'    
        call cpu_time(start_time)
        call compute_C_and_S(ver,fac,sigma,rr,Ma,nmax,C,S)
        call cpu_time(end_time)
        write(*,*) 'Done in ', end_time-start_time, ' [s].'
            if (flag_denorm) then
                write(*,*) 'Denormalizing SH coefficients.'
                call denormalize_C_and_S(C,S)
                write(*,*) 'Done.'
            end if
        write(*,*) 'Writing SH coefficients to file.'
        call write_C_and_S_and_rr(C,S,rr,Ma,output_address,flag_denorm,flag_no_sigma)
        write(*,*) 'Done!'
        stop
        else
        end if
    
         end
    