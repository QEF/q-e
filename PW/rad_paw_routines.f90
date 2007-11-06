!
! Copyright (C) 2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Notes for following table:
! 1. 'PAW_' prefix omitted, when present
! 2. '~' means "only for gradient correction".
! 3. calls to 'simpson' (flib/simpsn), 'infomsg' and 'errore' omitted
!
!> potential +
!            |
!            +- rho_lm [2]
!            +- h_potential +
!            |              +- hartree (radial_grids)
!            |
!            +- xc_potential +
!> newd +                    +- lm2rad
!       |                    +- vxc_t (Modules/vxc_t)
!       +- rho_lm [2]        +- exc_t (Modules/exc_t)
!                            +- rad2lm
!> integrate +               +~ gcxc_potential +
!            |                                 +- lm2rad
!            +- rho_lm [2]                     +- gradient +
!                                              |           + radial_gradient
!> ddot +                                      |             (Modules/vxcgc)
!       |                                      |
!       +- rho_lm [2]                          | +- gcxc (Modules/functionals)
!       +- h_potential +                       +{
!                      +- hartree              | +- gcc_spin (Modules/functionals)
!> init +                                      | +- gcx_spin (Modules/functionals)
!       |                                      |
!       + rad_init +                           +- rad2lm [2]
!                  +- weights                  +- divergence +
!                  +- ylmr2                                  +- rad2lm
!                  +~ dylmr2                                 +- radial_gradient
!                                                               (Modules/vxcgc)
!
! NOTE ON PARALLELIZATION:
! this code is parallelized on atoms, this means that each node computes potential,
! energy, newd coefficients, ddots and \int v \times n on a reduced number of atoms.
! The implementation assumes that divisions of atoms among the nodes is done always
! in the same way! Doing so we can avoid to allocate the potential for all the atoms
! on all the nodes, and (most important) we don't need to distribute the potential
! among the nodes after computing it.
!
MODULE rad_paw_routines
    !
    USE kinds,      ONLY : DP
    USe parameters, ONLY : ntypx
    !
    IMPLICIT NONE

    ! entry points:
    PUBLIC :: PAW_potential ! prepare paw potential and store it,
                            ! also computes energy if required
    PUBLIC :: PAW_integrate ! computes \int v(r) \times n(r) dr
    PUBLIC :: PAW_ddot      ! error estimate for mix_rho
    PUBLIC :: PAW_init      ! initialize
    PUBLIC :: PAW_newd      ! computes descreening coefficients
    !
    PRIVATE
    SAVE
    !
    ! Set to true after initialization, to prevent double allocs:
    LOGICAL              :: is_init = .false.

    ! We need a place to store the radial AE and pseudo potential,
    ! as different atoms may have different max_lm, and different max(|r|)
    ! using a derived type is the way to go
    TYPE paw_saved_potential
        REAL(DP),ALLOCATABLE :: &
            v(:,:,:,:)  ! indexes: |r|, lm, spin, {AE|PS}
    END TYPE
    TYPE(paw_saved_potential),ALLOCATABLE :: &
         saved(:) ! allocated in PAW_rad_init

    ! the following variables are used to convert spherical harmonics expansion
    ! to radial sampling, they are initialized for an angular momentum up to
    ! l = l_max and (l+1)**2 = lm_max
    ! see function PAW_rad_init for details
    INTEGER              :: l_max  = 0
    INTEGER              :: lm_max = 0
    INTEGER              :: nx     = 0
    REAL(DP),ALLOCATABLE :: ww(:)
    REAL(DP),ALLOCATABLE :: ylm(:,:) ! Y_lm(nx,lm_max)

    ! additional variables for gradient correction
    INTEGER,PARAMETER    :: xlm = 2     ! Additional angular momentum to
                                        ! integrate to have a good GC
    REAL(DP),ALLOCATABLE :: dylmt(:,:),&! |d(ylm)/dtheta|**2
                            dylmp(:,:)  ! |d(ylm)/dphi|**2
    REAL(DP),ALLOCATABLE :: cos_th(:),& ! cos(theta) (for divergence)
                            sin_th(:)   ! sin(theta) (for divergence)

    ! This type contains some useful data that has to be passed to all
    ! the functions, but cannot stay in global variables for parallel:
    TYPE paw_info
        INTEGER :: a ! atom index
        INTEGER :: t ! atom type index
        INTEGER :: m ! atom mesh = g(nt)%mesh
        INTEGER :: w ! w=1 --> all electron, w=2 --> pseudo
                     ! (used only for gradient correction)
    END TYPE

    ! the following macro controls the use of several fine-grained clocks
    ! set it to '! CALL' (without quotes) in order to disable them
#define OPTIONAL_CALL CALL

 CONTAINS

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! Computes V_h and V_xc using the "density" becsum provided, and stores the
!!! total potential in the global static save%v variable
!!
SUBROUTINE PAW_potential(becsum, energy, e_cmp)
    USE kinds,                  ONLY : DP
    USE atom,                   ONLY : g => rgrid
    USE ions_base,              ONLY : nat, ityp
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, lmaxq
    USE grid_paw_variables,     ONLY : pfunc, ptfunc, tpawp, &
                                       aerho_atc, psrho_atc, aug

    REAL(DP), INTENT(IN)    :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations
    REAL(DP),INTENT(OUT),OPTIONAL :: energy          ! if present compute E[rho]
    REAL(DP),INTENT(OUT),OPTIONAL :: e_cmp(nat,2,2)
    !
    INTEGER, PARAMETER      :: AE = 1, PS = 2,&      ! All-Electron and Pseudo
                               XC = 1, H  = 2        ! XC and Hartree
    REAL(DP), POINTER       :: rho_core(:,:)         ! pointer to AE/PS core charge density 
    TYPE(paw_info)          :: i                     ! minimal info on atoms
    INTEGER                 :: i_what                ! counter on AE and PS
    INTEGER                 :: is                    ! counter on AE and PS
    INTEGER                 :: na,first_nat,last_nat ! atoms counters and indexes
    !
    REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:) ! density expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: v_lm(:,:,:)   ! workspace: potential
    REAL(DP)                :: energy_xc, energy_h, energy_tot
    REAL(DP)                :: sgn

    CALL start_clock('PAW_pot')

    IF(present(energy)) energy_tot = 0._dp

    CALL divide (nat, first_nat, last_nat)
    !
    atoms: DO na = first_nat, last_nat
    !
    i%a = na         ! the index of the atom
    i%t = ityp(na)   ! the type of atom na
    i%m = g(i%t)%mesh! radial mesh size for atom na
    !
    ifpaw: IF (tpawp(i%t)) THEN
        !
        ! Arrays are allocated inside the cycle to allow reduced
        ! memory usage as differnt atoms have different meshes (they
        ! also have different lmax, I will fix this sooner or later)
        ALLOCATE(rho_lm(i%m,lmaxq**2,nspin))
        ALLOCATE(v_lm(i%m,lmaxq**2,nspin))
        !
        whattodo: DO i_what = AE, PS
            i%w = i_what     ! spherical_gradient likes to know
            ! STEP: 1 [ build rho_lm (PAW_rho_lm) ]
            NULLIFY(rho_core)
            IF (i_what == AE) THEN
                ! to pass the atom indes is dirtyer but faster and
                ! uses less memory than to pass only a hyperslice of the array
                CALL PAW_rho_lm(i, becsum, pfunc, rho_lm)
                ! used later for xc potential:
                rho_core => aerho_atc
                ! sign to sum up the enrgy
                sgn = +1._dp
            ELSE
                CALL PAW_rho_lm(i, becsum, ptfunc, rho_lm, aug)
                !    optional argument for pseudo part --> ^^^
                rho_core => psrho_atc ! as before
                sgn = -1._dp          ! as before
            ENDIF

        ! cleanup previously stored potentials
        saved(i%a)%v(:,:,:,i%w) = 0._dp
#ifdef __NO_HARTREE
        write(0,*) "########### skipping hartree paw potential"
#else
        ! First compute the Hartree potential (it does not depend on spin...):
#ifdef __SPHERICAL_TERM_ONLY
        write(0,*) "########### RADIAL GRID TEST: SPHERICALY AVERAGED RHO_LM"
        rho_lm(:,2:lmaxq**2,:) = 0.d0
#endif
        CALL PAW_h_potential(i, rho_lm, v_lm(:,:,1), energy)
        ! using "energy" as the in/out parameter I save a double call, but I have to do this:
        IF (present(energy)) energy_h = energy
        DO is = 1,nspin ! ... so it has to be copied to all spin components
            saved(i%a)%v(:,:,is,i%w) = v_lm(:,:,1)
        ENDDO
#endif

#ifdef __NO_XC
        write(0,*) "########### skipping XC paw potential"
#else
        ! Than the XC one:
        CALL PAW_xc_potential(i, rho_lm, rho_core, v_lm, energy)
        IF (present(energy)) energy_xc = energy
        saved(i%a)%v(:,:,:,i%w) = saved(i%a)%v(:,:,:,i_what) &
                                + v_lm(:,:,:)
#endif
        IF (present(energy)) energy_tot = energy_tot + sgn*(energy_xc + energy_h)
        IF (present(e_cmp)) THEN
            e_cmp(na, 1, i%w) = energy_xc
            e_cmp(na, 2, i%w) = energy_h
        ENDIF
        ENDDO whattodo
    DEALLOCATE(rho_lm, v_lm)
    !
    ENDIF ifpaw
    ENDDO atoms

#ifdef __PARA
    IF ( present(energy) ) &
        CALL reduce (1, energy_tot )
    ! note: potential doesn't need to be recollected,
    ! see note at the top of the file for details.
#endif

    ! put energy back in the output variable
    IF (present(energy)) energy = energy_tot

    CALL stop_clock('PAW_pot')

END SUBROUTINE PAW_potential


!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! Update the descreening coefficients:
!!! D_ij = \int (v_loc_scf + v_loc_at) p_ij
!!
!! This is subroutine does NOT cycle on atoms because newd likes it this way
SUBROUTINE PAW_newd(d_ae, d_ps)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : sqrtpi, e2
    USE lsda_mod,               ONLY : nspin
    USE ions_base,              ONLY : nat, ityp
    USE atom,                   ONLY : g => rgrid
    USE grid_paw_variables,     ONLY : pfunc, ptfunc, tpawp, aug, &
                                       aevloc_at, psvloc_at, ra=>nraug
    USE uspp_param,             ONLY : nh, nhm, lmaxq

    REAL(DP), TARGET, INTENT(INOUT) :: d_ae( nhm, nhm, nat, nspin)
    REAL(DP), TARGET, INTENT(INOUT) :: d_ps( nhm, nhm, nat, nspin)

    INTEGER                 :: i_what
    INTEGER                 :: na, first_nat, last_nat
    ! fake cross band occupations to select only one pfunc at a time:
    REAL(DP)                :: becfake(nhm*(nhm+1)/2,nat,nspin)

    INTEGER, PARAMETER      :: AE = 1, PS = 2  ! All-Electron and Pseudo
    INTEGER                 :: lm,k            ! counters on angmom and radial grid
    INTEGER                 :: nb, mb, nmb, is
    !
    REAL(DP), POINTER       :: v_at(:,:)       ! point to aevloc_at or psvloc_at
    REAL(DP), POINTER       :: d(:,:,:,:)      ! point to d_ae or d_ps
    REAL(DP), ALLOCATABLE   :: pfunc_lm(:,:,:) ! aux charge density
    REAL(DP)                :: integral        ! workspace
    !
    TYPE(paw_info)          :: i

    !
    ! At the moment this part is done using a fake becsum with all the states
    ! empty but one. It's inefficient, but safe.
    CALL start_clock ('PAW_newd')

    ! Some initialization
    becfake(:,:,:) = 0._dp
    !
    d_ae(:,:,:,:)  = 0._dp
    d_ps(:,:,:,:)  = 0._dp

    CALL divide (nat, first_nat, last_nat)

    atoms: DO na = first_nat, last_nat
    !
    i%a = na         ! the index of the atom
    i%t = ityp(na)   ! the type of atom na
    ! skip non-paw atoms
    IF ( .not. tpawp(i%t) ) CYCLE
    i%m = g(i%t)%mesh! radial mesh size for atom na
        !
        ! Different atoms may have different mesh sizes:
        ALLOCATE(pfunc_lm(i%m, lmaxq**2, nspin))
        !
        whattodo: DO i_what = AE, PS
            i%w = i_what
            !
            spins: DO is = 1, nspin
                nmb = 0
                ! loop on all pfunc for this kind of pseudo
                DO nb = 1, nh(i%t)
                DO mb = nb, nh(i%t)
                nmb = nmb+1 ! nmb = 1, nh*(nh+1)/2
                    !
                    ! compute the density from a single pfunc
                    becfake(nmb,na,is) = 1._dp
                    IF (i%w == AE) THEN
                        CALL PAW_rho_lm(i, becfake, pfunc, pfunc_lm)
                        v_at => aevloc_at
                        d    => d_ae
                    ELSE
                        CALL PAW_rho_lm(i, becfake, ptfunc, pfunc_lm, aug)
                        v_at => psvloc_at
                        d    => d_ps
                    ENDIF
                    !
                    ! Now I multiply the pfunc and the potential, I can use
                    ! pfunc_lm itself as workspace
                    DO lm = 1,lmaxq**2
#ifdef __NO_LOCAL
                    pfunc_lm(1:i%m,lm,is) = pfunc_lm(1:i%m,lm,is) * &
                            saved(i%a)%v(1:i%m,lm,is,i%w)
#else
                        IF ( lm == 1 ) THEN
                            pfunc_lm(1:i%m,lm,is) = pfunc_lm(1:i%m,lm,is) * &
                                ( saved(i%a)%v(1:i%m,lm,is,i%w) + e2*sqrtpi*v_at(1:i%m,i%t) )
                        ELSE
                            pfunc_lm(1:i%m,lm,is) = pfunc_lm(1:i%m,lm,is) * &
                                saved(i%a)%v(1:i%m,lm,is,i%w)
                        ENDIF
#endif
                        !
                        ! Integrate!
                        CALL simpson (ra(i%t),pfunc_lm(:,lm,is),g(i%t)%rab,integral)
                        d(nb,mb,i%a,is) = d(nb,mb,i%a,is) + integral
                    ENDDO
                    ! Symmetrize:
                    d(mb,nb, i%a,is) = d(nb,mb, i%a,is)
                    !
                    becfake(nmb,na,is) = 0._dp
                ENDDO ! mb
                ENDDO ! nb
            ENDDO spins
!#ifdef __PARA
!            CALL reduce (nhm*nhm*nat*nspin, d )
!#endif
        ENDDO whattodo

        ! cleanup
        DEALLOCATE(pfunc_lm)
        !
    ENDDO atoms

#ifdef __PARA
    CALL reduce (nhm*nhm*nat*nspin, d_ae )
    CALL reduce (nhm*nhm*nat*nspin, d_ps )
#endif

    CALL stop_clock ('PAW_newd')


END SUBROUTINE PAW_newd


!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! Compute a new density from becsum, and integrate it with the potential
!!! previously computed and stored in global static variable "saved(na)"
!!
FUNCTION PAW_integrate(becsum)
    USE kinds,                  ONLY : DP
    USE lsda_mod,               ONLY : nspin
    USE ions_base,              ONLY : nat, ityp
    USE atom,                   ONLY : g => rgrid
    USE grid_paw_variables,     ONLY : pfunc, ptfunc, tpawp, aug
    USE uspp_param,             ONLY : nhm, lmaxq

    REAL(DP)                :: PAW_integrate

    REAL(DP), INTENT(IN)    :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations

    INTEGER, PARAMETER      :: AE = 1, PS = 2        ! All-Electron and Pseudo
    INTEGER                 :: na,first_nat,last_nat ! atoms counters and indexes
    INTEGER                 :: lm,k                  ! counters on angmom and radial grid
    INTEGER                 :: is                    ! counter on spin
    !
    REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:) ! radial density expanded on Y_lm
    REAL(DP)                :: integral      ! workspace
    !
    TYPE(paw_info)          :: i
    INTEGER                 :: i_what,i_sign ! =1 => AE; =2 => PS

    CALL start_clock ('PAW_int')

    PAW_integrate = 0._dp

    CALL divide (nat, first_nat, last_nat)

    atoms: DO na = first_nat, last_nat
    !
    i%a = na         ! the index of the atom
    i%t = ityp(na)   ! the type of atom na
    i%m = g(i%t)%mesh! radial mesh size for atom na
    !
    ifpaw: IF (tpawp(i%t)) THEN
        !
        ALLOCATE(rho_lm(i%m,lmaxq**2,nspin))
        !
        whattodo: DO i_what = AE, PS
            i%w = i_what
            !
            IF (i%w == AE) THEN
                CALL PAW_rho_lm(i, becsum, pfunc, rho_lm)
                i_sign = +1._dp
            ELSE
                CALL PAW_rho_lm(i, becsum, ptfunc, rho_lm, aug)
                i_sign = -1._dp
            ENDIF
            !
            ! Compute the integral
            DO is = 1,nspin
            DO lm = 1, lmaxq**2
                ! I can use rho_lm itself as workspace
                rho_lm(1:i%m,lm,is) = rho_lm(1:i%m,lm,is) &
                                     * saved(i%a)%v(1:i%m,lm,is,i%w)
                CALL simpson (i%m,rho_lm(:,lm,is),g(i%t)%rab,integral)
                !WRITE(6,"(3i3,10f20.10)") i_what, is, lm, integral
                PAW_integrate = PAW_integrate + i_sign*integral
                !
            ENDDO
            ENDDO
        ENDDO whattodo
        !
        DEALLOCATE(rho_lm)
        !
    ENDIF ifpaw
    ENDDO atoms

#ifdef __PARA
    CALL reduce (1, PAW_integrate )
#endif

    CALL stop_clock ('PAW_int')


END FUNCTION PAW_integrate

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! As rho_ddot in mix_rho for radial grids
!!
FUNCTION PAW_ddot(bec1,bec2)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : pi
    USE lsda_mod,               ONLY : nspin
    USE ions_base,              ONLY : nat, ityp
    USE atom,                   ONLY : g => rgrid

    USE grid_paw_variables,     ONLY : pfunc, ptfunc, tpawp, aug
    USE uspp_param,             ONLY : nhm, lmaxq

    REAL(DP)                :: PAW_ddot

    REAL(DP), INTENT(IN) :: &
             bec1(nhm*(nhm+1)/2,nat,nspin), &! cross band occupations (previous step)
             bec2(nhm*(nhm+1)/2,nat,nspin)   ! cross band occupations (next step)

    INTEGER, PARAMETER      :: AE = 1, PS = 2        ! All-Electron and Pseudo
    INTEGER                 :: i_what                ! counter on AE and PS
    INTEGER                 :: na,first_nat,last_nat ! atoms counters and indexes
    INTEGER                 :: lm,k                  ! counters on angmom and radial grid

    ! hartree energy scalar fields expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:) ! radial density expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: v_lm(:,:)     ! hartree potential, summed on spins (from bec1)
    !
    REAL(DP)                :: i_sign        ! +1 for AE, -1 for PS
    REAL(DP)                :: integral      ! workspace
    TYPE(paw_info)          :: i

    CALL start_clock ('PAW_ddot')

    ! initialize for integration on angular momentum and gradient
    PAW_ddot = 0._dp

! !$OMP PARALLEL DO default(shared) reduction(+:paw_ddot) &
! !$OMP private(rho_lm, v_lm, i_sign, i_what, na)

    CALL divide (nat, first_nat, last_nat)
    !
    atoms: DO na = first_nat, last_nat
    !
    i%a = na         ! the index of the atom
    i%t = ityp(na)   ! the type of atom na
    i%m = g(i%t)%mesh ! radial mesh size for atom na
    !
    ifpaw: IF (tpawp(i%t)) THEN
        !
        ALLOCATE(rho_lm(i%m,lmaxq**2,nspin))
        ALLOCATE(v_lm(i%m,lmaxq**2))
        !
        whattodo: DO i_what = AE, PS
            ! Build rho from the occupations in bec1
            IF (i_what == AE) THEN
                CALL PAW_rho_lm(i, bec1, pfunc, rho_lm)
                i_sign = +1._dp
            ELSE
                CALL PAW_rho_lm(i, bec1, ptfunc, rho_lm, aug) !fun)
                i_sign = -1._dp
            ENDIF
            !
            ! Compute the hartree potential from bec1
            CALL PAW_h_potential(i, rho_lm, v_lm)
            !
            ! Now a new rho is computed, this time from bec2
            IF (i_what == AE) THEN
                CALL PAW_rho_lm(i, bec2, pfunc, rho_lm)
            ELSE
                CALL PAW_rho_lm(i, bec2, ptfunc, rho_lm, aug) !fun)
            ENDIF
            !
            ! Finally compute the integral
            DO lm = 1, lmaxq**2
                ! I can use v_lm as workspace
                DO k = 1, i%m
                    v_lm(k,lm) = v_lm(k,lm) * SUM(rho_lm(k,lm,1:nspin))
                ENDDO
                CALL simpson (i%m,v_lm(:,lm),g(i%t)%rab,integral)
                !
                ! Sum all the energies in PAW_ddot
                PAW_ddot = PAW_ddot + i_sign * integral
                !
            ENDDO
        ENDDO whattodo
        !
        DEALLOCATE(v_lm)
        DEALLOCATE(rho_lm)
        !
    ENDIF ifpaw
    ENDDO atoms
!!$OMP END PARALLEL DO
#ifdef __PARA
    CALL reduce (1, PAW_ddot )
#endif

    CALL stop_clock ('PAW_ddot')


END FUNCTION PAW_ddot


SUBROUTINE PAW_init()
    USE ions_base,              ONLY : nat, ityp
    USE grid_paw_variables,     ONLY : tpawp
    USE atom,                   ONLY : g => rgrid
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE funct,                  ONLY : dft_is_gradient

    INTEGER :: na, nt, first_nat, last_nat

    ! First a bit of generic initialization:
    ALLOCATE(saved(nat)) ! allocate space to store the potentials
    !
    ! Parallelizing this loop every node only allocs the potential
    ! for the atoms that it will actually use later.
    CALL divide (nat, first_nat, last_nat)
    DO na = first_nat, last_nat
        nt = ityp(na)
        ! note that if the atom is not paw it is left unallocated
        IF ( tpawp(nt) ) THEN
            IF (allocated(saved(na)%v)) DEALLOCATE(saved(na)%v)
            ALLOCATE( saved(na)%v(g(nt)%mesh, lmaxq**2, nspin, 2 ) )
                      !                                     {AE|PS}
        ENDIF
    ENDDO

    ! initialize for integration on angular momentum and gradient, integrating
    ! up to 2*lmaxq (twice the maximum angular momentum of rho) is enough for
    ! H energy and for XC energy. If I have gradient correction I have to go a bit higher
    IF ( dft_is_gradient() ) THEN
        CALL PAW_rad_init(2*lmaxq+xlm)
    ELSE
        CALL PAW_rad_init(2*lmaxq)
    ENDIF

END SUBROUTINE PAW_init

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! initialize several quantities related to radial integration: spherical harmonics and their 
!!! gradients along a few (depending on lmaxq) directions, weights for spherical integration
!!
SUBROUTINE PAW_rad_init(l)
    USE constants,              ONLY : pi, fpi, eps8
    USE funct,                  ONLY : dft_is_gradient
    INTEGER,INTENT(IN)          :: l ! max angular momentum component that will be
                                     ! integrated exactly (to numerical precision)

    REAL(DP),ALLOCATABLE        :: x(:),&       ! nx versors in smart directions
                                   w(:),&       ! temporary integration weights
                                   r(:,:),&     ! integration directions
                                   r2(:),&      ! square modulus of r
                                   ath(:),aph(:)! angles in sph coords for r

    INTEGER                     :: i,ii,n       ! counters
    INTEGER                     :: lm,lm2,m     ! indexes for ang.mom
    REAL(DP)                    :: phi,dphi,rho ! spherical coordinates
    REAL(DP)                    :: z            ! cartesian coordinates
    ! for gradient corrections:
    INTEGER                     :: ipol
    REAL(DP),ALLOCATABLE        :: aux(:,:),&   ! workspace
                                   s(:,:),&     ! integration directions + delta
                                   s2(:)        ! square modulus of s
    REAL(DP)                    :: vth(3), vph(3) !versors for theta and phi
    !
    CHARACTER(len=100)          :: message

    ! reinit if necessary
    IF( is_init ) THEN
        IF ( l /= l_max ) THEN
            CALL infomsg('PAW_rad_init',&
             'PAW radial integration already initialized but for a different l: reinitializing.')
            DEALLOCATE(ww, ylm)
            IF (ALLOCATEd(dylmt))  DEALLOCATE(dylmt)
            IF (ALLOCATEd(dylmp))  DEALLOCATE(dylmp)
            IF (ALLOCATEd(cos_th)) DEALLOCATE(cos_th)
            IF (ALLOCATEd(sin_th)) DEALLOCATE(sin_th)
        ELSE
            ! if already initialized correctly nothing to be done
            RETURN
        ENDIF
    ENDIF

    OPTIONAL_CALL start_clock ('PAW_rad_init')

    ! maximum value of l correctly integrated
    l_max = l
    ! volume element for angle phi
    dphi = 2._dp*pi/(l_max+1)
    ! number of samples for theta angle
    n = (l_max+2)/2
    ALLOCATE(x(n),w(n))
    ! compute weights for theta integration
    CALL weights(x,w,n)

    ! number of integration directions
    nx = n*(l_max+1)
    WRITE(message,"(a,i3,a,i2)") "Setup to integrate on ",nx," directions; integration exact up to l = ",l
    CALL infomsg('PAW_rad_init', message)
    ALLOCATE(r(3,nx),r2(nx), ww(nx), ath(nx), aph(nx))

    ! compute real weights multiplying theta and phi weights
    ii = 0
    do i=1,n
        z = x(i)
        rho=sqrt(1._dp-z**2)
        do m=0,l_max
            ii= ii+1
            phi = dphi*m
            r(1,ii) = rho*cos(phi)
            r(2,ii) = rho*sin(phi)
            r(3,ii) = z
            ww(ii) = w(i)*2._dp*pi/(l_max+1)
            r2(ii) = r(1,ii)**2+r(2,ii)**2+r(3,ii)**2
            ! these will be used later:
            ath(ii) = acos(z/sqrt(r2(ii)))
            aph(ii) = phi
        end do
    end do
    ! cleanup
    DEALLOCATE (x,w)

    ! initialize spherical harmonics that will be used
    ! to convert rho_lm to radial grid
    lm_max = (l_max+1)**2
    ALLOCATE(ylm(nx,lm_max))
    CALL ylmr2(lm_max, nx, r,r2,ylm)

    ! if gradient corrections will be used than we need
    ! to initialize the gradient of ylm, as we are working in spherical
    ! coordinates the formula involves \hat{theta} and \hat{phi}
    gradient: IF (dft_is_gradient()) THEN
        ALLOCATE(s(3,nx),s2(nx))
        ALLOCATE(dylmt(nx,lm_max),dylmp(nx,lm_max),aux(nx,lm_max))
        ALLOCATE(cos_th(nx), sin_th(nx))
        dylmt(:,:) = 0._dp
        dylmp(:,:) = 0._dp

        ! compute derivative along x, y and z => gradient, then compute the
        ! scalar products with \hat{theta} and \hat{phi} and store them in
        ! dylmt and dylmp respectively
        DO ipol = 1,3 !x,y,z
            CALL dylmr2(lm_max, nx, r,r2, aux, ipol)
            DO lm = 1, lm_max
            DO i = 1,nx
                vph = (/-sin(aph(i)), cos(aph(i)), 0._dp/)
                ! this is the explicit form, but the cross product trick (below) is much faster:
                ! vth = (/cos(aph(i))*cos(ath(i)), sin(aph(i))*cos(ath(i)), -sin(ath(i))/)
                vth = (/vph(2)*r(3,i)-vph(3)*r(2,i),&
                        vph(3)*r(1,i)-vph(1)*r(3,i),&
                        vph(1)*r(2,i)-vph(2)*r(1,i)/)
                dylmt(i,lm) = dylmt(i,lm) + aux(i,lm)*vth(ipol)
                ! CHECK: the 1/sin(th) factor should be correct, but deals wrong result, why?
                dylmp(i,lm) = dylmp(i,lm) + aux(i,lm)*vph(ipol) !/sin(ath(i))
                cos_th(i) = cos(ath(i))
                sin_th(i) = sin(ath(i))
            ENDDO
            ENDDO
        ENDDO
        DEALLOCATE(aux)
    ENDIF gradient
    ! cleanup
    DEALLOCATE (r,r2)

    ! success
    is_init = .true.

    OPTIONAL_CALL stop_clock ('PAW_rad_init')

 CONTAINS
    ! Computes weights for gaussian integrals,
    ! from numerical recipes
    SUBROUTINE weights(x,w,n)
    implicit none
    integer :: n, i,j,m
    real(8), parameter :: eps=1.d-14
    real(8) :: x(n),w(n), z,z1, p1,p2,p3,pp,pi
    
    pi = 4._dp*atan(1._dp)
    m=(n+1)/2
    do i=1,m
        z1 = 2._dp
        z=cos(pi*(i-0.25_dp)/(n+0.5_dp))
        do while (abs(z-z1).gt.eps)
        p1=1._dp
        p2=0._dp
        do j=1,n
            p3=p2
            p2=p1
            p1=((2._dp*j-1._dp)*z*p2-(j-1._dp)*p3)/j
        end do
        pp = n*(z*p1-p2)/(z*z-1._dp)
        z1=z
        z=z1-p1/pp
        end do
        x(i) = -z
        x(n+1-i) = z
        w(i) = 2._dp/((1._dp-z*z)*pp*pp)
        w(n+1-i) = w(i)
    end do

    END SUBROUTINE weights
END SUBROUTINE PAW_rad_init 



!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute xc potential and energy, as
!!! xc functional is not diagonal on angular momentum numerical integration is performed
SUBROUTINE PAW_xc_potential(i, rho_lm, rho_core, v_lm, energy)
    USE kinds,                  ONLY : DP
    USE parameters,             ONLY : ntypx
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : lmaxq
    USE ions_base,              ONLY : ityp
    USE radial_grids,           ONLY : ndmx
    USE atom,                   ONLY : g => rgrid
    USE funct,                  ONLY : dft_is_gradient
    USE constants,              ONLY : fpi ! REMOVE

    TYPE(paw_info)  :: i                               ! atom's minimal info
    REAL(DP), INTENT(IN)  :: rho_lm(i%m,lmaxq**2,nspin)! charge density as lm components
    REAL(DP), INTENT(IN)  :: rho_core(ndmx,ntypx)      ! core charge, radial and spherical
    REAL(DP), INTENT(OUT) :: v_lm(i%m,lmaxq**2,nspin)  ! potential density as lm components
    REAL(DP),OPTIONAL,INTENT(OUT) :: energy            ! XC energy (if required)
    !
    REAL(DP)              :: rho_loc(2)         ! local density (workspace), up and down
    REAL(DP)              :: v_loc(2)           ! local density (workspace), up and down
    REAL(DP)              :: v_rad(i%m,nx,nspin)! radial potential (to be integrated)
    REAL(DP)              :: rho_rad(i%m,nspin) ! workspace (only one radial slice of rho)
    !
    REAL(DP), ALLOCATABLE :: e_rad(:)           ! aux, used to store radial slices of energy
    REAL(DP)              :: e                  ! aux, used to integrate energy
    !
    INTEGER               :: ix,k               ! counters on directions and radial grid
    INTEGER               :: lm                 ! counter on angular momentum
    INTEGER               :: lsd                ! switch for local spin density

    REAL(DP), EXTERNAL    :: exc_t              ! computes XC energy

    REAL(DP) :: &
         vgc(ndmx,2),   & ! exchange-correlation potential (GGA only)
         egc(ndmx)        ! exchange correlation energy density (GGA only)

!#define __SURRENDER_OR_DIE
#ifdef __SURRENDER_OR_DIE
    ! REMOVE:
    REAL(DP) :: aux(i%m),aux2(i%m),auxc(i%m) ! charge density as lm components
#endif

    OPTIONAL_CALL start_clock ('PAW_xc_pot')
    !
    ! true if using spin
    lsd = nspin-1
    ! This will hold the "true" charge density, without r**2 or other factors
    rho_loc(:) = 0._dp
    !
    ! ALLOCATE(rho_rad(i%m,nspin))
    IF (present(energy)) THEN
        energy = 0._dp
        ALLOCATE(e_rad(i%m))
    ENDIF
    !
    DO ix = 1, nx
        ! *** LDA (and LSDA) part (no gradient correction) ***
        ! convert _lm density to real density along ix
        CALL PAW_lm2rad(i, ix, rho_lm, rho_rad)
        !
        ! compute the potential along ix
        DO k = 1,i%m
            rho_loc(1:nspin) = rho_rad(k,1:nspin)*g(i%t)%rm2(k)
            CALL vxc_t(rho_loc, rho_core(k,i%t), lsd, v_loc)
            v_rad(k,ix,1:nspin) = v_loc(1:nspin)
            IF (present(energy)) &
                e_rad(k) = exc_t(rho_loc, rho_core(k,i%t), lsd) &
                          * ( SUM(rho_rad(k,1:nspin)) + rho_core(k,i%t)*g(i%t)%r2(k) )
        ENDDO
        IF (present(energy)) THEN
            CALL simpson(i%m, e_rad, g(i%t)%rab, e)
            energy = energy + e * ww(ix)
        ENDIF
    ENDDO
    IF(present(energy)) DEALLOCATE(e_rad)

#ifdef __ONLY_GCXC
    write(0,*) "########### only GC for XC part"
    v_rad(:,:,:) = 0._dp
#endif

    ! Recompose the sph. harm. expansion
    CALL PAW_rad2lm(i, v_rad, v_lm, lmaxq)

#ifdef __SURRENDER_OR_DIE
    write(*,*) "__SURRENDER_OR_DIE"

    OPEN (6354, FILE='ld1.extracted', FORM='FORMATTED' )
    DO k = 1,i%m
        read(6354,'(4f20.10)') aux(k),aux2(k), auxc(k)
    ENDDO
    CLOSE(6354)
#endif

#ifdef __NO_GCXC
#else
!
!  Passing this charge : fpi*rho_lm(:,1,1) and this rho core : fpi*rho_core(:,1)*g(i%t)%r2(:) 
!  to the original vxcgc routine yeld almost the same ddd coefficients as using the charge loaded
!  from file, produced by the atomic code.
!  What the heck is happening in between????
!

    ! Add gradient correction, if necessary
    IF( dft_is_gradient() ) THEN
!           CALL vxcgc(ndmx,g(i%t)%mesh,nspin,g(i%t)%r,g(i%t)%r2,&
!                       sqrt(fpi)*rho_lm(:,1,1),fpi*rho_core(:,1)*g(i%t)%r2(:),vgc,egc,1)

!         v_lm(:,1,1) = 0._dp
!         rho_lm(1:i%m,1,1) = aux2(1:i%m)/sqrt(fpi)
        !v_lm(:,1,1) = 0._dp
        CALL PAW_gcxc_potential(i, rho_lm,rho_core, v_lm,energy)
!      v_lm(1:i%m,1,1:nspin) =  v_lm(1:i%m,1,1:nspin) + vgc(1:i%m,1:nspin)
    ENDIF

!     IF (i%w == 1) THEN
! #ifdef __SURRENDER_OR_DIE
!     write(12000,'(4f20.10)') (g(i%t)%r(k),rho_lm(k,1,1),aux2(k),v_lm(k,1,1),k=1,i%m)
! #else
!     write(12000,'(3f20.10)') (g(i%t)%r(k),rho_lm(k,1,1),v_lm(k,1,1),k=1,i%m)
! #endif
!     ENDIF
#endif

    ! cleanup
    !DEALLOCATE(rho_rad)

    OPTIONAL_CALL stop_clock ('PAW_xc_pot')

END SUBROUTINE PAW_xc_potential
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! add gradient correction to v_xc, code mostly adapted from ../atomic/vxcgc.f90
!!! in order to support non-spherical charges (as Y_lm expansion)
!!! Note that the first derivative in vxcgc becames a gradient, while the second is a divergence.
!!! We also have to temporary store some additional Y_lm components in order not to loose
!!! precision during teh calculation, even if only the ones up to lmaxq (the maximum in the
!!! density of charge) matter when computing \int v * rho 
SUBROUTINE PAW_gcxc_potential(i, rho_lm,rho_core, v_lm, energy)
    USE kinds,                  ONLY : DP
    USE ions_base,              ONLY : ityp
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    USE uspp_param,             ONLY : lmaxq
    USE radial_grids,           ONLY : ndmx
    USE parameters,             ONLY : ntypx
    USE constants,              ONLY : sqrtpi, fpi,pi,e2, eps => eps12, eps2 => eps24
    USE funct,                  ONLY : gcxc, gcx_spin, gcc_spin
    !
    TYPE(paw_info)  :: i                                  ! atom's minimal info
    REAL(DP), INTENT(IN)    :: rho_lm(i%m,lmaxq**2,nspin) ! charge density as lm components
    REAL(DP), INTENT(IN)    :: rho_core(ndmx,ntypx)       ! core charge, radial and spherical
    REAL(DP), INTENT(INOUT) :: v_lm(i%m,lmaxq**2,nspin)   ! potential to be updated
    REAL(DP),OPTIONAL,INTENT(INOUT) :: energy             ! if present, add GC to energy

    REAL(DP)                :: rho_rad(i%m,nx,nspin)      ! charge density sampled
    REAL(DP)                :: grad(i%m,3,nx,nspin)       ! gradient
    REAL(DP)                :: grad2(i%m,nx,nspin)        ! square modulus of gradient
                                                          ! (first of charge, than of hamiltonian)
    REAL(DP)                :: gc_rad(i%m,nx,nspin)       ! GC correction to V
    REAL(DP)                :: gc_lm(i%m,lmaxq**2,nspin)  ! GC correction to V
    REAL(DP)                :: h_rad(i%m,3,nx,nspin)      ! hamiltonian (vector field)
    REAL(DP)                :: h_lm(i%m,3,(lmaxq+xlm)**2,nspin)! hamiltonian (vector field)
                                                             !!! expanded to higher lm than rho !!!
    REAL(DP)                :: div_h(i%m,lmaxq**2,nspin)  ! div(hamiltonian)

    REAL(DP),ALLOCATABLE    :: e_rad(:)                   ! aux, used to store energy
    REAL(DP)                :: e                          ! aux, used to integrate energy

    INTEGER  :: k, ix, is, a, lm                          ! counters on spin and mesh
    REAL(DP) :: sx,sc,v1x,v2x,v1c,v2c                     ! workspace
    REAL(DP) :: v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw  ! workspace
    REAL(DP) :: sgn, arho                                 ! workspace
    REAL(DP) :: rup, rdw, gup, gdw, co2                   ! workspace
    REAL(DP) :: rh, zeta, grh2, grho2(2)

    OPTIONAL_CALL start_clock ('PAW_gcxc_v')

    IF (present(energy)) ALLOCATE(e_rad(i%m))

    ! Compute the gradient of rho along each direction
    DO ix = 1,nx
        ! FIXME: To pass non-consecutive slices is bad (not so bad if nspin == 2)
        !
        ! Here I need rho_rad for all the nx directions, not only one at a time:
        CALL PAW_lm2rad(i, ix, rho_lm, rho_rad(:,ix,:))
        CALL PAW_gradient(i, ix, rho_lm, rho_rad(:,ix,:), rho_core, &
                          grad2(:,ix,:), grad(:,:,ix,:))
    ENDDO

    spin:&
    IF ( nspin == 1 ) THEN
        !
        !     GGA case
        !
        DO ix = 1,nx
        DO k = 1, i%m
            ! arho is the absolute value of real charge, sgn is its sign
            arho = rho_rad(k,ix,1)*g(i%t)%rm2(k) + rho_core(k,i%t)
            sgn  = SIGN(1._dp,arho)
            arho = ABS(arho)

            ! I am using grad(rho)**2 here, so its eps has to be eps**2
            IF ( (arho>eps) .and. (grad2(k,ix,1)>eps2) ) THEN
                CALL gcxc(arho,grad2(k,ix,1), sx,sc,v1x,v2x,v1c,v2c)
                IF (present(energy)) &
                    e_rad(k)    = sgn *e2* (sx+sc) * g(i%t)%r2(k)
                gc_rad(k,ix,1)  = (v1x+v1c)*g(i%t)%rm2(k)
                h_rad(k,:,ix,1) = (v2x+v2c)*grad(k,:,ix,1)*g(i%t)%r2(k)
            ELSE
                IF (present(energy)) &
                    e_rad(k)    = 0._dp
                gc_rad(k,ix,1)  = 0._dp
                h_rad(k,:,ix,1) = 0._dp
            ENDIF
        ENDDO
        ! integrate energy (if required)
        IF (present(energy)) THEN
            CALL simpson(i%m, e_rad, g(i%t)%rab, e)
            energy = energy + e * ww(ix)
        ENDIF
        ENDDO
    ELSEIF ( nspin == 2 ) THEN
        !
        !   this is the \sigma-GGA case
        DO ix = 1,nx
        DO k = 1,i%m
            !
            ! Prepare the necessary quantities
            ! rho_core is considered half spin up and half spin down:
            co2 = rho_core(k,i%t)/2._dp
            ! than I build the real charge dividing by r**2
            rup = rho_rad(k,ix,1)*g(i%t)%rm2(k) + co2
            rdw = rho_rad(k,ix,2)*g(i%t)%rm2(k) + co2
            ! bang!
            CALL gcx_spin (rup, rdw, grad2(k,ix,1), grad2(k,ix,2), &
                sx, v1xup, v1xdw, v2xup, v2xdw)

            rh = rup + rdw ! total charge
            IF ( rh > eps ) THEN
                zeta = (rup - rdw ) / rh
                !grh2 = (grho (i, 1) + grho (i, 2) ) **2 
                grh2 =  (grad(k,1,ix,1) + grad(k,1,ix,2))**2 &
                      + (grad(k,2,ix,1) + grad(k,2,ix,2))**2 &
                      + (grad(k,3,ix,1) + grad(k,3,ix,2))**2
                CALL gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
            ELSE
                sc    = 0._dp
                v1cup = 0._dp
                v1cdw = 0._dp
                v2c   = 0._dp
            ENDIF
            IF (present(energy)) &
                e_rad(k)    = e2*(sx+sc)* g(i%t)%r2(k)
            gc_rad(k,ix,1)  = (v1xup+v1cup)*g(i%t)%rm2(k)
            gc_rad(k,ix,2)  = (v1xdw+v1cdw)*g(i%t)%rm2(k)
            !
            h_rad(k,:,ix,1) =( (v2xup+v2c)*grad(k,:,ix,1)+v2c*grad(k,:,ix,2) )*g(i%t)%r2(k)
            h_rad(k,:,ix,2) =( (v2xdw+v2c)*grad(k,:,ix,2)+v2c*grad(k,:,ix,1) )*g(i%t)%r2(k)
        ENDDO ! k
        ! integrate energy (if required)
        ! NOTE: this integration is duplicated for every spin, FIXME!
        IF (present(energy)) THEN
            CALL simpson(i%m, e_rad, g(i%t)%rab, e)
            energy = energy + e * ww(ix)
        ENDIF
        ENDDO ! ix
    ELSEIF ( nspin == 4 ) THEN
        CALL errore('PAW_gcxc_v', 'non-collinear not yet implemented!', -1)
    ELSE spin
        CALL errore('PAW_gcxc_v', 'unknown spin number', -2)
    ENDIF spin
    !
    ! convert the first part of the GC correction back to spherical harmonics
    CALL PAW_rad2lm(i, gc_rad, gc_lm(:,:,:), lmaxq)
    !
    ! We need the gradient of h to calculate the last part of the exchange
    ! and correlation potential. First we have to convert H to its Y_lm expansion
    DO a = 1,3
        CALL PAW_rad2lm(i, h_rad(:,a,:,:), h_lm(:,a,:,:), lmaxq+xlm)
    ENDDO
    !
    ! Compute div(H)
    CALL PAW_divergence(i, h_lm, div_h, lmaxq+xlm, lmaxq)
    !                         input max lm --^     ^-- output max lm

#ifdef __PAW_ONLYHAM
    write(0,*) "##### only HAM"
#else
#ifdef __PAW_NONHAM
    write(0,*) "##### NO HAM"
#else
#endif
#endif

    ! Finally sum it back into v_xc
    DO is = 1,nspin
    DO lm = 1,lmaxq**2
#ifdef __PAW_ONLYHAM
            v_lm(:,lm,is) = v_lm(:,lm,is) - e2*div_h(:,lm,is)*g(i%t)%rm2(:)
#else
#ifdef __PAW_NONHAM
            v_lm(:,lm,is) = v_lm(:,lm,is) + e2*gc_lm(:,lm,is) 
#else
            v_lm(1:i%m,lm,is) = v_lm(1:i%m,lm,is) + e2*(gc_lm(1:i%m,lm,is)*g(i%t)%r2(1:i%m)-div_h(1:i%m,lm,is))
#endif
#endif
    ENDDO
    ENDDO

!   DO k = 1,i%m
!      write(9000,'(100f20.10)') g(i%t)%r(k), rho_lm(k,1,1), rho_core(k,i%t),&
!                                gc_lm(k,1,1),div_h(k,1,1),e2*(gc_lm(k,1,1)*g(i%t)%r2(k)-div_h(k,1,1))
!   ENDDO


    OPTIONAL_CALL stop_clock ('PAW_gcxc_v')

END SUBROUTINE PAW_gcxc_potential
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! compute divergence of a vector field (actutally the hamiltonian)
!!! it is assumed that: 1. the input function is multiplied by r**2; 
!!! 2. the output function will be divided by r**2 outside
SUBROUTINE PAW_divergence(i, F_lm, div_F_lm, lmaxq_in, lmaxq_out)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : sqrtpi, fpi, eps12, e2
    !USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE ions_base,              ONLY : ityp
    USE parameters,             ONLY : ntypx
    USE atom,                   ONLY : g => rgrid

    TYPE(paw_info)  :: i              ! atom's minimal info
    INTEGER, INTENT(IN)  :: lmaxq_in  ! max angular momentum to derive
                                      ! (divergence is reliable up to lmaxq_in-2)
    INTEGER, INTENT(IN)  :: lmaxq_out ! max angular momentum to reconstruct for output
!     INTEGER, INTENT(IN)  :: ix ! line of the dylm2 matrix to use actually it is
!                                ! one of the nx spherical integration directions
    REAL(DP), INTENT(IN) :: F_lm(i%m,3,lmaxq_in**2,nspin)  ! Y_lm expansion of F
    !REAL(DP), INTENT(IN) :: F_rad(i%m,3,nspin)             ! F along direction ix
    REAL(DP), INTENT(OUT):: div_F_lm(i%m,lmaxq_out**2,nspin)! div(F) 
    !
    REAL(DP)             :: div_F_rad(i%m,nx,nspin)         ! div(F) on rad. grid
    REAL(DP)             :: aux(i%m),aux2(i%m)              ! workspace
    ! counters on: spin, angular momentum, radial grid point:
    INTEGER              :: is, lm, k, ix

    OPTIONAL_CALL start_clock ('PAW_div')

! This is the divergence in spherical coordinates:
!     {1 \over r^2}{\partial ( r^2 A_r ) \over \partial r} 
!   + {1 \over r\sin\theta}{\partial \over \partial \theta} (  A_\theta\sin\theta )
!   + {1 \over r\sin\theta}{\partial A_\phi \over \partial \phi}
!
! The derivative sum_LM d(Y_LM sin(theta) )/dtheta will be expanded as:
! sum_LM ( Y_lm cos(theta) + sin(theta) dY_lm/dtheta )

    ! The radial component of the divergence is computed last, for practical reasons

    CALL errore('PAW_divergence', 'More angular momentum components are needed (in input)'//&
                ' to provide the number you have requested (in output)', lmaxq_out-lmaxq_in+2)

    ! initialize
    div_F_rad(:,:,:) = 0._dp

!#ifdef __DONT_DO_THAT_THEN
    ! phi component
    DO is = 1,nspin
    DO ix = 1,nx
    aux(:) = 0._dp
        ! this derivative has no spherical component, so lm strarts from 2
        DO lm = 2,lmaxq_in**2
            aux(1:i%m) = aux(1:i%m) + dylmp(ix,lm)* (F_lm(1:i%m,2,lm,is)) &
                                    * g(i%t)%rm1(1:i%m) !/sin_th(ix) 
        ! as for PAW_gradient this is already present in dylmp --^
        ENDDO
        div_F_rad(1:i%m,ix,is) = div_F_rad(1:i%m,ix,is)+aux(1:i%m)
    ENDDO
    ENDDO
!#endif

!#ifdef __DONT_DO_THAT_THEN
    ! theta component
    DO is = 1,nspin
    DO ix = 1,nx
    aux(:) = 0._dp
        ! this derivative has a spherical component too!
        DO lm = 1,lmaxq_in**2
            aux(1:i%m) = aux(1:i%m) + F_lm(1:i%m,3,lm,is) &
                                    *( dylmt(ix,lm)*sin_th(ix) + ylm(ix,lm)*cos_th(ix) ) &
                                    * g(i%t)%rm1(1:i%m) /sin_th(ix)
        ENDDO
        div_F_rad(1:i%m,ix,is) = div_F_rad(1:i%m,ix,is)+aux(1:i%m)
    ENDDO
    ENDDO

    ! Convert what I have done so forth to Y_lm
    CALL PAW_rad2lm(i, div_F_rad, div_F_lm, lmaxq_out)
!#endif

!#ifdef __DONT_DO_THAT_THEN
    ! 1. compute the partial radial derivative d/dr
    !div_F_lm(:,:,:) = 0._dp !!! REMOVE: this kills radial components
    DO is = 1,nspin
    DO lm = 1,lmaxq_out**2
        ! Derive along \hat{r} (F should already contains a r**2 factor, otherwise
        ! it may be better to expand (1/r**2) d(A*r**2)/dr = dA/dr + 2A/r)
        CALL radial_gradient(F_lm(1:i%m,1,lm,is), aux, g(i%t)%r, i%m, 1)
        ! Sum it in the divergence: it is already in the right Y_lm form
        ! multiply by the r**-2 factor (notice that for the angular parts the reason is different
        ! than for the radial one!
        div_F_lm(1:i%m,lm,is) = (div_F_lm(1:i%m,lm,is) + aux(1:i%m))*g(i%t)%rm2(1:i%m)
    ENDDO
    ENDDO
!#endif

    OPTIONAL_CALL stop_clock ('PAW_div')

END SUBROUTINE PAW_divergence
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! build gradient of radial charge distribution from its spherical harmonics expansion
!!! uses pre-computed rho_rad
SUBROUTINE PAW_gradient(i, ix, rho_lm, rho_rad, rho_core, grho_rad2, grho_rad)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi
    USE radial_grids,           ONLY : ndmx
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE ions_base,              ONLY : ityp
    USE parameters,             ONLY : ntypx
    USE atom,                   ONLY : g => rgrid

    INTEGER, INTENT(IN)  :: ix ! line of the dylm2 matrix to use actually it is
                               ! one of the nx spherical integration directions
    TYPE(paw_info)  :: i                              ! atom's minimal info
    REAL(DP), INTENT(IN) :: rho_lm(i%m,lmaxq**2,nspin)! Y_lm expansion of rho
    REAL(DP), INTENT(IN) :: rho_rad(i%m,nspin)        ! radial density along direction ix
    REAL(DP), INTENT(IN) :: rho_core(ndmx,ntypx)      ! core density
    REAL(DP), INTENT(OUT):: grho_rad2(i%m,nspin)      ! |grad(rho)|^2 on rad. grid
    REAL(DP), OPTIONAL,INTENT(OUT):: grho_rad(i%m,3,nspin) ! vector gradient (only for gcxc)
    !              r, theta and phi components ---^
    !
    REAL(DP)             :: aux(i%m),aux2(i%m)       ! workspace
    ! counters on: spin, angular momentum, atom type, radial grid point:
    INTEGER              :: is, lm, k

    OPTIONAL_CALL start_clock ('PAW_grad')
    ! 1. build real charge density = rho/r**2 + rho_core
    ! 2. compute the partial derivative of rho_rad
    grho_rad2(:,:) = 0._dp
    DO is = 1,nspin
        ! build real charge density
        aux(1:i%m) = rho_rad(1:i%m,is)*g(i%t)%rm2(1:i%m) &
                          + rho_core(1:i%m,i%t)/nspin
        CALL radial_gradient(aux, aux2, g(i%t)%r, i%m, 1)
        ! compute the square
        grho_rad2(:,is) = aux2(:)**2
        ! store in vector gradient, if present:
        IF (present(grho_rad)) grho_rad(:,1,is) = aux2(:)
    ENDDO

    spin: &
    DO is = 1,nspin
        aux(:)  = 0._dp
        aux2(:) = 0._dp
        ! Spherical (lm=1) component (that would also include core correction) can be omitted
        ! as its contribution to non-radial derivative is zero
        DO lm = 2,lmaxq**2
            ! 5. [ \sum_{lm} rho(r) (dY_{lm}/dphi /cos(theta))  ]**2
            aux(1:i%m) = aux(1:i%m) + dylmp(ix,lm)* (rho_lm(1:i%m,lm,is))
            ! 6. [ \sum_{lm} rho(r) (dY_{lm}/dtheta)  ]**2
            aux2(1:i%m) = aux2(1:i%m) + dylmt(ix,lm)* (rho_lm(1:i%m,lm,is))
        ENDDO
        ! Square and sum up these 2 components, the (1/r**2)**3 factor come from:
        !  a. 1/r**2 from the derivative in spherical coordinates
        !  b. (1/r**2)**2 from rho_lm being multiplied by r**2 
        !     (as the derivative is orthogonal to r you can multiply after deriving)
        grho_rad2(1:i%m,is) = grho_rad2(1:i%m,is)&
                                + (aux(1:i%m)**2 + aux2(1:i%m)**2)&
                                    * g(i%t)%rm2(1:i%m)**3
        ! Store vector components:
        IF (present(grho_rad)) THEN
            grho_rad(1:i%m,2,is) = aux(1:i%m)  *g(i%t)%rm3(1:i%m) ! phi 
            grho_rad(1:i%m,3,is) = aux2(1:i%m) *g(i%t)%rm3(1:i%m) ! theta
        ENDIF
    ENDDO spin

    OPTIONAL_CALL stop_clock ('PAW_grad')

END SUBROUTINE PAW_gradient

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! computes H  potential from rho, used by PAW_h_energy and PAW_ddot
SUBROUTINE PAW_h_potential(i, rho_lm, v_lm, energy)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi, e2
    USE radial_grids,           ONLY : hartree
    USE uspp_param,             ONLY : lmaxq
    USE ions_base,              ONLY : ityp
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid

    TYPE(paw_info)  :: i                          ! atom's minimal info
    ! charge density as lm components already summed on spin:
    REAL(DP), INTENT(IN)  :: rho_lm(i%m,lmaxq**2,nspin)
    REAL(DP), INTENT(OUT) :: v_lm  (i%m,lmaxq**2) ! potential as lm components
    REAL(DP),INTENT(OUT),OPTIONAL :: energy       ! if present, compute energy
    !
    REAL(DP)              :: aux(i%m) ! workspace
    REAL(DP)              :: pref     ! workspace

    INTEGER               :: lm,l     ! counter on composite angmom lm = l**2 +m
    INTEGER               :: k        ! counter on radial grid (only for energy) 
    REAL(DP)              :: e        ! workspace

    OPTIONAL_CALL start_clock ('PAW_h_pot')

    ! this loop computes the hartree potential using the following formula:
    !               l is the first argument in hartree subroutine
    !               r1 = min(r,r'); r2 = MAX(r,r')
    ! V_h(r) = \sum{lm} Y_{lm}(\hat{r})/(2L+1) \int dr' 4\pi r'^2 \rho^{lm}(r') (r1^l/r2^{l+1})
    !     done here --> ^^^^^^^^^^^^^^^^^^^^^           ^^^^^^^^^^^^^^^^^^^^^^ <-- input to the hartree subroutine
    !                 output from the h.s. --> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    DO lm = 1, lmaxq**2
        l = INT(sqrt(DBLE(lm-1))) ! l has to start from *zero*
            pref = e2*fpi/REAL(2*l+1)
            DO k = 1, i%m
                aux(k) = pref * SUM(rho_lm(k,lm,1:nspin))
            ENDDO
            CALL hartree(l, 2*l+2, i%m, g(i%t), aux(:), v_lm(:,lm))
    ENDDO

    ! compute energy if required:
    ! E_h = \sum_lm \int v_lm(r) (rho_lm(r) r^2) dr
    IF(present(energy)) THEN
    energy = 0._dp
    DO lm = 1, lmaxq**2
        ! I can use v_lm as workspace
        DO k = 1, i%m
            aux(k) = v_lm(k,lm) * SUM(rho_lm(k,lm,1:nspin))
        ENDDO
        CALL simpson (i%m, aux, g(i%t)%rab, e)
        !
        ! Sum all the energies in PAW_ddot
        energy = energy + e
        !
    ENDDO
    ! fix double counting
    energy = energy/2._dp
    ENDIF

    OPTIONAL_CALL stop_clock ('PAW_h_pot')

END SUBROUTINE PAW_h_potential

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! sum up pfuncs x occupation to build radial density's angular momentum components
SUBROUTINE PAW_rho_lm(i, becsum, pfunc, rho_lm, aug)
    USE kinds,                  ONLY : DP
    USE ions_base,              ONLY : ntyp => nsp, nat 
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nh, lmaxq, nhm
    USE uspp,                   ONLY : indv, ap, nhtolm,lpl,lpx
    USE parameters,             ONLY : lqmax
    USE constants,              ONLY : eps12
    USE radial_grids,           ONLY : ndmx
    USE grid_paw_variables,     ONLY : augfun_t, nbrx
    USE atom,                   ONLY : g => rgrid

    TYPE(paw_info)  :: i                                    ! atom's minimal info
    REAL(DP), INTENT(IN)  :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupation
    REAL(DP), INTENT(IN)  :: pfunc(ndmx,nbrx,nbrx,ntyp)     ! psi_i * psi_j
    REAL(DP), INTENT(OUT) :: rho_lm(i%m,lmaxq**2,nspin)     ! AE charge density on rad. grid
    TYPE(augfun_t), OPTIONAL,INTENT(IN) :: &
                             aug(ntyp) ! augmentation functions (only for PS part)

    REAL(DP)                :: pref ! workspace (ap*becsum)

    INTEGER                 :: nb,mb, &     ! counters for pfunc nb,mb = 1, nh
                               nmb, &       ! composite "triangular" index for pfunc nmb = 1,nh*(nh+1)/2
                               lm,lp,l, &   ! counters for angular momentum lm = l**2+m
                               ispin        ! counter for spin (FIXME: may be unnecessary)

    ! This subroutine computes the angular momentum components of rho
    ! using the following formula:
    !   rho(\vec{r}) = \sum_{LM} Y_{LM} \sum_{i,j} (\hat{r}) a_{LM}^{(lm)_i(lm)_j} becsum_ij pfunc_ij(r)
    !
    ! actually different angular momentum components are stored separately:
    !   rho^{LM}(\vec{r}) = \sum_{i,j} (\hat{r}) a_{LM}^{(lm)_i(lm)_j} becsum_ij pfunc_ij(r)
    !
    ! notice that pfunc's are already multiplied by r^2 and they are indexed on the atom
    ! (they only depends on l, not on m), the augmentation charge depend only on l
    ! but the becsum depend on both l and m

    OPTIONAL_CALL start_clock ('PAW_rho_lm')

    ! initialize density
    rho_lm(:,:,:) = 0._dp

    spins: DO ispin = 1, nspin
    nmb = 0
        ! loop on all pfunc for this kind of pseudo
        DO nb = 1, nh(i%t)
        DO mb = nb, nh(i%t)
            nmb = nmb+1 ! nmb = 1, nh*(nh+1)/2
            IF (ABS(becsum(nmb,i%a,ispin)) < eps12) CYCLE
            !
            angular_momentum: &
            DO lp = 1, lpx (nhtolm(mb,i%t), nhtolm(nb,i%t)) !lmaxq**2
                ! the lpl array contains the possible combination of LM,lm_j,lm_j that
                ! have non-zero a_{LM}^{(lm)_i(lm)_j} (it saves some loops)
                lm = lpl (nhtolm(mb,i%t), nhtolm(nb,i%t), lp)
                ! becsum already contains a factor 2 for off-diagonal pfuncs
                pref = becsum(nmb,i%a,ispin) * ap(lm, nhtolm(nb,i%t), nhtolm(mb,i%t))
                !
                rho_lm(1:i%m,lm,ispin) = rho_lm(1:i%m,lm,ispin) +&
                                pref * pfunc(1:i%m, indv(nb,i%t), indv(mb,i%t), i%t)
                IF (present(aug)) THEN
                    ! if I'm doing the pseudo part I have to add the augmentation charge
                    l = INT(SQRT(DBLE(lm-1))) ! l has to start from zero
                    rho_lm(1:i%m,lm,ispin) = rho_lm(1:i%m,lm,ispin) +&
                                pref * aug(i%t)%fun(1:i%m, indv(nb,i%t), indv(mb,i%t), l)
                ENDIF ! augfun
            ENDDO angular_momentum 
        ENDDO !mb
        ENDDO !nb
    ENDDO spins

    OPTIONAL_CALL stop_clock ('PAW_rho_lm')

END SUBROUTINE PAW_rho_lm

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! build radial charge distribution from its spherical harmonics expansion
SUBROUTINE PAW_lm2rad(i, ix, F_lm, F_rad)
    USE kinds,                  ONLY : DP
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin

    TYPE(paw_info)              :: i  ! atom's minimal info
    INTEGER                     :: ix ! line of the ylm matrix to use
                                      ! actually it is one of the nx directions
    REAL(DP), INTENT(IN)        :: F_lm(i%m,lmaxq**2,nspin)! Y_lm expansion of rho
    REAL(DP), INTENT(OUT)       :: F_rad(i%m,nspin)        ! charge density on rad. grid
    !
    INTEGER                     :: ispin, lm ! counters on angmom and spin

    OPTIONAL_CALL start_clock ('PAW_lm2rad')
    F_rad(:,:) = 0._dp
    ! cycling on spin is a bit less general...
    spins: DO ispin = 1,nspin
        DO lm = 1, lmaxq**2
            F_rad(:,ispin) = F_rad(:,ispin) +&
                    ylm(ix,lm)*F_lm(:,lm,ispin)
        ENDDO ! lm
    ENDDO spins

    OPTIONAL_CALL stop_clock ('PAW_lm2rad')

END SUBROUTINE PAW_lm2rad

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! computes F_lm(r) = \int d \Omega F(r,th,ph) Y_lm(th,ph)
SUBROUTINE PAW_rad2lm(i, F_rad, F_lm, lmaxq_loc)
    USE kinds,                  ONLY : DP
    !USE uspp_param,             ONLY : lmaxq ! <-- now passed as parameter
    USE lsda_mod,               ONLY : nspin

    TYPE(paw_info)       :: i         ! atom's minimal info
    REAL(DP), INTENT(OUT):: F_lm(i%m, lmaxq_loc**2, nspin)! lm component of F up to lmaxq_loc
    REAL(DP), INTENT(IN) :: F_rad(i%m, nx, nspin)         ! radial samples of F
    INTEGER,  INTENT(IN) :: lmaxq_loc ! in some cases I have to keep higher angular components
                                      ! than the default ones (=lmaxq =the ones present in rho)
    !
    INTEGER                     :: ix    ! counter for integration
    INTEGER                     :: lm    ! counter for angmom
    INTEGER                     :: ispin ! counter for spin

    OPTIONAL_CALL start_clock ('PAW_rad2lm')
    F_lm(:,:,:) = 0._dp

    DO ispin = 1,nspin
    DO lm = 1,lmaxq_loc**2
    DO ix = 1, nx
        F_lm(1:i%m, lm, ispin) = F_lm(1:i%m, lm, ispin) &
                                + ww(ix) * ylm(ix,lm) * F_rad(1:i%m,ix,ispin)
                                ! |         ^- spherical harmonic
                                ! +- integration weight
    ENDDO
    ENDDO
    ENDDO

    OPTIONAL_CALL stop_clock ('PAW_rad2lm')

END SUBROUTINE PAW_rad2lm


END MODULE rad_paw_routines
