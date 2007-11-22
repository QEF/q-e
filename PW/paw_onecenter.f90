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
#include "f_defs.h"
MODULE paw_onecenter
    !
    USE kinds,          ONLY : DP
    USE paw_variables,  ONLY : paw_info, saved, xlm, rad
    !
    IMPLICIT NONE

    ! entry points:
    PUBLIC :: PAW_potential  ! prepare paw potential and store it,
                             ! also computes energy if required
    PUBLIC :: PAW_integrate  ! computes \int v(r) \times n(r) dr
    PUBLIC :: PAW_ddot       ! error estimate for mix_rho
    PUBLIC :: PAW_newd       ! computes descreening coefficients
    PUBLIC :: PAW_symmetrize ! symmetrize becsums
    !
    PRIVATE

    ! the following macro controls the use of several fine-grained clocks
    ! set it to 'if(.false.) CALL' (without quotes) in order to disable them,
    ! set it to 'CALL' to enable them.
#define OPTIONAL_CALL CALL

 CONTAINS

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! Computes V_h and V_xc using the "density" becsum provided, and stores the
!!! total potential in the global static save%v variable
!!
SUBROUTINE PAW_potential(becsum, energy, e_cmp)
    USE atom,              ONLY : g => rgrid
    USE ions_base,         ONLY : nat, ityp
    USE lsda_mod,          ONLY : nspin
    USE uspp_param,        ONLY : nhm, upf

    REAL(DP), INTENT(IN)    :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations
    REAL(DP),INTENT(OUT),OPTIONAL :: energy          ! if present compute E[rho]
    REAL(DP),INTENT(OUT),OPTIONAL :: e_cmp(nat,2,2)
    !
    INTEGER, PARAMETER      :: AE = 1, PS = 2,&      ! All-Electron and Pseudo
                               XC = 1, H  = 2        ! XC and Hartree
    REAL(DP), POINTER       :: rho_core(:)           ! pointer to AE/PS core charge density 
    TYPE(paw_info)          :: i                     ! minimal info on atoms
    INTEGER                 :: i_what                ! counter on AE and PS
    INTEGER                 :: is                    ! spin index
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
    i%b = upf(i%t)%nbeta
    i%l = upf(i%t)%paw%lmax_rho+1
    !
    ifpaw: IF (upf(i%t)%tpawp) THEN
        !
        ! Arrays are allocated inside the cycle to allow reduced
        ! memory usage as differnt atoms have different meshes (they
        ! also have different lmax, I will fix this sooner or later)
        ALLOCATE(rho_lm(i%m,i%l**2,nspin))
        ALLOCATE(v_lm(i%m,i%l**2,nspin))
        !
        whattodo: DO i_what = AE, PS
            ! STEP: 1 [ build rho_lm (PAW_rho_lm) ]
            NULLIFY(rho_core)
            IF (i_what == AE) THEN
                ! to pass the atom indes is dirtyer but faster and
                ! uses less memory than to pass only a hyperslice of the array
                CALL PAW_rho_lm(i, becsum, upf(i%t)%paw%pfunc, rho_lm)
                ! used later for xc potential:
                rho_core => upf(i%t)%paw%ae_rho_atc
                ! sign to sum up the enrgy
                sgn = +1._dp
            ELSE
                CALL PAW_rho_lm(i, becsum, upf(i%t)%paw%ptfunc, rho_lm, upf(i%t)%paw%aug)
                !    optional argument for pseudo part --> ^^^
                rho_core => upf(i%t)%rho_atc ! as before
                sgn = -1._dp                 ! as before
            ENDIF

        ! cleanup previously stored potentials
        saved(i%a)%v(:,:,:,i_what) = 0._dp
        ! First compute the Hartree potential (it does not depend on spin...):
        CALL PAW_h_potential(i, rho_lm, v_lm(:,:,1), energy)
        ! using "energy" as the in/out parameter I save a double call, but I have to do this:
        IF (present(energy)) energy_h = energy
        DO is = 1,nspin ! ... so it has to be copied to all spin components
            saved(i%a)%v(:,:,is,i_what) = v_lm(:,:,1)
        ENDDO

        ! Than the XC one:
        CALL PAW_xc_potential(i, rho_lm, rho_core, v_lm, energy)
        IF (present(energy)) energy_xc = energy
        saved(i%a)%v(:,:,:,i_what) = saved(i%a)%v(:,:,:,i_what) &
                                    + v_lm(:,:,:)
        IF (present(energy)) energy_tot = energy_tot + sgn*(energy_xc + energy_h)
        IF (present(e_cmp)) THEN
            e_cmp(na, 1, i_what) = energy_xc
            e_cmp(na, 2, i_what) = energy_h
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

SUBROUTINE PAW_symmetrize(becsum)
    USE lsda_mod,          ONLY : nspin
    USE uspp_param,        ONLY : nhm
    USE ions_base,         ONLY : nat, ityp
    USE symme,             ONLY : nsym, irt, d1, d2, d3
    USE uspp,              ONLY : nhtolm,nhtol
    USE uspp_param,        ONLY : nh, upf
    USE control_flags,     ONLY : nosym, gamma_only

    REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations

    REAL(DP)                :: becsym(nhm*(nhm+1)/2,nat,nspin)! symmetrized becsum

    INTEGER :: is, na, nt   ! counters on spin, atom, atom-type
    INTEGER :: ma           ! atom symmetric to na
    INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
               l_i, l_j, m_i, m_j
    INTEGER :: m_o, m_u     ! counters for sums on m
    INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    INTEGER :: isym         ! counter for symmetry operation
    INTEGER :: first_nat, last_nat ! to parallelize on atoms

    ! The following mess is necessary because the symmetrization operation
    ! in LDA+U code is simpler than in PAW, so the required quantities are
    ! represented in a simple but not general way.
    ! I will fix this when everything works.
    REAL(DP), TARGET :: d0(1,1,48)
    TYPE symmetryzation_tensor
        REAL(DP),POINTER :: d(:,:,:)
    END TYPE symmetryzation_tensor
    TYPE(symmetryzation_tensor) :: D(0:3)
    d0(1,1,:) = 1._dp
    D(0)%d => d0 ! d0(1,1,48)
    D(1)%d => d1 ! d1(3,3,48)
    D(2)%d => d2 ! d2(5,5,48)
    D(3)%d => d3 ! d3(7,7,48)

    IF( gamma_only .or. nosym ) RETURN

! => lm = l**2 + m
! => ih = lm + (l+proj)**2  <-- if the projector index starts from zero!
!       = lm + proj**2 + 2*l*proj
!       = m + l**2 + proj**2 + 2*l*proj
!        ^^^
! Known ih and m_i I can compute the index oh of a different m = m_o but
! the same augmentation channel (l_i = l_o, proj_i = proj_o):
!  oh = ih - m_i + m_o
! this expression should be general inside pwscf.

!#define __DEBUG_PAW_SYM
#ifdef __DEBUG_PAW_SYM
        nt = 1
        na = 1
        is = 1
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            IF (jh > ih) THEN
                ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
            ELSE
                ijh = nh(nt)*(jh-1) - jh*(jh-1)/2 + ih
            ENDIF
            write(*,"(1f10.5)", advance='no') becsum(ijh,na,is)
        ENDDO
            write(*,*)
        ENDDO
#endif

    CALL start_clock('PAW_symme')

    becsym(:,:,:) = 0._dp

    CALL divide (nat, first_nat, last_nat)
    DO is = 1, nspin
    !
    atoms: DO na = first_nat, last_nat
        nt = ityp(na)
        ! No need to symmetrize non-PAW atoms
        IF ( .not. upf(nt)%tpawp ) CYCLE
        !
        DO ih = 1, nh(nt)
        DO jh = ih, nh(nt) ! note: jh >= ih
            ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
            !
            lm_i  = nhtolm(ih,nt)
            lm_j  = nhtolm(jh,nt)
            !
            l_i   = nhtol(ih,nt)
            l_j   = nhtol(jh,nt)
            !
            m_i   = lm_i - l_i**2
            m_j   = lm_j - l_j**2
            !
            DO isym = 1,nsym
                ma = irt(isym,na)
                DO m_o = 1, 2*l_i +1
                DO m_u = 1, 2*l_j +1
                    oh = ih - m_i + m_o
                    uh = jh - m_j + m_u
                    IF ( oh >= uh ) THEN
                        ouh = nh(nt)*(uh-1) - uh*(uh-1)/2 + oh
                    ELSE
                        ouh = nh(nt)*(oh-1) - oh*(oh-1)/2 + uh
                    ENDIF
                    !
                    becsym(ijh, na, is) = becsym(ijh, na, is) &
                        + D(l_i)%d(m_o,m_i, isym) * D(l_j)%d(m_u,m_j, isym) &
                          * becsum(ouh, ma, is) / nsym
                ENDDO ! m_o
                ENDDO ! m_u
            ENDDO ! isym
        ENDDO ! ih
        ENDDO ! jh
    ENDDO atoms ! nat
    ENDDO ! nspin
#ifdef __PARA
    CALL reduce( (nhm*(nhm+1)/2) *nat*nspin, becsym)
#endif

#ifdef __DEBUG_PAW_SYM
        nt = 1
        na = 1
        is = 1
            write(*,*)
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            IF (jh > ih) THEN
                ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
            ELSE
                ijh = nh(nt)*(jh-1) - jh*(jh-1)/2 + ih
            ENDIF
            write(*,"(1f10.5)", advance='no') becsym(ijh,na,is)
        ENDDO
            write(*,*)
        ENDDO
            write(*,*)
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            IF (jh > ih) THEN
                ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
            ELSE
                ijh = nh(nt)*(jh-1) - jh*(jh-1)/2 + ih
            ENDIF
            write(*,"(1f10.5)", advance='no') becsym(ijh,na,is)-becsum(ijh,na,is)
        ENDDO
            write(*,*)
        ENDDO
        write(*,*) "________________________________________________"
#endif

    ! Apply symmetrization:
    becsum(:,:,:) = becsym(:,:,:)

    CALL stop_clock('PAW_symme')

END SUBROUTINE PAW_symmetrize

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! Update the descreening coefficients:
!!! D_ij = \int (v_loc_scf + v_loc_at) p_ij
!!
!! This is subroutine does NOT cycle on atoms because newd likes it this way
SUBROUTINE PAW_newd(d)
    USE constants,         ONLY : sqrtpi, e2
    USE lsda_mod,          ONLY : nspin
    USE ions_base,         ONLY : nat, ityp
    USE atom,              ONLY : g => rgrid
    USE uspp_param,        ONLY : nh, nhm, upf

    REAL(DP), INTENT(INOUT) :: d(:,:,:,:) ! descreening coefficients (AE - pseudo)

    INTEGER                 :: i_what
    INTEGER                 :: na, first_nat, last_nat
    ! fake cross band occupations to select only one pfunc at a time:
    REAL(DP)                :: becfake(nhm*(nhm+1)/2,nat,nspin)

    INTEGER, PARAMETER      :: AE = 1, PS = 2  ! All-Electron and Pseudo
    INTEGER                 :: lm              ! counters on angmom and radial grid
    INTEGER                 :: nb, mb, nmb, is
    !
    REAL(DP), POINTER       :: v_at(:)         ! point to aevloc_at or psvloc_at
    REAL(DP), ALLOCATABLE   :: pfunc_lm(:,:,:) ! aux charge density
    REAL(DP)                :: integral        ! workspace
    REAL(DP)                :: sgn             ! +1 for AE, -1 for pseudo
    !
    TYPE(paw_info)          :: i

    !
    ! At the moment this part is done using a fake becsum with all the states
    ! empty but one. It's inefficient, but safe.
    CALL start_clock ('PAW_newd')

    ! Some initialization
    becfake(:,:,:) = 0._dp
    !
    d(:,:,:,:) = 0._dp

    CALL divide (nat, first_nat, last_nat)

    atoms: DO na = first_nat, last_nat
    !
    i%a = na         ! the index of the atom
    i%t = ityp(na)   ! the type of atom na
    ! skip non-paw atoms
    IF ( .not. upf(i%t)%tpawp ) CYCLE
    i%m = g(i%t)%mesh! radial mesh size for atom na
    i%b = upf(i%t)%nbeta
    i%l = upf(i%t)%paw%lmax_rho+1
        !
        ! Different atoms may have different mesh sizes:
        ALLOCATE(pfunc_lm(i%m, i%l**2, nspin))
        !
        whattodo: DO i_what = AE, PS
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
                    IF (i_what == AE) THEN
                        CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%pfunc, pfunc_lm)
                        v_at => upf(i%t)%paw%ae_vloc
                        sgn = +1._dp
                    ELSE
                        CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%ptfunc, pfunc_lm, upf(i%t)%paw%aug)
                        v_at => upf(i%t)%vloc
                        sgn = -1._dp
                    ENDIF
                    !
                    ! Now I multiply the pfunc and the potential, I can use
                    ! pfunc_lm itself as workspace
                    DO lm = 1,i%l**2
                        IF ( lm == 1 ) THEN
                            pfunc_lm(1:i%m,lm,is) = pfunc_lm(1:i%m,lm,is) * &
                                ( saved(i%a)%v(1:i%m,lm,is,i_what) + e2*sqrtpi*v_at(1:i%m) )
                        ELSE
                            pfunc_lm(1:i%m,lm,is) = pfunc_lm(1:i%m,lm,is) * &
                                saved(i%a)%v(1:i%m,lm,is,i_what)
                        ENDIF
                        !
                        ! Integrate!
                        CALL simpson (upf(i%t)%paw%irmax,pfunc_lm(:,lm,is),g(i%t)%rab,integral)
                        d(nb,mb,i%a,is) = d(nb,mb,i%a,is) + sgn * integral
                    ENDDO
                    ! Symmetrize:
                    d(mb,nb, i%a,is) = d(nb,mb, i%a,is)
                    ! restore becfake to zero
                    becfake(nmb,na,is) = 0._dp
                ENDDO ! mb
                ENDDO ! nb
            ENDDO spins
        ENDDO whattodo

        ! cleanup
        DEALLOCATE(pfunc_lm)
        !
    ENDDO atoms

#ifdef __PARA
    CALL reduce (nhm*nhm*nat*nspin, d )
#endif

    CALL stop_clock ('PAW_newd')


END SUBROUTINE PAW_newd


!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! Compute a new density from becsum, and integrate it with the potential
!!! previously computed and stored in global static variable "saved(na)"
!!
FUNCTION PAW_integrate(becsum)
    USE lsda_mod,          ONLY : nspin
    USE ions_base,         ONLY : nat, ityp
    USE atom,              ONLY : g => rgrid
    USE uspp_param,        ONLY : nhm, upf

    REAL(DP)                :: PAW_integrate

    REAL(DP), INTENT(IN)    :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations

    INTEGER, PARAMETER      :: AE = 1, PS = 2        ! All-Electron and Pseudo
    INTEGER                 :: na,first_nat,last_nat ! atoms counters and indexes
    INTEGER                 :: lm                    ! counters on angmom and radial grid
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
    i%b = upf(i%t)%nbeta
    i%l = upf(i%t)%paw%lmax_rho+1
    !
    ifpaw: IF (upf(i%t)%tpawp) THEN
        !
        ALLOCATE(rho_lm(i%m,i%l**2,nspin))
        !
        whattodo: DO i_what = AE, PS
            !
            IF (i_what == AE) THEN
                CALL PAW_rho_lm(i, becsum, upf(i%t)%paw%pfunc, rho_lm)
                i_sign = +1._dp
            ELSE
                CALL PAW_rho_lm(i, becsum, upf(i%t)%paw%ptfunc, rho_lm, upf(i%t)%paw%aug)
                i_sign = -1._dp
            ENDIF
            !
            ! Compute the integral
            DO is = 1,nspin
            DO lm = 1, i%l**2
                ! I can use rho_lm itself as workspace
                rho_lm(1:i%m,lm,is) = rho_lm(1:i%m,lm,is) &
                                     * saved(i%a)%v(1:i%m,lm,is,i_what)
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
    USE constants,         ONLY : pi
    USE lsda_mod,          ONLY : nspin
    USE ions_base,         ONLY : nat, ityp
    USE atom,              ONLY : g => rgrid
    USE uspp_param,        ONLY : nhm, upf

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

    CALL divide (nat, first_nat, last_nat)
    !
    atoms: DO na = first_nat, last_nat
    !
    i%a = na           ! the index of the atom
    i%t = ityp(na)    ! the type of atom na
    i%m = g(i%t)%mesh ! radial mesh size for atom na
    i%b = upf(i%t)%nbeta
    i%l = upf(i%t)%paw%lmax_rho+1
    !
    ifpaw: IF (upf(i%t)%tpawp) THEN
        !
        ALLOCATE(rho_lm(i%m,i%l**2,nspin))
        ALLOCATE(v_lm(i%m,i%l**2))
        !
        whattodo: DO i_what = AE, PS
            ! Build rho from the occupations in bec1
            IF (i_what == AE) THEN
                CALL PAW_rho_lm(i, bec1, upf(i%t)%paw%pfunc, rho_lm)
                i_sign = +1._dp
            ELSE
                CALL PAW_rho_lm(i, bec1, upf(i%t)%paw%ptfunc, rho_lm, upf(i%t)%paw%aug)
                i_sign = -1._dp
            ENDIF
            !
            ! Compute the hartree potential from bec1
            CALL PAW_h_potential(i, rho_lm, v_lm)
            !
            ! Now a new rho is computed, this time from bec2
            IF (i_what == AE) THEN
                CALL PAW_rho_lm(i, bec2, upf(i%t)%paw%pfunc, rho_lm)
            ELSE
                CALL PAW_rho_lm(i, bec2, upf(i%t)%paw%ptfunc, rho_lm, upf(i%t)%paw%aug)
            ENDIF
            !
            ! Finally compute the integral
            DO lm = 1, i%l**2
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

#ifdef __PARA
    CALL reduce (1, PAW_ddot )
#endif

    CALL stop_clock ('PAW_ddot')


END FUNCTION PAW_ddot
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute xc potential and energy, as
!!! xc functional is not diagonal on angular momentum numerical integration is performed
SUBROUTINE PAW_xc_potential(i, rho_lm, rho_core, v_lm, energy)
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    USE funct,                  ONLY : dft_is_gradient
    USE constants,              ONLY : fpi ! REMOVE

    TYPE(paw_info)  :: i                               ! atom's minimal info
    REAL(DP), INTENT(IN)  :: rho_lm(i%m,i%l**2,nspin)! charge density as lm components
    REAL(DP), INTENT(IN)  :: rho_core(i%m)      ! core charge, radial and spherical
    REAL(DP), INTENT(OUT) :: v_lm(i%m,i%l**2,nspin)  ! potential density as lm components
    REAL(DP),OPTIONAL,INTENT(OUT) :: energy            ! XC energy (if required)
    !
    REAL(DP)              :: rho_loc(2)         ! local density (workspace), up and down
    REAL(DP)              :: v_loc(2)           ! local density (workspace), up and down
    REAL(DP)              :: v_rad(i%m,rad(i%t)%nx,nspin)! radial potential (to be integrated)
    REAL(DP)              :: rho_rad(i%m,nspin) ! workspace (only one radial slice of rho)
    !
    REAL(DP), ALLOCATABLE :: e_rad(:)           ! aux, used to store radial slices of energy
    REAL(DP)              :: e                  ! aux, used to integrate energy
    !
    INTEGER               :: ix,k               ! counters on directions and radial grid
    INTEGER               :: lsd                ! switch for local spin density

    REAL(DP), EXTERNAL    :: exc_t              ! computes XC energy

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
    DO ix = 1, rad(i%t)%nx
        ! *** LDA (and LSDA) part (no gradient correction) ***
        ! convert _lm density to real density along ix
        CALL PAW_lm2rad(i, ix, rho_lm, rho_rad)
        !
        ! compute the potential along ix
        DO k = 1,i%m
            rho_loc(1:nspin) = rho_rad(k,1:nspin)*g(i%t)%rm2(k)
            CALL vxc_t(rho_loc, rho_core(k), lsd, v_loc)
            v_rad(k,ix,1:nspin) = v_loc(1:nspin)
            IF (present(energy)) &
                e_rad(k) = exc_t(rho_loc, rho_core(k), lsd) &
                          * ( SUM(rho_rad(k,1:nspin)) + rho_core(k)*g(i%t)%r2(k) )
        ENDDO
        IF (present(energy)) THEN
            CALL simpson(i%m, e_rad, g(i%t)%rab, e)
            energy = energy + e * rad(i%t)%ww(ix)
        ENDIF
    ENDDO
    IF(present(energy)) DEALLOCATE(e_rad)

    ! Recompose the sph. harm. expansion
    CALL PAW_rad2lm(i, v_rad, v_lm, i%l)

!
!  Passing this charge : fpi*rho_lm(:,1,1) and this rho core : fpi*rho_core(:,1)*g(i%t)%r2(:) 
!  to the original vxcgc routine yeld almost the same ddd coefficients as using the charge loaded
!  from file, produced by the atomic code.
!  What the heck is happening in between????
!
    ! Add gradient correction, if necessary
    IF( dft_is_gradient() ) &
        CALL PAW_gcxc_potential(i, rho_lm,rho_core, v_lm,energy)


    OPTIONAL_CALL stop_clock ('PAW_xc_pot')

END SUBROUTINE PAW_xc_potential
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! add gradient correction to v_xc, code mostly adapted from ../atomic/vxcgc.f90
!!! in order to support non-spherical charges (as Y_lm expansion)
!!! Note that the first derivative in vxcgc becames a gradient, while the second is a divergence.
!!! We also have to temporary store some additional Y_lm components in order not to loose
!!! precision during teh calculation, even if only the ones up to lmax_rho (the maximum in the
!!! density of charge) matter when computing \int v * rho 
SUBROUTINE PAW_gcxc_potential(i, rho_lm,rho_core, v_lm, energy)
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    USE constants,              ONLY : sqrtpi, fpi,pi,e2, eps => eps12, eps2 => eps24
    USE funct,                  ONLY : gcxc, gcx_spin, gcc_spin
    !
    TYPE(paw_info)  :: i                                  ! atom's minimal info
    REAL(DP), INTENT(IN)    :: rho_lm(i%m,i%l**2,nspin) ! charge density as lm components
    REAL(DP), INTENT(IN)    :: rho_core(i%m)              ! core charge, radial and spherical
    REAL(DP), INTENT(INOUT) :: v_lm(i%m,i%l**2,nspin)   ! potential to be updated
    REAL(DP),OPTIONAL,INTENT(INOUT) :: energy             ! if present, add GC to energy

    REAL(DP)                :: rho_rad(i%m,nspin)! charge density sampled
    REAL(DP)                :: grad(i%m,3,nspin) ! gradient
    REAL(DP)                :: grad2(i%m,nspin)  ! square modulus of gradient
                                                             ! (first of charge, than of hamiltonian)
    REAL(DP)                :: gc_rad(i%m,rad(i%t)%nx,nspin) ! GC correction to V
    REAL(DP)                :: gc_lm(i%m,i%l**2,nspin)     ! GC correction to V
    REAL(DP)                :: h_rad(i%m,3,rad(i%t)%nx,nspin)! hamiltonian (vector field)
    REAL(DP)                :: h_lm(i%m,3,(i%l+xlm)**2,nspin)! hamiltonian (vector field)
                                                             !!! expanded to higher lm than rho !!!
    REAL(DP)                :: div_h(i%m,i%l**2,nspin)  ! div(hamiltonian)

    REAL(DP),ALLOCATABLE    :: e_rad(:)                   ! aux, used to store energy
    REAL(DP)                :: e                          ! aux, used to integrate energy

    INTEGER  :: k, ix, is, lm                             ! counters on spin and mesh
    REAL(DP) :: sx,sc,v1x,v2x,v1c,v2c                     ! workspace
    REAL(DP) :: v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw  ! workspace
    REAL(DP) :: sgn, arho                                 ! workspace
    REAL(DP) :: rup, rdw, co2                             ! workspace
    REAL(DP) :: rh, zeta, grh2

    OPTIONAL_CALL start_clock ('PAW_gcxc_v')

    IF (present(energy)) ALLOCATE(e_rad(i%m))

    spin:&
    IF ( nspin == 1 ) THEN
        !
        !     GGA case
        !
        DO ix = 1,rad(i%t)%nx
        !
        !  WARNING: the next 2 calls are duplicated for spin==2
        CALL PAW_lm2rad(i, ix, rho_lm, rho_rad)
        CALL PAW_gradient(i, ix, rho_lm, rho_rad, rho_core, &
                          grad2, grad)
        DO k = 1, i%m
            ! arho is the absolute value of real charge, sgn is its sign
            arho = rho_rad(k,1)*g(i%t)%rm2(k) + rho_core(k)
            sgn  = SIGN(1._dp,arho)
            arho = ABS(arho)

            ! I am using grad(rho)**2 here, so its eps has to be eps**2
            IF ( (arho>eps) .and. (grad2(k,1)>eps2) ) THEN
                CALL gcxc(arho,grad2(k,1), sx,sc,v1x,v2x,v1c,v2c)
                IF (present(energy)) &
                    e_rad(k)    = sgn *e2* (sx+sc) * g(i%t)%r2(k)
                gc_rad(k,ix,1)  = (v1x+v1c)!*g(i%t)%rm2(k)
                h_rad(k,:,ix,1) = (v2x+v2c)*grad(k,:,1)*g(i%t)%r2(k)
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
            energy = energy + e * rad(i%t)%ww(ix)
        ENDIF
        ENDDO
    ELSEIF ( nspin == 2 ) THEN
        !
        !   this is the \sigma-GGA case
        DO ix = 1,rad(i%t)%nx
        !
        CALL PAW_lm2rad(i, ix, rho_lm, rho_rad)
        CALL PAW_gradient(i, ix, rho_lm, rho_rad, rho_core, &
                          grad2, grad)
        !
        DO k = 1,i%m
            !
            ! Prepare the necessary quantities
            ! rho_core is considered half spin up and half spin down:
            co2 = rho_core(k)/2._dp
            ! than I build the real charge dividing by r**2
            rup = rho_rad(k,1)*g(i%t)%rm2(k) + co2
            rdw = rho_rad(k,2)*g(i%t)%rm2(k) + co2
            ! bang!
            CALL gcx_spin (rup, rdw, grad2(k,1), grad2(k,2), &
                           sx, v1xup, v1xdw, v2xup, v2xdw)

            rh = rup + rdw ! total charge
            IF ( rh > eps ) THEN
                zeta = (rup - rdw ) / rh ! polarization
                !grh2 = (grho (i, 1) + grho (i, 2) ) **2 
                grh2 =  (grad(k,1,1) + grad(k,1,2))**2 &
                      + (grad(k,2,1) + grad(k,2,2))**2 &
                      + (grad(k,3,1) + grad(k,3,2))**2
                CALL gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
            ELSE
                sc    = 0._dp
                v1cup = 0._dp
                v1cdw = 0._dp
                v2c   = 0._dp
            ENDIF
            IF (present(energy)) &
                e_rad(k)    = e2*(sx+sc)* g(i%t)%r2(k)
            gc_rad(k,ix,1)  = (v1xup+v1cup)!*g(i%t)%rm2(k)
            gc_rad(k,ix,2)  = (v1xdw+v1cdw)!*g(i%t)%rm2(k)
            !
            h_rad(k,:,ix,1) =( (v2xup+v2c)*grad(k,:,1)+v2c*grad(k,:,2) )*g(i%t)%r2(k)
            h_rad(k,:,ix,2) =( (v2xdw+v2c)*grad(k,:,2)+v2c*grad(k,:,1) )*g(i%t)%r2(k)
        ENDDO ! k
        ! integrate energy (if required)
        ! NOTE: this integration is duplicated for every spin, FIXME!
        IF (present(energy)) THEN
            CALL simpson(i%m, e_rad, g(i%t)%rab, e)
            energy = energy + e * rad(i%t)%ww(ix)
        ENDIF
        ENDDO ! ix
    ELSEIF ( nspin == 4 ) THEN
        CALL errore('PAW_gcxc_v', 'non-collinear not yet implemented!', 1)
    ELSE spin
        CALL errore('PAW_gcxc_v', 'unknown spin number', 2)
    ENDIF spin
    !
    ! convert the first part of the GC correction back to spherical harmonics
    CALL PAW_rad2lm(i, gc_rad, gc_lm, i%l)
    !
    ! We need the gradient of h to calculate the last part of the exchange
    ! and correlation potential. First we have to convert H to its Y_lm expansion
    CALL PAW_rad2lm3(i, h_rad, h_lm, i%l+xlm)
    !
    ! Compute div(H)
    CALL PAW_divergence(i, h_lm, div_h, i%l+xlm, i%l)
    !                         input max lm --^     ^-- output max lm

    ! Finally sum it back into v_xc
            OPTIONAL_CALL start_clock ('PAW_gcxc_3c')
    DO is = 1,nspin
    DO lm = 1,i%l**2
            !v_lm(1:i%m,lm,is) = v_lm(1:i%m,lm,is) + e2*(gc_lm(1:i%m,lm,is)*g(i%t)%r2(1:i%m)-div_h(1:i%m,lm,is))
            v_lm(1:i%m,lm,is) = v_lm(1:i%m,lm,is) + e2*(gc_lm(1:i%m,lm,is)-div_h(1:i%m,lm,is))
    ENDDO
    ENDDO

    OPTIONAL_CALL stop_clock ('PAW_gcxc_v')

END SUBROUTINE PAW_gcxc_potential
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! compute divergence of a vector field (actutally the hamiltonian)
!!! it is assumed that: 1. the input function is multiplied by r**2; 
!!! 2. the output function is multiplied by r**2 too
SUBROUTINE PAW_divergence(i, F_lm, div_F_lm, lmaxq_in, lmaxq_out)
    USE constants,              ONLY : sqrtpi, fpi, eps12, e2
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid

    TYPE(paw_info)  :: i              ! atom's minimal info
    INTEGER, INTENT(IN)  :: lmaxq_in  ! max angular momentum to derive
                                      ! (divergence is reliable up to lmaxq_in-2)
    INTEGER, INTENT(IN)  :: lmaxq_out ! max angular momentum to reconstruct for output
    REAL(DP), INTENT(IN) :: F_lm(i%m,3,lmaxq_in**2,nspin)   ! Y_lm expansion of F
    REAL(DP), INTENT(OUT):: div_F_lm(i%m,lmaxq_out**2,nspin)! div(F) 
    !
    REAL(DP)             :: div_F_rad(i%m,rad(i%t)%nx,nspin)! div(F) on rad. grid
    REAL(DP)             :: aux(i%m)!,aux2(i%m)              ! workspace
    ! counters on: spin, angular momentum, radial grid point:
    INTEGER              :: is, lm, ix

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

    ! phi component
    DO is = 1,nspin
    DO ix = 1,rad(i%t)%nx
    aux(:) = 0._dp
        ! this derivative has no spherical component, so lm strarts from 2
        DO lm = 2,lmaxq_in**2
            aux(1:i%m) = aux(1:i%m) + rad(i%t)%dylmp(ix,lm)* (F_lm(1:i%m,2,lm,is))! &
                                    !* g(i%t)%rm1(1:i%m) !/sin_th(ix) 
        ! as for PAW_gradient this is already present in dylmp --^
        ENDDO
        div_F_rad(1:i%m,ix,is) = aux(1:i%m)
    ENDDO
    ENDDO

    ! theta component
    DO is = 1,nspin
    DO ix = 1,rad(i%t)%nx
    aux(:) = 0._dp
        ! this derivative has a spherical component too!
        DO lm = 1,lmaxq_in**2
            aux(1:i%m) = aux(1:i%m) + F_lm(1:i%m,3,lm,is) &
                                    *( rad(i%t)%dylmt(ix,lm)  &
                                     + rad(i%t)%ylm(ix,lm) * rad(i%t)%cotg_th(ix) )
        ENDDO
        div_F_rad(1:i%m,ix,is) = div_F_rad(1:i%m,ix,is)+aux(1:i%m)
    ENDDO
    ENDDO

    ! Convert what I have done so forth to Y_lm
    CALL PAW_rad2lm(i, div_F_rad, div_F_lm, lmaxq_out)
    ! Multiply by 1/r**3: 1/r is for theta and phi componente only
    ! 1/r**2 is common to all the three components.
    DO is = 1,nspin
    DO lm = 1,lmaxq_out**2
        div_F_lm(1:i%m,lm,is) = div_F_lm(1:i%m,lm,is) * g(i%t)%rm3(1:i%m)
    ENDDO
    ENDDO

    ! Compute partial radial derivative d/dr
    DO is = 1,nspin
    DO lm = 1,lmaxq_out**2
        ! Derive along \hat{r} (F already contains a r**2 factor, otherwise
        ! it may be better to expand (1/r**2) d(A*r**2)/dr = dA/dr + 2A/r)
        CALL radial_gradient(F_lm(1:i%m,1,lm,is), aux, g(i%t)%r, i%m, 1)
        ! Sum it in the divergence: it is already in the right Y_lm form
        aux(1:i%m) = aux(1:i%m)*g(i%t)%rm2(1:i%m)
        !
        div_F_lm(1:i%m,lm,is) = div_F_lm(1:i%m,lm,is) + aux(1:i%m)
    ENDDO
    ENDDO

    OPTIONAL_CALL stop_clock ('PAW_div')

END SUBROUTINE PAW_divergence
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! build gradient of radial charge distribution from its spherical harmonics expansion
SUBROUTINE PAW_gradient(i, ix, rho_lm, rho_rad, rho_core, grho_rad2, grho_rad)
    USE constants,              ONLY : fpi
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid

    INTEGER, INTENT(IN)  :: ix ! line of the dylm2 matrix to use actually it is
                               ! one of the nx spherical integration directions
    TYPE(paw_info)  :: i                              ! atom's minimal info
    REAL(DP), INTENT(IN) :: rho_lm(i%m,i%l**2,nspin)! Y_lm expansion of rho
    REAL(DP), INTENT(IN) :: rho_rad(i%m,nspin)        ! radial density along direction ix
    REAL(DP), INTENT(IN) :: rho_core(i%m)             ! core density
    REAL(DP), INTENT(OUT):: grho_rad2(i%m,nspin)      ! |grad(rho)|^2 on rad. grid
    REAL(DP), OPTIONAL,INTENT(OUT):: grho_rad(i%m,3,nspin) ! vector gradient (only for gcxc)
    !              r, theta and phi components ---^
    !
    REAL(DP)             :: aux(i%m),aux2(i%m)       ! workspace
    INTEGER              :: is, lm ! counters on: spin, angular momentum

    OPTIONAL_CALL start_clock ('PAW_grad')
    ! 1. build real charge density = rho/r**2 + rho_core
    ! 2. compute the partial derivative of rho_rad
    grho_rad2(:,:) = 0._dp
    DO is = 1,nspin
        ! build real charge density
        aux(1:i%m) = rho_rad(1:i%m,is)*g(i%t)%rm2(1:i%m) &
                          + rho_core(1:i%m)/nspin
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
        DO lm = 2,i%l**2
            ! 5. [ \sum_{lm} rho(r) (dY_{lm}/dphi /cos(theta))  ]**2
            aux(1:i%m) = aux(1:i%m) + rad(i%t)%dylmp(ix,lm)* (rho_lm(1:i%m,lm,is))
            ! 6. [ \sum_{lm} rho(r) (dY_{lm}/dtheta)  ]**2
            aux2(1:i%m) = aux2(1:i%m) + rad(i%t)%dylmt(ix,lm)* (rho_lm(1:i%m,lm,is))
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
    USE constants,              ONLY : fpi, e2
    USE radial_grids,           ONLY : hartree
    USE ions_base,              ONLY : ityp
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid

    TYPE(paw_info)  :: i                          ! atom's minimal info
    ! charge density as lm components already summed on spin:
    REAL(DP), INTENT(IN)  :: rho_lm(i%m,i%l**2,nspin)
    REAL(DP), INTENT(OUT) :: v_lm  (i%m,i%l**2) ! potential as lm components
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
    DO lm = 1, i%l**2
        l = INT(sqrt(DBLE(lm-1))) ! l has to start from *zero*
            pref = e2*fpi/REAL(2*l+1)
            DO k = 1, i%m
                aux(k) = pref * SUM(rho_lm(k,lm,1:nspin))
            ENDDO
            !write(*,*) "yadda--_", l, 2*l+2, i%m
            CALL hartree(l, 2*l+2, i%m, g(i%t), aux(:), v_lm(:,lm))
    ENDDO

    ! compute energy if required:
    ! E_h = \sum_lm \int v_lm(r) (rho_lm(r) r^2) dr
    IF(present(energy)) THEN
    energy = 0._dp
    DO lm = 1, i%l**2
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
    USE ions_base,         ONLY : nat 
    USE lsda_mod,          ONLY : nspin
    USE uspp_param,        ONLY : nh, nhm
    USE uspp,              ONLY : indv, ap, nhtolm,lpl,lpx
    USE constants,         ONLY : eps12
    USE atom,              ONLY : g => rgrid

    TYPE(paw_info)  :: i                                    ! atom's minimal info
    REAL(DP), INTENT(IN)  :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupation
    REAL(DP), INTENT(IN)  :: pfunc(i%m,i%b,i%b)             ! psi_i * psi_j
    REAL(DP), INTENT(OUT) :: rho_lm(i%m,i%l**2,nspin)       ! AE charge density on rad. grid
    REAL(DP), OPTIONAL,INTENT(IN) :: &
                             aug(i%m,i%b,i%b,0:2*i%l) ! augmentation functions (only for PS part)

    REAL(DP)                :: pref ! workspace (ap*becsum)

    INTEGER                 :: nb,mb, &     ! counters for pfunc nb,mb = 1, nh
                               nmb, &       ! composite "triangular" index for pfunc nmb = 1,nh*(nh+1)/2
                               lm,lp,l, &   ! counters for angular momentum lm = l**2+m
                               ispin        ! counter for spin (FIXME: may be unnecessary)

    ! This subroutine computes the angular momentum components of rho
    ! using the following formula:
    !   rho(\vec{r}) = \sum_{LM} Y_{LM} \sum_{i,j} (\hat{r}) a_{LM}^{(lm)_i(lm)_j} becsum_ij pfunc_ij(r)
    ! where a_{LM}^{(lm)_i(lm)_j} are the Clebsh-Gordan coefficients.
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

!     write(*,*) i
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
!                 write(*,"(2i4,i6,i4,i6,i4,f12.6)") nb,mb,indv(nb,i%t), indv(mb,i%t),lp,lm,&
!                                        ap(lm, nhtolm(nb,i%t), nhtolm(mb,i%t))
                ! becsum already contains a factor 2 for off-diagonal pfuncs
                pref = becsum(nmb,i%a,ispin) * ap(lm, nhtolm(nb,i%t), nhtolm(mb,i%t))
                !
                rho_lm(1:i%m,lm,ispin) = rho_lm(1:i%m,lm,ispin) +&
                                pref * pfunc(1:i%m, indv(nb,i%t), indv(mb,i%t))
                IF (present(aug)) THEN
                    ! if I'm doing the pseudo part I have to add the augmentation charge
                    l = INT(SQRT(DBLE(lm-1))) ! l has to start from zero
                    rho_lm(1:i%m,lm,ispin) = rho_lm(1:i%m,lm,ispin) +&
                                pref * aug(1:i%m, indv(nb,i%t), indv(mb,i%t), l)
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
    USE lsda_mod,               ONLY : nspin

    TYPE(paw_info)              :: i  ! atom's minimal info
    INTEGER                     :: ix ! line of the ylm matrix to use
                                      ! actually it is one of the nx directions
    REAL(DP), INTENT(IN)        :: F_lm(i%m,i%l**2,nspin)! Y_lm expansion of rho
    REAL(DP), INTENT(OUT)       :: F_rad(i%m,nspin)      ! charge density on rad. grid
    !
    INTEGER                     :: ispin, lm ! counters on angmom and spin

    OPTIONAL_CALL start_clock ('PAW_lm2rad')
    F_rad(:,:) = 0._dp
    ! cycling on spin is a bit less general...
    spins: DO ispin = 1,nspin
        DO lm = 1, i%l**2
            F_rad(:,ispin) = F_rad(:,ispin) +&
                    rad(i%t)%ylm(ix,lm)*F_lm(:,lm,ispin)
        ENDDO ! lm
    ENDDO spins

    OPTIONAL_CALL stop_clock ('PAW_lm2rad')

END SUBROUTINE PAW_lm2rad

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! computes F_lm(r) = \int d \Omega F(r,th,ph) Y_lm(th,ph)
SUBROUTINE PAW_rad2lm(i, F_rad, F_lm, lmax_loc)
    USE lsda_mod,               ONLY : nspin

    TYPE(paw_info)       :: i        ! atom's minimal info
    REAL(DP), INTENT(OUT):: F_lm(i%m, lmax_loc**2, nspin) ! lm component of F up to lmax_loc
    REAL(DP), INTENT(IN) :: F_rad(i%m, rad(i%t)%nx, nspin)! radial samples of F
    INTEGER,  INTENT(IN) :: lmax_loc ! in some cases I have to keep higher angular components
                                     ! than the default ones (=lmaxq =the ones present in rho)
    !
    INTEGER                     :: ix    ! counter for integration
    INTEGER                     :: lm    ! counter for angmom
    INTEGER                     :: ispin ! counter for spin

    OPTIONAL_CALL start_clock ('PAW_rad2lm')
    F_lm(:,:,:) = 0._dp

    DO ispin = 1,nspin
    DO lm = 1,lmax_loc**2
    DO ix = 1, rad(i%t)%nx
        F_lm(1:i%m, lm, ispin) = F_lm(1:i%m, lm, ispin) &
                                + F_rad(1:i%m,ix,ispin)* rad(i%t)%wwylm(ix,lm)
    ENDDO
    ENDDO
    ENDDO

    OPTIONAL_CALL stop_clock ('PAW_rad2lm')

END SUBROUTINE PAW_rad2lm
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! computes F_lm(r) = \int d \Omega F(r,th,ph) Y_lm(th,ph)
!!! duplicated version to work on vector fields, necessary for performance reasons
SUBROUTINE PAW_rad2lm3(i, F_rad, F_lm, lmax_loc)
    USE lsda_mod,               ONLY : nspin

    TYPE(paw_info)       :: i        ! atom's minimal info
    REAL(DP), INTENT(OUT):: F_lm(i%m, 3, lmax_loc**2, nspin) ! lm component of F up to lmax_loc
    REAL(DP), INTENT(IN) :: F_rad(i%m, 3, rad(i%t)%nx, nspin)! radial samples of F
    INTEGER,  INTENT(IN) :: lmax_loc ! in some cases I have to keep higher angular components
                                     ! than the default ones (=lmaxq =the ones present in rho)
    !
    REAL(DP)             :: aux(i%m) ! optimization
    INTEGER              :: ix    ! counter for integration
    INTEGER              :: lm    ! counter for angmom
    INTEGER              :: ispin ! counter for spin

    OPTIONAL_CALL start_clock ('PAW_rad2lm3')

    ! Third try: 50% faster than blind implementation (60% with prefetch)
    DO ispin = 1,nspin
    DO lm = 1,lmax_loc**2
      aux(:) = 0._dp
      DO ix = 1, rad(i%t)%nx
            aux(1:i%m) = aux(1:i%m) + F_rad(1:i%m,1,ix,ispin) * rad(i%t)%wwylm(ix,lm)
            !CALL MM_PREFETCH( F_rad(1:i%m,1,MIN(ix+1,rad(i%t)%nx),ispin), 1 )
      ENDDO
      F_lm(1:i%m, 1, lm, ispin) = aux(1:i%m)
      !
      aux(:) = 0._dp
      DO ix = 1, rad(i%t)%nx
            aux(1:i%m) = aux(1:i%m) + F_rad(1:i%m,2,ix,ispin) * rad(i%t)%wwylm(ix,lm)
            !CALL MM_PREFETCH( F_rad(1:i%m,2,MIN(ix+1,rad(i%t)%nx),ispin), 1 )
      ENDDO
      F_lm(1:i%m, 2, lm, ispin) = aux(1:i%m)
      !
      aux(:) = 0._dp
      DO ix = 1, rad(i%t)%nx
            aux(1:i%m) = aux(1:i%m) + F_rad(1:i%m,3,ix,ispin) * rad(i%t)%wwylm(ix,lm)
            !CALL MM_PREFETCH( F_rad(1:i%m,3,MIN(ix+1,rad(i%t)%nx),ispin), 1 )
      ENDDO
      F_lm(1:i%m, 3, lm, ispin) = aux(1:i%m)

    ENDDO
    ENDDO

    OPTIONAL_CALL stop_clock ('PAW_rad2lm3')

END SUBROUTINE PAW_rad2lm3



END MODULE paw_onecenter
