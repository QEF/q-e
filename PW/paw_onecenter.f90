!
! Copyright (C) 2007-2009 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! NOTE ON PARALLELIZATION:
! this code is parallelized on atoms, i.e. each node computes potential, energy,
! newd coefficients, ddots and \int v \times n on a reduced number of atoms.
! The implementation assumes that divisions of atoms among the nodes is always
! done in the same way! By doing so we can avoid to allocate the potential for 
! all the atoms on all the nodes, and (most important) we don't need to
!  distribute the potential among the nodes after computing it.
!
#include "f_defs.h"
MODULE paw_onecenter
    !
    USE kinds,          ONLY : DP
    USE paw_variables,  ONLY : paw_info, rad, radial_grad_style
    USE mp_global,      ONLY : nproc_image, me_image, intra_image_comm
    USE mp,             ONLY : mp_sum
    !
    IMPLICIT NONE

    ! entry points:
    PUBLIC :: PAW_potential  ! prepare paw potential and store it,
                             ! also computes energy if required
    PUBLIC :: PAW_ddot       ! error estimate for mix_rho
    PUBLIC :: PAW_symmetrize ! symmetrize becsums
    PUBLIC :: PAW_desymmetrize! symmetrize dbecsums for electric field
    PUBLIC :: PAW_dusymmetrize! symmetrize dbecsums for phonon modes
    PUBLIC :: PAW_dumqsymmetrize! symmetrize dbecsums for phonon modes
                             ! with respect to minus_q
    PUBLIC :: PAW_dpotential ! calculate change of the paw potential 
                             ! and derivatives of D^1-~D^1 coefficients
    PUBLIC :: PAW_rho_lm     ! uses becsum to generate one-center charges
                             ! (all-electron and pseudo) on radial grid
    !
    INTEGER, SAVE :: paw_comm, me_paw, nproc_paw
    !
    PRIVATE
    !
    ! the following macro controls the use of several fine-grained clocks
    ! set it to 'if(.false.) CALL' (without quotes) in order to disable them,
    ! set it to 'CALL' to enable them.
    !
#define OPTIONAL_CALL if(.false.) CALL
    !
    INTEGER, EXTERNAL :: ldim_block
    INTEGER, EXTERNAL :: gind_block

 CONTAINS

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! Computes V_h and V_xc using the "density" becsum provided and then
!!! 
!!! Update the descreening coefficients:
!!! D_ij = \int v_Hxc p_ij - \int vt_Hxc (pt_ij + augfun_ij)
!!!
!!! calculate the onecenter contribution to the energy
!!!
SUBROUTINE PAW_potential(becsum, d, energy, e_cmp)
   USE atom,              ONLY : g => rgrid
   USE ions_base,         ONLY : nat, ityp
   USE lsda_mod,          ONLY : nspin
   USE uspp_param,        ONLY : nh, nhm, upf
   USE mp,                ONLY : mp_barrier, mp_comm_split, mp_comm_free, mp_size, mp_rank

   REAL(DP), INTENT(IN)  :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations
   REAL(DP), INTENT(OUT) :: d(nhm*(nhm+1)/2,nat,nspin) ! descreening coefficients (AE - PS)
   REAL(DP), INTENT(OUT), OPTIONAL :: energy           ! if present compute E[rho]
   REAL(DP), INTENT(OUT), OPTIONAL :: e_cmp(nat, 2, 2) ! components of the energy
   !                                          {AE!PS}
   INTEGER, PARAMETER      :: AE = 1, PS = 2,&      ! All-Electron and Pseudo
                              XC = 1, H  = 2        ! XC and Hartree
   REAL(DP), POINTER       :: rho_core(:)           ! pointer to AE/PS core charge density 
   TYPE(paw_info)          :: i                     ! minimal info on atoms
   INTEGER                 :: i_what                ! counter on AE and PS
   INTEGER                 :: is                    ! spin index
   INTEGER                 :: lm                    ! counters on angmom and radial grid
   INTEGER                 :: nb, mb, nmb           ! augfun indexes
   INTEGER                 :: ia,ia_s,ia_e          ! atoms counters and indexes
   INTEGER                 :: mykey                 ! my index in the atom group
   INTEGER                 :: j, l2, kkbeta, imesh
   !
   REAL(DP), ALLOCATABLE   :: v_lm(:,:,:)   ! workspace: potential
   REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:) ! density expanded on Y_lm
   REAL(DP), ALLOCATABLE   :: savedv_lm(:,:,:)   ! workspace: potential
   ! fake cross band occupations to select only one pfunc at a time:
   REAL(DP)                :: becfake(nhm*(nhm+1)/2,nat,nspin)
   REAL(DP)                :: integral           ! workspace
   REAL(DP)                :: energy_xc, energy_h, energy_tot
   REAL(DP)                :: sgn                ! +1 for AE -1 for PS

   CALL start_clock('PAW_pot')
   ! Some initialization
   becfake(:,:,:) = 0._dp
   d(:,:,:) = 0._dp
   energy_tot = 0._dp
   !
   !
   ! Parallel: divide tasks among all the processor for this image
   ! (i.e. all the processors except for NEB and similar)
   !
   CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )
   !
   ! build the group of all the procs associated with the same atom
   !
   CALL mp_comm_split( intra_image_comm, ia_s - 1, me_image, paw_comm )
   !
   me_paw    = mp_rank( paw_comm )
   nproc_paw = mp_size( paw_comm )
   !
   atoms: DO ia = ia_s, ia_e
      !
      i%a = ia                      ! atom's index
      i%t = ityp(ia)                ! type of atom ia
      i%m = g(i%t)%mesh             ! radial mesh size for atom i%t
      i%b = upf(i%t)%nbeta          ! number of beta functions for i%t
      i%l = upf(i%t)%lmax_rho+1 ! max ang.mom. in augmentation for ia
      l2  = i%l**2
      kkbeta = upf(i%t)%kkbeta
      imesh  = i%m
      !
      ifpaw: IF (upf(i%t)%tpawp) THEN
         !
         ! Arrays are allocated inside the cycle to allow reduced
         ! memory usage as differnt atoms have different meshes
         ALLOCATE(v_lm(i%m,l2,nspin))
         ALLOCATE(savedv_lm(i%m,l2,nspin))
         ALLOCATE(rho_lm(i%m,l2,nspin))
         !
         whattodo: DO i_what = AE, PS
            ! STEP: 1 [ build rho_lm (PAW_rho_lm) ]
            NULLIFY(rho_core)
            IF (i_what == AE) THEN
               ! Compute rho spherical harmonics expansion from becsum and pfunc
               CALL PAW_rho_lm(i, becsum, upf(i%t)%paw%pfunc, rho_lm)
               ! used later for xc potential:
               rho_core => upf(i%t)%paw%ae_rho_atc
               ! sign to sum up the enrgy
               sgn = +1._dp
            ELSE
               CALL PAW_rho_lm(i, becsum, upf(i%t)%paw%ptfunc, rho_lm, upf(i%t)%qfuncl)
               !          optional argument for pseudo part (aug. charge) --> ^^^
               rho_core => upf(i%t)%rho_atc ! as before
               sgn = -1._dp                 ! as before
            ENDIF
            ! cleanup auxiliary potentials
            savedv_lm(:,:,:) = 0._dp

            ! First compute the Hartree potential (it does not depend on spin...):
            CALL PAW_h_potential(i, rho_lm, v_lm(:,:,1), energy)
      ! NOTE: optional variables works recursively: e.g. if energy is not present here
            ! it will not be present in PAW_h_potential too!
            !IF (present(energy)) write(*,*) 'H',i%a,i_what,sgn*energy
            IF (present(energy) .AND. mykey == 0 ) energy_tot = energy_tot + sgn*energy
            IF (present(e_cmp) .AND. mykey == 0 ) e_cmp(ia, H, i_what) = energy
            DO is = 1,nspin ! ... v_H has to be copied to all spin components
               savedv_lm(:,:,is) = v_lm(:,:,1)
            ENDDO

            ! Then the XC one:
            CALL PAW_xc_potential(i, rho_lm, rho_core, v_lm, energy)
            !IF (present(energy)) write(*,*) 'X',i%a,i_what,sgn*energy
            IF (present(energy) .AND. mykey == 0 ) energy_tot = energy_tot + sgn*energy
            IF (present(e_cmp) .AND. mykey == 0 )  e_cmp(ia, XC, i_what) = energy
            savedv_lm(:,:,:) = savedv_lm(:,:,:) + v_lm(:,:,:)
            !
            spins: DO is = 1, nspin
               nmb = 0
               ! loop on all pfunc for this kind of pseudo
               DO nb = 1, nh(i%t)
                  DO mb = nb, nh(i%t)
                     nmb = nmb+1 ! nmb = 1, nh*(nh+1)/2
                     !
                     ! compute the density from a single pfunc
                     becfake(nmb,ia,is) = 1._dp
                     IF (i_what == AE) THEN
                        CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%pfunc, rho_lm)
                     ELSE
                        CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%ptfunc, rho_lm, upf(i%t)%qfuncl)
                        !                  optional argument for pseudo part --> ^^^
                     ENDIF
                     !
                     ! Now I multiply the rho_lm and the potential, I can use
                     ! rho_lm itself as workspace
                     DO lm = 1, l2
                        DO j = 1, imesh
                           rho_lm(j,lm,is) = rho_lm(j,lm,is) * savedv_lm(j,lm,is)
                        END DO
                        ! Integrate!
                        CALL simpson ( kkbeta, rho_lm(1,lm,is), g(i%t)%rab(1), integral)
                        d(nmb,i%a,is) = d(nmb,i%a,is) + sgn * integral
                     ENDDO
                     ! restore becfake to zero
                     becfake(nmb,ia,is) = 0._dp
                  ENDDO ! mb
               ENDDO ! nb
            ENDDO spins
         ENDDO whattodo
         ! cleanup
         DEALLOCATE(rho_lm)
         DEALLOCATE(savedv_lm)
         DEALLOCATE(v_lm)
         !
      ENDIF ifpaw
   ENDDO atoms

#ifdef __PARA
   ! recollect D coeffs and total one-center energy
   IF( mykey /= 0 ) energy_tot = 0.0d0
   CALL mp_sum(energy_tot, intra_image_comm)
   IF( mykey /= 0 ) d = 0.0d0
   CALL mp_sum(d, intra_image_comm)
#endif
   ! put energy back in the output variable
   IF ( present(energy) ) energy = energy_tot
   !
   CALL mp_comm_free( paw_comm )
   !
   CALL stop_clock('PAW_pot')

END SUBROUTINE PAW_potential

SUBROUTINE PAW_symmetrize(becsum)
    USE lsda_mod,          ONLY : nspin
    USE uspp_param,        ONLY : nhm
    USE ions_base,         ONLY : nat, ityp
    USE symme,             ONLY : nsym, irt, d1, d2, d3
    USE uspp,              ONLY : nhtolm,nhtol,ijtoh
    USE uspp_param,        ONLY : nh, upf
    USE io_global,         ONLY : stdout, ionode

    REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations

    REAL(DP)                :: becsym(nhm*(nhm+1)/2,nat,nspin)! symmetrized becsum
    REAL(DP) :: pref, usym

    INTEGER :: ia,mykey,ia_s,ia_e 
                            ! atoms counters and indexes
    INTEGER :: is, nt       ! counters on spin, atom-type
    INTEGER :: ma           ! atom symmetric to na
    INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
               l_i, l_j, m_i, m_j
    INTEGER :: m_o, m_u     ! counters for sums on m
    INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    INTEGER :: isym         ! counter for symmetry operation

    ! The following mess is necessary because the symmetrization operation
    ! in LDA+U code is simpler than in PAW, so the required quantities are
    ! represented in a simple but not general way.
    ! I will fix this when everything works.
    REAL(DP), TARGET :: d0(1,1,48)
    TYPE symmetryzation_tensor
        REAL(DP),POINTER :: d(:,:,:)
    END TYPE symmetryzation_tensor
    TYPE(symmetryzation_tensor) :: D(0:3)

    IF( nsym==1 ) RETURN
    d0(1,1,:) = 1._dp
    D(0)%d => d0 ! d0(1,1,48)
    D(1)%d => d1 ! d1(3,3,48)
    D(2)%d => d2 ! d2(5,5,48)
    D(3)%d => d3 ! d3(7,7,48)

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


    CALL start_clock('PAW_symme')

    becsym(:,:,:) = 0._dp
    usym = 1._dp / REAL(nsym)

    ! Parallel: divide among processors for the same image
    CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )

    DO is = 1, nspin
    !
    atoms: DO ia = ia_s, ia_e
        nt = ityp(ia)
        ! No need to symmetrize non-PAW atoms
        IF ( .not. upf(nt)%tpawp ) CYCLE
        !
        DO ih = 1, nh(nt)
        DO jh = ih, nh(nt) ! note: jh >= ih
            !ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
            ijh = ijtoh(ih,jh,nt)
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
                ma = irt(isym,ia)
                DO m_o = 1, 2*l_i +1
                DO m_u = 1, 2*l_j +1
                    oh = ih - m_i + m_o
                    uh = jh - m_j + m_u
                    ouh = ijtoh(oh,uh,nt)
                    ! In becsum off-diagonal terms are multiplied by 2, I have
                    ! to neutralize this factor and restore it later
                    IF ( oh == uh ) THEN
                        pref = 2._dp * usym
                    ELSE
                        pref = usym
                    ENDIF
                    !
                    becsym(ijh, ia, is) = becsym(ijh, ia, is) &
                        + D(l_i)%d(m_o,m_i, isym) * D(l_j)%d(m_u,m_j, isym) &
                          * pref * becsum(ouh, ma, is)
                ENDDO ! m_o
                ENDDO ! m_u
            ENDDO ! isym
            !
            ! Put the prefactor back in:
            IF ( ih == jh ) becsym(ijh,ia,is) = .5_dp * becsym(ijh,ia,is)
        ENDDO ! ih
        ENDDO ! jh
    ENDDO atoms ! nat
    ENDDO ! nspin
#ifdef __PARA
    IF( mykey /= 0 ) becsym = 0.0_dp
    CALL mp_sum(becsym, intra_image_comm)
#endif

#ifdef __DEBUG_PAW_SYM
   write(stdout,*) "------------"
    if(ionode) then
        ia = 1
        nt = ityp(ia)
        DO is = 1, nspin
            write(*,*) is
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            ijh = ijtoh(ih,jh,nt)
            write(stdout,"(1f10.3)", advance='no') becsym(ijh,ia,is)
        ENDDO
            write(stdout,*)
        ENDDO
            write(stdout,*)
        ENDDO
    endif
   write(stdout,*) "------------"
#endif

    ! Apply symmetrization:
    becsum(:,:,:) = becsym(:,:,:)

    CALL stop_clock('PAW_symme')

END SUBROUTINE PAW_symmetrize

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
    INTEGER                 :: ia,mykey,ia_s,ia_e  
                                                     ! atoms counters and indexes
    INTEGER                 :: lm,k                  ! counters on angmom and radial grid

    ! hartree energy scalar fields expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:) ! radial density expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: v_lm(:,:)     ! hartree potential, summed on spins (from bec1)
    !
    REAL(DP)                :: i_sign        ! +1 for AE, -1 for PS
    REAL(DP)                :: integral      ! workspace
    TYPE(paw_info)          :: i

    CALL start_clock ('PAW_ddot')
    ! initialize 
    PAW_ddot = 0._dp

    ! Parallel: divide among processors for the same image
    CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )
    !
    atoms: DO ia = ia_s, ia_e
    !
    i%a = ia           ! the index of the atom
    i%t = ityp(ia)    ! the type of atom ia
    i%m = g(i%t)%mesh ! radial mesh size for atom ia
    i%b = upf(i%t)%nbeta
    i%l = upf(i%t)%lmax_rho+1
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
                CALL PAW_rho_lm(i, bec1, upf(i%t)%paw%ptfunc, rho_lm, upf(i%t)%qfuncl)
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
                CALL PAW_rho_lm(i, bec2, upf(i%t)%paw%ptfunc, rho_lm, upf(i%t)%qfuncl)
            ENDIF
            !
            ! Finally compute the integral
            DO lm = 1, i%l**2
                ! I can use v_lm as workspace
                DO k = 1, i%m
                    v_lm(k,lm) = v_lm(k,lm) * SUM(rho_lm(k,lm,1:nspin))
                ENDDO
                CALL simpson (upf(i%t)%kkbeta,v_lm(:,lm),g(i%t)%rab,integral)
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
    IF( mykey /= 0 ) PAW_ddot = 0.0_dp
    CALL mp_sum(PAW_ddot, intra_image_comm)
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
    USE funct,                  ONLY : dft_is_gradient, evxc_t_vec
    USE constants,              ONLY : fpi ! REMOVE

    TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
    REAL(DP), INTENT(IN)  :: rho_lm(i%m,i%l**2,nspin)! charge density as lm components
    REAL(DP), INTENT(IN)  :: rho_core(i%m)           ! core charge, radial and spherical
    REAL(DP), INTENT(OUT) :: v_lm(i%m,i%l**2,nspin)  ! potential density as lm components
    REAL(DP),OPTIONAL,INTENT(OUT) :: energy          ! XC energy (if required)
    !
    REAL(DP), ALLOCATABLE :: rho_loc(:,:)         ! local density (workspace), up and down
    REAL(DP)              :: v_rad(i%m,rad(i%t)%nx,nspin)! radial potential (to be integrated)
    REAL(DP), ALLOCATABLE :: rho_rad(:,:)       ! workspace (only one radial slice of rho)
    !
    REAL(DP), ALLOCATABLE :: e_rad(:)           ! aux, used to store radial slices of energy
    REAL(DP), ALLOCATABLE :: e_of_tid(:)        ! aux, for openmp parallel reduce
    REAL(DP)              :: e                  ! aux, used to integrate energy
    !
    INTEGER               :: ix,k               ! counters on directions and radial grid
    INTEGER               :: lsd                ! switch for local spin density

    REAL(DP)              :: exc_ret, stmp
    !
    INTEGER               :: nx_loc, ix_s, ix_e
    INTEGER               :: mytid, ntids
#ifdef __OPENMP
    INTEGER, EXTERNAL     :: omp_get_thread_num, omp_get_num_threads
#endif

    OPTIONAL_CALL start_clock ('PAW_xc_pot')
    !
    ! true if using spin
    lsd = nspin-1
    !
    nx_loc = ldim_block( rad(i%t)%nx, nproc_paw, me_paw )
    ix_s   = gind_block( 1, rad(i%t)%nx, nproc_paw, me_paw )
    ix_e   = ix_s + nx_loc - 1
    !
!$omp parallel default(private), &
!$omp shared(i,rad,v_lm,rho_lm,rho_core,v_rad,ix_s,ix_e,energy,e_of_tid,nspin,g,lsd)
#ifdef __OPENMP
    mytid = omp_get_thread_num()+1 ! take the thread ID
    ntids = omp_get_num_threads()  ! take the number of threads
#else
    mytid = 1
    ntids = 1
#endif
    ! This will hold the "true" charge density, without r**2 or other factors
    ALLOCATE( rho_loc(i%m,2) ) 
    rho_loc = 0._dp
    !
    ALLOCATE( rho_rad(i%m,nspin) ) 
    !
    IF (present(energy)) THEN
!$omp single
        energy = 0._dp
        ALLOCATE(e_of_tid(ntids))
!$omp end single
        ALLOCATE(e_rad(i%m))
        e_of_tid(mytid) = 0._dp
    ENDIF

!$omp workshare
    v_rad = 0.0_dp
!$omp end workshare

!$omp do
    DO ix = ix_s, ix_e
        !
        ! *** LDA (and LSDA) part (no gradient correction) ***
        ! convert _lm density to real density along ix
        !
        CALL PAW_lm2rad(i, ix, rho_lm, rho_rad)
        !
        ! compute the potential along ix
        !
        IF( nspin < 2 ) THEN
           DO k = 1,i%m
              rho_loc(k,1) = rho_rad(k,1)*g(i%t)%rm2(k)
           ENDDO
        ELSE
           DO k = 1,i%m
              rho_loc(k,1) = rho_rad(k,1)*g(i%t)%rm2(k)
              rho_loc(k,2) = rho_rad(k,2)*g(i%t)%rm2(k)
           ENDDO
        END IF
        !
        ! Integrate to obtain the energy
        !
        IF (present(energy)) THEN
           CALL evxc_t_vec(rho_loc, rho_core, lsd, i%m, v_rad(:,ix,:), e_rad)
           IF( nspin < 2 ) THEN
              e_rad = e_rad * ( rho_rad(:,1) + rho_core*g(i%t)%r2 )
           ELSE
              e_rad = e_rad * ( rho_rad(:,1) + rho_rad(:,2) + rho_core*g(i%t)%r2 )
           END IF
           ! Integrate to obtain the energy
           CALL simpson(i%m, e_rad, g(i%t)%rab, e)
           e_of_tid(mytid) = e_of_tid(mytid) + e * rad(i%t)%ww(ix)
        ELSE
           CALL evxc_t_vec(rho_loc, rho_core, lsd, i%m, v_rad(:,ix,:))
        ENDIF
    ENDDO
!$omp end do nowait

    IF(present(energy)) THEN
       DEALLOCATE(e_rad)
    END IF

    DEALLOCATE( rho_rad ) 
    DEALLOCATE( rho_loc ) 

!$omp end parallel

    CALL mp_sum( v_rad, paw_comm )

    IF(present(energy)) THEN
       energy = sum(e_of_tid)
       DEALLOCATE(e_of_tid)
       CALL mp_sum( energy, paw_comm )
    END IF

    ! Recompose the sph. harm. expansion
    CALL PAW_rad2lm(i, v_rad, v_lm, i%l)

    ! Add gradient correction, if necessary
    IF( dft_is_gradient() ) &
        CALL PAW_gcxc_potential( i, rho_lm, rho_core, v_lm, energy )


    OPTIONAL_CALL stop_clock ('PAW_xc_pot')

END SUBROUTINE PAW_xc_potential
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! add gradient correction to v_xc, code mostly adapted from ../atomic/vxcgc.f90
!!! in order to support non-spherical charges (as Y_lm expansion)
!!! Note that the first derivative in vxcgc becames a gradient, while the second is a divergence.
!!! We also have to temporary store some additional Y_lm components in order not to loose
!!! precision during the calculation, even if only the ones up to lmax_rho (the maximum in the
!!! density of charge) matter when computing \int v * rho 
SUBROUTINE PAW_gcxc_potential(i, rho_lm,rho_core, v_lm, energy)
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    USE constants,              ONLY : sqrtpi, fpi,pi,e2, eps => eps12, eps2 => eps24
    USE funct,                  ONLY : gcxc, gcx_spin_vec, gcc_spin, gcx_spin
    USE mp,                     ONLY : mp_sum
    !
    TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
    REAL(DP), INTENT(IN)    :: rho_lm(i%m,i%l**2,nspin) ! charge density as lm components
    REAL(DP), INTENT(IN)    :: rho_core(i%m)            ! core charge, radial and spherical
    REAL(DP), INTENT(INOUT) :: v_lm(i%m,i%l**2,nspin)   ! potential to be updated
    REAL(DP),OPTIONAL,INTENT(INOUT) :: energy           ! if present, add GC to energy

    REAL(DP),ALLOCATABLE    :: rho_rad(:,:)! charge density sampled
    REAL(DP),ALLOCATABLE    :: grad(:,:,:) ! gradient
    REAL(DP),ALLOCATABLE    :: grad2(:,:)  ! square modulus of gradient
                                                             ! (first of charge, than of hamiltonian)
    REAL(DP),ALLOCATABLE    :: gc_rad(:,:,:) ! GC correction to V (radial samples)
    REAL(DP),ALLOCATABLE    :: gc_lm(:,:,:)       ! GC correction to V (Y_lm expansion)
    REAL(DP),ALLOCATABLE    :: h_rad(:,:,:,:)! hamiltonian (vector field)
    REAL(DP),ALLOCATABLE    :: h_lm(:,:,:,:)! hamiltonian (vector field)
                                       !!! ^^^^^^^^^^^^^^^^^^ expanded to higher lm than rho !!!
    REAL(DP),ALLOCATABLE    :: div_h(:,:,:)  ! div(hamiltonian)

    REAL(DP),ALLOCATABLE    :: e_rad(:)                   ! aux, used to store energy
    REAL(DP)                :: e, e_gcxc                  ! aux, used to integrate energy

    INTEGER  :: k, ix, is, lm                             ! counters on spin and mesh
    REAL(DP) :: sx,sc,v1x,v2x,v1c,v2c                     ! workspace
    REAL(DP) :: v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw  ! workspace
    REAL(DP) :: sgn, arho                                 ! workspace
    REAL(DP) :: rup, rdw, co2                             ! workspace
    REAL(DP) :: rh, zeta, grh2
    REAL(DP), ALLOCATABLE :: rup_vec(:), rdw_vec(:)
    REAL(DP), ALLOCATABLE :: sx_vec(:)
    REAL(DP), ALLOCATABLE :: v1xup_vec(:), v1xdw_vec(:)
    REAL(DP), ALLOCATABLE :: v2xup_vec(:), v2xdw_vec(:)

    
    INTEGER :: nx_loc, ix_s, ix_e
    INTEGER :: mytid, ntids
#ifdef __OPENMP
    INTEGER, EXTERNAL :: omp_get_thread_num, omp_get_num_threads
#endif
    REAL(DP),ALLOCATABLE :: egcxc_of_tid(:)


    OPTIONAL_CALL start_clock ('PAW_gcxc_v')


    nx_loc = ldim_block( rad(i%t)%nx, nproc_paw, me_paw )
    ix_s   = gind_block( 1, rad(i%t)%nx, nproc_paw, me_paw )
    ix_e   = ix_s + nx_loc - 1
    e_gcxc = 0._dp

    ALLOCATE( gc_rad(i%m,rad(i%t)%nx,nspin) )! GC correction to V (radial samples)
    ALLOCATE( gc_lm(i%m,i%l**2,nspin)       )! GC correction to V (Y_lm expansion)
    ALLOCATE( h_rad(i%m,3,rad(i%t)%nx,nspin))! hamiltonian (vector field)
    ALLOCATE( h_lm(i%m,3,(i%l+rad(i%t)%ladd)**2,nspin) ) 
                                       !!! ^^^^^^^^^^^^^^^^^^ expanded to higher lm than rho !!!
    ALLOCATE( div_h(i%m,i%l**2,nspin) )

!$omp parallel default(private) &
!$omp& shared(i,g,nspin,rad,e_gcxc,egcxc_of_tid,gc_rad,h_rad,rho_lm,rho_core,energy,ix_s,ix_e)
    mytid = 1
    ntids = 1
#ifdef __OPENMP
    mytid = omp_get_thread_num()+1 ! take the thread ID
    ntids = omp_get_num_threads()  ! take the number of threads
#endif
    ALLOCATE( rho_rad(i%m,nspin))! charge density sampled
    ALLOCATE( grad(i%m,3,nspin) )! gradient
    ALLOCATE( grad2(i%m,nspin)  )! square modulus of gradient
                                                             ! (first of charge, than of hamiltonian)
!$omp workshare
    gc_rad = 0.0d0
    h_rad  = 0.0d0
!$omp end workshare nowait

    IF (present(energy)) THEN
!$omp single
        allocate(egcxc_of_tid(ntids))
!$omp end single
        egcxc_of_tid(mytid) = 0.0_dp
        ALLOCATE(e_rad(i%m))
    ENDIF

    spin:&
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    IF ( nspin == 1 ) THEN
        !
        !     GGA case
        !
!omp do
        DO ix = ix_s, ix_e
           !
           !  WARNING: the next 2 calls are duplicated for spin==2
           CALL PAW_lm2rad(i, ix, rho_lm, rho_rad)
           CALL PAW_gradient(i, ix, rho_lm, rho_rad, rho_core, grad2, grad)

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
           !
           ! integrate energy (if required)
           IF (present(energy)) THEN
               CALL simpson(i%m, e_rad, g(i%t)%rab, e)
               egcxc_of_tid(mytid) = egcxc_of_tid(mytid) + e * rad(i%t)%ww(ix)
           ENDIF
        ENDDO
!omp end do
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ELSEIF ( nspin == 2 ) THEN
        ALLOCATE( rup_vec(i%m) )
        ALLOCATE( rdw_vec(i%m) )
        ALLOCATE( sx_vec(i%m) )
        ALLOCATE( v1xup_vec(i%m) )
        ALLOCATE( v1xdw_vec(i%m) )
        ALLOCATE( v2xup_vec(i%m) )
        ALLOCATE( v2xdw_vec(i%m) )
        !
        !   this is the \sigma-GGA case
        !
!$omp do
        DO ix = ix_s, ix_e
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
            rup_vec(k) = rho_rad(k,1)*g(i%t)%rm2(k) + co2
            rdw_vec(k) = rho_rad(k,2)*g(i%t)%rm2(k) + co2
        END DO
        ! bang!
        CALL gcx_spin_vec (rup_vec, rdw_vec, grad2(:,1), grad2(:,2), &
                sx_vec, v1xup_vec, v1xdw_vec, v2xup_vec, v2xdw_vec, i%m)
        DO k = 1,i%m
            rh = rup_vec(k) + rdw_vec(k) ! total charge
            IF ( rh > eps ) THEN
                zeta = (rup_vec(k) - rdw_vec(k) ) / rh ! polarization
                !
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
               e_rad(k)    = e2*(sx_vec(k)+sc)* g(i%t)%r2(k)
            gc_rad(k,ix,1)  = (v1xup_vec(k)+v1cup)!*g(i%t)%rm2(k)
            gc_rad(k,ix,2)  = (v1xdw_vec(k)+v1cdw)!*g(i%t)%rm2(k)
            !
            h_rad(k,:,ix,1) =( (v2xup_vec(k)+v2c)*grad(k,:,1)+v2c*grad(k,:,2) )*g(i%t)%r2(k)
            h_rad(k,:,ix,2) =( (v2xdw_vec(k)+v2c)*grad(k,:,2)+v2c*grad(k,:,1) )*g(i%t)%r2(k)
        ENDDO ! k
        ! integrate energy (if required)
        ! NOTE: this integration is duplicated for every spin, FIXME!
        IF (present(energy)) THEN
            CALL simpson(i%m, e_rad, g(i%t)%rab, e)
            egcxc_of_tid(mytid) = egcxc_of_tid(mytid) + e * rad(i%t)%ww(ix)
        ENDIF
        ENDDO ! ix
!$omp end do nowait
        DEALLOCATE( rup_vec )
        DEALLOCATE( rdw_vec )
        DEALLOCATE( sx_vec )
        DEALLOCATE( v1xup_vec )
        DEALLOCATE( v1xdw_vec )
        DEALLOCATE( v2xup_vec )
        DEALLOCATE( v2xdw_vec )
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ELSEIF ( nspin == 4 ) THEN
!$omp master
        CALL errore('PAW_gcxc_v', 'non-collinear not yet implemented!', 1)
!$omp end master
    ELSE spin
!$omp master
        CALL errore('PAW_gcxc_v', 'unknown spin number', 2)
!$omp end master
    ENDIF spin
    !
    IF (present(energy)) THEN
       DEALLOCATE(e_rad)
    ENDIF

    DEALLOCATE( rho_rad )
    DEALLOCATE( grad )
    DEALLOCATE( grad2 )
!$omp end parallel
	!
    CALL mp_sum( gc_rad, paw_comm )
    CALL mp_sum( h_rad, paw_comm )
    !
    IF (present(energy)) THEN
       e_gcxc = sum(egcxc_of_tid)
       CALL mp_sum( e_gcxc, paw_comm )
       energy = energy + e_gcxc
    ENDIF

    !
    ! convert the first part of the GC correction back to spherical harmonics
    CALL PAW_rad2lm(i, gc_rad, gc_lm, i%l)
    !
    ! We need the gradient of h to calculate the last part of the exchange
    ! and correlation potential. First we have to convert H to its Y_lm expansion
    CALL PAW_rad2lm3(i, h_rad, h_lm, i%l+rad(i%t)%ladd)
    !
    ! Compute div(H)
    CALL PAW_divergence(i, h_lm, div_h, i%l+rad(i%t)%ladd, i%l)
    !                         input max lm --^     ^-- output max lm

    ! Finally sum it back into v_xc
    DO is = 1,nspin
    DO lm = 1,i%l**2
            !v_lm(1:i%m,lm,is) = v_lm(1:i%m,lm,is) + e2*(gc_lm(1:i%m,lm,is)*g(i%t)%r2(1:i%m)-div_h(1:i%m,lm,is))
            v_lm(1:i%m,lm,is) = v_lm(1:i%m,lm,is) + e2*(gc_lm(1:i%m,lm,is)-div_h(1:i%m,lm,is))
    ENDDO
    ENDDO

    DEALLOCATE( gc_rad )
    DEALLOCATE( gc_lm )
    DEALLOCATE( h_rad )
    DEALLOCATE( h_lm )
    DEALLOCATE( div_h )


    !if(present(energy)) write(*,*) "gcxc -->", e_gcxc
    OPTIONAL_CALL stop_clock ('PAW_gcxc_v')

END SUBROUTINE PAW_gcxc_potential
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! compute divergence of a vector field (actutally the hamiltonian)
!!! it is assumed that: 1. the input function is multiplied by r**2; 
!!! 2. the output function is multiplied by r**2 too
SUBROUTINE PAW_divergence(i, F_lm, div_F_lm, lmaxq_in, lmaxq_out)
    USE constants,              ONLY : sqrtpi, fpi, e2
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid

    TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
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

!     CALL errore('PAW_divergence', 'More angular momentum components are needed (in input)'//&
!                 ' to provide the number you have requested (in output)', lmaxq_out-lmaxq_in+2)

    ! phi component
    DO is = 1,nspin
    DO ix = 1,rad(i%t)%nx
    aux(:) = 0._dp
        ! this derivative has no spherical component, so lm starts from 2
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
        CALL radial_gradient(F_lm(1:i%m,1,lm,is), aux, g(i%t)%r, i%m, radial_grad_style)
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
    TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
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
        CALL radial_gradient(aux, aux2, g(i%t)%r, i%m, radial_grad_style)
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
            aux(1:i%m) = aux(1:i%m) + rad(i%t)%dylmp(ix,lm)* rho_lm(1:i%m,lm,is)
            ! 6. [ \sum_{lm} rho(r) (dY_{lm}/dtheta)  ]**2
            aux2(1:i%m) = aux2(1:i%m) + rad(i%t)%dylmt(ix,lm)* rho_lm(1:i%m,lm,is)
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

    TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
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
    ! V_h(r) = \sum{lm} Y_{lm}(\hat{r})/(2l+1) \int dr' 4\pi r'^2 \rho^{lm}(r') (r1^l/r2^{l+1})
    !     done here --> ^^^^^^^^^^^^^^^^^^^^^           ^^^^^^^^^^^^^^^^^^^^^^ <-- input to the hartree subroutine
    !                 output from the h.s. --> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    DO lm = 1, i%l**2
        l = INT(sqrt(DBLE(lm-1))) ! l has to start from *zero*
            pref = e2*fpi/REAL(2*l+1)
            DO k = 1, i%m
                aux(k) = pref * SUM(rho_lm(k,lm,1:nspin))
            ENDDO
            !
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

    TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
    REAL(DP), INTENT(IN)  :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupation
    REAL(DP), INTENT(IN)  :: pfunc(i%m,i%b,i%b)             ! psi_i * psi_j
    REAL(DP), INTENT(OUT) :: rho_lm(i%m,i%l**2,nspin)       ! AE charge density on rad. grid
    REAL(DP), OPTIONAL,INTENT(IN) :: &
                             aug(i%m,i%b*(i%b+1)/2,0:2*i%l) ! augmentation functions (only for PS part)

    REAL(DP)                :: pref ! workspace (ap*becsum)

    INTEGER                 :: ih,jh, &     ! counters for pfunc ih,jh = 1, nh (CRYSTAL index)
                               nb,mb, &     ! counters for pfunc nb,mb = 1, nbeta (ATOMIC index)
                               ijh,nmb, &   ! composite "triangular" index for pfunc nmb = 1,nh*(nh+1)/2
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
    ! but the becsum depend on both l and m.

    OPTIONAL_CALL start_clock ('PAW_rho_lm')

    ! initialize density
    rho_lm(:,:,:) = 0._dp

    spins: DO ispin = 1, nspin
    ijh = 0
        ! loop on all pfunc for this kind of pseudo
        DO ih = 1, nh(i%t)
        DO jh = ih, nh(i%t)
            ijh = ijh+1
            nb = indv(ih,i%t)
            mb = indv(jh,i%t)
            nmb = mb * (mb-1)/2 + nb  ! mb has to be .ge. nb
            !write(*,'(99i4)') nb,mb,nmb
            IF (ABS(becsum(ijh,i%a,ispin)) < eps12) CYCLE
            !
            angular_momentum: &
            DO lp = 1, lpx (nhtolm(jh,i%t), nhtolm(ih,i%t)) !lmaxq**2
                ! the lpl array contains the possible combination of LM,lm_j,lm_j that
                ! have non-zero a_{LM}^{(lm)_i(lm)_j} (it saves some loops)
                lm = lpl (nhtolm(jh,i%t), nhtolm(ih,i%t), lp)
                !
                ! becsum already contains a factor 2 for off-diagonal pfuncs
                pref = becsum(ijh,i%a,ispin) * ap(lm, nhtolm(ih,i%t), nhtolm(jh,i%t))
                !
                rho_lm(1:i%m,lm,ispin) = rho_lm(1:i%m,lm,ispin) &
                                        +pref * pfunc(1:i%m, nb, mb)
                IF (present(aug)) THEN
                    ! if I'm doing the pseudo part I have to add the augmentation charge
                    l = INT(SQRT(DBLE(lm-1))) ! l has to start from zero, lm = l**2 +m
                    rho_lm(1:i%m,lm,ispin) = rho_lm(1:i%m,lm,ispin) &
                                            +pref * aug(1:i%m, nmb, l)
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

    TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
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

    TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
    INTEGER,  INTENT(IN) :: lmax_loc ! in some cases I have to keep higher angular components
                                     ! than the default ones (=lmaxq =the ones present in rho)
    REAL(DP), INTENT(OUT):: F_lm(i%m, lmax_loc**2, nspin) ! lm component of F up to lmax_loc
    REAL(DP), INTENT(IN) :: F_rad(i%m, rad(i%t)%nx, nspin)! radial samples of F
    !
    INTEGER                     :: ix    ! counter for integration
    INTEGER                     :: lm    ! counter for angmom
    INTEGER                     :: ispin ! counter for spin
    INTEGER                     :: j

    OPTIONAL_CALL start_clock ('PAW_rad2lm')

!$omp parallel default(shared), private(ispin,lm,ix,j)
    DO ispin = 1,nspin
!$omp do
    DO lm = 1,lmax_loc**2
    F_lm(:,lm,ispin) = 0._dp
    DO ix = 1, rad(i%t)%nx
    DO j  = 1, i%m
        F_lm(j, lm, ispin) = F_lm(j, lm, ispin) + F_rad(j,ix,ispin)* rad(i%t)%wwylm(ix,lm)
    ENDDO
    ENDDO
    ENDDO
!$omp end do
    ENDDO
!$omp end parallel

    OPTIONAL_CALL stop_clock ('PAW_rad2lm')

END SUBROUTINE PAW_rad2lm
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! computes F_lm(r) = \int d \Omega F(r,th,ph) Y_lm(th,ph)
!!! duplicated version to work on vector fields, necessary for performance reasons
SUBROUTINE PAW_rad2lm3(i, F_rad, F_lm, lmax_loc)
    USE lsda_mod,               ONLY : nspin

    TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
    INTEGER,  INTENT(IN) :: lmax_loc  ! in some cases I have to keep higher angular components
                                      ! than the default ones (=lmaxq =the ones present in rho)
    REAL(DP), INTENT(OUT):: F_lm(i%m, 3, lmax_loc**2, nspin) ! lm component of F up to lmax_loc
    REAL(DP), INTENT(IN) :: F_rad(i%m, 3, rad(i%t)%nx, nspin)! radial samples of F
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
!
! Computes dV_h and dV_xc using the "change of density" dbecsum provided 
! Update the change of the descreening coefficients:
! D_ij = \int dv_Hxc p_ij - \int dvt_Hxc (pt_ij + augfun_ij)
!
!
SUBROUTINE PAW_dpotential(dbecsum, becsum, int3, npe, max_irr_dim)
   USE atom,              ONLY : g => rgrid
   USE ions_base,         ONLY : nat, ityp
   USE lsda_mod,          ONLY : nspin
   USE uspp_param,        ONLY : nh, nhm, upf

   INTEGER, INTENT(IN) :: npe     ! number of perturbations
  
   INTEGER, INTENT(IN) :: max_irr_dim   ! maximum number of perturbations

   REAL(DP), INTENT(IN) :: becsum(nhm*(nhm+1)/2,nat,nspin) ! cross band 
                                                           ! occupations 
   COMPLEX(DP), INTENT(IN) :: dbecsum(nhm*(nhm+1)/2,nat,nspin,npe)! 
   
   COMPLEX(DP), INTENT(OUT) :: int3(nhm,nhm,max_irr_dim,nat,nspin) ! change of 
                                           !descreening coefficients (AE - PS)
   INTEGER, PARAMETER      :: AE = 1, PS = 2,&      ! All-Electron and Pseudo
                              XC = 1, H  = 2        ! XC and Hartree
   REAL(DP), POINTER       :: rho_core(:)           ! pointer to AE/PS core charge density 
   TYPE(paw_info)          :: i                     ! minimal info on atoms
   INTEGER                 :: i_what                ! counter on AE and PS
   INTEGER                 :: is                    ! spin index
   INTEGER                 :: lm                    ! counters on angmom and radial grid
   INTEGER                 :: nb, mb, nmb           ! augfun indexes
   INTEGER                 :: ia,mykey,ia_s,ia_e    ! atoms counters and indexes
   !
   REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:) ! density expanded on Y_lm
   REAL(DP), ALLOCATABLE   :: dv_lm(:,:,:) ! workspace: change of potential
   REAL(DP), ALLOCATABLE   :: drhor_lm(:,:,:,:) ! change of density expanded 
                                              ! on Y_lm (real part)
   REAL(DP), ALLOCATABLE   :: drhoi_lm(:,:,:,:) ! change of density expanded 
                                              ! on Y_lm (imaginary part)
   REAL(DP), ALLOCATABLE   :: savedvr_lm(:,:,:,:)   ! workspace: potential
   REAL(DP), ALLOCATABLE   :: savedvi_lm(:,:,:,:)   ! workspace: potential
   REAL(DP), ALLOCATABLE   :: aux_lm(:) ! auxiliary radial function
   ! fake cross band occupations to select only one pfunc at a time:
   REAL(DP)                :: becfake(nhm*(nhm+1)/2,nat,nspin)
   REAL(DP)                :: integral_r           ! workspace
   REAL(DP)                :: integral_i           ! workspace
   REAL(DP)                :: sgn                ! +1 for AE -1 for PS
   COMPLEX(DP)             :: sumd
   INTEGER  :: ipert

   CALL start_clock('PAW_dpot')
   ! Some initialization
   becfake(:,:,:) = 0._dp
   int3 = (0.0_DP, 0.0_DP)
   !
   ! Parallel: divide tasks among all the processor for this image
   ! (i.e. all the processors except for NEB and similar)
   CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )
   !
   atoms: DO ia = ia_s, ia_e
      !
      i%a = ia                      ! atom's index
      i%t = ityp(ia)                ! type of atom ia
      i%m = g(i%t)%mesh             ! radial mesh size for atom i%t
      i%b = upf(i%t)%nbeta          ! number of beta functions for i%t
      i%l = upf(i%t)%lmax_rho+1 ! max ang.mom. in augmentation for ia
      !
      ifpaw: IF (upf(i%t)%tpawp) THEN
         !
         ! Arrays are allocated inside the cycle to allow reduced
         ! memory usage as differnt atoms have different meshes
         !
         ALLOCATE(dv_lm(i%m,i%l**2,nspin))
         ALLOCATE(savedvr_lm(i%m,i%l**2,nspin,npe))
         ALLOCATE(savedvi_lm(i%m,i%l**2,nspin,npe))
         ALLOCATE(rho_lm(i%m,i%l**2,nspin))
         ALLOCATE(drhor_lm(i%m,i%l**2,nspin,npe))
         ALLOCATE(drhoi_lm(i%m,i%l**2,nspin,npe))
         ALLOCATE(aux_lm(i%m))
         !
         whattodo: DO i_what = AE, PS
            NULLIFY(rho_core)
            IF (i_what == AE) THEN
               CALL PAW_rho_lm(i, becsum, upf(i%t)%paw%pfunc, rho_lm)
               rho_core => upf(i%t)%paw%ae_rho_atc
               sgn = +1._dp
            ELSE
               CALL PAW_rho_lm(i, becsum, upf(i%t)%paw%ptfunc, &
                                          rho_lm, upf(i%t)%qfuncl)
               rho_core => upf(i%t)%rho_atc 
               sgn = -1._dp                 
            ENDIF
!
!           Compute the change of the charge density. Complex because the
!           displacements might be complex
!
            DO ipert=1,npe
               IF (i_what == AE) THEN
                  becfake(:,ia,:)=DBLE(dbecsum(:,ia,:,ipert))
                  CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%pfunc, &
                                     drhor_lm(1,1,1,ipert))
                  becfake(:,ia,:)=AIMAG(dbecsum(:,ia,:,ipert))
                  CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%pfunc, &
                                     drhoi_lm(1,1,1,ipert))
               ELSE
                  becfake(:,ia,:)=DBLE(dbecsum(:,ia,:,ipert))
                  CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%ptfunc, &
                                     drhor_lm(1,1,1,ipert), upf(i%t)%qfuncl)
                  becfake(:,ia,:)=AIMAG(dbecsum(:,ia,:,ipert))
                  CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%ptfunc, &
                                     drhoi_lm(1,1,1,ipert), upf(i%t)%qfuncl)
               END IF
            END DO

            savedvr_lm(:,:,:,:) = 0._dp
            savedvi_lm(:,:,:,:) = 0._dp

            DO ipert=1,npe
               !
               ! Change of Hartree potential
               !
               CALL PAW_h_potential(i, drhor_lm(1,1,1,ipert), dv_lm(:,:,1))
               DO is = 1,nspin 
                  savedvr_lm(:,:,is,ipert) = dv_lm(:,:,1)
               ENDDO
               CALL PAW_h_potential(i, drhoi_lm(1,1,1,ipert), dv_lm(:,:,1))
               DO is = 1,nspin 
                  savedvi_lm(:,:,is,ipert) = dv_lm(:,:,1)
               ENDDO 
               !
               ! Change of Exchange-correlation potential
               !
               CALL PAW_dxc_potential(i, drhor_lm(1,1,1,ipert), &
                                         rho_lm, rho_core, dv_lm)
               savedvr_lm(:,:,:,ipert) = savedvr_lm(:,:,:,ipert)+dv_lm(:,:,:)

               CALL PAW_dxc_potential(i, drhoi_lm(1,1,1,ipert), &
                                         rho_lm, rho_core, dv_lm)
               savedvi_lm(:,:,:,ipert) = savedvi_lm(:,:,:,ipert)+dv_lm(:,:,:)
            END DO
            !
            spins: DO is = 1, nspin
               nmb = 0
               ! loop on all pfunc for this kind of pseudo
               becfake=0.0_DP
               DO nb = 1, nh(i%t)
                  DO mb = nb, nh(i%t)
                     nmb = nmb+1 
                     becfake(nmb,ia,is) = 1._dp
                     IF (i_what == AE) THEN
                        CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%pfunc, rho_lm)
                     ELSE
                        CALL PAW_rho_lm(i, becfake, upf(i%t)%paw%ptfunc, &
                                           rho_lm, upf(i%t)%qfuncl)
                     ENDIF
!
!                 Integrate the change of Hxc potential and the partial waves
!                 to find the change of the D coefficients: D^1-~D^1
!
                     DO ipert=1,npe
                        DO lm = 1,i%l**2
                           aux_lm(1:i%m)=rho_lm(1:i%m,lm,is)* &
                                               savedvr_lm(1:i%m,lm,is,ipert) 
                           CALL simpson (upf(i%t)%kkbeta,aux_lm, &
                                                      g(i%t)%rab,integral_r)
                           aux_lm(1:i%m)=rho_lm(1:i%m,lm,is)* &
                                            savedvi_lm(1:i%m,lm,is,ipert) 
                           CALL simpson (upf(i%t)%kkbeta,aux_lm, &
                                                      g(i%t)%rab,integral_i)
                           int3(nb,mb,ipert,i%a,is) = &
                                        int3(nb,mb,ipert,i%a,is) &
                                       + sgn * CMPLX(integral_r, integral_i)
                        ENDDO
                        IF (nb /= mb)  int3(mb,nb,ipert,i%a,is) = &
                                                    int3(nb,mb,ipert,i%a,is) 
                     ENDDO
                     becfake(nmb,ia,is) = 0._dp
                  ENDDO ! mb
               ENDDO ! nb
            ENDDO spins
         ENDDO whattodo
         ! cleanup
         DEALLOCATE(rho_lm)
         DEALLOCATE(drhor_lm)
         DEALLOCATE(drhoi_lm)
         DEALLOCATE(savedvr_lm)
         DEALLOCATE(savedvi_lm)
         DEALLOCATE(dv_lm)
         DEALLOCATE(aux_lm)
         !
      ENDIF ifpaw
   ENDDO atoms

#ifdef __PARA
    IF( mykey /= 0 ) int3 = 0.0_dp
    CALL mp_sum(int3, intra_image_comm)
#endif
   CALL stop_clock('PAW_dpot')

END SUBROUTINE PAW_dpotential

SUBROUTINE PAW_dxc_potential(i, drho_lm, rho_lm, rho_core, v_lm)
!
!  This routine computes the change of the exchange and correlation 
!  potential in the spherical basis. It receives as input the charge
!  density and its variation.
!
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    USE funct,                  ONLY : dmxc, dmxc_spin, dmxc_nc, &
                                       dft_is_gradient

    TYPE(paw_info), INTENT(IN) :: i                   ! atom's minimal info
    REAL(DP), INTENT(IN)  :: rho_lm(i%m,i%l**2,nspin) ! charge density as 
                                                      ! lm components
    REAL(DP), INTENT(IN)  :: drho_lm(i%m,i%l**2,nspin)! change of charge 
                                                      ! density as lm components
    REAL(DP), INTENT(IN)  :: rho_core(i%m)            ! core charge, radial 
                                                      ! and spherical
    REAL(DP), INTENT(OUT) :: v_lm(i%m,i%l**2,nspin)   ! potential density 
                                                      ! as lm components
    REAL(DP), ALLOCATABLE  :: dmuxc(:,:,:)            ! fxc in the lsda case
    REAL(DP), ALLOCATABLE  :: v_rad(:,:,:)            ! radial potential 
                                                      ! (to be integrated)
    REAL(DP), ALLOCATABLE  :: rho_rad(:,:)            ! workspace (only one 
                                                      ! radial slice of rho)
    REAL(DP)              :: rho_loc(nspin)           ! workspace 
    
    REAL(DP) :: rhotot, rhoup, rhodw                  ! auxiliary
    
    INTEGER               :: ix,k                     ! counters on directions 
                                                      ! and radial grid

    CALL start_clock ('PAW_dxc_pot')
    ALLOCATE(dmuxc(i%m,nspin,nspin))
    ALLOCATE(v_rad(i%m,rad(i%t)%nx,nspin))
    ALLOCATE(rho_rad(i%m,nspin))
    !
    DO ix = 1, rad(i%t)%nx
!
! *** LDA (and LSDA) part (no gradient correction) ***
! convert _lm density to real density along ix
!
       CALL PAW_lm2rad(i, ix, rho_lm, rho_rad)
!
!      Compute the fxc function on the radial mesh along ix
!
       DO k = 1,i%m
          rho_loc(1:nspin) = rho_rad(k,1:nspin)*g(i%t)%rm2(k)
          IF (nspin==2) THEN
             rhoup = rho_loc(1)  + 0.5_DP * rho_core (k)
             rhodw = rho_loc(2)  + 0.5_DP * rho_core (k)
             CALL dmxc_spin (rhoup, rhodw, dmuxc(k,1,1), dmuxc(k,2,1),  &
                                           dmuxc(k,1,2), dmuxc(k,2,2) )
          ELSE
             rhotot = rho_loc(1) + rho_core (k)
             IF (rhotot.GT.1.d-30) v_rad (k,ix,1) = dmxc (rhotot)
             IF (rhotot.LT. - 1.d-30) v_rad(k, ix, 1) = - dmxc ( - rhotot)
             IF (rhotot.LT.1.d-30.AND.rhotot.GT.-1.d-30) v_rad(k,ix,1)=0.0_DP
          ENDIF
       ENDDO
!
!   Compute the change of the charge on the radial mesh along ix
!
       CALL PAW_lm2rad(i, ix, drho_lm, rho_rad)
!
!   fxc * dn
!
       IF (nspin==2) THEN
          DO k = 1,i%m
             v_rad(k,ix,1)= dmuxc(k,1,1)*rho_rad(k,1)*g(i%t)%rm2(k) &
                          + dmuxc(k,1,2)*rho_rad(k,2)*g(i%t)%rm2(k) 
             v_rad(k,ix,2)= dmuxc(k,2,1)*rho_rad(k,1)*g(i%t)%rm2(k) &
                          + dmuxc(k,2,2)*rho_rad(k,2)*g(i%t)%rm2(k) 
          ENDDO
       ELSE
          DO k = 1,i%m
             v_rad(k,ix,1)=v_rad(k,ix,1)*rho_rad(k,1)*g(i%t)%rm2(k) 
          ENDDO
       ENDIF
    ENDDO
!
! Recompose the sph. harm. expansion
!
    CALL PAW_rad2lm(i, v_rad, v_lm, i%l)
!
! Add gradient correction, if necessary
!
    IF( dft_is_gradient() ) &
       CALL PAW_dgcxc_potential(i,rho_lm,rho_core,drho_lm,v_lm)

    DEALLOCATE(rho_rad)
    DEALLOCATE(v_rad)
    DEALLOCATE(dmuxc)

    CALL stop_clock ('PAW_dxc_pot')

    RETURN
END SUBROUTINE PAW_dxc_potential

SUBROUTINE PAW_desymmetrize(dbecsum)
!
! This routine similar to PAW_symmetrize, symmetrize the change of 
! dbecsum due to an electric field perturbation. 
!
    USE lsda_mod,          ONLY : nspin
    USE uspp_param,        ONLY : nhm
    USE ions_base,         ONLY : nat, ityp
    USE symme,             ONLY : nsym, irt, d1, d2, d3, s
    USE uspp,              ONLY : nhtolm,nhtol,ijtoh
    USE uspp_param,        ONLY : nh, upf
    USE io_global,         ONLY : stdout, ionode

    COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2,nat,nspin,3)! cross band occupations

    COMPLEX(DP)                :: becsym(nhm*(nhm+1)/2,nat,nspin,3)! symmetrized becsum
    REAL(DP) :: pref, usym

    INTEGER :: ia, mykey,ia_s,ia_e   ! atoms counters and indexes
    INTEGER :: is, nt       ! counters on spin, atom-type
    INTEGER :: ma           ! atom symmetric to na
    INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
               l_i, l_j, m_i, m_j
    INTEGER :: m_o, m_u     ! counters for sums on m
    INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    INTEGER :: isym         ! counter for symmetry operation
    INTEGER :: ipol, jpol

    ! The following mess is necessary because the symmetrization operation
    ! in LDA+U code is simpler than in PAW, so the required quantities are
    ! represented in a simple but not general way.
    ! I will fix this when everything works.
    REAL(DP), TARGET :: d0(1,1,48)
    TYPE symmetryzation_tensor
        REAL(DP),POINTER :: d(:,:,:)
    END TYPE symmetryzation_tensor
    TYPE(symmetryzation_tensor) :: D(0:3)

    IF( nsym == 1 ) RETURN
    d0(1,1,:) = 1._dp
    D(0)%d => d0 ! d0(1,1,48)
    D(1)%d => d1 ! d1(3,3,48)
    D(2)%d => d2 ! d2(5,5,48)
    D(3)%d => d3 ! d3(7,7,48)

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


    CALL start_clock('PAW_dsymme')

    becsym(:,:,:,:) = (0.0_DP,0.0_DP)
    usym = 1._dp / REAL(nsym)

    ! Parallel: divide among processors for the same image
    CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )

    DO is = 1, nspin
    !
    atoms: DO ia = ia_s, ia_e
        nt = ityp(ia)
        ! No need to symmetrize non-PAW atoms
        IF ( .not. upf(nt)%tpawp ) CYCLE
        !
        DO ih = 1, nh(nt)
        DO jh = ih, nh(nt) ! note: jh >= ih
            !ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
            ijh = ijtoh(ih,jh,nt)
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
                ma = irt(isym,ia)
                DO m_o = 1, 2*l_i +1
                DO m_u = 1, 2*l_j +1
                    oh = ih - m_i + m_o
                    uh = jh - m_j + m_u
                    ouh = ijtoh(oh,uh,nt)
                    ! In becsum off-diagonal terms are multiplied by 2, I have
                    ! to neutralize this factor and restore it later
                    IF ( oh == uh ) THEN
                        pref = 2._dp * usym
                    ELSE
                        pref = usym
                    ENDIF
                    !
                    DO ipol=1,3
                       DO jpol=1,3
                          becsym(ijh, ia, is, ipol) = becsym(ijh, ia, is,ipol) &
                        + D(l_i)%d(m_o,m_i, isym) * D(l_j)%d(m_u,m_j, isym) &
                          * pref * dbecsum(ouh, ma, is, jpol) * s(ipol,jpol,isym)
                       ENDDO
                    ENDDO
                ENDDO ! m_o
                ENDDO ! m_u
            ENDDO ! isym
            !
            ! Put the prefactor back in:
            IF ( ih == jh ) becsym(ijh,ia,is,:) = .5_dp * becsym(ijh,ia,is,:)
        ENDDO ! ih
        ENDDO ! jh
    ENDDO atoms ! nat
    ENDDO ! nspin
#ifdef __PARA
    IF( mykey /= 0 ) becsym = 0.0_dp
    CALL mp_sum(becsym, intra_image_comm)
#endif

#ifdef __DEBUG_PAW_SYM
   write(stdout,*) "------------"
    if(ionode) then
        ia = 1
        nt = ityp(ia)
        DO is = 1, nspin
            write(*,*) is
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            ijh = ijtoh(ih,jh,nt)
            DO ipol=1,3
               write(stdout,"(1f10.3)", advance='no') becsym(ijh,ia,is,ipol)
            ENDDO
        ENDDO
            write(stdout,*)
        ENDDO
            write(stdout,*)
        ENDDO
    endif
   write(stdout,*) "------------"
#endif

    ! Apply symmetrization:
    dbecsum(:,:,:,:) = becsym(:,:,:,:)

    CALL stop_clock('PAW_dsymme')

END SUBROUTINE PAW_desymmetrize

SUBROUTINE PAW_dusymmetrize(dbecsum,npe,irr,max_irr_dim,nsymq,irgq,rtau,xq,t)
!
! This routine similar to PAW_symmetrize, symmetrize the change of 
! dbecsum due to an electric field perturbation. 
!
    USE lsda_mod,          ONLY : nspin
    USE uspp_param,        ONLY : nhm
    USE ions_base,         ONLY : nat, ityp
    USE symme,             ONLY : irt, d1, d2, d3
    USE constants,         ONLY : tpi
    USE uspp,              ONLY : nhtolm,nhtol,ijtoh
    USE uspp_param,        ONLY : nh, upf
    USE io_global,         ONLY : stdout, ionode

    COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2,nat,nspin,npe)! cross band occupations

    COMPLEX(DP)                :: becsym(nhm*(nhm+1)/2,nat,nspin,npe)! symmetrized becsum
    REAL(DP) :: pref, usym

    INTEGER, INTENT(IN) :: npe, irr, max_irr_dim, nsymq, irgq(48)
    REAL(DP), INTENT(IN) :: rtau(3,48,nat), xq(3)
    COMPLEX(DP), INTENT(IN) :: t(max_irr_dim, max_irr_dim, 48, 3*nat)
    INTEGER :: ia, mykey,ia_s,ia_e   ! atoms counters and indexes
    INTEGER :: is, nt       ! counters on spin, atom-type
    INTEGER :: ma           ! atom symmetric to na
    INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
               l_i, l_j, m_i, m_j
    INTEGER :: m_o, m_u     ! counters for sums on m
    INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    INTEGER :: isym, irot   ! counter for symmetry operation
    INTEGER :: ipol, jpol
    COMPLEX(DP) :: fase(48,nat)
    REAL(DP) :: arg, ft(3)

    ! The following mess is necessary because the symmetrization operation
    ! in LDA+U code is simpler than in PAW, so the required quantities are
    ! represented in a simple but not general way.
    ! I will fix this when everything works.
    REAL(DP), TARGET :: d0(1,1,48)
    TYPE symmetryzation_tensor
        REAL(DP),POINTER :: d(:,:,:)
    END TYPE symmetryzation_tensor
    TYPE(symmetryzation_tensor) :: D(0:3)

    IF( nsymq==1 ) RETURN
    d0(1,1,:) = 1._dp
    D(0)%d => d0 ! d0(1,1,48)
    D(1)%d => d1 ! d1(3,3,48)
    D(2)%d => d2 ! d2(5,5,48)
    D(3)%d => d3 ! d3(7,7,48)

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

    CALL start_clock('PAW_dsymme')

    becsym(:,:,:,:) = (0.0_DP,0.0_DP)
    usym = 1._dp / REAL(nsymq)

    do ia=1,nat
       do isym=1,nsymq
          irot = irgq (isym)
          arg = 0.0_DP
          do ipol = 1, 3
             arg = arg + xq (ipol) *  rtau(ipol,irot,ia)
          enddo
          arg = arg * tpi
          fase(irot,ia) = CMPLX (cos (arg),  sin (arg) )
       enddo
    enddo

    ! Parallel: divide among processors for the same image
    CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )

    DO is = 1, nspin
    !
    atoms: DO ia = ia_s, ia_e
        nt = ityp(ia)
        ! No need to symmetrize non-PAW atoms
        IF ( .not. upf(nt)%tpawp ) CYCLE
        !
        DO ih = 1, nh(nt)
        DO jh = ih, nh(nt) ! note: jh >= ih
            !ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
            ijh = ijtoh(ih,jh,nt)
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
            DO isym = 1,nsymq
                irot = irgq (isym)
                ma = irt(irot,ia)
                DO m_o = 1, 2*l_i +1
                DO m_u = 1, 2*l_j +1
                    oh = ih - m_i + m_o
                    uh = jh - m_j + m_u
                    ouh = ijtoh(oh,uh,nt)
                    ! In becsum off-diagonal terms are multiplied by 2, I have
                    ! to neutralize this factor and restore it later
                    IF ( oh == uh ) THEN
                        pref = 2._dp * usym
                    ELSE
                        pref = usym
                    ENDIF
                    !
                    DO ipol=1,npe
                       DO jpol=1,npe
                          becsym(ijh, ia, is, ipol) = becsym(ijh, ia, is,ipol) &
                        + D(l_i)%d(m_o,m_i, irot) * D(l_j)%d(m_u,m_j, irot) &
                          * pref * dbecsum(ouh, ma, is, jpol) * &
                          t(jpol,ipol,irot,irr) * fase(irot,ia)
                       ENDDO
                    ENDDO
                ENDDO ! m_o
                ENDDO ! m_u
            ENDDO ! isym
            !
            ! Put the prefactor back in:
            IF ( ih == jh ) becsym(ijh,ia,is,:) = .5_dp * becsym(ijh,ia,is,:)
        ENDDO ! ih
        ENDDO ! jh
    ENDDO atoms ! nat
    ENDDO ! nspin
#ifdef __PARA
    IF( mykey /= 0 ) becsym = 0.0_dp
    CALL mp_sum(becsym, intra_image_comm)
#endif

#ifdef __DEBUG_PAW_SYM
   write(stdout,*) "------------"
    if(ionode) then
        ia = 1
        nt = ityp(ia)
        DO is = 1, nspin
            write(*,*) is
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            ijh = ijtoh(ih,jh,nt)
            DO ipol=1,npe
               write(stdout,"(1f10.3)", advance='no') becsym(ijh,ia,is,ipol)
            ENDDO
        ENDDO
            write(stdout,*)
        ENDDO
            write(stdout,*)
        ENDDO
    endif
   write(stdout,*) "------------"
#endif

    ! Apply symmetrization:
    dbecsum(:,:,:,:) = becsym(:,:,:,:)

    CALL stop_clock('PAW_dsymme')

END SUBROUTINE PAW_dusymmetrize

SUBROUTINE PAW_dumqsymmetrize(dbecsum,npe,irr,max_irr_dim,isymq,rtau,xq,tmq)
!
! This routine similar to PAW_symmetrize, symmetrize the change of 
! dbecsum due to an electric field perturbation. 
!
    USE lsda_mod,          ONLY : nspin
    USE uspp_param,        ONLY : nhm
    USE ions_base,         ONLY : nat, ityp
    USE constants,         ONLY : tpi
    USE symme,             ONLY : nsym, irt, d1, d2, d3
    USE uspp,              ONLY : nhtolm,nhtol,ijtoh
    USE uspp_param,        ONLY : nh, upf
    USE io_global,         ONLY : stdout, ionode

    COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2,nat,nspin,npe)! cross band occupations

    COMPLEX(DP)                :: becsym(nhm*(nhm+1)/2,nat,nspin,npe)! symmetrized becsum
    REAL(DP), INTENT(IN) :: rtau(3,48,nat), xq(3)
    REAL(DP) :: pref

    INTEGER, INTENT(IN) :: npe, irr, max_irr_dim
    INTEGER, INTENT(IN) :: isymq         ! counter for symmetry operation
    COMPLEX(DP), INTENT(IN) :: tmq(max_irr_dim, max_irr_dim, 3*nat)
    INTEGER :: ia, mykey,ia_s,ia_e   ! atoms counters and indexes
    INTEGER :: is, nt       ! counters on spin, atom-type
    INTEGER :: ma           ! atom symmetric to na
    INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
               l_i, l_j, m_i, m_j
    INTEGER :: m_o, m_u     ! counters for sums on m
    INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    INTEGER :: ipol, jpol
    REAL(DP) :: arg
    COMPLEX(DP) :: fase(nat)

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

    CALL start_clock('PAW_dsymme')

    becsym(:,:,:,:) = (0.0_DP,0.0_DP)
    do ia=1,nat
       arg = 0.0_DP
       do ipol = 1, 3
          arg = arg + xq (ipol) *  rtau(ipol,isymq,ia)
       enddo
       arg = arg * tpi
       fase(ia) = CMPLX (cos (arg),  sin (arg) )
    enddo


    ! Parallel: divide among processors for the same image
    CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )

    DO is = 1, nspin
    !
    atoms: DO ia = ia_s, ia_e
        nt = ityp(ia)
        ! No need to symmetrize non-PAW atoms
        IF ( .not. upf(nt)%tpawp ) CYCLE
        !
        DO ih = 1, nh(nt)
        DO jh = ih, nh(nt) ! note: jh >= ih
            !ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
            ijh = ijtoh(ih,jh,nt)
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
            ma = irt(isymq,ia)
            DO m_o = 1, 2*l_i +1
            DO m_u = 1, 2*l_j +1
               oh = ih - m_i + m_o
               uh = jh - m_j + m_u
               ouh = ijtoh(oh,uh,nt)
               ! In becsum off-diagonal terms are multiplied by 2, I have
               ! to neutralize this factor and restore it later
               IF ( oh == uh ) THEN
                   pref = 2._dp 
               ELSE
                   pref = 1._DP
               ENDIF
               !
               DO ipol=1,npe
                  DO jpol=1,npe
                     becsym(ijh, ia, is, ipol) = becsym(ijh, ia, is,ipol) &
                        + D(l_i)%d(m_o,m_i, isymq) * D(l_j)%d(m_u,m_j, isymq) &
                          * pref * dbecsum(ouh, ma, is, jpol) * &
                          tmq(jpol,ipol,irr)*fase(ia)
                  ENDDO
               ENDDO
            ENDDO ! m_o
            ENDDO ! m_u
            !
            ! Put the prefactor back in:
            IF ( ih == jh ) becsym(ijh,ia,is,:) = .5_dp * becsym(ijh,ia,is,:)
            becsym(ijh, ia, is,:)=(CONJG(becsym(ijh, ia, is, :))+ &
                                            dbecsum(ijh, ia, is, :))*0.5_DP
        ENDDO ! ih
        ENDDO ! jh
    ENDDO atoms ! nat
    ENDDO ! nspin
#ifdef __PARA
    IF( mykey /= 0 ) becsym = 0.0_dp
    CALL mp_sum(becsym, intra_image_comm)
#endif

#ifdef __DEBUG_PAW_SYM
   write(stdout,*) "------------"
    if(ionode) then
        ia = 1
        nt = ityp(ia)
        DO is = 1, nspin
            write(*,*) is
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            ijh = ijtoh(ih,jh,nt)
            DO ipol=1,npe
               write(stdout,"(1f10.3)", advance='no') becsym(ijh,ia,is,ipol)
            ENDDO
        ENDDO
            write(stdout,*)
        ENDDO
            write(stdout,*)
        ENDDO
    endif
   write(stdout,*) "------------"
#endif

    ! Apply symmetrization:
    dbecsum(:,:,:,:) = becsym(:,:,:,:)

    CALL stop_clock('PAW_dsymme')

END SUBROUTINE PAW_dumqsymmetrize
!
! add gradient correction to dvxc. Both unpolarized and
! spin polarized cases are supported. 
!
SUBROUTINE PAW_dgcxc_potential(i,rho_lm,rho_core, drho_lm, v_lm)
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    USE constants,              ONLY : pi,e2, eps => eps12, eps2 => eps24
    USE funct,                  ONLY : gcxc, gcx_spin, gcc_spin, dgcxc, &
                                       dgcxc_spin
    !
    TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
    REAL(DP), INTENT(IN)    :: rho_lm(i%m,i%l**2,nspin) ! charge density as lm components
    REAL(DP), INTENT(IN)    :: drho_lm(i%m,i%l**2,nspin) ! change of charge density as lm components
    REAL(DP), INTENT(IN)    :: rho_core(i%m)            ! core charge, radial and spherical
    REAL(DP), INTENT(INOUT) :: v_lm(i%m,i%l**2,nspin)   ! potential to be updated
    REAL(DP)                :: zero(i%m)            ! dcore charge, not used
    REAL(DP)                :: rho_rad(i%m,nspin)! charge density sampled
    REAL(DP)                :: drho_rad(i%m,nspin)! charge density sampled
    REAL(DP)                :: grad(i%m,3,nspin) ! gradient
    REAL(DP)                :: grad2(i%m,nspin)  ! square modulus of gradient
                                                             ! (first of charge, than of hamiltonian)
    REAL(DP)                :: dgrad(i%m,3,nspin) ! gradient
    REAL(DP)                :: dgrad2(i%m,nspin)  ! square modulus of gradient
                                                  ! of dcharge
    REAL(DP)                :: gc_rad(i%m,rad(i%t)%nx,nspin) ! GC correction to V (radial samples)
    REAL(DP)                :: gc_lm(i%m,i%l**2,nspin)       ! GC correction to V (Y_lm expansion)
    REAL(DP)                :: h_rad(i%m,3,rad(i%t)%nx,nspin)! hamiltonian (vector field)
    REAL(DP)                :: h_lm(i%m,3,(i%l+rad(i%t)%ladd)**2,nspin)! hamiltonian (vector field)
                                       !!! ^^^^^^^^^^^^^^^^^^ expanded to higher lm than rho !!!
    REAL(DP)                :: div_h(i%m,i%l**2,nspin)  ! div(hamiltonian)


    INTEGER  :: k, ix, is, lm                             ! counters on spin and mesh
    REAL(DP) :: sx,sc,v1x,v2x,v1c,v2c                     ! workspace
    REAL(DP) :: v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw  ! workspace
    REAL(DP) :: vrrx,vsrx,vssx,vrrc,vsrc,vssc             ! workspace
    REAL(DP) :: dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s         ! workspace
    REAL(DP) :: vrrxup, vrrxdw, vrsxup, vrsxdw, vssxup, vssxdw, &
                vrrcup, vrrcdw, vrscup, vrscdw, vrzcup, vrzcdw
    REAL(DP) :: dsvxc_rr(2,2), dsvxc_sr(2,2), &
                dsvxc_ss(2,2), dsvxc_s(2,2) ! workspace
    REAL(DP) :: a(2,2,2), b(2,2,2,2), c(2,2,2)
    REAL(DP) :: arho, s1                                  ! workspace
    REAL(DP) :: rup, rdw, co2                             ! workspace
    REAL(DP) :: rh, zeta, grh2
    REAL(DP) :: grho(3,2), ps(2,2), ps1(3,2,2), ps2(3,2,2,2)
    INTEGER :: js, ls, ks, ipol

    OPTIONAL_CALL start_clock ('PAW_dgcxc_v')

    zero=0.0_DP
    gc_rad=0.0_DP
    h_rad=0.0_DP
    IF ( nspin == 1 ) THEN
       !
       !     GGA case - no spin polarization
       !
       DO ix = 1,rad(i%t)%nx
       !
          CALL PAW_lm2rad(i, ix, rho_lm, rho_rad)
          CALL PAW_gradient(i, ix, rho_lm, rho_rad, rho_core, grad2, grad)
          CALL PAW_lm2rad(i, ix, drho_lm, drho_rad)
          CALL PAW_gradient(i, ix, drho_lm, drho_rad, zero, dgrad2, dgrad)
          DO k = 1, i%m
             ! arho is the absolute value of real charge, sgn is its sign
             arho = rho_rad(k,1)*g(i%t)%rm2(k) + rho_core(k)
             arho = ABS(arho)
             s1 = grad (k, 1, 1) * dgrad(k, 1, 1) + &
                  grad (k, 2, 1) * dgrad(k, 2, 1) + &
                  grad (k, 3, 1) * dgrad(k, 3, 1)

           ! I am using grad(rho)**2 here, so its eps has to be eps**2
             IF ( (arho>eps) .and. (grad2(k,1)>eps2) ) THEN
                CALL gcxc(arho,grad2(k,1),sx,sc,v1x,v2x,v1c,v2c)
                CALL dgcxc(arho,grad2(k,1),vrrx,vsrx,vssx,vrrc,vsrc,vssc)
                dvxc_rr = vrrx + vrrc
                dvxc_sr = vsrx + vsrc
                dvxc_ss = vssx + vssc
                dvxc_s  = v2x + v2c
                gc_rad(k,ix,1)  = dvxc_rr*drho_rad(k, 1)*g(i%t)%rm2(k) &
                                + dvxc_sr*s1
                h_rad(k,:,ix,1) = ((dvxc_sr*drho_rad(k, 1)*g(i%t)%rm2(k) + &
                                   dvxc_ss*s1)*grad(k,:, 1) + &
                                   dvxc_s*dgrad(k,:,1))*g(i%t)%r2(k)
             ELSE
                gc_rad(k,ix,1)  = 0._dp
                h_rad(k,:,ix,1) = 0._dp
             ENDIF
          ENDDO
       ENDDO
    ELSEIF ( nspin == 2 ) THEN
       !
       !    \sigma-GGA case - spin polarization
       !
       DO ix = 1,rad(i%t)%nx
       !
          CALL PAW_lm2rad(i, ix, rho_lm, rho_rad)
          CALL PAW_gradient(i, ix, rho_lm, rho_rad, rho_core, &
                         grad2, grad)
          CALL PAW_lm2rad(i, ix, drho_lm, drho_rad)
          CALL PAW_gradient(i, ix, drho_lm, drho_rad, zero, dgrad2, dgrad)
       !
       DO k = 1,i%m
          !
          ! Prepare the necessary quantities
          ! rho_core is considered half spin up and half spin down:
          co2 = rho_core(k)/2._dp
          rup = rho_rad(k,1)*g(i%t)%rm2(k) + co2
          rdw = rho_rad(k,2)*g(i%t)%rm2(k) + co2
          CALL gcx_spin (rup, rdw, grad2(k,1), grad2(k,2), &
                          sx, v1xup, v1xdw, v2xup, v2xdw)
          grho(:,:)=grad(k,:,:)
          CALL dgcxc_spin (rup, rdw, grho (1,1), grho (1, 2), vrrxup, &
             vrrxdw, vrsxup, vrsxdw, vssxup, vssxdw, &
             vrrcup, vrrcdw, vrscup, vrscdw, vssc, vrzcup, vrzcdw)

          rh = rup + rdw ! total charge
          IF ( rh > eps ) THEN
             zeta = (rup - rdw ) / rh ! polarization
             !
             grh2 =  (grad(k,1,1) + grad(k,1,2))**2 &
                   + (grad(k,2,1) + grad(k,2,2))**2 &
                   + (grad(k,3,1) + grad(k,3,2))**2
             CALL gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
             dsvxc_rr (1, 1) = vrrxup + vrrcup + vrzcup *(1.d0 - zeta) / rh
             dsvxc_rr (1, 2) = vrrcup - vrzcup * (1.d0 + zeta) / rh
             dsvxc_rr (2, 1) = vrrcdw + vrzcdw * (1.d0 - zeta) / rh
             dsvxc_rr (2, 2) = vrrxdw + vrrcdw - vrzcdw *(1.d0 + zeta) / rh
             dsvxc_s (1, 1) = v2xup + v2c
             dsvxc_s (1, 2) = v2c
             dsvxc_s (2, 1) = v2c
             dsvxc_s (2, 2) = v2xdw + v2c
          ELSE
             sc    = 0._DP
             v1cup = 0._DP
             v1cdw = 0._DP
             v2c   = 0._DP
             dsvxc_rr = 0._DP
             dsvxc_s = 0._DP
          ENDIF
          dsvxc_sr (1, 1) = vrsxup + vrscup
          dsvxc_sr (1, 2) = vrscup
          dsvxc_sr (2, 1) = vrscdw
          dsvxc_sr (2, 2) = vrsxdw + vrscdw
          dsvxc_ss (1, 1) = vssxup + vssc
          dsvxc_ss (1, 2) = vssc
          dsvxc_ss (2, 1) = vssc
          dsvxc_ss (2, 2) = vssxdw + vssc
          ps (:,:) = (0._DP, 0._DP)
          DO is = 1, nspin
             DO js = 1, nspin
                ps1(:, is, js)=drho_rad(k,is)*g(i%t)%rm2(k)*grad(k,:,js)
                DO ipol=1,3
                   ps(is, js)=ps(is,js)+grad(k,ipol,is)*dgrad(k,ipol,js)
                ENDDO
                DO ks = 1, nspin
                   IF (is == js .AND. js == ks) THEN
                      a (is, js, ks) = dsvxc_sr (is, is)
                      c (is, js, ks) = dsvxc_sr (is, is)
                   ELSE
                      IF (is == 1) THEN
                         a (is, js, ks) = dsvxc_sr (1, 2)
                      ELSE
                         a (is, js, ks) = dsvxc_sr (2, 1)
                      ENDIF
                      IF (js == 1) THEN
                         c (is, js, ks) = dsvxc_sr (1, 2)
                      ELSE
                         c (is, js, ks) = dsvxc_sr (2, 1)
                      ENDIF
                   ENDIF
                   ps2 (:, is, js, ks) = ps (is, js) * grad (k,:,ks)
                   DO ls = 1, nspin
                      IF (is == js .AND. js == ks .AND. ks == ls) THEN
                         b (is, js, ks, ls) = dsvxc_ss (is, is)
                      ELSE
                         IF (is == 1) THEN
                            b (is, js, ks, ls) = dsvxc_ss (1, 2)
                         ELSE
                            b (is, js, ks, ls) = dsvxc_ss (2, 1)
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          DO is = 1, nspin
             DO js = 1, nspin
                gc_rad(k,ix,is)  = gc_rad(k,ix,is)+ dsvxc_rr (is,js) &
                                           *drho_rad(k, js)*g(i%t)%rm2(k)
                h_rad(k,:,ix,is) = h_rad(k,:,ix,is) + &
                                    dsvxc_s (is,js) * dgrad(k,:,js)
                DO ks = 1, nspin
                   gc_rad(k,ix,is) = gc_rad(k,ix,is)+a(is,js,ks)*ps(js,ks)
                   h_rad(k,:,ix,is) = h_rad(k,:,ix,is) + &
                         c (is, js, ks) * ps1 (:, js, ks)
                   DO ls = 1, nspin
                      h_rad(k,:,ix,is) = h_rad(k,:,ix,is) + &
                            b (is, js, ks, ls) * ps2 (:, js, ks, ls)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          h_rad(k,:,ix,:)=h_rad(k,:,ix,:)*g(i%t)%r2(k)
        ENDDO ! k
     ENDDO ! ix
    ELSEIF ( nspin == 4 ) THEN
        CALL errore('PAW_gcxc_v', 'non-collinear not yet implemented!', 1)
    ELSE 
        CALL errore('PAW_gcxc_v', 'unknown spin number', 2)
    ENDIF 
    !
    ! convert the first part of the GC correction back to spherical harmonics
    CALL PAW_rad2lm(i, gc_rad, gc_lm, i%l)
    !
    ! We need the divergence of h to calculate the last part of the exchange
    ! and correlation potential. First we have to convert H to its Y_lm expansion
    CALL PAW_rad2lm3(i, h_rad, h_lm, i%l+rad(i%t)%ladd)
    !
    ! Compute div(H)
    CALL PAW_divergence(i, h_lm, div_h, i%l+rad(i%t)%ladd, i%l)
    !                         input max lm --^     ^-- output max lm
    ! Finally sum it back into v_xc
    DO is = 1,nspin
       DO lm = 1,i%l**2
            v_lm(1:i%m,lm,is) = v_lm(1:i%m,lm,is) + &
                   e2*(gc_lm(1:i%m,lm,is)-div_h(1:i%m,lm,is))
       ENDDO
    ENDDO

    OPTIONAL_CALL stop_clock ('PAW_dgcxc_v')

END SUBROUTINE PAW_dgcxc_potential

#undef OPTIONAL_CALL
END MODULE paw_onecenter
