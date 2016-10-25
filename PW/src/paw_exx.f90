!
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Lorenzo Paulatto <paulatz@gmail.com> October 2012
!
! [1] J. Chem. Phys. 122, 234102 (2005)
!=----------------------------------------------------------------------------=!
MODULE paw_exx
  !=----------------------------------------------------------------------------=!
  USE kinds, ONLY : DP
  TYPE paw_keeq_type
    REAL(DP),POINTER :: k(:,:,:,:)
  END TYPE paw_keeq_type
  TYPE(paw_keeq_type),ALLOCATABLE :: ke(:)
  LOGICAL,PRIVATE :: paw_has_init_keeq = .false.

  CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE PAW_newdxx(weight, becphi, becpsi, deexx)
    !-----------------------------------------------------------------------
    ! This subroutine computes some sort of EXX contribution to the non-local 
    ! part of the hamiltonian. PAW one-center terms are computed here.
    USE ions_base,      ONLY : nat, ntyp => nsp, ityp
    USE uspp_param,     ONLY : upf, nh
    USE uspp,           ONLY : nkb
    USE paw_variables,  ONLY : okpaw
    USE mp_images,      ONLY : me_image
    USE uspp,           ONLY : indv_ijkb0
    IMPLICIT NONE
    !
    ! In input I get a slice of <beta|left> and <beta|right> only for this kpoint and this band
    COMPLEX(DP),INTENT(in)    :: becphi(nkb)
    COMPLEX(DP),INTENT(in)    :: becpsi(nkb) 
    COMPLEX(DP),INTENT(inout) :: deexx(nkb)
    REAL(DP)                  :: weight
    !
    ! ... local variables
    INTEGER :: ijkb0, ih, jh, na, np, ikb
    !
    IF(.not.paw_has_init_keeq) &
      CALL errore("PAW_deexx", "you have to initialize paw keeq before", 1)
    !
    CALL start_clock( 'PAW_newdxx' )
    !
    IF(okpaw) RETURN
    ! Worst possible parallelisation:
    IF(me_image/=0) RETURN
    !
    DO np = 1, ntyp
      ONLY_FOR_PAW : &
      IF ( upf(np)%tpawp ) THEN
        !
        DO ih = 1, nh(np)
        DO jh = 1, nh(np)
            !
            ATOMS_LOOP : &
            DO na = 1, nat
            IF (ityp(na)==np) THEN
                !
                ! NOTE: see addusxx_g for the next line:
                ijkb0 = indv_ijkb0(na)
                ikb = ijkb0 + ih
                deexx(ikb) = deexx(ikb) &
                        - weight*PAW_deexx(na, ih, jh, ijkb0, becphi, becpsi)
                !
            END IF
            ENDDO ATOMS_LOOP ! nat
        ENDDO ! jh
        ENDDO ! ih
      END IF &
      ONLY_FOR_PAW
    ENDDO
    !
    CALL stop_clock( 'PAW_newdxx' )
    !
    RETURN
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE PAW_newdxx
  !-----------------------------------------------------------------------
  !
  !=----------------------------------------------------------------------------=!
  FUNCTION PAW_deexx(na, ih, jh, ijkb0, becphi, becpsi)
    !=----------------------------------------------------------------------------=!
    ! Compute the 2-electron 4-wavefunctions integral 
    ! Integral over bands and kpoints is done outside (doing it here does not fit properly with exx.f90)
    USE ions_base,          ONLY : nat, ityp
    USE uspp_param,         ONLY : nh, upf
    USE uspp,               ONLY : nkb
    IMPLICIT NONE
    INTEGER,INTENT(in) :: na, ih, jh, ijkb0
    COMPLEX(DP),INTENT(in) :: becphi(nkb), becpsi(nkb)
    !
    COMPLEX(DP) :: PAW_deexx
    !
    INTEGER :: np
    INTEGER :: oh, uh
    INTEGER :: ikb, jkb, okb, ukb 
    !
    PAW_deexx = 0._dp
    np = ityp(na)
    IF(.not.upf(np)%tpawp) RETURN
    !
    ! CALL start_clock("PAW_deexx")
    ikb = ijkb0 + ih
    jkb = ijkb0 + jh

    DO oh = 1, nh(np)
      okb = ijkb0 + oh
      DO uh = 1, nh(np)
        ukb = ijkb0 + uh
        ! Eq. 35 + 32 Ref. 1, the 1/2 factor comes from 32
        PAW_deexx = PAW_deexx  &
                  +  0.5_DP *ke(np)%k(ih,jh,oh,uh) * becphi(jkb) &
                            * CONJG(becphi(ukb)) * becpsi(okb)
        !
      ENDDO !uh, ukb
    ENDDO !oh, okb
    !
    ! CALL stop_clock("PAW_deexx")
    RETURN
    !=----------------------------------------------------------------------------=!
  END FUNCTION PAW_deexx
  !=----------------------------------------------------------------------------=!
  !
  !=----------------------------------------------------------------------------=!
  FUNCTION PAW_xx_energy(becphi, becpsi)
    !=----------------------------------------------------------------------------=!
    ! Compute the energy: 2-electron 4-wavefunctions integral and sum with weights and <beta|psi>
    ! Integral over bands and kpoints is done outside (doing it here would not fit properly with exx.f90)
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    USE uspp_param,         ONLY : nh, upf
    USE uspp,               ONLY : nkb
    USE mp_images,          ONLY : me_image
    IMPLICIT NONE
    COMPLEX(DP),INTENT(in) :: becphi(nkb), becpsi(nkb)
    !
    REAL(DP) :: PAW_xx_energy
    !
    INTEGER :: np, na
    INTEGER :: ih, jh, oh, uh
    INTEGER :: ikb, jkb, okb, ukb, ijkb0
    IF(.not.paw_has_init_keeq) &
        CALL errore("PAW_xx_energy", "you have to initialize paw keeq before", 1)
    !
    PAW_xx_energy = 0._dp
    IF(me_image/=0) RETURN
    !
    CALL start_clock("PAW_xx_nrg")
    !
    ijkb0 = 0
    DO np = 1, ntyp
      ONLY_FOR_PAW : &
      IF ( upf(np)%tpawp ) THEN
          DO na = 1, nat
          IF (ityp(na)==np) THEN
              !
              DO ih = 1, nh(np)
                ikb = ijkb0 + ih
                DO jh = 1, nh(np)
                  jkb = ijkb0 + jh
                  DO oh = 1, nh(np)
                    okb = ijkb0 + oh
                    DO uh = 1, nh(np)
                      ukb = ijkb0 + uh
                      ! Eq. 32 and 42 Ref. 1 :
                      PAW_xx_energy = PAW_xx_energy - 0.5_dp * ke(np)%k(ih,jh,oh,uh) &
                                    * CONJG(becpsi(ikb)) * becpsi(okb) & ! \rho_ik eq. 31 ref. 1
                                    * becphi(jkb) * CONJG(becphi(ukb))   ! \rho_lj eq. 31 ref. 1
                      !
                    ENDDO !uh, ukb
                  ENDDO !oh, okb
                ENDDO !jh, jkb
              ENDDO !ih, ikb
              !
              ijkb0 = ijkb0 + nh(np)
          END IF
          ENDDO ! nat
      ELSE ONLY_FOR_PAW 
          DO na = 1, nat
            IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
          ENDDO
      END IF &
      ONLY_FOR_PAW 
    ENDDO
    !
    CALL stop_clock("PAW_xx_nrg")
    RETURN
    !=----------------------------------------------------------------------------=!
  END FUNCTION PAW_xx_energy
  !=----------------------------------------------------------------------------=!
  !
  !=----------------------------------------------------------------------------=!
  SUBROUTINE PAW_init_keeq()
    !=----------------------------------------------------------------------------=!
    ! Driver to compute the 2-electron 4-wavefunctions integrals
    USE kinds,             ONLY : DP
    USE ions_base,         ONLY : ntyp => nsp
    USE uspp_param,        ONLY : nh
    IMPLICIT NONE
    INTEGER :: ns, ih,jh,oh,uh
    REAL(DP),ALLOCATABLE :: k_ae(:,:,:,:), k_ps(:,:,:,:)

    IF(paw_has_init_keeq) RETURN !CALL errore("PAW_init_keeq", "already init paw keeq", 1)
    paw_has_init_keeq = .true.
    !
    ! We have one matrix for the all electron and one for the pseudo part for each atomic specie
    ALLOCATE(ke(ntyp))
    CALL allocate_keeq(ntyp, nh, ke)

    DO ns = 1,ntyp
      !
      ALLOCATE(k_ae(nh(ns),nh(ns),nh(ns),nh(ns)))
      CALL PAW_keeq('AE', ns, k_ae)
      !
      ALLOCATE(k_ps(nh(ns),nh(ns),nh(ns),nh(ns)))
      CALL PAW_keeq('PS', ns, k_ps)
      !
      ! Symmetrize wrt the on-site wavefunctions indexes as the hartree kernel is not 
      ! perfectly symmetrical: the asymmetry accumulates and causes S matrix to be non-positive
      ! definite (especially with many k-points)
      DO ih = 1, nh(ns)
      DO jh = 1, nh(ns)
      DO oh = 1, nh(ns)
      DO uh = 1, nh(ns)
        !
        ke(ns)%k(ih,jh,oh,uh) = 0.25_dp * ( &
                k_ae(ih,jh,oh,uh)-k_ps(ih,jh,oh,uh) &
              + k_ae(oh,uh,ih,jh)-k_ps(oh,uh,ih,jh) &
              + k_ae(jh,ih,uh,oh)-k_ps(jh,ih,uh,oh) &
              + k_ae(uh,oh,jh,ih)-k_ps(uh,oh,jh,ih) )
        !
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !
      DEALLOCATE(k_ae, k_ps)
      !
    ENDDO
    !=----------------------------------------------------------------------------=!
  END SUBROUTINE PAW_init_keeq
  !=----------------------------------------------------------------------------=!
  !
  !=----------------------------------------------------------------------------=!
  SUBROUTINE PAW_destroy_keeq()
    !=----------------------------------------------------------------------------=!
    ! ke_ae and ke_ps for later use
    USE ions_base,         ONLY : ityp, ntyp => nsp
    IMPLICIT NONE

    IF(.not.paw_has_init_keeq) RETURN !CALL errore("PAW_destroy_keeq", "nothing to destroy :(", 1)
    paw_has_init_keeq = .false.
    !
    ! We have one matrix for the all electron and one for the pseudo part for each atomic specie
    CALL deallocate_keeq(ntyp, ke)
    DEALLOCATE(ke)
    !
    !=----------------------------------------------------------------------------=!
  END SUBROUTINE PAW_destroy_keeq
  !=----------------------------------------------------------------------------=!
  !
  !=----------------------------------------------------------------------------=!
  SUBROUTINE allocate_keeq(ntp, nh, keeq)
    !=----------------------------------------------------------------------------=!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ntp
    INTEGER,INTENT(in) :: nh(ntp)
    TYPE(paw_keeq_type),INTENT(inout) :: keeq(ntp)
    INTEGER :: i
    !
    DO i = 1,ntp
      ALLOCATE(keeq(i)%k(nh(i),nh(i),nh(i),nh(i)))
    ENDDO
    RETURN
    !=----------------------------------------------------------------------------=!
  END SUBROUTINE allocate_keeq
  !=----------------------------------------------------------------------------=!
  !=----------------------------------------------------------------------------=!
  SUBROUTINE deallocate_keeq(ntp, keeq)
    !=----------------------------------------------------------------------------=!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ntp
    TYPE(paw_keeq_type),INTENT(inout) :: keeq(ntp)
    INTEGER :: i
    !
    DO i = 1,ntp
      DEALLOCATE(keeq(i)%k)
    ENDDO
    !
    RETURN
    !=----------------------------------------------------------------------------=!
  END SUBROUTINE deallocate_keeq
  !=----------------------------------------------------------------------------=!
  !
  !=----------------------------------------------------------------------------=!
  SUBROUTINE PAW_keeq(what, np, keeq)
    !=----------------------------------------------------------------------------=!
    ! Compute the 2-electron 4-wavefunctions integrals and i.e. the exchange integral
    ! between two one-center wavefunctions. Includes augmentation in the pseudo case.
    ! Store it in global variables ke_ae and ke_ps for later use. 
    USE constants,         ONLY : e2
    USE atom,              ONLY : g => rgrid
    USE ions_base,         ONLY : nat, ityp, ntyp => nsp
    USE uspp_param,        ONLY : nh, nhm, upf
    USE paw_variables,     ONLY : paw_info
    USE paw_onecenter,     ONLY : PAW_h_potential
    USE lsda_mod,  ONLY : nspin
    IMPLICIT NONE
    INTEGER, INTENT(in)  :: np ! atomic type
    REAL(DP),INTENT(out) :: keeq(nh(np),nh(np),nh(np),nh(np))
    CHARACTER(len=2),INTENT(in) :: what ! "AE"= all-electron or "PS"=pseudo
    !
    TYPE(paw_info) :: i
    REAL(DP), ALLOCATABLE   :: v_lm(:,:)        ! workspace: potential
    REAL(DP), ALLOCATABLE   :: rho_lm_ij(:,:,:) ! density expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: rho_lm_ou(:,:)   ! density expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: aux(:)   ! density expanded on Y_lm
    INTEGER  :: ih,jh,oh,uh,k,lm
    REAL(DP) :: kexx, e

    IF(what/="AE" .and. what /="PS") CALL errore("PAW_keeq", "can only do all-electron or pseudo", 1)

    ! Only wavefunctions on the same atom exchange, and the result only depends on the atom type
    IF (.not.upf(np)%tpawp) THEN
        keeq = 0._dp
        RETURN
    ENDIF
    !
    CALL start_clock('PAW_keeq')
    !
    i%a = -1                  ! atom's index (UNUSED HERE)
    i%t = np                  ! type of atom
    i%m = g(i%t)%mesh         ! radial mesh size for atom i%t
    i%b = upf(i%t)%nbeta      ! number of beta functions for i%t
    i%l = upf(i%t)%lmax_rho+1 ! max ang.mom. in augmentation for ia
    !
    ! Note: PAW_h_potential computes the V_H for every spin component, 
    !       here they will be all zero except the 1st
    ALLOCATE(rho_lm_ij(i%m,i%l**2,nspin))
    ALLOCATE(rho_lm_ou(i%m,i%l**2))
    ALLOCATE(v_lm(i%m,i%l**2))
    ALLOCATE(aux(i%m))
    !
    DO ih = 1, nh(i%t)
    DO jh = 1, nh(i%t)
      !
      rho_lm_ij = 0._dp
      IF(what=="PS")THEN
        CALL PAW_rho_lm_ij(i, ih, jh, upf(i%t)%paw%ptfunc, rho_lm_ij(:,:,1), upf(i%t)%qfuncl)
      ELSE
        CALL PAW_rho_lm_ij(i, ih, jh, upf(i%t)%paw%pfunc, rho_lm_ij(:,:,1))
      ENDIF
      CALL PAW_h_potential(i, rho_lm_ij, v_lm) 
      !
      DO oh = 1, nh(i%t)
      DO uh = 1, nh(i%t)
        !
        rho_lm_ou = 0._dp
        IF(what=="PS")THEN
          CALL PAW_rho_lm_ij(i, oh, uh, upf(i%t)%paw%ptfunc, rho_lm_ou, upf(i%t)%qfuncl)
        ELSE
          CALL PAW_rho_lm_ij(i, oh, uh, upf(i%t)%paw%pfunc, rho_lm_ou)
        ENDIF
        !
        ! Now I have rho_ij and rho_ou, I have to compute the 4-wfcs hartree kernel keeq=K_ijou
        kexx = 0._dp
        DO lm = 1, i%l**2
            DO k = 1, i%m
                aux(k) = v_lm(k,lm) * rho_lm_ou(k,lm)
            ENDDO
            CALL simpson(upf(i%t)%kkbeta, aux, g(i%t)%rab, e)
            kexx = kexx + e
            !
        ENDDO
        keeq(ih,jh,oh,uh) = e2*kexx ! = K_ijlk : Eq. 33 Ref. 1
        !
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    !
    DEALLOCATE(aux, v_lm, rho_lm_ij, rho_lm_ou)
    !
    CALL stop_clock('PAW_keeq')
    !
    !=----------------------------------------------------------------------------=!
  END SUBROUTINE PAW_keeq
  !=----------------------------------------------------------------------------=!
  !
  !=----------------------------------------------------------------------------=!
  SUBROUTINE PAW_rho_lm_ij(i, ih, jh, pfunc, rho_lm, aug)
    !=----------------------------------------------------------------------------=!
    ! Computes the fake two-wavefunctions density i.e. phi_i(r)phi_j(r)^*, 
    ! Represent it as spherical harmonics. Details: this is a generalized version of PAW_rho_lm.
    USE uspp_param,        ONLY : upf
    USE uspp,              ONLY : indv, ap, nhtolm,lpl,lpx
    USE atom,              ONLY : g => rgrid
    USE paw_variables,     ONLY : paw_info
    IMPLICIT NONE
    TYPE(paw_info), INTENT(IN) :: i   ! atom's minimal info
    INTEGER,INTENT(in) :: ih, jh
    REAL(DP), INTENT(IN)  :: pfunc(i%m,i%b,i%b)             ! psi_i * psi_j
    REAL(DP), INTENT(OUT) :: rho_lm(i%m,i%l**2)       ! AE charge density on rad. grid
    REAL(DP), OPTIONAL,INTENT(IN) :: &
                            aug(i%m,i%b*(i%b+1)/2,0:2*upf(i%t)%lmax) ! augmentation functions (only for PS part)

    INTEGER                 :: nb, mb, &
                               nmb, &    ! composite "triangular" index for pfunc nmb = 1,nh*(nh+1)/2
                               lm,lp,l   ! counters for angular momentum lm = l**2+m
    REAL(DP) :: pref
    ! initialize density
    rho_lm(:,:) = 0._dp
    ! loop on all pfunc for this kind of pseudo
    nb = indv(ih,i%t)
    mb = indv(jh,i%t)
    nmb = mb * (mb-1)/2 + nb  ! mb has to be .ge. nb
    !
    angular_momentum: &
    DO lp = 1, lpx (nhtolm(jh,i%t), nhtolm(ih,i%t)) !lmaxq**2
        ! the lpl array contains the possible combination of LM,lm_j,lm_j that
        ! have non-zero Clebch-Goordon coefficient a_{LM}^{(lm)_i(lm)_j} and
        ! lpx is the number of non-zero CG coeffs... (this way we save some loops)
        lm = lpl (nhtolm(jh,i%t), nhtolm(ih,i%t), lp)
        !
        ! pref is "just" the Clebsh-Gordon coefficients
        pref = ap(lm, nhtolm(ih,i%t), nhtolm(jh,i%t))
        !
        rho_lm(1:i%m,lm) = rho_lm(1:i%m,lm) + pref * pfunc(1:i%m, nb, mb)
        IF (present(aug)) THEN
            ! if I'm doing the pseudo part I have to add the augmentation charge
            l = INT(SQRT(DBLE(lm-1))) ! l has to start from zero, lm = l**2 +m
            rho_lm(1:i%m,lm) = rho_lm(1:i%m,lm) +pref * aug(1:i%m, nmb, l)
        ENDIF ! augfun
    ENDDO angular_momentum 
    !
    !=----------------------------------------------------------------------------=!
  END SUBROUTINE PAW_rho_lm_ij
  !=----------------------------------------------------------------------------=!
  !
!=----------------------------------------------------------------------------=!
END MODULE paw_exx
!=----------------------------------------------------------------------------=!
