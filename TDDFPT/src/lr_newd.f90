!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! OBM
!  050608 gamma_only correction
!        tvanp --> upf%tvanp
!  160608 reduce --> mp_sum

#include "f_defs.h"
!
!----------------------------------------------------------------------
SUBROUTINE lr_newd( vr, deeq )
  !----------------------------------------------------------------------
  !
  !   This routine computes the integral of the perturbed potential with
  !   the Q function 
  !
  USE kinds,                   ONLY : DP
  USE realus,                  ONLY : real_space, real_space_debug
  USE ions_base,               ONLY : ityp, nat
  USE cell_base,               ONLY : omega
  USE gvect,                   ONLY : nr1, nr2, nr3, nrxx 
  USE lsda_mod,                ONLY : nspin
  USE uspp,                    ONLY : okvan
  USE uspp_param,              ONLY : nh, nhm, upf
  use mp,                      only : mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE lr_variables,   ONLY : lr_verbosity
   USE io_global,      ONLY : stdout
 !
  IMPLICIT NONE
  !
  REAL(kind=dp), intent(in)  :: vr(nrxx)
  REAL(kind=dp), intent(out) :: deeq( nhm, nhm, nat, nspin )
  !
  IF ( .NOT. okvan ) RETURN
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_newd>")')
  endif
  CALL start_clock ('lr_newd')
  !
  !
  IF ( real_space_debug>6 ) THEN
     !
     call lr_newd_r()
     !
  ELSE
     !
     call lr_newd_g()
     !
  END IF
  !
  CALL stop_clock ('lr_newd')
  !
  RETURN
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE lr_newd_g()
    !------------------------------------------------------------------------
    !
    ! ... this subroutine is the version of lr_newd in G space
    !
    USE ions_base,               ONLY : ntyp=>nsp
    USE gvect,                   ONLY : nrx1, nrx2, nrx3
    USE gvect,                   ONLY : g, gg, ngm, gstart, ig1, ig2, ig3, nl
    USE gvect,                   ONLY : eigts1, eigts2, eigts3
    USE uspp_param,              ONLY : lmaxq
    USE wavefunctions_module,    ONLY : psic
    USE control_flags,                   ONLY : gamma_only
    !
    IMPLICIT NONE
    !
    INTEGER :: ig, nt, ih, jh, na, is
    !
    ! counters on g vectors, atom type, beta functions x 2, atoms, spin
    !
    COMPLEX(kind=DP), allocatable :: aux (:,:), qgm_na (:)
    !
    ! work space
    !
    REAL(kind=DP), allocatable  :: ylmk0 (:,:), qmod (:)
    !
    ! spherical harmonics, modulus of G
    !
    REAL(kind=DP) :: fact, DDOT
    !
    ! passed in self consistent potential
    !
    COMPLEX(KIND=DP), ALLOCATABLE :: qgm(:)  ! complete fourier transform of Q
    !
    ALLOCATE (qgm( ngm))
    !
    IF ( gamma_only ) THEN
       !
       fact = 2.d0
       !
    ELSE
       !
       fact = 1.d0
       !
    END IF
    !
    ALLOCATE ( aux( ngm, nspin ) )
    ALLOCATE ( qgm_na( ngm ) )
    ALLOCATE ( qmod( ngm ) )
    ALLOCATE ( ylmk0( ngm, lmaxq*lmaxq ) )
    !
    deeq(:,:,:,:) = 0.0d0
    !
    CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
    !
    DO ig = 1, ngm
       !
       qmod (ig) = sqrt (gg (ig) )
       !
    END DO
    !
    ! fourier transform of the total effective potential
    !
#ifdef DEBUG_NEWD
    CALL start_clock ('newd:fftvg')
#endif
    !
    psic (:) = vr (:)
    !
    CALL cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
    !
    DO ig = 1, ngm
       !
       aux (ig, 1) = psic (nl (ig) )
       !
    END DO
    !
#ifdef DEBUG_NEWD
    CALL stop_clock ('newd:fftvg')
#endif
    !
    ! here we compute the integral Q*V for each atom,
    !       I = sum_G exp(-iR.G) Q_nm v^*
    !
    DO nt = 1, ntyp 
       !
       IF ( upf(nt)%tvanp ) then
          !
          DO ih = 1, nh (nt)
             !
             DO jh = ih, nh (nt)
                ! 
#ifdef DEBUG_NEWD
                CALL start_clock ('newd:qvan2')
#endif
                !
                ! The Q(r) for this atomic species without structure factor
                !
                CALL qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
                !
#ifdef DEBUG_NEWD
                CALL stop_clock ('newd:qvan2')
#endif
                !
                DO na = 1, nat
                   !
                   IF ( ityp (na) .eq. nt ) THEN
                      !
#ifdef DEBUG_NEWD
                      CALL start_clock ('newd:int1')
#endif
                      !
                      ! The Q(r) for this specific atom
                      !
                      DO ig = 1, ngm
                         !
                         qgm_na (ig) = qgm(ig) * eigts1 (ig1(ig), na) &
                                     * eigts2 (ig2(ig), na) &
                                     * eigts3 (ig3(ig), na) 
                         !
                      END DO
                      !
#ifdef DEBUG_NEWD
                      CALL stop_clock ('newd:int1')
#endif
                      !
#ifdef DEBUG_NEWD
                      CALL start_clock ('newd:int2')
#endif
                      !
                      !  and the product with the Q functions
                      !
                      DO is = 1, nspin
                         !
                         deeq (ih, jh, na, is) = fact * omega * &
                                                 DDOT (2 * ngm, aux(1,is), 1, qgm_na, 1)
                         !
                         IF ( gamma_only .and. gstart==2 ) &
                              deeq (ih, jh, na, is) = deeq (ih, jh, na, is) - &
                              omega*real(aux(1,is)*qgm_na(1),dp)
                         !
                         deeq (jh, ih, na, is) = deeq (ih, jh, na, is)
                         !
                      END DO
                      !
#ifdef DEBUG_NEWD
                      CALL stop_clock ('newd:int2')
#endif
                      !
                   END IF
                   !
                END DO
                !
             END DO
             !
          END DO
          !
       END IF
       !
    END DO
    !
#ifdef __PARA
    !CALL reduce (nhm * nhm * nat * nspin, deeq)
    call mp_sum(deeq, intra_pool_comm)
#endif
    !
    DEALLOCATE ( aux, qgm_na, qmod, ylmk0 )
    DEALLOCATE ( qgm )
    !
    RETURN
    !
  END SUBROUTINE lr_newd_g
  !
  !------------------------------------------------------------------------
  SUBROUTINE lr_newd_r()
    !------------------------------------------------------------------------
    !
    ! ... this subroutine is the version of lr_newd in real space
    ! what is the difference between newd_r ?
    !
    USE constants,        ONLY : pi, fpi
    USE uspp,             ONLY : okvan
    USE realus,           ONLY : box, maxbox, qsave
    !
    IMPLICIT NONE
    !
    REAL(DP), ALLOCATABLE :: aux(:)
    INTEGER               :: ia, ih, jh, ijh, is, ir, nt
    INTEGER               :: mbia, iqs
    !
    deeq(:,:,:,:) = 0.D0
    !
    ALLOCATE( aux( nrxx ) )
    !
    DO is = 1, nspin
       !
       aux(:) = vr(:)
       !
       iqs = 0
       !
       DO ia = 1, nat
          !
          mbia = maxbox(ia)
          !
          IF ( mbia == 0 ) CYCLE
          !
          nt = ityp(ia)
          !
          IF ( .NOT. upf(nt)%tvanp ) CYCLE
          !
          DO ih = 1, nh(nt)
             !
             DO jh = ih, nh(nt)
                !
                DO ir = 1, mbia
                   !
                   iqs = iqs + 1
                   deeq(ih,jh,ia,is)= deeq(ih,jh,ia,is) + &
                                      qsave(iqs)*aux(box(ir,ia))
                   !
                END DO
                !
                deeq(jh,ih,ia,is) = deeq(ih,jh,ia,is)
                !
             END DO
             !
          END DO
          !
       END DO
       !
    END DO
    !
    deeq(:,:,:,:) = deeq(:,:,:,:)*omega/(nr1*nr2*nr3)
    !
    DEALLOCATE( aux )
    !
#ifdef __PARA
    !CALL reduce( nhm*nhm*nat*nspin, deeq )
    call mp_sum(deeq, intra_pool_comm)
#endif
    !
    RETURN
    !
    END SUBROUTINE lr_newd_r
    !
END SUBROUTINE lr_newd
