!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
SUBROUTINE lr_addusdens(rho_1)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the perturbed charge density the part which is due to
  !  the US augmentation.
  !
  ! OBM
  ! 050608 gamma_point correction
  !        tvanp --> upf%tvanp
  USE kinds,          ONLY : DP
  USE realus,   ONLY : real_space, real_space_debug
  USE ions_base,      ONLY : ityp,nat
  USE gvect,          ONLY : nr1, nr2, nr3, nrxx
  USE lsda_mod,       ONLY : nspin
  USE uspp,           ONLY : okvan, becsum
  USE uspp_param,     ONLY : nh, upf
  USE lr_variables,   ONLY : lr_verbosity
  USE io_global,      ONLY : stdout

  !
  IMPLICIT NONE
  !
  ! perturbed charge density
  !
  REAL(kind=dp), intent(inout) :: rho_1(nrxx)
  !
  IF ( .NOT. okvan ) RETURN
  !
  call start_clock ('lr_addusdens')
  !
  if (lr_verbosity > 5) THEN
   WRITE(stdout,'("<lr_addusdens>")')
  ENDIF
  IF ( real_space_debug > 6  ) THEN
     !
     call lr_addusdens_r()
     !
  ELSE
     !
     call lr_addusdens_g()
     !
  END IF
  !
  if (lr_verbosity > 5) THEN
   WRITE(stdout,'("<end of lr_addusdens>")')
  ENDIF
  !
  call stop_clock ('lr_addusdens')
  !
  RETURN
  !
CONTAINS
  !----------------------------------------------------------------------
  SUBROUTINE lr_addusdens_g()
    !----------------------------------------------------------------------
    !
    ! ... G space version
    !
    USE ions_base,               ONLY : ntyp=>nsp
    USE gvect,                   ONLY : nrx1, nrx2, nrx3, ngm, nl, nlm, gg, g
    USE gvect,                   ONLY : eigts1, eigts2, eigts3, ig1, ig2, ig3
    USE uspp_param,              ONLY : lmaxq
    USE control_flags,           ONLY : gamma_only
    USE wavefunctions_module,    ONLY : psic
    !
    IMPLICIT NONE
    !
    INTEGER :: ig, na, nt, ih, jh, ijh, is
    !
    ! counters
    !
    REAL(kind=DP), allocatable :: qmod (:), ylmk0 (:,:)
    !
    ! the modulus of G
    ! the spherical harmonics
    !
    COMPLEX(kind=DP) :: skk
    COMPLEX(kind=DP), allocatable ::  aux (:,:)
    !
    ! work space for FFT
    ! work space for rho(G,nspin)
    !
    COMPLEX(KIND=DP), ALLOCATABLE :: qgm(:)  ! complete fourier transform of Q
    !
    IF ( .NOT.okvan ) RETURN
    !
    ALLOCATE (qgm( ngm))
    ALLOCATE ( aux ( ngm, nspin) )    
    ALLOCATE ( qmod( ngm) )    
    ALLOCATE ( ylmk0( ngm, lmaxq * lmaxq ) )    
    !
    aux (:,:) = (0.d0, 0.d0)
    !
    CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
    !
    DO ig = 1, ngm
       !
       qmod (ig) = sqrt (gg (ig) )
       !
    END DO
    !
    DO nt = 1, ntyp
       !
       IF ( upf(nt)%tvanp ) then
          !
          ijh = 0
          ! 
          DO ih = 1, nh (nt)
             !
             DO jh = ih, nh (nt)
                !
#ifdef DEBUG_ADDUSDENS
                CALL start_clock ('addus:qvan2')
#endif
                !
                CALL qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
                !
#ifdef DEBUG_ADDUSDENS
                CALL stop_clock ('addus:qvan2')
#endif
                !
                ijh = ijh + 1
                !
                DO na = 1, nat
                   !
                   IF (ityp (na) .eq.nt) then
                      !
                      !  Multiply becsum and qg with the correct structure factor
                      !
#ifdef DEBUG_ADDUSDENS
                      call start_clock ('addus:aux')
#endif
                      !
                      DO is = 1, nspin
                         !
                         DO ig = 1, ngm
                            !
                            skk = eigts1 (ig1 (ig), na) * &
                                  eigts2 (ig2 (ig), na) * &
                                  eigts3 (ig3 (ig), na)
                            aux(ig,is)=aux(ig,is) + qgm(ig)*skk*becsum(ijh,na,is)
                            !
                         END DO
                         !
                      END DO
                      !
#ifdef DEBUG_ADDUSDENS
                      CALL stop_clock ('addus:aux')
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
    DEALLOCATE (ylmk0)
    DEALLOCATE (qmod)
    !
    !     convert aux to real space and add to the charge density
    !
    psic(:) = (0.d0, 0.d0)
    psic( nl(:) ) = aux(:,1)
    !
    IF (gamma_only) psic( nlm(:) ) = CONJG(aux(:,1))
    !
    CALL cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
    !
    rho_1 (:) = rho_1 (:) + REAL(psic(:),dp)
    !
    DEALLOCATE (aux)
    DEALLOCATE (qgm)
    !
    RETURN
    !
  END SUBROUTINE lr_addusdens_g
  !
  !------------------------------------------------------------------------
  SUBROUTINE lr_addusdens_r()
    !------------------------------------------------------------------------
    !
    ! ... Real space version
    !
    USE cell_base,        ONLY : omega
    USE uspp_param,       ONLY : upf, nh
    USE realus,           ONLY : box, maxbox, qsave
    !
    IMPLICIT NONE
    !
    INTEGER  :: ia, nt, ir, irb, ih, jh, ijh, is, mbia, iqs
    REAL(DP) :: charge
    !
    !
    DO is = 1, nspin
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
          ijh = 0
          !
          DO ih = 1, nh(nt)
             DO jh = ih, nh(nt)
                !
                ijh = ijh + 1
                !
                DO ir = 1, mbia
                   !
                   irb = box(ir,ia)
                   iqs = iqs + 1
                   !
                   rho_1(irb) = rho_1(irb) + qsave(iqs)*becsum(ijh,ia,is)
                END DO
             END DO
          END DO
       END DO
       !
    END DO
    !
    RETURN
    !
   END SUBROUTINE lr_addusdens_r
   !

END SUBROUTINE lr_addusdens
