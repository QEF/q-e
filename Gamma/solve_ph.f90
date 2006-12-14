!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE solve_ph ( )
  !-----------------------------------------------------------------------
  !
  USE io_global,             ONLY : stdout, ionode
  USE io_files,              ONLY : iunres
  USE pwcom
  USE wavefunctions_module,  ONLY : evc
  USE becmod,                ONLY : rbecp
  USE cgcom

  IMPLICIT NONE

  INTEGER :: nu, i, ibnd, jbnd, info, iter, mode_done, kpoint
  REAL(DP), ALLOCATABLE ::  diag(:)
  COMPLEX(DP), ALLOCATABLE :: gr(:,:), h(:,:), work(:,:)
  REAL(DP), ALLOCATABLE :: overlap(:,:)
  LOGICAL :: orthonormal, precondition, startwith0, exst
  EXTERNAL A_h
  !
  CALL start_clock('solve_ph')
  !
  ALLOCATE ( rbecp( nkb,nbnd) )
  ALLOCATE ( diag( npwx) )
  ALLOCATE ( overlap( nbnd, nbnd) )
  ALLOCATE ( work( npwx, nbnd) )
  ALLOCATE ( gr  ( npwx, nbnd) )
  ALLOCATE ( h   ( npwx, nbnd) )
  !
  kpoint = 1
  DO i = 1,npw
     g2kin(i) = ( (xk(1,kpoint)+g(1,igk(i)))**2 +                   &
                  (xk(2,kpoint)+g(2,igk(i)))**2 +                   &
                  (xk(3,kpoint)+g(3,igk(i)))**2 ) * tpiba2
  END DO
  !
  orthonormal = .FALSE.
  precondition= .TRUE.
  !
  IF (precondition) THEN
     DO i = 1,npw
        diag(i) = 1.0d0/MAX(1.d0,g2kin(i))
     END DO
     CALL zvscal(npw,npwx,nbnd,diag,evc,work)
     CALL pw_gemm ('Y',nbnd, nbnd, npw, work, npwx, evc, npwx, overlap, nbnd)
     CALL DPOTRF('U',nbnd,overlap,nbnd,info)
     IF (info.NE.0) CALL errore('solve_ph','cannot factorize',info)
  END IF
  !
  WRITE( stdout,'(/" ***  Starting Conjugate Gradient minimization",   &
       &            9x,"***")')
  !
  !  check if a restart file exists
  !
  if (recover) THEN
     CALL seqopn( iunres, 'restartph', 'FORMATTED', exst )
     IF (.NOT. exst) GO TO 1
     READ (iunres,*,err=1,END=1) mode_done
     READ (iunres,*,err=1,END=1) dyn
     CLOSE(unit=iunres)
     PRINT '("  Phonon: modes up to mode ",i3," already done")', mode_done
     go to 2
1    CLOSE(unit=iunres)
  END IF
  !  initialisation if not restarting from previous calculation
  CALL  dynmat_init
  mode_done=0
2 CONTINUE
  !
  DO nu = 1, nmodes
     IF ( has_equivalent((nu-1)/3+1).EQ.1) THEN
        ! calculate only independent modes
        WRITE( stdout,'(" ***  mode # ",i3," : using symmetry")') nu
        GOTO 10
     END IF
     IF ( nu.LE.mode_done) THEN
        ! do not recalculate modes already done
        WRITE( stdout,'(" ***  mode # ",i3," : using previous run")') nu
        GOTO 10
     END IF
     IF ( asr .AND. (nu-1)/3+1.EQ.nasr ) THEN
        ! impose ASR on last atom instead of calculating mode
        WRITE( stdout,'(" ***  mode # ",i3," : using asr")') nu
        GOTO 10
     END IF
     ! calculate |b> = dV/dtau*psi
     CALL dvpsi_kb(kpoint,nu)
     ! initialize delta psi
     startwith0=.TRUE.
     dpsi(:,:) = (0.d0, 0.d0)
     ! solve the linear system
     ! NB: dvpsi is used also as work space and is destroyed by cgsolve
     CALL cgsolve (A_h,npw,evc,npwx,nbnd,overlap,nbnd, &
                   orthonormal,precondition,diag,      &
                   startwith0,et(1,kpoint),dvpsi,gr,h, &
                   dvpsi,work,niter_ph,tr2_ph,iter,dpsi)
     ! < DeltaPsi | DeltaV | Psi > contribution to the dynamical matrix
     CALL drhodv(nu)
     ! save partial result
     !
     IF ( ionode ) THEN
        !
        CALL seqopn( iunres, 'restartph', 'FORMATTED', exst )
        WRITE(iunres,*) nu
        WRITE(iunres,*) dyn
        CLOSE(unit=iunres)
        !
     END IF
     !
     WRITE( stdout,'(" ***  mode # ",i3," : ",i3," iterations")')  &
          &          nu, iter
10   CONTINUE
  END DO
  !
  DEALLOCATE(h)
  DEALLOCATE(gr)
  DEALLOCATE(overlap)
  DEALLOCATE(work)
  DEALLOCATE(diag)
  DEALLOCATE(rbecp)
  !
  CALL stop_clock('solve_ph')
  !
  RETURN
END SUBROUTINE solve_ph
!
!---------------------------------------------------------------------------
SUBROUTINE set_asr(nat,nasr,dyn)
  !---------------------------------------------------------------------------
  !
  ! Impose Acoustic Sum Rule on the dynamical matrix
  ! We assume that (3*nat-1) columns have been calculated
  ! and that the missing column corresponds to atom nasr
  !
  IMPLICIT NONE
  INTEGER nat, nasr
  REAL(8) :: dyn(3*nat,3*nat)
  !
  INTEGER na, nb, i,j
  REAL(8) :: sum

  IF (nasr.LE.0 .OR. nasr.GT.nat) RETURN
  DO j=1,3
     DO i=1,3
        DO nb=1,nat
           sum=0.d0
           DO na=1,nat
              IF (na.NE.nasr) sum = sum + dyn(3*(na-1)+i,3*(nb-1)+j)
           END DO
           dyn(3*(nasr-1)+i,3*(nb-1)+j)= -sum
        END DO
     END DO
  END DO

  RETURN
END SUBROUTINE set_asr
