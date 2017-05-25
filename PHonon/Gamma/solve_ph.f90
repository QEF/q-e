!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE solve_ph ( )
  !-----------------------------------------------------------------------
  !
  USE io_global,             ONLY : stdout, ionode,ionode_id
  USE io_files,              ONLY : iunres, seqopn
  USE mp_world,              ONLY : world_comm
  USE mp,                    ONLY : mp_bcast
  USE uspp,                  ONLY : nkb
  USE wavefunctions_module,  ONLY : evc
  USE becmod,                ONLY : bec_type, becp, calbec, &
                                    allocate_bec_type, deallocate_bec_type
  USE klist,                 ONLY : xk, ngk
  USE wvfct,                 ONLY : nbnd, npwx, g2kin, et
  USE gvect,                 ONLY : g
  USE cell_base,             ONLY : tpiba2
  USE cgcom

  IMPLICIT NONE

  INTEGER :: npw, nu, i, ibnd, jbnd, info, iter, mode_done, ik
  REAL(DP), ALLOCATABLE ::  diag(:)
  COMPLEX(DP), ALLOCATABLE :: gr(:,:), h(:,:), work(:,:)
  REAL(DP), ALLOCATABLE :: overlap(:,:)
  LOGICAL :: orthonormal, precondition, startwith0, exst
  EXTERNAL A_h
  !
  CALL start_clock('solve_ph')
  !
  CALL allocate_bec_type ( nkb, nbnd, becp )
  ALLOCATE ( diag( npwx) )
  ALLOCATE ( overlap( nbnd, nbnd) )
  ALLOCATE ( work( npwx, nbnd) )
  ALLOCATE ( gr  ( npwx, nbnd) )
  ALLOCATE ( h   ( npwx, nbnd) )
  !
  ik = 1
  npw = ngk(ik)
  DO i = 1,npw
     g2kin(i) = ( g(1,i)**2 + g(2,i)**2 + g(3,i)**2 ) * tpiba2
  ENDDO
  !
  orthonormal = .false.
  precondition= .true.
  !
  IF (precondition) THEN
     DO i = 1,npw
        diag(i) = 1.0d0/max(1.d0,g2kin(i))
     ENDDO
     CALL zvscal(npw,npwx,nbnd,diag,evc,work)
     CALL calbec (npw, work, evc, overlap)
     CALL DPOTRF('U',nbnd,overlap,nbnd,info)
     IF (info/=0) CALL errore('solve_ph','cannot factorize',info)
  ENDIF
  !
  WRITE( stdout,'(/" ***  Starting Conjugate Gradient minimization",   &
       &            9x,"***")')
  !
  !  check if a restart file exists
  !
  IF (recover) THEN
     IF (ionode) CALL seqopn( iunres, 'restartph', 'FORMATTED', exst )
     CALL mp_bcast(exst,ionode_id,world_comm)
     IF (.not. exst) GOTO 1
     IF (ionode) THEN
        READ (iunres,*,err=1,END=1) mode_done
        READ (iunres,*,err=1,END=1) dyn
        CLOSE(unit=iunres)
     END IF
     CALL mp_bcast(mode_done,ionode_id,world_comm)
     CALL mp_bcast(dyn,ionode_id,world_comm)
     PRINT '("  Phonon: modes up to mode ",i3," already done")', mode_done
     GOTO 2
1    CLOSE(unit=iunres)
  ENDIF
  !  initialisation if not restarting from previous calculation
  CALL  dynmat_init
  mode_done=0
2 CONTINUE
  !
  DO nu = 1, nmodes
     IF ( has_equivalent((nu-1)/3+1)==1) THEN
        ! calculate only independent modes
        WRITE( stdout,'(" ***  mode # ",i3," : using symmetry")') nu
        GOTO 10
     ENDIF
     IF ( nu<=mode_done) THEN
        ! do not recalculate modes already done
        WRITE( stdout,'(" ***  mode # ",i3," : using previous run")') nu
        GOTO 10
     ENDIF
     IF ( asr .and. (nu-1)/3+1==nasr ) THEN
        ! impose ASR on last atom instead of calculating mode
        WRITE( stdout,'(" ***  mode # ",i3," : using asr")') nu
        GOTO 10
     ENDIF
     ! calculate |b> = dV/dtau*psi
     CALL dvpsi_kb(ik,nu)
     ! initialize delta psi
     startwith0=.true.
     dpsi(:,:) = (0.d0, 0.d0)
     ! solve the linear system
     ! NB: dvpsi is used also as work space and is destroyed by cgsolve
     CALL cgsolve (A_h,npw,evc,npwx,nbnd,overlap,nbnd, &
                   orthonormal,precondition,diag,      &
                   startwith0,et(1,ik),dvpsi,gr,h, &
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
     ENDIF
     !
     WRITE( stdout,'(" ***  mode # ",i3," : ",i3," iterations")')  &
          &          nu, iter
10   CONTINUE
  ENDDO
  !
  DEALLOCATE(h)
  DEALLOCATE(gr)
  DEALLOCATE(overlap)
  DEALLOCATE(work)
  DEALLOCATE(diag)
  CALL deallocate_bec_type (becp)
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

  IF (nasr<=0 .or. nasr>nat) RETURN
  DO j=1,3
     DO i=1,3
        DO nb=1,nat
           sum=0.d0
           DO na=1,nat
              IF (na/=nasr) sum = sum + dyn(3*(na-1)+i,3*(nb-1)+j)
           ENDDO
           dyn(3*(nasr-1)+i,3*(nb-1)+j)= -sum
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE set_asr
