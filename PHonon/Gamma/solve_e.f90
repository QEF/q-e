!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE solve_e
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  USE io_files,  ONLY : seqopn
  USE uspp,      ONLY : nkb
  USE becmod,    ONLY : bec_type, becp, calbec, allocate_bec_type, &
                        deallocate_bec_type
  USE cell_base, ONLY : tpiba2
  USE gvect,     ONLY : g
  USE klist,     ONLY : xk
  USE wvfct,     ONLY : nbnd, npwx, npw, g2kin, et
  USE wavefunctions_module,  ONLY: evc
  USE cgcom
  !
  IMPLICIT NONE
  !
  INTEGER :: ipol, nrec, i, ibnd, jbnd, info, iter, ik
  real(DP), ALLOCATABLE ::diag(:)
  COMPLEX(DP), ALLOCATABLE :: gr(:,:), h(:,:), work(:,:)
  real(DP), ALLOCATABLE :: overlap(:,:)
  LOGICAL :: orthonormal, precondition,startwith0,here
  CHARACTER(len=7) :: fildwf, filbar
  EXTERNAL A_h
  !
  CALL start_clock('solve_e')
  !
  CALL allocate_bec_type ( nkb, nbnd, becp)
  ALLOCATE ( diag( npwx) )
  ALLOCATE ( overlap( nbnd, nbnd) )
  ALLOCATE ( work( npwx, nbnd) )
  ALLOCATE ( gr  ( npwx, nbnd) )
  ALLOCATE ( h   ( npwx, nbnd) )
  !
  ik = 1
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
     IF (info/=0) CALL errore('solve_e','cannot factorize',info)
  ENDIF
  !
  WRITE( stdout,'(/" ***  Starting Conjugate Gradient minimization",   &
       &            9x,"***")')
  nrec=0
  !
  DO ipol = 1,3
     !  read |b> = dV/dtau*psi
     iubar=ipol
     WRITE(filbar,'("filbar",i1)') ipol
     CALL seqopn (iubar,filbar,'unformatted',here)
     IF (.not.here) CALL errore('solve_e','file '//filbar//          &
          &        'mysteriously vanished',ipol)
     READ (iubar) dvpsi
     CLOSE(unit=iubar,status='keep')
     !
     iudwf=10+ipol
     WRITE(fildwf,'("fildwx",i1)') ipol
     CALL  seqopn (iudwf,fildwf,'unformatted',here)
!!!         if (.not.here) then
     !  calculate Delta*psi  (if not already done)
     dpsi(:,:) = (0.d0, 0.d0)
     startwith0= .true.
!!!         else
     !  otherwise restart from Delta*psi that is found on file
!!!            read(iudwf) dpsi
!!!         end if
     CALL cgsolve (A_h,npw,evc,npwx,nbnd,overlap,nbnd, &
                   orthonormal,precondition,diag,      &
                   startwith0,et(1,ik),dvpsi,gr,h, &
                   dvpsi,work,niter_ph,tr2_ph,iter,dpsi)
     !  write Delta*psi for an electric field
     REWIND (iudwf)
     WRITE (iudwf) dpsi
     CLOSE(unit=iudwf)
     !
     WRITE( stdout,'(" ***  pol. # ",i3," : ",i3," iterations")')  &
          &              ipol, iter
  ENDDO
  !
  DEALLOCATE(h)
  DEALLOCATE(gr)
  DEALLOCATE(overlap)
  DEALLOCATE(work)
  DEALLOCATE(diag)
  CALL deallocate_bec_type (becp)
  !
  CALL stop_clock('solve_e')
  !
  RETURN
END SUBROUTINE solve_e
