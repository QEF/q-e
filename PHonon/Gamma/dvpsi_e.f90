!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE dvpsi_e(ik,ipol)
  !----------------------------------------------------------------------
  !
  ! Calculates x * psi_k  for each k-point and for the 3 polarizations
  ! Requires on input: vkb, evc
  !
  USE kinds, ONLY: DP
  USE ions_base, ONLY : ntyp => nsp, nat, ityp
  USE uspp, ONLY: nkb, vkb, dvan
  USE uspp_param, ONLY: nh
  USE wavefunctions_module,  ONLY: evc
  USE becmod, ONLY: bec_type, becp, calbec, allocate_bec_type, &
      deallocate_bec_type
  USE cell_base, ONLY : tpiba
  USE gvect,     ONLY : g
  USE klist,     ONLY : xk, ngk
  USE wvfct,     ONLY : nbnd, npwx,g2kin, et
  USE cgcom
  !
  IMPLICIT NONE
  INTEGER :: ik, ipol
  INTEGER :: npw, i,l, na,nt, ibnd,jbnd, info, ih,jkb, iter
  real(DP) :: upol(3,3)
  real(DP), ALLOCATABLE :: gk(:,:), q(:), overlap(:,:), &
       becp_(:,:), dbec(:,:), dbec_(:,:)
  COMPLEX(DP), ALLOCATABLE :: dvkb(:,:), dvkb1(:,:), work(:,:), &
       &           gr(:,:), h(:,:)
  LOGICAL:: precondition, orthonormal,startwith0
  EXTERNAL H_h
  data upol /1.0d0,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 0.0d0,0.0d0,1.0d0/
  !
  CALL start_clock('dvpsi_e')
  !
  !   becp contains <beta|psi> - used in H_h
  !
  CALL allocate_bec_type ( nkb, nbnd, becp )
  ALLOCATE ( gk   ( 3, npwx) )
  ALLOCATE ( dvkb ( npwx, nkb) )
  ALLOCATE ( dvkb1( npwx, nkb) )
  ALLOCATE ( becp_(nkb,nbnd), dbec ( nkb, nbnd), dbec_(nkb, nbnd) )
  !
  !  g2kin is used in H_h, called by cgsolve below
  !
  npw = ngk(ik)
  DO i = 1,npw
     gk(1,i) = g(1,i)*tpiba
     gk(2,i) = g(2,i)*tpiba
     gk(3,i) = g(3,i)*tpiba
     g2kin(i)= gk(1,i)**2 + gk(2,i)**2 + gk(3,i)**2
  ENDDO
  !
  !  this is  the kinetic contribution to [H,x]:  -2i (k+G)_ipol * psi
  !
  DO ibnd = 1,nbnd
     DO i = 1,npw
        dpsi(i,ibnd) = gk(ipol,i)*(0.0d0,-2.0d0) * evc(i,ibnd)
     ENDDO
  ENDDO
  !
  DO i = 1,npw
     IF (g2kin(i)>1.0d-10) THEN
        gk(1,i) = gk(1,i)/sqrt(g2kin(i))
        gk(2,i) = gk(2,i)/sqrt(g2kin(i))
        gk(3,i) = gk(3,i)/sqrt(g2kin(i))
     ENDIF
  ENDDO
  !
  ! and these are the contributions from nonlocal pseudopotentials
  ! ( upol(3,3) are the three unit vectors along x,y,z)
  !
  CALL gen_us_dj(ik,dvkb)
  CALL gen_us_dy(ik,upol(1,ipol),dvkb1)
  !
  DO jkb = 1, nkb
     DO i = 1,npw
        dvkb(i,jkb) =(0.d0,-1.d0)*(dvkb1(i,jkb) + dvkb(i,jkb)*gk(ipol,i))
     ENDDO
  ENDDO
  !
  CALL calbec ( npw,  vkb, evc,  becp )
  CALL calbec ( npw, dvkb, evc,  dbec )
  !
  jkb = 0
  DO nt=1, ntyp
     DO na = 1,nat
        IF (nt==ityp(na)) THEN
           DO ih=1,nh(nt)
              jkb=jkb+1
              DO ibnd = 1,nbnd
                 dbec_(jkb,ibnd) = dbec(jkb,ibnd)*dvan(ih,ih,nt)
                 becp_(jkb,ibnd) =becp%r(jkb,ibnd)*dvan(ih,ih,nt)
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  IF (jkb/=nkb) CALL errore('dvpsi_e','unexpected error',1)
  !
  CALL dgemm ('N', 'N', 2*npw, nbnd, nkb,-1.d0, vkb, &
       2*npwx, dbec_, nkb, 1.d0, dpsi, 2*npwx)
  CALL dgemm ('N', 'N', 2*npw, nbnd, nkb, 1.d0,dvkb, &
       2*npwx, becp_, nkb, 1.d0, dpsi, 2*npwx)
  !
  DEALLOCATE(dbec, dbec_, becp_)
  DEALLOCATE(dvkb1)
  DEALLOCATE(dvkb)
  DEALLOCATE(gk)
  !
  !   dpsi contains now [H,x] psi_v  for the three cartesian polarizations.
  !   Now solve the linear systems (H-e_v)*(x*psi_v) = [H,x]*psi_v
  !
  ALLOCATE  ( overlap( nbnd, nbnd))
  ALLOCATE  ( work( npwx, nbnd))
  ALLOCATE  ( gr( npwx, nbnd))
  ALLOCATE  ( h ( npwx, nbnd))
  ALLOCATE  ( q ( npwx))
  !
  orthonormal = .false.
  precondition= .true.
  !
  IF (precondition) THEN
     DO i = 1,npw
        q(i) = 1.0d0/max(1.d0,g2kin(i))
     ENDDO
     CALL zvscal(npw,npwx,nbnd,q,evc,work)
     CALL calbec ( npw, work, evc, overlap )
     CALL DPOTRF('U',nbnd,overlap,nbnd,info)
     IF (info/=0) CALL errore('solve_ph','cannot factorize',info)
  ENDIF
  !
  startwith0= .true.
  dvpsi(:,:) = (0.d0, 0.d0)
  !
  CALL cgsolve (H_h,npw,evc,npwx,nbnd,overlap,nbnd,   &
       orthonormal,precondition,q,startwith0,et(1,ik),&
       dpsi,gr,h,dpsi,work,niter_ph,tr2_ph,iter,dvpsi)
  !
  DEALLOCATE(q)
  DEALLOCATE(h)
  DEALLOCATE(gr)
  DEALLOCATE(work)
  DEALLOCATE(overlap)
  CALL deallocate_bec_type ( becp )
  !
  CALL stop_clock('dvpsi_e')
  !
  RETURN
END SUBROUTINE dvpsi_e
