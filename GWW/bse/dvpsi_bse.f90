!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE dvpsi_e(kpoint,ipol,dvpsi2,l_lr)
  !----------------------------------------------------------------------
  ! MARGHE: DA RIPULIRE E OTTIMIZZARE
  ! Calculates x * psi_k  for each k-points and for the 3 polarizations
  ! Requires on input: vkb, evc, igk
  !
  USE ions_base, ONLY : ntyp => nsp, nat, ityp
  USE kinds, ONLY: DP
  USE pwcom
  USE cell_base, ONLY : tpiba, tpiba2
  USE uspp, ONLY: nkb, vkb, dvan
  USE uspp_param, ONLY: nh
  USE wavefunctions_module,  ONLY: evc
  USE wvfct,    ONLY : nbnd, npwx,et
  USE becmod, ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files, ONLY: nwordwfc,iunwfc
  USE lsda_mod,              ONLY : nspin
  USE control_flags,    ONLY : gamma_only
  USE io_global,   ONLY : ionode,stdout
  USE gvect,       ONLY : gstart, ngm, g
  USE mp, ONLY : mp_sum, mp_barrier
  USE mp_world,             ONLY : world_comm
  use bse_wannier, ONLY:num_nbndv
  USE gvecw,              ONLY :  ecutwfc
  USE klist, ONLY : igk_k
!  USE cgcom
  !
  IMPLICIT NONE

  LOGICAL :: l_lr!if true calculates first-order valence wave-functions

  INTEGER :: kpoint, ipol,is,niter_ph
  INTEGER :: i,l, na,nt, ibnd,jbnd, info, ih,jkb, iter
  real(DP) :: upol(3,3),tr2_ph
  real(DP), ALLOCATABLE :: gk(:,:), q(:), overlap(:,:), &
       becp_(:,:), dbec(:,:), dbec_(:,:)
  real(DP), ALLOCATABLE :: becp2_(:,:), dbec2(:,:), dbec2_(:,:)
  COMPLEX(DP), ALLOCATABLE :: dvkb(:,:), dvkb1(:,:), work(:,:), &
       &           gr(:,:), h(:,:)
  COMPLEX(DP), ALLOCATABLE:: evc2(:,:),dpsi2(:,:)
  COMPLEX(DP) :: dvpsi(npwx , nbnd)
  COMPLEX(DP) :: dvpsi2(npwx , num_nbndv(1))
  
  LOGICAL:: precondition, orthonormal,startwith0,debug
  EXTERNAL H_h
  data upol /1.0d0,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 0.0d0,0.0d0,1.0d0/
  real(kind=DP), allocatable :: omat(:,:)
  integer j

  !
  CALL start_clock('dvpsi_e')

  CALL gk_sort (xk(1,kpoint),ngm,g,ecutwfc/tpiba2,npw,igk_k(1,1),g2kin)
  CALL init_us_2 (npw, igk_k(1,1), xk(1,kpoint), vkb)
  
  gamma_only=.true.
  debug=.true.


  allocate(omat(nbnd,nbnd))

  allocate( evc( npwx, nbnd ) )
  allocate( evc2( npwx, num_nbndv(1) ) )
  do is=1,nspin
     call davcio(evc,2*nwordwfc,iunwfc,is,-1) 
  enddo
  evc2(1:npwx,1:num_nbndv(1))=evc(1:npwx, 1:num_nbndv(1))

  !
  !   becp contains <beta|psi> - used in H_h
  !
  CALL allocate_bec_type ( nkb, nbnd, becp )
  ALLOCATE ( gk   ( 3, npwx) )
  ALLOCATE ( dvkb ( npwx, nkb) )
  ALLOCATE ( dvkb1( npwx, nkb) )
  ALLOCATE ( becp_(nkb,nbnd), dbec ( nkb, nbnd), dbec_(nkb, nbnd) )
  ALLOCATE ( becp2_(nkb,num_nbndv(1)), dbec2 ( nkb, num_nbndv(1)), dbec2_(nkb, num_nbndv(1)) )
  ALLOCATE ( dpsi2 (npwx , num_nbndv(1)))
  !
  DO i = 1,npw
     gk(1,i) = (xk(1,kpoint)+g(1,igk_k(i,1)))*tpiba
     gk(2,i) = (xk(2,kpoint)+g(2,igk_k(i,1)))*tpiba
     gk(3,i) = (xk(3,kpoint)+g(3,igk_k(i,1)))*tpiba
     g2kin(i)= gk(1,i)**2 + gk(2,i)**2 + gk(3,i)**2
  ENDDO
  !
  !  this is  the kinetic contribution to [H,x]:  -2i (k+G)_ipol * psi
  !

  dpsi2(1:npwx,1:num_nbndv(1))=(0.d0,0.d0)  
  
  DO ibnd = 1,num_nbndv(1)
     DO i = 1,npw
        dpsi2(i,ibnd) = gk(ipol,i)*(0.0d0,-2.0d0) * evc2(i,ibnd)
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


  CALL gen_us_dj(kpoint,dvkb)
  CALL gen_us_dy(kpoint,upol(1,ipol),dvkb1)

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
!                     if(ionode) write(*,*) 'dbec(j,ib)',  dbec(jkb,ibnd)
!                     if(ionode) write(*,*) 'dvan(ih,ih,nt)', dvan(ih,ih,nt)
                 becp_(jkb,ibnd) =becp%r(jkb,ibnd)*dvan(ih,ih,nt)
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO

 dbec2_(1:nkb,1:num_nbndv(1))=dbec_(1:nkb,1:num_nbndv(1)) 
 becp2_(1:nkb,1:num_nbndv(1))=becp_(1:nkb,1:num_nbndv(1)) 



  !
  IF (jkb/=nkb) CALL errore('dvpsi_e','unexpected error',1)
  !


  CALL dgemm ('N', 'N', 2*npw, num_nbndv(1), nkb,-1.d0, vkb, &
       2*npwx, dbec2_, nkb, 1.d0, dpsi2, 2*npwx)

  FLUSH( stdout )


  CALL dgemm ('N', 'N', 2*npw, num_nbndv(1), nkb, 1.d0,dvkb, &
       2*npwx, becp2_, nkb, 1.d0, dpsi2, 2*npwx)
  !
  DEALLOCATE(dbec, dbec_, becp_)
  DEALLOCATE(dbec2, dbec2_, becp2_)
  DEALLOCATE(dvkb1)
  DEALLOCATE(dvkb)
  DEALLOCATE(gk)

  !
  !   dpsi contains now [H,x] psi_v  for the three cartesian polarizations.
  !   Now solve the linear systems (H-e_v)*(x*psi_v) = [H,x]*psi_v
  !
!  ALLOCATE  ( overlap( nbnd, nbnd))
!  ALLOCATE  ( work(npwx, nbnd))
!  ALLOCATE  ( gr( npwx, nbnd))
!  ALLOCATE  ( h ( npwx, nbnd))
!  ALLOCATE  ( q ( npwx))
  !
  ALLOCATE  ( overlap( num_nbndv(1), num_nbndv(1)))
  ALLOCATE  ( work(npwx, num_nbndv(1)))
  ALLOCATE  ( gr( npwx, num_nbndv(1)))
  ALLOCATE  ( h ( npwx, num_nbndv(1)))
  ALLOCATE  ( q ( npwx))
!
  orthonormal = .false.
  precondition= .true.
  !
  IF (precondition) THEN
     DO i = 1,npw
        q(i) = 1.0d0/max(1.d0,g2kin(i))
     ENDDO
     CALL zvscal(npw,npwx,num_nbndv(1),q,evc2,work)
     CALL calbec ( npw, work, evc2, overlap)
     CALL DPOTRF('U',num_nbndv(1),overlap,num_nbndv(1),info)
     IF (info/=0) CALL errore('solve_ph','cannot factorize',info)
  ENDIF
  !
  startwith0= .true.
  dvpsi2(:,:) = (0.d0, 0.d0)
  niter_ph = 50
  tr2_ph = 1.0d-12
  !
  CALL cgsolve (npw,evc2,npwx,num_nbndv(1),overlap,num_nbndv(1),   &
       orthonormal,precondition,q,startwith0,et(1,1),&
       dpsi2,gr,h,dpsi2,work,niter_ph,tr2_ph,iter,dvpsi2)

  if(l_lr) then
     dpsi2(1:npwx,1:num_nbndv(1))=dvpsi2(1:npwx,1:num_nbndv(1))
     IF (precondition) THEN
        DO i = 1,npw
           q(i) = 1.0d0/max(1.d0,g2kin(i))
        ENDDO
        CALL zvscal(npw,npwx,num_nbndv(1),q,evc2,work)
        CALL calbec ( npw, work, evc2, overlap)
        CALL DPOTRF('U',num_nbndv(1),overlap,num_nbndv(1),info)
        IF (info/=0) CALL errore('solve_ph','cannot factorize',info)
     ENDIF
     startwith0= .true.
     dvpsi2(:,:) = (0.d0, 0.d0)
     niter_ph = 50
     tr2_ph = 1.0d-12
     CALL cgsolve (npw,evc2,npwx,num_nbndv(1),overlap,num_nbndv(1),   &
          orthonormal,precondition,q,startwith0,et(1,1),&
          dpsi2,gr,h,dpsi2,work,niter_ph,tr2_ph,iter,dvpsi2)

  endif

  DEALLOCATE(q)
  DEALLOCATE(h)
  DEALLOCATE(gr)
  DEALLOCATE(work)
  DEALLOCATE(overlap)
  DEALLOCATE(evc)
  DEALLOCATE(dpsi2)
  DEALLOCATE(evc2)
  DEALLOCATE(omat)

  CALL deallocate_bec_type ( becp )
  
  !
  CALL stop_clock('dvpsi_e')
  !
  RETURN
END SUBROUTINE dvpsi_e
