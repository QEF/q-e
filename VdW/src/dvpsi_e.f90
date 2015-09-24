!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE dvpsi_e_vdw (kpoint, ipol)
  !----------------------------------------------------------------------
  !
  ! On output: dvpsi contains P_c^+ x | psi_kpoint > in crystal axis
  !            (projected on at(*,ipol) )
  !
  ! dvpsi is READ from file if this_pcxpsi_is_on_file(kpoint,ipol)=.true.
  ! otherwise dvpsi is COMPUTED and WRITTEN on file (vkb,evc,igk must be set)
  !
  !
  USE ions_base,             ONLY : nat, ityp, ntyp => nsp
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : DP
  USE becmod,                ONLY : becp, calbec
  USE uspp,                  ONLY : okvan, nkb, vkb, qq, deeq
  USE uspp_param,            ONLY : nh
  USE eff_v,                 ONLY : dvext, evc => evc_veff
  USE phcom
  USE pwcom
  USE mp_global,             ONLY : intra_pool_comm
  USE mp,                    ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ipol, kpoint
  !
  ! Local variables
  !
  INTEGER :: ig, na, ibnd, jbnd, ikb, jkb, nt, lter, ih, jh, ijkb0, nrec
  ! counters

  real(kind=DP), ALLOCATABLE  :: gk (:,:), h_diag (:,:),  eprec1 (:)
  ! the derivative of |k+G|
  real(kind=DP) ::   anorm, thresh
  ! preconditioning cut-off
  ! the desired convergence of linter

  LOGICAL :: conv_root
  ! true if convergence has been achieved

  COMPLEX(kind=DP), ALLOCATABLE :: ps (:,:), dvkb (:,:), dvkb1 (:,:), &
       work (:,:), becp2(:,:), spsi(:,:)
  COMPLEX(kind=DP), EXTERNAL :: zdotc
  ! the scalar products
  EXTERNAL ch_psi_all, cg_psi
  !
  CALL start_clock ('dvpsi_e')
  IF (this_pcxpsi_is_on_file(kpoint,ipol)) THEN
     !
     CALL zcopy (npwx,  dvext (1, ipol, 1), 1, dvpsi (1, 1), 1)
     CALL stop_clock ('dvpsi_e')
     !
     RETURN
     !
  ENDIF
  !
  nh = 0   ! working for norm-conseving pp only
!  nkb = 0
  !
  IF (nkb > 0) THEN
     ALLOCATE (work ( npwx, nkb))
  ELSE
     ALLOCATE (work ( npwx, 1))
  ENDIF
  ALLOCATE (gk ( 3, npwx))
  ALLOCATE (h_diag( npwx , nbnd))
  ALLOCATE (ps ( 2 , nbnd))
  ALLOCATE (spsi ( npwx, nbnd))
  ALLOCATE (eprec1 ( nbnd))
  IF (nkb > 0) THEN
     ALLOCATE (becp2 (nkb, nbnd), dvkb (npwx, nkb), dvkb1(npwx, nkb))
     dvkb (:,:) = (0.d0, 0.d0)
     dvkb1(:,:) = (0.d0, 0.d0)
     dvpsi(:,:) = (0.d0, 0.d0)
  ELSE
     ALLOCATE (dvkb1(0,0))
  ENDIF
  DO ig = 1, npw
     gk (1, ig) = (xk (1, kpoint) + g (1, igk (ig) ) ) * tpiba
     gk (2, ig) = (xk (2, kpoint) + g (2, igk (ig) ) ) * tpiba
     gk (3, ig) = (xk (3, kpoint) + g (3, igk (ig) ) ) * tpiba
     g2kin (ig) = gk (1, ig) **2 + gk (2, ig) **2 + gk (3, ig) **2
  ENDDO
  !
  ! this is  the kinetic contribution to [H,x]:  -2i (k+G)_ipol * psi
  !
  DO ibnd = 1, nbnd_occ (kpoint)
     DO ig = 1, npw
        dpsi (ig, ibnd) =  (at(1, ipol) * gk(1, ig) + &
                            at(2, ipol) * gk(2, ig) + &
                            at(3, ipol) * gk(3, ig) ) &
                          *(0.d0,-2.d0)*evc (ig, ibnd)
     ENDDO
  ENDDO
  !
  ! and this is the contribution from nonlocal pseudopotentials
  !
!  call gen_us_dj (kpoint, dvkb)        !this call is not needed bcs nkb=0
  CALL gen_us_dy (kpoint, at (1, ipol), dvkb1)
  DO ig = 1, npw
     IF (g2kin (ig) < 1.0d-10) THEN
        gk (1, ig) = 0.d0
        gk (2, ig) = 0.d0
        gk (3, ig) = 0.d0
     ELSE
        gk (1, ig) = gk (1, ig) / sqrt (g2kin (ig) )
        gk (2, ig) = gk (2, ig) / sqrt (g2kin (ig) )
        gk (3, ig) = gk (3, ig) / sqrt (g2kin (ig) )
     ENDIF
  ENDDO

  jkb = 0
  DO nt = 1, ntyp
     DO na = 1, nat
        IF (nt == ityp (na)) THEN
           DO ikb = 1, nh (nt)
              jkb = jkb + 1
              DO ig = 1, npw
                 work (ig,jkb) = dvkb1 (ig, jkb) + dvkb (ig, jkb) * &
                      (at (1, ipol) * gk (1, ig) + &
                      at (2, ipol) * gk (2, ig) + &
                      at (3, ipol) * gk (3, ig) )
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  IF (nkb > 0) CALL calbec (npw, work, evc, becp2)

  ijkb0 = 0
  DO nt = 1, ntyp
     DO na = 1, nat
        IF (nt == ityp (na)) THEN
           DO ih = 1, nh (nt)
              ikb = ijkb0 + ih
              ps(:,:)=(0.d0,0.d0)
              DO jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 DO ibnd = 1, nbnd_occ (kpoint)
                    ps (1, ibnd) = ps(1,ibnd)+ becp2(jkb,ibnd)* &
                         (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin) &
                         -et(ibnd,kpoint)*qq(ih,jh,nt))
                    ps (2, ibnd) = ps(2,ibnd) +becp1(kpoint)%k(jkb,ibnd) * &
                         (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin)&
                         -et(ibnd,kpoint)*qq(ih,jh,nt))
                 ENDDO
              ENDDO
              DO ibnd = 1, nbnd_occ (kpoint)
                 CALL zaxpy(npw,ps(1,ibnd),vkb(1,ikb),1,dpsi(1,ibnd),1)
                 CALL zaxpy(npw,ps(2,ibnd),work(1,ikb),1,dpsi(1,ibnd),1)
              ENDDO
           ENDDO
           ijkb0=ijkb0+nh(nt)
        ENDIF
     ENDDO
  ENDDO
  IF (jkb /= nkb) CALL errore ('dvpsi_e', 'unexpected error', 1)
  !
  !    orthogonalize dpsi to the valence subspace
  !
  DO ibnd = 1, nbnd_occ (kpoint)
     work (:,1) = (0.d0, 0.d0)
     DO jbnd = 1, nbnd_occ (kpoint)
        ps (1, jbnd) = - zdotc(npw,evc(1,jbnd),1,dpsi(1, ibnd),1)
     ENDDO
#ifdef __MPI
     CALL mp_sum( ps, intra_pool_comm )
#endif
     DO jbnd = 1, nbnd_occ (kpoint)
        CALL zaxpy (npw, ps (1, jbnd), evc (1, jbnd), 1, work, 1)
     ENDDO
     IF ( nkb > 0 ) CALL calbec (npw, vkb, work, becp, 1)
     CALL s_psi (npwx, npw, 1, work, spsi)
     CALL daxpy (2 * npw, 1.0d0, spsi, 1, dpsi (1, ibnd), 1)
  ENDDO
  !
  !   dpsi contains now P^+_c [H-eS,x] psi_v  for the three crystal
  !   polarizations
  !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
  !
  thresh = 1.d-5
  DO ibnd = 1, nbnd_occ (kpoint)
     conv_root = .true.
     DO ig = 1, npwq
        work (ig,1) = g2kin (ig) * evc (ig, ibnd)
     ENDDO
     eprec1 (ibnd) = 1.35d0 * zdotc (npwq, evc (1, ibnd), 1, work, 1)
  ENDDO
#ifdef __MPI
  CALL mp_sum( eprec1( 1 : nbnd_occ(kpoint) ), intra_pool_comm )
#endif
  DO ibnd = 1, nbnd_occ (kpoint)
     DO ig = 1, npwq
        h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec1 (ibnd) )
     ENDDO
  ENDDO
  !
  dvpsi = (0.d0,0.d0)
  !
  CALL cgsolve_all (ch_psi_all, cg_psi, et (1, kpoint), dpsi, dvpsi, &
       h_diag, npwx, npw, thresh, kpoint, lter, conv_root, anorm, &
       nbnd_occ (kpoint),1 )
  !
  IF (.not.conv_root) WRITE( stdout, '(5x,"kpoint",i4," ibnd",i4, &
       & " linter: root not converged ",e10.3)') &
       kpoint, ibnd, anorm
  !
#ifdef FLUSH
  FLUSH (6)
#endif
  !
  ! we have now obtained P_c x |psi>.
  ! In the case of USPP this quantity is needed for the Born
  ! effective charges, so we save it to disc
  !
  ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
  ! therefore we apply S again, and then subtract the additional term
  ! furthermore we add the term due to dipole of the augmentation charges.
  !
  IF (okvan) THEN
     !
     ! for effective charges
     !
!     nrec = (ipol - 1) * nksq + kpoint
!     call davcio (dvpsi, lrcom, iucom, nrec, 1)
     !
     IF (nkb > 0) CALL calbec (npw, vkb, dvpsi, becp)
     CALL s_psi(npwx,npw,nbnd,dvpsi,spsi)
     CALL dcopy(2*npwx*nbnd,spsi,1,dvpsi,1)
     CALL adddvepsi_us(becp1(kpoint),becp2,ipol,kpoint,dvpsi)
  ENDIF

  IF (nkb > 0) DEALLOCATE (dvkb1, dvkb, becp2)
  DEALLOCATE (eprec1)
  DEALLOCATE (spsi)
  DEALLOCATE (ps)
  DEALLOCATE (h_diag)
  DEALLOCATE (gk)
  DEALLOCATE (work)

  CALL zcopy (npwx, dvpsi (1, 1), 1, dvext (1, ipol, 1), 1)
  this_pcxpsi_is_on_file(kpoint,ipol) = .true.
  !
  CALL stop_clock ('dvpsi_e')
  !
  RETURN
  !
END SUBROUTINE dvpsi_e_vdw
