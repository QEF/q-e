!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE compute_ppsi (ppsi, ppsi_us, ik, ipol, nbnd_occ, current_spin)
  !----------------------------------------------------------------------
  !
  ! On output: ppsi contains P_c^+ p | psi_ik > for the ipol cartesian
  !            coordinate
  !            ppsi_us contains the additional term required for US PP.
  !            See J. Chem. Phys. 120, 9935 (2004) Eq. 10.
  !
  ! (important: vkb and evc must have been initialized for this k-point)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE io_global,            ONLY : stdout
  USE wvfct,                ONLY : nbnd, npwx
  USE uspp,                 ONLY : nkb, vkb, qq_nt, qq_so, okvan
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin, npol, lspinorb
  USE klist,                ONLY : ngk
  USE becmod,               ONLY : bec_type, becp, calbec, allocate_bec_type, &
                                   deallocate_bec_type
  USE uspp_param,           ONLY : nh, nhm
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ipol, ik, nbnd_occ, current_spin
  !
  COMPLEX(DP) :: ppsi(npwx,npol,nbnd_occ), ppsi_us(npwx,npol,nbnd_occ)
  ! Local variables
  !
  INTEGER :: npw, na, ibnd, ikb, jkb, nt, ih, jh, ip, ijkb0
  ! counters

  REAL(DP)  :: vers(3)

  COMPLEX(DP), ALLOCATABLE :: ps(:), ps_nc(:,:), dpqq_so(:,:,:,:,:)

  REAL(DP), ALLOCATABLE :: dpqq(:,:,:,:)
  !
  TYPE(bec_type) :: becp2
  !
  ! polarization vector in Cartesian coordinates
  vers = 0.d0
  vers(ipol) = 1.d0
  !
  CALL allocate_bec_type(nkb, nbnd, becp2)
  !
  CALL commutator_Hx_psi (ik, nbnd_occ, vers, becp, becp2, ppsi)
  !
  ! commutator_Hx_psi calculates [H, x], here we need i/2 [H, x]
  ppsi = ppsi * (0.0_DP, 0.5_DP)
  !
  !   ppsi contains p - i/2 [x, V_{nl}-eS] psi_v for the ipol polarization
  !
  !   In the US case there is another term in the matrix element.
  !   This term has to be multiplied by the difference of the eigenvalues,
  !   so it is calculated separately here and multiplied in the calling
  !   routine.
  !
  ! Note: The remaining part of this routine does a similar thing as
  !       LR_Modules/adddvepsi_us.f90.
  !
  npw = ngk(ik)
  !
  IF (okvan) THEN
     ppsi_us=(0.d0,0.d0)
     ALLOCATE (dpqq( nhm, nhm, 3, ntyp))
     CALL compute_qdipol(dpqq)
     IF (noncolin) THEN
        ALLOCATE (ps_nc(nbnd_occ,npol))
        IF (lspinorb) THEN
           ALLOCATE (dpqq_so( nhm, nhm, nspin, 3, ntyp))
           CALL compute_qdipol_so(dpqq, dpqq_so)
        ENDIF
     ELSE
        ALLOCATE (ps(nbnd_occ))
     ENDIF
     ijkb0 = 0
     DO nt = 1, ntyp
        DO na = 1, nat
           IF (ityp(na)==nt) THEN
              DO ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 IF (noncolin) THEN
                    ps_nc = (0.d0,0.d0)
                 ELSE
                    ps = (0.d0,0.d0)
                 ENDIF
                 DO jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    DO ibnd=1, nbnd_occ
                       IF (noncolin) THEN
                          DO ip=1,npol
                             IF (lspinorb) THEN
                                ps_nc(ibnd,ip)=ps_nc(ibnd,ip) +          &
                                    (0.d0,1.d0)*(becp2%nc(jkb,1,ibnd)*   &
                                    qq_so(ih,jh,1+(ip-1)*2,nt) +         &
                                    becp2%nc(jkb,2,ibnd) *               &
                                    qq_so(ih,jh,2+(ip-1)*2,nt) )         &
                                  + becp%nc(jkb,1,ibnd)*                 &
                                    dpqq_so(ih,jh,1+(ip-1)*2,ipol,nt)    &
                                  + becp%nc(jkb,2,ibnd)*                 &
                                    dpqq_so(ih,jh,2+(ip-1)*2,ipol,nt)
                             ELSE
                                ps_nc(ibnd,ip)=ps_nc(ibnd,ip)+           &
                                    becp2%nc(jkb,ip,ibnd)*(0.d0,1.d0)*   &
                                    qq_nt(ih,jh,nt)+becp%nc(jkb,ip,ibnd)    &
                                                   *dpqq(ih,jh,ipol,nt)
                             ENDIF
                          ENDDO
                       ELSE
                          ps(ibnd) = ps(ibnd) + becp2%k(jkb,ibnd) *  &
                                (0.d0,1.d0) * qq_nt(ih,jh,nt)   +  &
                                becp%k(jkb,ibnd) * dpqq(ih,jh,ipol,nt)
                       ENDIF
                    ENDDO
                 ENDDO
                 DO ibnd = 1, nbnd_occ
                    IF (noncolin) THEN
                       DO ip=1,npol
                          CALL zaxpy(npw,ps_nc(ibnd,ip),vkb(1,ikb),1,&
                                     ppsi_us(1,ip,ibnd),1)
                       ENDDO
                    ELSE
                       CALL zaxpy(npw,ps(ibnd),vkb(1,ikb),1,ppsi_us(1,1,ibnd),1)
                    ENDIF
                 ENDDO
              ENDDO
              ijkb0=ijkb0+nh(nt)
           ENDIF
        ENDDO
     ENDDO
     !
     IF (jkb/=nkb) CALL errore ('compute_ppsi', 'unexpected error', 1)
     !
     DEALLOCATE(dpqq)
     IF (noncolin) THEN
        DEALLOCATE(ps_nc)
        IF (lspinorb) DEALLOCATE(dpqq_so)
     ELSE
        DEALLOCATE(ps)
     ENDIF
  ENDIF
  !
  CALL deallocate_bec_type(becp2)
  !
  RETURN
END SUBROUTINE compute_ppsi
