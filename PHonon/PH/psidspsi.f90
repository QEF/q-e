!
! Copyright (C) 2003-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine psidspsi (ik, uact, pdsp)
  !----------========----------------------------------------------------
  !! This routine calculates \(\langle \psi_v'|ds/du|\psi_v\rangle\)
  !! at \(q=0\). The displacements are described by the vector \(\text{uact}\).
  !! The result is stored in \(\text{pdsp}\). The routine is called for each
  !! k-point and for each pattern u. It computes simultaneously all the bands.
  !
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : tpiba
  USE gvect,            ONLY : g
  USE klist,            ONLY : xk, ngk, igk_k
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE lsda_mod,         ONLY : lsda, current_spin, isk
  USE noncollin_module, ONLY : noncolin, npol, lspinorb
  USE wavefunctions,    ONLY : evc
  USE wvfct,            ONLY : nbnd, npwx
  USE uspp,             ONLY : nkb, vkb, qq_nt, qq_so, ofsbeta
  USE uspp_param,       ONLY : nh, nhm
  USE phus,             ONLY : alphap
  USE lrus,             ONLY : becp1
  USE control_lr,       ONLY : lgamma
  USE qpoint,           ONLY : ikks

  implicit none
  !
  integer, intent(in) :: ik
  !! input: the k point
  complex(DP) :: uact (3*nat)
  !! input: the pattern of displacements
  complex(DP) :: pdsp(nbnd,nbnd)
  !! output: \(\langle \psi_v'|ds/du|\psi_v\rangle\)
  !
  ! ... local variables
  !
  integer :: na, mu, ikk, ikq, ig, igg, nt, ibnd, ijkb0, &
       ikb, jkb, ih, jh, ipol, is, js, ijs, npw

  real(DP), parameter :: eps = 1.d-12

  complex(DP), ALLOCATABLE :: ps1 (:,:), ps2 (:,:,:), aux (:), dspsi(:,:)
  complex(DP), ALLOCATABLE :: ps1_nc(:,:,:), ps2_nc(:,:,:,:)
  ! temporary arrays for optimization
  complex(DP) :: temp_ps1, temp_ps2
  complex(DP) :: temp_ps1_nc(npol), temp_ps2_nc(npol)
  complex(DP) :: ps2_col(nbnd)
  ! small buffers for loop optimization (allocated once with max size)
  complex(DP) :: alphap_buf_nc(nhm, npol), becp1_buf_nc(nhm, npol)
  complex(DP) :: alphap_buf_k(nhm), becp1_buf_k(nhm)
  complex(DP) :: qq_so_buf(nhm, npol*npol), qq_nt_buf(nhm)

  logical :: ok

  if (noncolin) then
     allocate (ps1_nc ( nkb, npol, nbnd ))
     allocate (ps2_nc ( nkb, npol, nbnd, 3))
  else
     allocate (ps1 ( nkb, nbnd ))
     allocate (ps2 ( nkb, nbnd, 3))
  endif
  allocate (dspsi (npwx*npol, nbnd))
  allocate (aux ( npwx ))

  if (lgamma) then
     ikk = ikks(ik)
     ikq = ikk
     npw = ngk(ikk)
  else
     call infomsg ('psidspsi', 'called for lgamma .eq. false')
  endif
  if (lsda) current_spin = isk (ikk)

  if (noncolin) then
     ps1_nc = (0.d0, 0.d0)
     ps2_nc = (0.d0, 0.d0)
  else
     ps1(:,:)   = (0.d0, 0.d0)
     ps2(:,:,:) = (0.d0, 0.d0)
  endif
  pdsp(:,:) = (0.d0, 0.d0)
  dspsi = (0.d0, 0.d0)
  !
  DO na = 1, nat
     nt = ityp(na)
     ijkb0 = ofsbeta(na)
     mu = 3 * (na - 1)
     IF ( ABS(uact(mu+1)) + ABS(uact(mu+2)) + ABS(uact(mu+3)) > eps ) THEN
        DO ih = 1, nh(nt)
           ikb = ijkb0 + ih
           DO ipol = 1, 3
              DO ibnd = 1, nbnd
                 ! Load small buffers for this (ih, ipol, ibnd)
                 qq_nt_buf(1:nh(nt)) = qq_nt(ih, 1:nh(nt), nt)
                 IF (noncolin) THEN
                    DO jh = 1, nh(nt)
                       jkb = ijkb0 + jh
                       alphap_buf_nc(jh, 1:npol) = alphap(ipol, ik)%nc(jkb, 1:npol, ibnd)
                       becp1_buf_nc(jh, 1:npol)  = becp1(ik)%nc(jkb, 1:npol, ibnd)
                    END DO
                    IF (lspinorb) THEN
                       DO jh = 1, nh(nt)
                          qq_so_buf(jh, 1:npol*npol) = qq_so(ih, jh, 1:npol*npol, nt)
                       END DO
                    END IF
                 ELSE
                    DO jh = 1, nh(nt)
                       jkb = ijkb0 + jh
                       alphap_buf_k(jh) = alphap(ipol, ik)%k(jkb, ibnd)
                       becp1_buf_k(jh)  = becp1(ik)%k(jkb, ibnd)
                    END DO
                 END IF
                 ! Initialize temp accumulators
                 IF (noncolin) THEN
                    temp_ps1_nc = (0.d0, 0.d0)
                    temp_ps2_nc = (0.d0, 0.d0)
                 ELSE
                    temp_ps1 = (0.d0, 0.d0)
                    temp_ps2 = (0.d0, 0.d0)
                 END IF
                 ! Accumulate over jh
                 DO jh = 1, nh(nt)
                    IF (noncolin) THEN
                       IF (lspinorb) THEN
                          ijs = 0
                          DO is = 1, npol
                             DO js = 1, npol
                                ijs = ijs + 1
                                temp_ps1_nc(is) = temp_ps1_nc(is) + &
                                   qq_so_buf(jh, ijs) * alphap_buf_nc(jh, js) * uact(mu + ipol)
                                temp_ps2_nc(is) = temp_ps2_nc(is) + &
                                   qq_so_buf(jh, ijs) * becp1_buf_nc(jh, js) * uact(mu + ipol)
                             END DO
                          END DO
                       ELSE
                          DO is = 1, npol
                             temp_ps1_nc(is) = temp_ps1_nc(is) + &
                                qq_nt_buf(jh) * alphap_buf_nc(jh, is) * uact(mu + ipol)
                             temp_ps2_nc(is) = temp_ps2_nc(is) + &
                                qq_nt_buf(jh) * becp1_buf_nc(jh, is) * uact(mu + ipol)
                          END DO
                       END IF
                    ELSE
                       temp_ps1 = temp_ps1 + qq_nt_buf(jh) * alphap_buf_k(jh) * uact(mu + ipol)
                       temp_ps2 = temp_ps2 + qq_nt_buf(jh) * becp1_buf_k(jh) * uact(mu + ipol)
                    END IF
                 END DO
                 ! ps1(ikb,ibnd)      = sum_{jh,ipol} q_{ih,jh} <d(beta_{jh})/d(tau)*u | psi_ibnd>
                 ! ps2(ikb,ibnd,ipol) = sum_{jh}      q_{ih,jh} <beta_{jh} | psi_ibnd> * uact(ipol)
                 IF (noncolin) THEN
                    DO is = 1, npol
                       ps1_nc(ikb, is, ibnd) = ps1_nc(ikb, is, ibnd) + temp_ps1_nc(is)
                       ps2_nc(ikb, is, ibnd, ipol) = temp_ps2_nc(is)
                    END DO
                 ELSE
                    ps1(ikb, ibnd) = ps1(ikb, ibnd) + temp_ps1
                    ps2(ikb, ibnd, ipol) = temp_ps2
                 END IF
              END DO
           END DO
        END DO
     END IF
  END DO
  !
  !      This term is proportional to beta(k+q+G)
  !
  if (nkb > 0) then
     if (noncolin) then
        call zgemm ('N', 'N', npw, nbnd*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1_nc, nkb, (1.d0, 0.d0) , dspsi, npwx)
     else
        call zgemm ('N', 'N', npw, nbnd, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1, nkb, (1.d0, 0.d0) , dspsi, npwx)
     endif
  endif
  IF (noncolin) THEN
     DEALLOCATE (ps1_nc)
  ELSE
     DEALLOCATE (ps1)
  END IF
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  DO ikb = 1, nkb
     DO ipol = 1, 3
        IF (noncolin) THEN
           ok = ANY(ABS(ps2_nc(ikb, 1, 1:nbnd, ipol)) > eps) .OR. &
                ANY(ABS(ps2_nc(ikb, 2, 1:nbnd, ipol)) > eps)
        ELSE
           ok = ANY(ABS(ps2(ikb, 1:nbnd, ipol)) > eps)
        END IF
        IF (ok) THEN
           DO ig = 1, npw
              igg = igk_k(ig, ikk)
              aux(ig) = vkb(ig, ikb) * (0.d0, -1.d0) * tpiba * &
                        (xk(ipol, ik) + g(ipol, igg))
           END DO
           IF (noncolin) THEN
              ps2_col(1:nbnd) = ps2_nc(ikb, 1, 1:nbnd, ipol)
              CALL zgemm('N', 'T', npw, nbnd, 1, (1.d0,0.d0), aux, npwx, &
                         ps2_col, nbnd, (1.d0,0.d0), dspsi(1,1), npwx*npol)
              ps2_col(1:nbnd) = ps2_nc(ikb, 2, 1:nbnd, ipol)
              CALL zgemm('N', 'T', npw, nbnd, 1, (1.d0,0.d0), aux, npwx, &
                         ps2_col, nbnd, (1.d0,0.d0), dspsi(npwx+1,1), npwx*npol)
           ELSE
              ps2_col(1:nbnd) = ps2(ikb, 1:nbnd, ipol)
              CALL zgemm('N', 'T', npw, nbnd, 1, (1.d0,0.d0), aux, npwx, &
                         ps2_col, nbnd, (1.d0,0.d0), dspsi, npwx*npol)
           END IF
        END IF
     END DO
  END DO
  DEALLOCATE (aux)
  IF (noncolin) THEN
     DEALLOCATE (ps2_nc)
  ELSE
     DEALLOCATE (ps2)
  END IF
  !
  !      pdsp(ibnd,jbnd) = <evc(ibnd)|dspsi(jbnd)>
  !
  CALL zgemm('C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0), &
             evc, npwx*npol, dspsi, npwx*npol, (0.d0,0.d0), pdsp, nbnd)
  DEALLOCATE (dspsi)

  return
end subroutine psidspsi
