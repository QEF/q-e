!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine incdrhous_nc (drhoscf, weight, ik, dbecsum, evcr, wgg, becq, &
     alpq, mode)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density
  !     and of the magnetization due
  !     to the displacement of the augmentation charge. Only the
  !     smooth part is computed here.
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : ntyp => nsp, nat, ityp
  USE fft_base,  ONLY : dffts, dfftp
  USE fft_interfaces, ONLY: invfft
  USE gvecs,     ONLY : nls
  USE lsda_mod,  ONLY : nspin
  USE spin_orb,  ONLY : lspinorb, domag
  USE noncollin_module, ONLY : npol, nspin_mag
  USE uspp,      ONLY : nkb, qq, qq_so
  USE uspp_param,ONLY : nhm, nh
  USE wvfct,     ONLY : nbnd, npwx
  USE modes,     ONLY : u
  USE phus,      ONLY : alphap
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE becmod,    ONLY : bec_type
  USE klist,     ONLY : ngk, igk_k
  USE qpoint,    ONLY : nksq, ikqs, ikks
  USE eqv,       ONLY : dpsi, evq
  USE control_lr,ONLY : nbnd_occ
  USE lrus,      ONLY : becp1

  implicit none

  integer :: ik, mode
  ! input: the k point
  ! input: the mode which is computed
  ! input: the quantity to compute (1 charge, 2-4 magnetization)

  real(DP) :: weight, wgg (nbnd, nbnd, nksq)
  ! input: the weight of the k point
  ! input: the weights

  complex(DP) :: evcr (dffts%nnr, npol, nbnd), drhoscf(dfftp%nnr,nspin_mag), &
       dbecsum(nhm, nhm, nat, nspin)
  ! input: the wavefunctions at k in real
  ! output: the change of the charge densi
  ! inp/out: the accumulated dbec
  type (bec_type) :: becq(nksq), & ! nkb, nbnd)
                     alpq (3, nksq)
  ! input: the becp with psi_{k+q}
  ! input: the alphap with psi_{k+q}
  !
  !   here the local variable
  !
  real(DP) :: wgt
  ! the effective weight of the k point

  complex(DP), allocatable :: ps1 (:,:), dpsir (:,:)
  ! auxiliary space
  ! the change of wavefunctions in real sp

  integer :: ibnd, jbnd, nt, na, mu, ih, jh, ikb, jkb, ijkb0, &
       startb, lastb, ipol, ikk, ir, ig, ijs, is1, is2, ikq, npwq
  ! counters

  call start_clock ('incdrhous')
  allocate (dpsir(dffts%nnr,npol))
  allocate (ps1 (nbnd, nbnd))

  call divide (intra_bgrp_comm, nbnd, startb, lastb)
  ps1 (:,:) = (0.d0, 0.d0)
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npwq= ngk(ikq)
  !
  !   Here we prepare the two terms
  !
  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) == nt) then
           mu = 3 * (na - 1)
           if (abs(u(mu+1,mode)) + abs(u(mu+2,mode)) &
                                 + abs(u(mu+3,mode)) > 1.0d-12) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ibnd = 1, nbnd
                       do jbnd = startb, lastb
                          do ipol = 1, 3
                             mu = 3 * (na - 1) + ipol
                             IF (lspinorb) THEN
                                ijs=0
                                DO is1=1,npol
                                   DO is2=1,npol
                                      ijs=ijs+1
                                      ps1(ibnd,jbnd)=ps1(ibnd,jbnd)-       &
                                          qq_so(ih,jh,ijs,nt) *            &
                                       (alphap(ipol,ik)%nc(jkb,is2,ibnd)*&
                                      CONJG(becq(ik)%nc(ikb,is1,jbnd)) +   &
                                            becp1(ik)%nc(jkb,is2,ibnd) *   &
                                      CONJG(alpq(ipol,ik)%nc(ikb,is1,jbnd)) )* &
                                      wgg (ibnd, jbnd, ik) * u (mu, mode)
                                   END DO
                                END DO
                             ELSE
                                ps1(ibnd,jbnd)=ps1(ibnd,jbnd)-qq(ih,jh,nt)*&
                                (alphap(ipol,ik)%nc(ikb,1,ibnd) *&
                                      CONJG(becq(ik)%nc(jkb,1,jbnd)) +     &
                                            becp1(ik)%nc(ikb,1,ibnd) *     &
                                      CONJG(alpq(ipol,ik)%nc(jkb,1,jbnd)) +    &
                                 alphap(ipol,ik)%nc(ikb,2,ibnd) *&
                                      CONJG(becq(ik)%nc(jkb,2,jbnd)) +     &
                                            becp1(ik)%nc(ikb,2,ibnd) *     &
                                      CONJG(alpq(ipol,ik)%nc(jkb,2,jbnd)) )*  &
                                      wgg (ibnd, jbnd, ik) * u (mu, mode)
                             END IF
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
           endif
           ijkb0 = ijkb0 + nh (nt)
        endif
     enddo
  enddo
#if defined(__MPI)
  call mp_sum (ps1, intra_bgrp_comm)
#endif
  dpsi (:,:) = (0.d0, 0.d0)
  wgt = 2.d0 * weight / omega
  do ibnd = 1, nbnd_occ (ikk)
     do jbnd = 1, nbnd
        call zaxpy (npwx*npol,ps1(ibnd,jbnd),evq(1,jbnd),1,dpsi(1,ibnd), 1)
     enddo
     dpsir = (0.d0, 0.d0)
     do ig = 1, npwq
        dpsir(nls(igk_k(ig,ikq)),1) = dpsi (ig, ibnd)
        dpsir(nls(igk_k(ig,ikq)),2) = dpsi (ig+npwx, ibnd)
     enddo
     CALL invfft ('Wave', dpsir(:,1), dffts)
     CALL invfft ('Wave', dpsir(:,2), dffts)
     do ir = 1, dffts%nnr
        drhoscf(ir,1)=drhoscf(ir,1)+wgt* &
                                   (dpsir(ir,1)*CONJG(evcr(ir,1,ibnd))+ &
                                    dpsir(ir,2)*CONJG(evcr(ir,2,ibnd)) )
        IF (domag) THEN
           drhoscf(ir,2)=drhoscf(ir,2)+ &
                                wgt*(dpsir(ir,1)*CONJG(evcr(ir,2,ibnd))+ &
                                     dpsir(ir,2)*CONJG(evcr(ir,1,ibnd)))
           drhoscf(ir,3)=drhoscf(ir,3)+ &
                       wgt*(dpsir(ir,2)*CONJG(evcr(ir,1,ibnd)) - &
                       dpsir(ir,1)*CONJG(evcr(ir,2,ibnd)) ) *(0.d0,-1.d0)
           drhoscf(ir,4)=drhoscf(ir,4)+wgt* &
                                    (dpsir(ir,1)*CONJG(evcr(ir,1,ibnd)) - &
                                     dpsir(ir,2)*CONJG(evcr(ir,2,ibnd)) )
        END IF
     enddo
  enddo

  call addusdbec_nc (ik, weight, dpsi, dbecsum)
  deallocate (ps1)
  deallocate (dpsir)

  call stop_clock ('incdrhous')
  return
end subroutine incdrhous_nc
