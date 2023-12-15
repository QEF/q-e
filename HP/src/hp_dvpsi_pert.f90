!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine hp_dvpsi_pert (ik, nrec)
  !----------------------------------------------------------------------
  !
  ! This routine performes an action of the perturbing DFT+U potential on to
  ! the unperturbed KS wavefunction.
  ! See Eq. (46) in Ref. [1]
  ! [1] Phys. Rev. B 98, 085127 (2018)
  !
  ! On output: dvpsi(ig, ibnd) = 
  !     \sum_m | S(k+q)*phi_(k+q,J,m)> * < S(k)*phi_(k,J,m)| psi(ibnd,k) >
  !
  ! dvpsi is for a given "k", "q" and "J"
  !
  ! (evc, swfcatomk, swfcatomkpq must be set)
  !
  USE kinds,                ONLY : DP
  USE io_files,             ONLY : nwordwfcU
  USE ions_base,            ONLY : nat, ityp
  USE wavefunctions,        ONLY : evc
  USE klist,                ONLY : ngk
  USE buffers,              ONLY : save_buffer, get_buffer
  USE wvfct,                ONLY : npwx, nbnd
  USE mp_pools,             ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
  USE eqv,                  ONLY : dvpsi
  USE qpoint,               ONLY : nksq, ikks, ikqs
  USE units_lr,             ONLY : iuatswfc
  USE control_lr,           ONLY : lgamma
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, offsetU, nwfcU
  USE ldaU_hp,              ONLY : nqsh, perturbed_atom, iudvwfc, lrdvwfc
  USE ldaU_lr,              ONLY : swfcatomk, swfcatomkpq
  USE noncollin_module,     ONLY : noncolin, npol, domag
  USE qpoint_aux,           ONLY : ikmks, ikmkmqs
  USE io_global,            ONLY : stdout
  USE mp_bands,             ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik, nrec
  ! ik:    given k point
  ! nrec:  record number for evc and dvpsi 
  !
  INTEGER :: ikk, ikq, na, nt, m, ihubst, ibnd, ig, counter, ldim
  INTEGER :: npw, npwq
  ! number of plane waves at k and k+q
  COMPLEX (DP), ALLOCATABLE :: proj(:,:) 
  !
  CALL start_clock ('hp_dvpsi_pert')
  !
  ! Check that there is one perturbed atom
  !
  counter = 0
  DO na = 1, nat
     IF (perturbed_atom(na)) counter = counter + 1
  ENDDO
  IF (counter /= 1) CALL errore( 'hp_dvpsi_pert', "One perturbed atom must be specified", 1)
  !
  dvpsi(:,:) = (0.0d0, 0.0d0)
  !
  ! Compute dvpsi for ik and write on buffer iudvwfc
  !
  ALLOCATE (proj(nwfcU,nbnd))
  proj(:,:) = (0.0d0, 0.0d0)
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  !
  ! Read atomic wfc's : S(k)*phi(k) and S(k+q)*phi(k+q)
  !  
  CALL get_buffer(swfcatomk, nwordwfcU, iuatswfc, ikk)
  IF (.NOT.lgamma)  CALL get_buffer(swfcatomkpq, nwordwfcU, iuatswfc, ikq)
  !
  IF ( nrec > nksq ) then
   !CALL apply_trev(swfcatomk, ikk, ikk, size(swfcatomk, 2))
   !CALL apply_trev(swfcatomkpq, ikq, ikq, size(swfcatomkpq, 2))
  ENDIF
  ! Calculate proj(ibnd, ihubst) = < S(k)*phi(k,I,m)| psi(ibnd,k) >
  !
  DO na = 1, nat
     IF ( perturbed_atom(na) ) THEN
        nt = ityp(na)
        ldim = (2 * Hubbard_l(nt) + 1) * npol 
        IF (noncolin) then
           ihubst = offsetU(na) + 1
           CALL ZGEMM('C','N', ldim, nbnd, npwx*npol, (1.d0,0.d0), &
                    swfcatomk(1,ihubst), npwx*npol, evc, npwx*npol,&
                    (0.d0,0.d0), proj(ihubst,1), nwfcU)
        ELSE
           DO m = 1, ldim
              ihubst = offsetU(na) + m   ! I m index
              DO ibnd = 1, nbnd
                 ! FIXME: use ZGEMM instead of dot_product
                 proj(ihubst, ibnd) = DOT_PRODUCT( swfcatomk(1:npwx*npol,ihubst ), evc(1:npwx*npol,ibnd) )
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ENDDO
  !
#if defined (__MPI)
  CALL mp_sum(proj, intra_pool_comm)  
#endif
  !
  DO na = 1, nat
     IF ( perturbed_atom(na) ) THEN
        nt = ityp(na)
        ldim = (2 * Hubbard_l(nt) + 1) * npol
        DO m = 1, ldim
           ihubst = offsetU(na) + m
           DO ibnd = 1, nbnd
              DO ig = 1, npwx
                 dvpsi(ig, ibnd) = dvpsi(ig, ibnd) + &
                        & swfcatomkpq(ig,ihubst) * proj(ihubst,ibnd)
                 IF (noncolin) &
                    dvpsi(ig+npwx, ibnd) = dvpsi(ig+npwx, ibnd) + &
                        & swfcatomkpq(ig+npwx,ihubst) * proj(ihubst,ibnd)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  ! Write dvpsi on file.
  !
  CALL save_buffer(dvpsi, lrdvwfc, iudvwfc, nrec)
  !
  DEALLOCATE (proj)
  !
  CALL stop_clock ('hp_dvpsi_pert')
  !
  RETURN
  !
end subroutine hp_dvpsi_pert
