!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine hp_dvpsi_pert (ik)
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
  ! dvpsi is READ from file if this_pert_is_on_file(ik) = .TRUE.
  ! otherwise dvpsi is COMPUTED and WRITTEN on file 
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
  USE ldaU_hp,              ONLY : nqsh, perturbed_atom, this_pert_is_on_file, &
                                   iudvwfc, lrdvwfc, swfcatomk, swfcatomkpq
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  ! given k point
  !
  INTEGER :: ikk, ikq, na, nt, m, ihubst, ibnd, ig, counter
  INTEGER :: npw, npwq
  ! number of plane waves at k and k+q
  COMPLEX (DP), ALLOCATABLE :: proj(:,:) 
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  CALL start_clock ('hp_dvpsi_pert')
  !
  ! Check that there is one perturbed atom
  !
  counter = 0
  DO na = 1, nat
     IF (perturbed_atom(na)) counter = counter + 1
  ENDDO
  IF (counter.NE.1) CALL errore( 'hp_dvpsi_pert', "One perturbed atom must be specified", 1) 
  !
  dvpsi(:,:) = (0.0d0, 0.0d0)
  !
  ! If this is not the first iteration, hence dvpsi was already
  ! computed before. So read it from file and exit.
  !
  IF (this_pert_is_on_file(ik)) THEN
     !
     CALL get_buffer(dvpsi, lrdvwfc, iudvwfc, ik)
     CALL stop_clock ('hp_dvpsi_pert')
     RETURN  
     !  
  ENDIF
  !
  ! If this is a first iteration, then dvpsi must be computed
  ! and written on file.
  !
  ALLOCATE (proj(nbnd,nwfcU))
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
  ! Calculate proj(ibnd, ihubst) = < S(k)*phi(k,I,m)| psi(ibnd,k) >
  !
  DO na = 1, nat
     nt = ityp(na)
     IF ( perturbed_atom(na) ) THEN
        DO m = 1, 2 * Hubbard_l(nt) + 1
           ihubst = offsetU(na) + m   ! I m index
           DO ibnd = 1, nbnd
              proj(ibnd, ihubst) = ZDOTC(npw, swfcatomk(:,ihubst ), 1, evc(:,ibnd), 1)
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
#if defined (__MPI)
  CALL mp_sum(proj, intra_pool_comm)
#endif
  !
  DO na = 1, nat
     nt = ityp(na)
     IF ( perturbed_atom(na) ) THEN
        DO m = 1, 2 * Hubbard_l(nt) + 1
           ihubst = offsetU(na) + m
           DO ibnd = 1, nbnd
              DO ig = 1, npwq
                 dvpsi(ig, ibnd) = dvpsi(ig, ibnd) + &
                        & swfcatomkpq(ig,ihubst) * proj(ibnd,ihubst)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  ! Write dvpsi on file.
  !
  CALL save_buffer(dvpsi, lrdvwfc, iudvwfc, ik)
  this_pert_is_on_file(ik) = .true.
  !
  DEALLOCATE (proj)
  !
  CALL stop_clock ('hp_dvpsi_pert')
  !
  RETURN
  !
end subroutine hp_dvpsi_pert
