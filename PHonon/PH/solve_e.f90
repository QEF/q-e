!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_e
  !-----------------------------------------------------------------------
  !! This routine is a driver for the solution of the linear system which
  !! defines the change of the wavefunction due to an electric field.
  !! It performs the following tasks:
  !! a) computes the bare potential term times \(|\psi\rangle \);
  !! b) adds to it the screening term \(\Delta V_\text{SCF}|psi\rangle\).
  !!    If \(\text{lda_plus_u}=\text{TRUE}\) compute also the SCF part
  !!    of the response Hubbard potential;
  !! c) applies \(P_c^+\) (orthogonalization to valence states);
  !! d) calls \(\texttt{cgsolve_all}\) to solve the linear system;
  !! e) computes \(\Delta \rho\), \(\Delta V_\text{SCF}|psi\rangle\) and
  !!    symmetrizes them;
  !! f) if \(\text{lda_plus_u}=\text{TRUE}\) compute also the response
  !!    occupation matrices dnsscf.
  !! Step b, c, d are done in \(\text{sternheimer_kernel}\).
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : ionode
  USE io_files,              ONLY : diropn
  USE mp,                    ONLY : mp_sum
  USE klist,                 ONLY : ltetra, lgauss, xk, ngk, igk_k
  USE gvecs,                 ONLY : doublegrid
  USE fft_base,              ONLY : dfftp, dffts
  USE lsda_mod,              ONLY : lsda, current_spin, isk
  USE check_stop,            ONLY : check_stop_now
  USE buffers,               ONLY : get_buffer
  USE wavefunctions,         ONLY : evc
  USE uspp,                  ONLY : vkb
  USE uspp_param,            ONLY : nhm
  USE noncollin_module,      ONLY : nspin_mag
  USE paw_variables,         ONLY : okpaw
  USE ldaU,                  ONLY : lda_plus_u
  USE units_ph,              ONLY : lrdrho, iudrho, lrebar, iuebar
  USE units_lr,              ONLY : iuwfc, lrwfc
  USE output,                ONLY : fildrho
  USE control_ph,            ONLY : ext_recover
  USE recover_mod,           ONLY : read_rec
  USE qpoint,                ONLY : nksq, ikks
  USE control_lr,            ONLY : lgamma, convt, rec_code_read, rec_code, where_rec
  USE uspp_init,             ONLY : init_us_2
  USE dfpt_kernels,          ONLY : dfpt_kernel
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst
  !!
  INTEGER :: ikk, npw, iter0, ipol, ik
  !! counters
  REAL(DP) :: dr2
  !! self-consistency error

  COMPLEX(DP), ALLOCATABLE, TARGET :: dvscfin (:,:,:)
  !! change of the scf potential (input)
  COMPLEX(DP), POINTER :: dvscfins (:,:,:)
  !! change of the scf potential (smooth)
  COMPLEX(DP), ALLOCATABLE :: drhos(:, :, :)
  !! change of the charge density (smooth part only, dffts)
  COMPLEX(DP), ALLOCATABLE :: drhop(:, :, :)
  !! change of the charge density (smooth and hard parts, dfftp)
  COMPLEX(DP), ALLOCATABLE :: dbecsum(:,:,:,:)
  !! the becsum with dpsi
  INTEGER :: nnr
  !
  call start_clock ('solve_e')
  !
  !  This routine is task group aware
  !
  allocate (dvscfin( dfftp%nnr, nspin_mag, 3))
  nnr = dfftp%nnr
  dvscfin=(0.0_DP,0.0_DP)
  if (doublegrid) then
     allocate (dvscfins(dffts%nnr, nspin_mag, 3))
     nnr = dffts%nnr
  else
     dvscfins => dvscfin
  endif
  !$acc enter data create(dvscfins(1:nnr, 1:nspin_mag, 1:3))
  ALLOCATE(drhos(dffts%nnr, nspin_mag, 3))
  ALLOCATE(drhop(dfftp%nnr, nspin_mag, 3))
  allocate (dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, 3))
  dbecsum = (0.d0, 0.d0)


  if (rec_code_read == -20.AND.ext_recover) then
     ! restarting in Electric field calculation
     IF (okpaw) THEN
        CALL read_rec(dr2, iter0, 3, dvscfin, dvscfins, drhop, dbecsum)
     ELSE
        CALL read_rec(dr2, iter0, 3, dvscfin, dvscfins)
     ENDIF
  else if (rec_code_read > -20 .AND. rec_code_read <= -10) then
     ! restarting in Raman: proceed
     convt = .true.
  else
     convt = .false.
     iter0 = 0
     dr2 = 0.d0
  endif
  !
  IF (rec_code_read > -20) convt=.TRUE.
  !
  if (convt) go to 155
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  if ( (lgauss .or. ltetra) .or..not.lgamma) call errore ('solve_e', &
       'called in the wrong case', 1)
  !
  ! Compute P_c^+ x psi for all polarization and k points and store in buffer
  !
  DO ik = 1, nksq
     DO ipol = 1, 3
        ikk = ikks(ik)
        npw = ngk(ikk)
        IF (lsda) current_spin = isk(ikk)
        !
        ! reads unperturbed wavefunctions psi_k in G_space, for all bands
        !
        IF (nksq > 1) THEN
           CALL get_buffer(evc, lrwfc, iuwfc, ikk)
           !$acc update device(evc)
        ENDIF
        !
        CALL init_us_2(npw, igk_k(1, ikk), xk(1, ikk), vkb, .true.)
        !$acc update host(vkb)
        !
        ! computes P_c^+ x psi_kpoint, written to buffer iuebar.
        !
        CALL dvpsi_e(ik, ipol)
        !
     ENDDO ! ipol
  ENDDO ! ik
  !
  ! Set records for restart
  !
  rec_code = -20  ! Electric field
  where_rec = 'solve_e...'
  !
  !   Solve DFPT fixed-point equation
  !
  CALL dfpt_kernel('PHONON', 3, iter0, lrebar, iuebar, dr2, drhos, drhop, dvscfins, dvscfin, dbecsum, 1, 0, 'efield')
  !
  IF (lda_plus_u) CALL dnsq_store(3, 0)
  !
  IF ( fildrho /= ' ') THEN
     IF ( ionode ) THEN
        INQUIRE (UNIT = iudrho, OPENED = exst)
        IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
        CALL diropn (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
     END IF
     DO ipol=1,3
        CALL davcio_drho(dvscfin(1,1,ipol),lrdrho, iudrho,ipol,+1)
     END DO
  END IF
  !
155 continue
  !
  deallocate (dbecsum)
  deallocate (drhos)
  deallocate (drhop)
  !$acc exit data delete(dvscfins)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)

  call stop_clock ('solve_e')
  return
end subroutine solve_e
