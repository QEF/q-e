!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------
SUBROUTINE lr_init_nfo()
  !-------------------------------------------------------------
  !
  !  This subroutine prepares several variables which are needed in the
  !  TDDFPT program:
  !  1) Optical case: initialization of igk_k and ngk. 
  !  2) Initialization of ikks, ikqs, and nksq.
  !  3) EELS: Calculate phases associated with a q vector.
  !  4) Compute the number of occupied bands for each k point.
  !  5) Computes alpha_pv (needed in PH/orthogonalize.f90 when lgauss=.true.)
  !
  ! Created by Osman Baris Malcioglu (2009)
  ! Modified by Iurii Timrov (2013)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, tau
  USE klist,                ONLY : nks,degauss,lgauss,ngauss,xk,wk,ngk,&
                                   igk_k,nelec, two_fermi_energies, nelup, neldw
  USE wvfct,                ONLY : nbnd, et, igk, npw, g2kin
  USE realus,               ONLY : real_space
  USE lr_variables,         ONLY : lr_verbosity, eels, nwordd0psi, &
                                  & nwordrestart, restart, size_evc
  USE io_global,            ONLY : stdout
  USE constants,            ONLY : pi, tpi, degspin, eps8
  USE noncollin_module,     ONLY : noncolin, npol
  USE mp,                   ONLY : mp_max, mp_min
  USE mp_global,            ONLY : inter_pool_comm
  USE gvect,                ONLY : ngm, g
  USE cell_base,            ONLY : at, bg, omega
  USE ener,                 ONLY : ef, ef_up, ef_dw
  USE ktetra,               ONLY : ltetra
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE control_ph,           ONLY : tmp_dir_phq
  USE wvfct,                ONLY : npwx, wg
  USE gvecw,                ONLY : gcutw
  USE io_files,             ONLY : iunigk, seqopn, tmp_dir, prefix, &
                                 & diropn, nwordwfc, wfc_dir
  USE gvecs,                ONLY : doublegrid
  USE units_ph,             ONLY : iuwfc, lrwfc
  USE fft_base,             ONLY : dfftp 
  USE uspp,                 ONLY : vkb, okvan, nkb
  USE wavefunctions_module, ONLY : evc
  USE becmod,               ONLY : calbec, allocate_bec_type

  USE lrus,                 ONLY : becp1
  USE control_lr,           ONLY : alpha_pv, nbnd_occ
  USE qpoint,               ONLY : xq, npwq, igkq, ikks, ikqs, nksq, eigqts
  USE eqv,                  ONLY : evq
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  REAL(kind=DP) :: small, emin, emax, xmax, fac, targ, arg
  INTEGER       :: i, ik, ibnd, ipol, ikk, ikq, ios, isym, na
  REAL(DP), ALLOCATABLE :: wg_up(:,:), wg_dw(:,:)
  LOGICAL       :: exst ! logical variable to check file existence
  !
  ! 1) Optical case: initialize igk_k and ngk
  !    Open shell related
  !
  IF (.NOT.eels) THEN
     !
     IF ( .not. allocated( igk_k ) )  ALLOCATE(igk_k(npwx,nks))
     IF ( .not. allocated( ngk ) )    ALLOCATE(ngk(nks))
     !
     IF (.not. real_space) THEN
        !
        DO ik = 1, nks
           !
           CALL gk_sort( xk(1,ik), ngm, g, gcutw, npw, igk, g2kin )
           !
           ngk(ik) = npw
           igk_k(:,ik) = igk(:)
           !
        ENDDO
        !
     ENDIF
     !
  ENDIF
  !
  ! 2) Initialization of ikks, ikqs, and nksq.
  !
  IF (eels) THEN
     !
     ! EELS
     !
     ! nksq is the number of k-points, NOT including k+q points
     ! The following block was copied from PH/initialize_ph.f90.
     !
     nksq = nks / 2
     !
     ALLOCATE(ikks(nksq), ikqs(nksq))
     !
     DO ik = 1, nksq
        ikks(ik) = 2 * ik - 1
        ikqs(ik) = 2 * ik
     ENDDO
     !
  ELSE
     !
     ! Optical case
     !
     nksq = nks
     !
  ENDIF
  !
  ! The length of the arrays d0psi, evc1 etc.
  !
  nwordd0psi   = 2 * nbnd * npwx * npol * nksq
  nwordrestart = 2 * nbnd * npwx * npol * nksq
  nwordwfc     =     nbnd * npwx * npol 
  !
  ! 3) EELS-specific operations
  !
  IF (eels) THEN
     !
     ! Open the file to read the wavefunctions at k and k+q 
     ! after the nscf calculation.
     !
     iuwfc = 21
     lrwfc = 2 * nbnd * npwx * npol
     IF (restart) wfc_dir = tmp_dir_phq
     CALL diropn (iuwfc, 'wfc', lrwfc, exst)
     IF (.NOT.exst) THEN
        CALL errore ('lr_init_nfo', 'file '//trim(prefix)//'.wfc not found', 1)
     ENDIF
     !
     size_evc = nksq * nbnd * npwx * npol 
     !
     ! If restart=.true. recalculate the small group of q.
     !
     IF (restart) CALL lr_smallgq (xq)
     !
     ! USPP-specific initializations
     !
     IF (okvan) THEN
        !
        ! Calculate phases associated with the q vector
        !
        ALLOCATE (eigqts(nat))
        !
        DO na = 1, nat
           arg = ( xq(1) * tau(1,na) + &
                   xq(2) * tau(2,na) + &
                   xq(3) * tau(3,na) ) * tpi
           eigqts(na) = cmplx(cos(arg), -sin(arg) ,kind=DP)
        ENDDO
        !
        ! Calculate becp1 = <vkb|evc>
        !
        ALLOCATE (becp1(nksq))
        !
        DO ik = 1, nksq
           !
           CALL allocate_bec_type (nkb,nbnd,becp1(ik))
           !
           ikk = ikks(ik)
           !
           ! Determination of npw and igk.
           CALL gk_sort( xk(1,ikk), ngm, g, gcutw, npw,  igk,  g2kin )
           !
           ! Read the wavefunction evc
           CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
           !
           ! Calculate beta-functions vkb at k point
           CALL init_us_2(npw, igk, xk(1,ikk), vkb)
           !
           ! Calculate becp1
           CALL calbec (npw, vkb, evc, becp1(ik))
           !
        ENDDO
        !
     ENDIF
     !
  ENDIF
  !
  ! 4) Compute the number of occupied bands for each k point (PH/phq_setup.f90)
  !
  !if (.not. allocated (nbnd_occ) allocate( nbnd_occ (nks) )
  !
  IF (lgauss) THEN
     !
     ! discard conduction bands such that w0gauss(x,n) < small
     !
     ! hint
     !   small = 1.0333492677046d-2  ! corresponds to 2 gaussian sigma
     !   small = 6.9626525973374d-5  ! corresponds to 3 gaussian sigma
     !   small = 6.3491173359333d-8  ! corresponds to 4 gaussian sigma
     !
     small = 6.9626525973374d-5
     !
     ! - appropriate limit for gaussian broadening (used for all ngauss)
     !
     xmax = sqrt ( - log (sqrt (pi) * small) )
     !
     ! - appropriate limit for Fermi-Dirac
     !
     IF (ngauss== - 99) THEN
        fac = 1.d0 / sqrt (small)
        xmax = 2.d0 * log (0.5d0 * (fac + sqrt (fac * fac - 4.d0) ) )
     ENDIF
     targ = ef + xmax * degauss
     DO ik = 1, nks
        DO ibnd = 1, nbnd
           IF (et (ibnd, ik) <targ) nbnd_occ (ik) = ibnd
        ENDDO
        IF (nbnd_occ (ik) ==nbnd) WRITE( stdout, '(5x,/,&
             &"Possibly too few bands at point ", i4,3f10.5)') &
             ik,  (xk (ipol, ik) , ipol = 1, 3)
     ENDDO
  ELSEIF (ltetra) THEN
     CALL errore('lr_init_nfo','tddfpt + tetrahedra not implemented', 1)
  ELSE
     IF (noncolin) THEN
        nbnd_occ = nint (nelec)
     ELSE
        IF ( two_fermi_energies ) THEN
           !
           ALLOCATE(wg_up(nbnd,nks))
           ALLOCATE(wg_dw(nbnd,nks))
           CALL iweights( nks, wk, nbnd, nelup, et, ef_up, wg_up, 1, isk )
           CALL iweights( nks, wk, nbnd, neldw, et, ef_dw, wg_dw, 2, isk )
           DO ik = 1, nks
              DO ibnd=1,nbnd
                 IF (isk(ik)==1) THEN
                    IF (wg_up(ibnd,ik) > 0.0_DP) nbnd_occ (ik) = nbnd_occ(ik)+1
                 ELSE
                    IF (wg_dw(ibnd,ik) > 0.0_DP) nbnd_occ (ik) = nbnd_occ(ik)+1
                 ENDIF
              ENDDO
           ENDDO
           !
           ! the following line to prevent NaN in Ef
           !
           ef = ( ef_up + ef_dw ) / 2.0_dp
           !
           DEALLOCATE(wg_up)
           DEALLOCATE(wg_dw)
        ELSE
           IF (lsda) call infomsg('lr_init_nfo', &
                                 'occupation numbers probably wrong')
           DO ik = 1, nks
              nbnd_occ (ik) = nint (nelec) / degspin
           ENDDO
        ENDIF
     ENDIF
  ENDIF
  !
  ! 5) Computes alpha_pv
  !
  IF (eels) THEN
     !
     alpha_pv = 0.0d0
     !
  ELSE
     !
     emin = et (1, 1)
     DO ik = 1, nks
        DO ibnd = 1, nbnd
           emin = min (emin, et (ibnd, ik) )
        ENDDO
     ENDDO
     !
#ifdef __MPI
  ! find the minimum across pools
  CALL mp_min( emin, inter_pool_comm )
#endif
     !
     IF (lgauss) THEN
        emax = targ
        alpha_pv = emax - emin
     ELSE
        emax = et (1, 1)
        DO ik = 1, nks
           DO ibnd = 1, nbnd_occ(ik)
              emax = max (emax, et (ibnd, ik) )
           ENDDO
        ENDDO
        !
#ifdef __MPI
     ! find the maximum across pools
     CALL mp_max( emax, inter_pool_comm )
#endif
        !
        alpha_pv = 2.d0 * (emax - emin)
     ENDIF
     ! avoid zero value for alpha_pv
     alpha_pv = max (alpha_pv, 1.0d-2)
     !
  ENDIF
  !
  ! Initialize npw, igk, npwq, igkq.
  ! Copied from PH/phq_init.f90
  !
  ! Open the file with npw, igk and npwq, igkq.
  ! The igk at a given k and k+q
  !
  !iunigk = 24
  !IF (nksq > 1) CALL seqopn (iunigk, 'igk', 'unformatted', exst)
  !
  !IF ( nksq > 1 ) REWIND( iunigk )
  !
  !allocate (igkq(npwx))
  !
  !DO ik = 1, nksq
     !
  !   ikk  = ikks(ik)
  !   ikq  = ikqs(ik)
     !
     ! ... g2kin is used here as work space
     !
  !   CALL gk_sort( xk(1,ikk), ngm, g, gcutw, npw, igk, g2kin )
     !
     ! ... if there is only one k-point evc, evq, npw, igk stay in memory
     !
  !   IF ( nksq > 1 ) WRITE( iunigk ) npw, igk
     !
  !   CALL gk_sort( xk(1,ikq), ngm, g, gcutw, npwq, igkq, g2kin )
     !
  !   IF ( nksq > 1 ) WRITE( iunigk ) npwq, igkq
     !
     !IF ( ABS( xq(1) - ( xk(1,ikq) - xk(1,ikk) ) ) > eps8 .OR. &
     !     ABS( xq(2) - ( xk(2,ikq) - xk(2,ikk) ) ) > eps8 .OR. &
     !     ABS( xq(3) - ( xk(3,ikq) - xk(3,ikk) ) ) > eps8 ) THEN
     !WRITE( stdout,'(/,5x,"k points #",i6," and ", &
     !       & i6,5x," total number ",i6)') ikk, ikq, nksq
     !WRITE( stdout, '(  5x,"Expected q ",3f10.7)')(xq(ipol), ipol=1,3)
     !WRITE( stdout, '(  5x,"Found      ",3f10.7)')((xk(ipol,ikq) &
     !                                   -xk(ipol,ikk)), ipol = 1, 3)
     !CALL errore( 'lr_init_nfo', 'wrong order of k points', 1 )
     !END IF
     !
  !ENDDO
  !
  RETURN
  !
END SUBROUTINE lr_init_nfo
