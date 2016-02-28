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
  !  5) Computes alpha_pv (needed by orthogonalize.f90 when lgauss=.true.)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, tau
  USE klist,                ONLY : nks,xk,ngk,igk_k
  USE wvfct,                ONLY : nbnd, igk, npw, g2kin
  USE realus,               ONLY : real_space
  USE lr_variables,         ONLY : lr_verbosity, eels, nwordd0psi, &
                                   nwordrestart, restart, size_evc, tmp_dir_lr
  USE io_global,            ONLY : stdout
  USE constants,            ONLY : tpi, eps8
  USE noncollin_module,     ONLY : npol
  USE gvect,                ONLY : ngm, g
  USE cell_base,            ONLY : at, bg, omega
  USE lsda_mod,             ONLY : current_spin, nspin
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
  USE control_lr,           ONLY : alpha_pv
  USE qpoint,               ONLY : xq, npwq, igkq, ikks, ikqs, nksq, eigqts
  USE eqv,                  ONLY : evq
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  REAL(kind=DP) :: arg
  INTEGER       :: i, ik, ibnd, ipol, ikk, ikq, ios, isym, na
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
     IF (restart) wfc_dir = tmp_dir_lr
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
  ! 4) Compute the number of occupied bands for each k point
  !
  CALL setup_nbnd_occ()
  !
  ! 5) Compute alpha_pv
  !
  IF (eels) THEN
     !
     alpha_pv = 0.0d0
     !
  ELSE
     !
     CALL setup_alpha_pv()
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
