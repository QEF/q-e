!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
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
  !  1) Initialization of ikks, ikqs, and nksq.
  !  2) EELS: Calculate phases associated with a q vector.
  !  3) Compute the number of occupied bands for each k point.
  !  4) Computes alpha_pv (needed by orthogonalize.f90 when lgauss=.true.)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, tau
  USE klist,                ONLY : nks,xk,ngk,igk_k
  USE wvfct,                ONLY : nbnd
  USE lr_variables,         ONLY : lr_verbosity, eels, size_evc, calculator, &
                                 & iund0psi, iudwf, iu1dwf,&
                                 & iundvpsi, magnons, iunTwfc, & 
                                 & restart
  USE io_global,            ONLY : stdout
  USE constants,            ONLY : tpi, eps8
  USE noncollin_module,     ONLY : npol, nspin_mag
  USE gvect,                ONLY : ngm, g
  USE cell_base,            ONLY : at, bg, omega, tpiba2
  USE lsda_mod,             ONLY : current_spin, nspin
  USE wvfct,                ONLY : npwx, wg
  USE gvecw,                ONLY : gcutw
  USE io_files,             ONLY : prefix, iunwfc, nwordwfc, wfc_dir
  USE gvecs,                ONLY : doublegrid
  USE fft_base,             ONLY : dfftp 
  USE uspp,                 ONLY : vkb, okvan, nkb
  USE wavefunctions,        ONLY : evc
  USE eqv,                  ONLY : evq
  USE becmod,               ONLY : calbec, allocate_bec_type
  USE lrus,                 ONLY : becp1
  USE control_lr,           ONLY : alpha_pv, alpha_mix, tr2_ph, nbnd_occ, &
                                   nbnd_occx
  USE qpoint,               ONLY : xq, ikks, ikqs, nksq, eigqts, npwq
  USE eqv,                  ONLY : evq
  USE buffers,              ONLY : open_buffer, get_buffer, save_buffer
  USE control_flags,        ONLY : io_level, use_gpu
  USE lr_magnons_routines,  ONLY : T_rev
  USE gvecw,                ONLY : ecutwfc
  USE mp_global,            ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_max
  use mp_world, only : mpime
  USE uspp_init,             ONLY : init_us_2
#if defined (__CUDA)
  USE cudafor
#endif
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  REAL(kind=DP) :: arg
  INTEGER :: na,  & ! dummy index which runs over atoms 
             ik,  & ! dummy index which runs over k points
             ikk, & ! index of the point k
             npw, & ! number of the plane-waves at point k
             it
  LOGICAL :: exst, exst_mem 
  !
  ! Magnons local variables 
  INTEGER :: ikq, imk, imkq, ibnd
  COMPLEX(DP), allocatable :: Tevc(:,:)
  INTEGER :: npw_mk, npw_mkq, istat
  !
  ! logical variable to check file exists
  ! logical variable to check file exists in memory
  !
  ! 1) Initialization of ikks, ikqs, and nksq.
  !    Inspired by PH/initialize_ph.f90.
  !
  IF (eels) THEN
     !
     ! EELS (q!=0)
     !
     ! nksq is the number of k-points, NOT including k+q points
     !
     nksq = nks / 2
     !
     ALLOCATE(ikks(nksq), ikqs(nksq))
     DO ik = 1, nksq
        ikks(ik) = 2 * ik - 1
        ikqs(ik) = 2 * ik
     ENDDO
     !
  ELSEIF (magnons) THEN
     !
     ! Magnons 
     !
     nksq = nks / 3
     !
     ALLOCATE(ikks(nksq), ikqs(nksq))
     DO ik = 1, nksq
        ! position of k
        ikks(ik) = 3 * ik - 2
        ! position of k+q
        ikqs(ik) = 3 * ik -1
     ENDDO
     !    
  ELSE
     !
     ! Optical case (q=0)
     !
     nksq = nks
     !
     ! Note: in this case the arrays ikks and ikqs 
     ! are needed only in h_prec and ch_psi_all.
     !
     ALLOCATE(ikks(nksq), ikqs(nksq))
     DO ik = 1, nksq
        ikks(ik) = ik
        ikqs(ik) = ik
     ENDDO
     !
  ENDIF
  !
  ! 2) EELS-specific operations
  !
  IF (eels) THEN
     !
     size_evc = nbnd * npwx * npol * nksq
     nwordwfc = nbnd * npwx * npol
     IF (trim(calculator)=='sternheimer') THEN
        CALL open_buffer ( iundvpsi, 'dvpsi.', nwordwfc, io_level, exst_mem, exst)
        CALL open_buffer ( iudwf, 'dwf', nwordwfc, io_level, exst_mem, exst)
        CALL open_buffer ( iu1dwf, 'mwf', nwordwfc, io_level, exst_mem, exst)
        !
     ENDIF
     !
     CALL open_buffer (iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst)
     ! 
     IF (.NOT.exst .AND. .NOT.exst_mem) THEN
        CALL errore ('lr_init_nfo', 'file '//trim(prefix)//'.wfc not found', 1)
     ENDIF
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
           npw = ngk(ikk)
           !
           ! Read the wavefunction evc
           CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
           !
           ! Calculate beta-functions vkb at point k
           CALL init_us_2(npw, igk_k(1,ikk), xk(1,ikk), vkb, .true.)
           !$acc update host(vkb)
           !
           ! Calculate becp1=<vkb|evc>
           CALL calbec (npw, vkb, evc, becp1(ik))
           !
        ENDDO
        !
     ENDIF
     !
  ENDIF
  !
  ! 2) MAGNONS-specific operations
  !
  ! The length of the arrays d0psi, evc1 etc.
  !
  IF (magnons) THEN
     !
     nwordwfc = nbnd * npwx * npol
     !
     ! additional 2 for the anti-resonant part
     !
     ! Move at teh end because we need nbnd_occx
     !
     ! Open the file to read the wavefunctions at k and k+Q and k-Q
     ! after the nscf calculation.
     !
     size_evc = nbnd * npwx * npol * nksq
     !
     CALL open_buffer (iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst)
     ! 
     IF (.NOT.exst .AND. .NOT.exst_mem) THEN
        CALL errore ('lr_init_nfo', 'file '//trim(prefix)//'.wfc not found', 1)
     ENDIF
     !
     ! This used to be in lr_alloc_init, moved here because I need these variables now
     IF (allocated(evc)) THEN
#if defined(__CUDA)
        !$acc exit data delete(evc)
        IF(use_gpu) istat = cudaHostUnregister(C_LOC(evc(1,1)))
#endif
        DEALLOCATE(evc)
        ALLOCATE(evc(npwx*npol,nbnd))
#if defined(__CUDA)
        IF(use_gpu) istat = cudaHostRegister(C_LOC(evc(1,1)), sizeof(evc), cudaHostRegisterMapped)
        !$acc enter data create(evc)
#endif
        evc(:,:) = (0.0d0, 0.0d0)
        !$acc update device(evc)
     ENDIF 
     !
     ALLOCATE(evq(npwx*npol,nbnd))
     evq(:,:) = (0.0d0, 0.0d0)     
     !$acc enter data copyin(evq)
     !
     ! Write T-rev wfc to disk and keep unit open, they are going to be used
     ! throughout the whole calculation
     ! 
     ! If not a restart, create T-rev wfc and write them to disk
     !
     IF ( .not. restart ) THEN
        !
        allocate( Tevc(npwx*npol,nbnd) )
        Tevc(:,:)   = (0.0d0, 0.0d0)
        !
        CALL open_buffer (iunTwfc, 'Twfc', nwordwfc, io_level, exst_mem, exst)
        !
        DO ik = 1, nksq
           !
           ikk = ikks(ik)        ! position of +k
           ikq = ikqs(ik)        ! position of +k+Q
           !
           if ( mod(ik,2) == 0) then
              imk  = ikk -3      ! position of -k
              imkq = ikk -1      ! position of -k-Q
           else
              imk  = ikk +3      ! position of -k
              imkq = ikk +5      ! position of -k-Q
           endif
           ! Order G vectors   
           npw     = ngk(ikk)
           npwq    = ngk(ikq)
           npw_mk  = ngk(imk)
           npw_mkq = ngk(imkq)
           !
           ! -k wfc
           ! 
           ! read u(-k)
           CALL get_buffer (evc, nwordwfc, iunwfc, imk)
           !
           ! form Tu(-k)
           do ibnd = 1, nbnd
              call T_rev( evc(:,ibnd), npw_mk, igk_k(:,imk), npw, igk_k(:,ikk), Tevc(:,ibnd), .true. )
           enddo
           !
           ! write it to disk
           CALL save_buffer (Tevc, nwordwfc, iunTwfc, 2*ik-1)
           !
           ! -k-q wfc
           !
           ! read u(-k-q)
           CALL get_buffer (evc, nwordwfc, iunwfc, imkq)
           !
           ! form Tu(-k-q)
           do ibnd = 1, nbnd
              call T_rev( evc(:,ibnd), npw_mkq, igk_k(:,imkq), npwq, igk_k(:,ikq), Tevc(:,ibnd), .true. )
           enddo
           !
           ! write it to disk
           CALL save_buffer (Tevc, nwordwfc, iunTwfc, 2*ik)
           !
         ENDDO 
         !
         CLOSE( unit = iunTwfc)
         !
      ENDIF
      ! Keep dir open throughout the calculation
      CALL open_buffer (iunTwfc, 'Twfc', nwordwfc, io_level, exst_mem, exst)
      ! 
      IF (.NOT.exst .AND. .NOT.exst_mem) THEN
         CALL errore ('lr_init_nfo', 'file '//trim(prefix)//'.Twfc not found', 1)
      ENDIF
      !
  ENDIF
  !
  ! 3) Compute the number of occupied bands for each k point
  !
  CALL setup_nbnd_occ()
  !
  if ( allocated(Tevc) ) deallocate( Tevc )
  !
  !IF (magnons) THEN
  ! do ik = 1, nks
  !   write(stdout,*) mpime, ik, nbnd_occ(ik), nbnd_occx
  ! enddo
  !ENDIF
  !
  ! 4) Compute alpha_pv
  !
  IF (eels .OR. magnons) THEN
     !
     alpha_pv = 0.0d0
     !
     IF (trim(calculator)=='sternheimer') THEN
        !
        ! Setup all gradient correction stuff
        !
        call setup_dgc()
        !
        ! 4) Computes the inverse of each matrix of the crystal symmetry group
        !
        CALL setup_alpha_pv()
        ! 
        DO it = 2, 100
           IF (alpha_mix (it) .eq. 0.0d0) alpha_mix(it) = alpha_mix (it - 1)
        ENDDO
        tr2_ph = 1.D-12
     ENDIF
     !
  ELSE
     !
     CALL setup_alpha_pv()
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE lr_init_nfo
