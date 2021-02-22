!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
!-----------------------------------------------------------------------
subroutine kc_setup
  !-----------------------------------------------------------------------
  !
  !!  This subroutine prepares several variables which are needed in the
  !!  KC calculation:
  !!  * computes the total local potential (external+scf) on the smooth
  !!    grid to be used in h_psi and similia
  !!  * set the inverse of every matrix invs
  !!  * for metals sets the occupied bands
  !!  * Open buffer for the KS and eventualy the minimizing WFCs
  !!  * Open buffer for the KS states in the WANNIER gauge
  !!  * Read the U matrix from Wannier
  !!  * Rotate the KS state to the localized gauge
  !!  * Compute the periodic part of the wannier orbital and store in the buffer 
  !!    a separate directory is generated for each q 
  !
  !
  USE kinds,             ONLY : DP
  USE ions_base,         ONLY : nat, ntyp => nsp, ityp
  USE io_files,          ONLY : tmp_dir
  USE lsda_mod,          ONLY : nspin, starting_magnetization, lsda, isk
  USE scf,               ONLY : v, vrs, vltot,  kedtau
  USE fft_base,          ONLY : dfftp, dffts
  USE gvect,             ONLY : ngm
  USE gvecs,             ONLY : doublegrid
  USE uspp_param,        ONLY : upf
  USE spin_orb,          ONLY : domag
  USE noncollin_module,  ONLY : noncolin, m_loc, angle1, angle2, ux!, nspin_mag, npol
  USE wvfct,             ONLY : nbnd
  USE nlcc_ph,           ONLY : drc
  USE uspp,              ONLY : nlcc_any
  USE funct,             ONLY : dft_is_gradient
  !
  USE units_lr,          ONLY : iuwfc
  USE wvfct,             ONLY : npwx
  USE control_flags,     ONLY : io_level
  USE io_files,          ONLY : prefix
  USE buffers,           ONLY : open_buffer, save_buffer, close_buffer, get_buffer
  USE control_kc_wann,   ONLY : evc0, iuwfc_wann, iuwfc_wann_allk, kc_iverbosity, &
                                spin_component, isq, read_unitary_matrix, & 
                                num_wann, num_wann_occ, occ_mat !, wq, nqstot
  USE io_global,         ONLY : stdout
  USE klist,             ONLY : nkstot, xk
  USE cell_base,         ONLY : at !, bg
  USE fft_base,          ONLY : dffts
  USE disp,              ONLY : x_q, lgamma_iq
  !
  USE control_lr,       ONLY : nbnd_occ, lgamma
  USE control_ph,       ONLY : tmp_dir_ph, tmp_dir_phq
  USE scf,              ONLY : rho
  USE save_ph,          ONLY : tmp_dir_save
  USE io_global,        ONLY : ionode, ionode_id
  USE mp_images,        ONLY : intra_image_comm
  USE mp,               ONLY : mp_bcast
  USE io_files,         ONLY : create_directory
  USE io_rho_xml,       ONLY : write_scf
  !
  USE mp_bands,         ONLY : inter_bgrp_comm
  USE io_kcwann,    ONLY : write_rhowann
  !
  USE mp,               ONLY : mp_sum
  !
  implicit none

  integer :: na, i, ik
  ! counters
  !
  INTEGER   :: lrwfc, iun_qlist
  LOGICAL   :: exst, exst_mem
  INTEGER :: iq, nqs
  REAL(DP) :: xq(3)
  COMPLEX(DP), ALLOCATABLE :: rhowann(:,:), rhowann_aux(:)
  CHARACTER (LEN=256) :: filename, file_base
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  !
  INTEGER :: &
       igk_k_all(npwx,nkstot),&    ! index of G corresponding to a given index of k+G
       ngk_all(nkstot)             ! number of plane waves for each k point
  !
  call start_clock ('kc_setup')
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  !
  ! 2) Set non linear core correction stuff
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  if (nlcc_any) allocate (drc( ngm, ntyp))
  !
  !  3) If necessary calculate the local magnetization. This information is
  !      needed in find_sym
  !
  IF (.not.ALLOCATED(m_loc)) ALLOCATE( m_loc( 3, nat ) )
  IF (noncolin.and.domag) THEN
     DO na = 1, nat
        !
        m_loc(1,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
        m_loc(2,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
        m_loc(3,na) = starting_magnetization(ityp(na)) * &
                      COS( angle1(ityp(na)) )
     END DO
     ux=0.0_DP
     if (dft_is_gradient()) call compute_ux(m_loc,ux,nat)
  ENDIF
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  call setup_nbnd_occ ( ) 
  !
  ! Open buffers for the KS and eventualy the minimizing WFCs 
  ! 
  iuwfc = 20
  lrwfc = nbnd * npwx
  io_level = 1
  CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
  IF (.NOT.exst.AND..NOT.exst_mem) THEN
     CALL errore ('kc_setup', 'file '//trim(prefix)//'.wfc not found', 1)
  END IF
  if (kc_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for KS wfcs, OPENED")')
  !
  ! 8) READ the U matrix from Wannier and set the toal number of WFs
  IF (read_unitary_matrix) THEN 
    CALL read_wannier ( )
    if (kc_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Unitary matrix, READ from file")')
  ELSE
    num_wann = nbnd
    num_wann_occ = nbnd_occ(1) ! Assuming insulating
  ENDIF
  !
  ! ... Open a file to store the KS states in the WANNIER gauge
  iuwfc_wann = 21
  io_level = 1
  lrwfc = num_wann * npwx 
  CALL open_buffer ( iuwfc_wann, 'wfc_wann', lrwfc, io_level, exst )
  if (kc_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for WFs, OPENED")')
  !
  ! ... Open an other buffer for the KS states in the WANNIER gauge which contains
  ! ... all the k points (not just the one in this pool). This is needed for each k-poit 
  ! ... to have access to all the other k-points. MEMORY INTENSE
  iuwfc_wann_allk = 210
  io_level = 1
  lrwfc = num_wann * npwx
  CALL open_buffer ( iuwfc_wann_allk, 'wfc_wann_allk', lrwfc, io_level, exst )
  if (kc_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for WFs ALL-k, OPENED")')
  !
  !
  ALLOCATE (rhowann ( dffts%nnr, num_wann), rhowann_aux(dffts%nnr) )
  ALLOCATE ( evc0(npwx, num_wann) )
  ALLOCATE ( occ_mat (num_wann, num_wann, nkstot) )
  occ_mat = 0.D0
  ! 
  ! ... Rotate the KS state to the localized gauge nd save on a buffer
  CALL rotate_ks () 
  !
  ! ... pass all the k points to all the pool
  CALL bcast_wfc ( igk_k_all, ngk_all )
  !
  DEALLOCATE ( nbnd_occ )  ! otherwise allocate_ph complains: FIXME
  !
  ! 8)
  !CALL compute_map_ikq ()
  !
  ! Compute the peridoic part of the wannier orbital and store it on file 
  !
  WRITE(stdout,'(/)')
  WRITE( stdout, '(5X,"INFO: PREPARING THE KCWANN CALCULATION ...")')
  !
  ! write the list of q points on a file
  ! and store the q coordinates
  !
  iun_qlist = 127
  OPEN (iun_qlist, file = TRIM(tmp_dir)//TRIM(prefix)//'.qlist')
  !
  nqs = nkstot/nspin
  ALLOCATE (x_q (3, nqs) )
  ALLOCATE ( isq(nqs) )
  iq=1
  IF (ionode) THEN 
     WRITE(iun_qlist,'(i5)') nkstot/nspin 
     DO ik = 1, nkstot
       !WRITE(*,*) ik, isk(ik)
       IF (lsda .AND. isk(ik) /= spin_component) CYCLE
       WRITE(iun_qlist, '(3f12.8)') xk(:,ik)
       x_q(:,iq) = xk(:,ik)  
       isq(iq) = isk(ik) 
       iq = iq + 1
     ENDDO 
  ENDIF
  CALL mp_bcast (x_q, ionode_id, intra_image_comm)
  CALL mp_bcast (isq, ionode_id, intra_image_comm)
  !
  ALLOCATE ( lgamma_iq(nqs) )

  lgamma_iq(:) = .FALSE.
  !
  WRITE( stdout, '(/, 5X,"INFO: Compute Wannier-orbital Densities ...")')
  !
  DO iq = 1, nqs
    !! For each q in the mesh 
    !
    xq = x_q(:,iq)
    !
    ! IF (ionode) WRITE(iun_qlist,'(3f12.8)') xq
    !
    lgamma_iq(iq)=(x_q(1,iq)==0.D0.AND.x_q(2,iq)==0.D0.AND.x_q(3,iq)==0.D0)
    CALL cryst_to_cart(1, xq, at, -1)
    WRITE( stdout, '(/,8X, 78("="))')
    WRITE( stdout, '(  8X, "iq = ", i5)') iq
    WRITE( stdout, '(  8X, "The Wannier density at  q = ",3F12.7, "  [Cart ]")') x_q(:,iq)
    WRITE( stdout, '(  8X, "The Wannier density at  q = ",3F12.7, "  [Cryst]")') xq(:)
    WRITE( stdout, '(  8X, 78("="),/)')
    !
    CALL compute_map_ikq_single (iq)
    ! The map to identify which k point in the 1BZ corresponds to k+q and the G vector that produce the mapping
    ! The results are stored in the global variable map_ikq and shift_1bz (used inside rho_of_q) 
    ! can (should) be moved inside rho_of_q ( )
    !
    rhowann(:,:)=ZERO
    !! Initialize the periodic part of the wannier orbtal density at this q
    !
    CALL rho_of_q (rhowann, ngk_all, igk_k_all)
    ! Compute the peridic part rho_q(r) of the wannier density rho(r)
    ! rho(r)   = \sum_q exp[iqr]rho_q(r)
    !
    WRITE( stdout, '(8X,"INFO: rho_q(r) DONE ",/)')
    !
    !
    ! ... each q /= gamma is saved on a different directory
    lgamma = lgamma_iq(iq)
    !
    tmp_dir_phq= TRIM (tmp_dir_ph) // TRIM(prefix) // '_q' &
                  & // TRIM(int_to_char(iq))//'/'
    filename=TRIM(tmp_dir_phq)//TRIM(prefix)//'.save/charge-density.dat'
    IF (ionode) inquire (file =TRIM(filename), exist = exst)
    !
    CALL mp_bcast( exst, ionode_id, intra_image_comm )
    !
    CALL create_directory( tmp_dir_phq )
    tmp_dir=tmp_dir_phq
    CALL write_scf( rho, nspin )
    ! write the periodic part of the wannier orbital density on file
    DO i = 1, num_wann
      rhowann_aux (:) = rhowann(:,i) 
      file_base=TRIM(tmp_dir_phq)//TRIM(prefix)//'.save/rhowann_iwann_'//TRIM(int_to_char(i))
      CALL write_rhowann( file_base, rhowann_aux, dffts, ionode, inter_bgrp_comm )
    ENDDO
    tmp_dir=tmp_dir_save
    !
  ENDDO
  !
  WRITE( stdout, '(5X,"INFO: PREPARING THE KCWANN CALCULATION ... DONE")')
  WRITE(stdout,'(/)')
  !
  CALL close_buffer  ( iuwfc, 'KEEP' )
  !
  CALL stop_clock ('kc_setup')
  !
  DEALLOCATE (rhowann, rhowann_aux )
  ! 
  RETURN
END SUBROUTINE kc_setup
