!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
!-----------------------------------------------------------------------
subroutine kcw_setup_ham
  !-----------------------------------------------------------------------
  !
  !! As kcw_setup.f90 plus hamiltonian specific setups
  !
  !
  USE kinds,             ONLY : DP
  USE ions_base,         ONLY : nat, ityp
  USE io_files,          ONLY : tmp_dir
  USE scf,               ONLY : v, vrs, vltot,  kedtau
  USE fft_base,          ONLY : dfftp, dffts
  USE gvecs,             ONLY : doublegrid
  USE uspp_init,         ONLY : init_us_2
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux!, nspin_mag, npol
  USE wvfct,             ONLY : nbnd
  USE uspp,              ONLY : nkb, vkb
  !USE funct,             ONLY : dft_is_gradient
  USE xc_lib,            ONLY : xclib_dft_is
  !
  USE units_lr,          ONLY : iuwfc
  USE wvfct,             ONLY : npwx, current_k, npw
  USE control_flags,     ONLY : io_level
  USE io_files,          ONLY : prefix
  USE buffers,           ONLY : open_buffer, save_buffer, close_buffer
  USE control_kcw,       ONLY : alpha_final, evc0, iuwfc_wann, iurho_wann, kcw_iverbosity, &
                                read_unitary_matrix, hamlt, alpha_corr_done, &
                                num_wann, num_wann_occ, num_wann_emp, i_orb, iorb_start, iorb_end, &
                                calculation, nqstot, occ_mat ,alpha_final_full, spin_component, &
                                tmp_dir_kcw, tmp_dir_kcwq, x_q, lgamma_iq !, wq
  USE io_global,         ONLY : stdout
  USE klist,             ONLY : nkstot, xk, nks, ngk, igk_k
  USE cell_base,         ONLY : at !, bg
  USE fft_base,          ONLY : dffts
  !
  USE control_lr,        ONLY : nbnd_occ
  USE mp,                ONLY : mp_bcast
  USE eqv,               ONLY : dmuxc
  !
  USE io_kcw,            ONLY : read_rhowann, read_mlwf
  USE lsda_mod,          ONLY : lsda, isk, nspin, current_spin, starting_magnetization
  !
  USE coulomb,           ONLY : setup_coulomb
  !
  !
  !USE symm_base,       ONLY : s, t_rev, irt, nrot, nsym, invsym, nosym, &
  !                            d1,d2,d3, time_reversal, set_sym_bl, &
  !                            find_sym, inverse_s, copy_sym
  !USE cell_base,         ONLY : bg
  !
  implicit none

  integer :: na, i
  ! counters
  !
  INTEGER   :: lrwfc, lrrho, iun_qlist
  LOGICAL   :: exst, exst_mem
  INTEGER :: iq, nqs
  REAL(DP) :: xq(3)
  COMPLEX(DP), ALLOCATABLE :: rhowann(:,:), rhowann_aux(:)
  CHARACTER (LEN=256) :: file_base
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  !
  INTEGER :: ik, ik_eff
  CHARACTER(LEN=256)  :: dirname
  INTEGER, EXTERNAL :: global_kpoint_index
  LOGICAL :: mlwf_from_u = .FALSE.
  !
  !LOGICAL :: skip_equivalence
  !INTEGER :: nk1, nk2, nk3, k1, k2, k3
  !REAL(DP) :: wq_aux(100), xq_aux(3,100)
  !
  call start_clock ('kcw_setup')
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
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
     !if (dft_is_gradient()) call compute_ux(m_loc,ux,nat)
     IF ( xclib_dft_is('gradient') ) CALL compute_ux(m_loc,ux,nat)
  ENDIF
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  call setup_nbnd_occ ( ) 
  !
  ALLOCATE ( alpha_final(nbnd) )
  alpha_final(:) = 1.D0
  IF (nkstot/nspin == 1 ) THEN 
    ALLOCATE (alpha_final_full(nbnd))
    alpha_final_full(:) = 1.D0
  ENDIF
  !
  ! 6) Computes the derivative of the XC potential
  !
  IF (calculation == 'ham') THEN 
    allocate (dmuxc ( dfftp%nnr , nspin , nspin))
    dmuxc = 0.d0
    call setup_dmuxc()
    CALL kcw_R_points ()
    CALL read_alpha ()
    IF (nkstot/nspin == 1 ) alpha_final_full = alpha_final
  ! CALL setup_dgc()
  ENDIF
  !
  ! Open buffers for the KS and eventualy the minimizing WFCs 
  ! 
  iuwfc = 20
  lrwfc = nbnd * npwx
  io_level = 1
  CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
  IF (.NOT.exst.AND..NOT.exst_mem) THEN
     CALL errore ('kcw_setup_ham', 'file '//trim(prefix)//'.wfc not found', 1)
  END IF
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for KS wfcs, OPENED")')
  !
  ! 8) READ the U matrix from Wannier and set the toal number of WFs
  IF (read_unitary_matrix) THEN 
    CALL read_wannier ( ) ! This is not needed if the Wannier ware written/read from file.
                          ! For now needed to set-up the number of wannier. To be removed 
                          ! as soon as a proper data-file will be written by wann2kcw: FIXME
    if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Unitary matrix, READ from file")')
  ELSE
    num_wann = nbnd
    num_wann_occ = nbnd_occ(1)
    num_wann_emp = num_wann-num_wann_occ
  ENDIF
  ! 
  ! Open a file to store the KS states in the WANNIER gauge
  !
  iuwfc_wann = 21
  io_level = 1
  lrwfc = num_wann * npwx 
  CALL open_buffer ( iuwfc_wann, 'wfc_wann', lrwfc, io_level, exst )
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for WFs, OPENED")')
  !
  ! Open a buffer for the wannier orbital densities. Those have been written by wann2kcw
  ! and must be in the outdir. If not STOP
  !
  iurho_wann = 22
  io_level = 1
  lrrho=num_wann*dffts%nnr
  CALL open_buffer ( iurho_wann, 'rho_wann', lrrho, io_level, exst )
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for WF rho, OPENED")')
  !
  ALLOCATE ( rhowann ( dffts%nnr, num_wann), rhowann_aux(dffts%nnr) )
  ALLOCATE ( evc0(npwx, num_wann) )
  ALLOCATE ( hamlt(nkstot, num_wann, num_wann) )
  ALLOCATE ( alpha_corr_done (num_wann) ) 
  ALLOCATE ( occ_mat (num_wann, num_wann, nkstot) )
  occ_mat = 0.D0
  alpha_corr_done = .FALSE.
  hamlt(:,:,:) = ZERO
  !
  mlwf_from_u = .FALSE. ! Set this to true to revert to the previous behaviour (March 2021)
  IF (mlwf_from_u) THEN 
    !! ... Rotate the KS state to the localized gauge nd save on a buffer
    WRITE(stdout,'(/,5X, "INFO: MLWF from U matrix: &
                  Reading collected, re-writing distributed wavefunctions")')
    CALL rotate_ks () 
    !
  ELSE 
    dirname = TRIM (tmp_dir_kcw) 
    WRITE(stdout,'(/,5X, "INFO: MLWF read from file: &
                  Reading collected, re-writing distributed wavefunctions")')
    DO ik = 1, nks
        !
        current_k = ik
        IF ( lsda ) current_spin = isk(ik)
        IF ( lsda .AND. isk(ik) /= spin_component) CYCLE
        npw = ngk(ik)
        IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
        CALL read_mlwf ( dirname, ik, evc0 )
        ik_eff = ik-(spin_component-1)*nkstot/nspin
        CALL save_buffer ( evc0, lrwfc, iuwfc_wann, ik_eff )
        CALL ks_hamiltonian(evc0, ik, num_wann) 
    END DO
  ENDIF
  !
  DEALLOCATE ( nbnd_occ )  ! otherwise allocate_ph complains: FIXME
  !
  call setup_coulomb()
  !call setup_coulomb_exx()
  !
  WRITE( stdout, '(/, 5X,"INFO: READING Wannier-orbital Densities ...")') 
  !
  ! ... Read the q-point grid written by wann2kcw
  iun_qlist = 127
  OPEN (iun_qlist, file = TRIM(tmp_dir_kcw)//'qlist.txt')
  !
  READ(iun_qlist,'(i5)') nqs
  nqstot = nqs 
  !
  ALLOCATE (x_q (3, nqs) ) 
  ALLOCATE ( lgamma_iq(nqs) )
  lgamma_iq(:) = .FALSE.
  !
  DO iq = 1, nqs
    !! For each q in the mesh 
    !
    READ (iun_qlist,'(3f12.8)') xq
    x_q(:,iq) = xq ! Store the q vectors
    lgamma_iq(iq)=(x_q(1,iq)==0.D0.AND.x_q(2,iq)==0.D0.AND.x_q(3,iq)==0.D0)
    CALL cryst_to_cart(1, xq, at, -1)
    !
    WRITE( stdout,'(/,8X, 78("="))')
    WRITE( stdout, '( 8X, "iq = ", i5)') iq
    WRITE( stdout, '( 8X,"The  Wannier density at  q = ",3F12.7, "  [Cart ]")')  xk(:,iq)
    WRITE( stdout, '( 8X,"The  Wannier density at  q = ",3F12.7, "  [Cryst]")')  xq(:)
    WRITE( stdout, '( 8X, 78("="),/)')
    ! 
    tmp_dir_kcwq= TRIM (tmp_dir_kcw) // 'q' &
                & // TRIM(int_to_char(iq))//'/'
    !
    DO i = 1, num_wann
      file_base=TRIM(tmp_dir_kcwq)//'rhowann_iwann_'//TRIM(int_to_char(i))
      CALL read_rhowann( file_base, dffts, rhowann_aux )
      rhowann(:,i) = rhowann_aux(:)
    ENDDO
    !
    ! ... Save the rho_q on file
    !
    lrrho=num_wann*dffts%nnr
    CALL save_buffer (rhowann, lrrho, iurho_wann, iq)
    !
  ENDDO
  !
  ! ... Which orbital to compute
  !
  IF (i_orb == -1 ) THEN 
    !
    iorb_start = 1
    iorb_end = num_wann
    WRITE( stdout, '(5X,"INFO: total number of wannier to compute", i5, " from ", i5, " to ", i5)')  &
    iorb_end-iorb_start+1, iorb_start, iorb_end
    !
  ELSE IF  (i_orb .gt. num_wann ) THEN 
    !
    CALL errore('kcw_setup','i_orb > num_wann' ,ABS(num_wann))
    !
  ELSE
    !
    iorb_start = i_orb
    iorb_end = i_orb
    WRITE( stdout, '(5X,"INFO: total number of wannier to compute", i5, " from ", i5, " to ", i5)') & 
    1, iorb_start, iorb_end
    !
  ENDIF
  !
  WRITE( stdout, '(5X,"INFO: PREPARING THE KCW CALCULATION ... DONE")')
  WRITE(stdout,'(/)')
  ! 
  ! ##### ALL the following was ment to restore the symmetry and compute only non equivalent k/q points. 
  ! ##### But it does not seem to properly work! NEED TO UNDERSTAND BETTER THE SYMMETRY
  !skip_equivalence = .FALSE.
  !nk1 = 2; nk2 = 2; nk3 = 2
  !k1  = 0; k2  = 0; k3  = 0
  !CALL kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev, &
  !                  bg, nk1*nk2*nk3, k1,k2,k3, nk1,nk2,nk3, nqstot, xq_aux, wq_aux)
  !WRITE(stdout,*) "NICOLA", nqstot
  !DEALLOCATE(x_q) 
  !ALLOCATE ( x_q(3,nqstot) )
  !ALLOCATE( wq(nqstot) )
  !wq(1:nqstot) = wq_aux(1:nqstot)
  !x_q(:,1:nqstot) = xq_aux(:,1:nqstot)
  ! 
  !WRITE(stdout,'("NICOLA wq =", f12.8, 3x, 3f12.8 )') (wq(iq), x_q(:,iq), iq=1,nqstot)
  !STOP

  !CALL close_buffer  ( iuwfc, 'KEEP' )
  !
  CALL stop_clock ('kcw_setup')
  !
  DEALLOCATE (rhowann, rhowann_aux)
  !
  RETURN
  !
END SUBROUTINE kcw_setup_ham
