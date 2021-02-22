!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
!-----------------------------------------------------------------------
subroutine kc_setup_screen
  !-----------------------------------------------------------------------
  !
  !! As kc_setup.f90 plus screening specific setups
  !
  !
  USE kinds,             ONLY : DP
  USE ions_base,         ONLY : nat, ntyp => nsp, ityp
  USE io_files,          ONLY : tmp_dir
  USE lsda_mod,          ONLY : nspin, starting_magnetization
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
  USE buffers,           ONLY : open_buffer, save_buffer, close_buffer
  USE control_kc_wann,   ONLY : alpha_final, iurho_wann, kc_iverbosity, &
                                read_unitary_matrix, num_wann, num_wann_occ, i_orb, iorb_start, &
                                iorb_end, nqstot, occ_mat, l_do_alpha, group_alpha!, wq
  USE io_global,         ONLY : stdout
  USE klist,             ONLY : xk, nkstot
  USE cell_base,         ONLY : at !, bg
  USE fft_base,          ONLY : dffts
  USE disp,              ONLY : x_q,  done_iq, lgamma_iq
  !
  USE control_lr,        ONLY : nbnd_occ
  USE control_ph,        ONLY : tmp_dir_ph, tmp_dir_phq
  USE mp,                ONLY : mp_bcast
  USE io_kcwann,     ONLY : read_rhowann
  !
!  USE xml_io_base,     ONLY : write_rho
!  USE symm_base,       ONLY : s, t_rev, irt, nrot, nsym, invsym, nosym, &
!                              d1,d2,d3, time_reversal, sname, set_sym_bl, &
!                              find_sym, inverse_s, no_t_rev, copy_sym
!  USE cell_base,         ONLY : bg
!  USE control_kc_wann,  ONLY : wq, nqstot
!  USE symm_base,         ONLY : nosym, nrot
  !
  implicit none
  !
  integer :: na, i
  ! counters
  !
  INTEGER   :: lrwfc, lrrho, iun_qlist
  ! Lenght record sor wfc and density, iunit for the q point file (written by wann2kc.x)
  !
  LOGICAL   :: exst, exst_mem
  ! Check on the existence of the buffers
  !
  INTEGER :: iq, nqs
  ! Counters on the q points, total number of q points
  !
  REAL(DP) :: xq(3)
  ! the q-point coordinatew
  !
  COMPLEX(DP), ALLOCATABLE :: rhowann(:,:), rhowann_aux(:)
  ! the periodic part of the wannier orbital density
  !
  CHARACTER (LEN=256) :: file_base
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  LOGICAL :: iorb_from_input
  !
!  LOGICAL :: skip_equivalence
!  INTEGER :: nk1, nk2, nk3, k1, k2, k3
!  REAL(DP) :: wq_aux(100), xq_aux(3,100)
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
  !
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
  ALLOCATE ( alpha_final(nbnd) )
  alpha_final(:) = 1.D0
  !
  ! ... Open buffers for the KS and eventualy the minimizing WFCs 
  !
  iuwfc = 20
  lrwfc = nbnd * npwx
  io_level = 1
  CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
  IF (.NOT.exst.AND..NOT.exst_mem) THEN
     CALL errore ('kc_setup_screen', 'file '//trim(prefix)//'.wfc not found', 1)
  END IF
  if (kc_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for KS wfcs, OPENED")')
  !
  ! 8) READ the U matrix from Wannier and set the toal number of WFs
  !
  IF (read_unitary_matrix) THEN 
    CALL read_wannier ( )
    if (kc_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Unitary matrix, READ from file")')
  ELSE
    num_wann = nbnd
    num_wann_occ = nbnd_occ(1)
  ENDIF
  ! 
  ! ... Open a buffer for the wannier orbital densities. Those have been written by wann2kc
  ! ... and must be in the outdir. If not STOP
  !
  iurho_wann = 22
  io_level = 1
  lrrho=num_wann*dffts%nnr
  CALL open_buffer ( iurho_wann, 'rho_wann', lrrho, io_level, exst )
  if (kc_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Buffer for WF rho, OPENED")')
  !
  ALLOCATE (rhowann ( dffts%nnr, num_wann), rhowann_aux(dffts%nnr) )
  ALLOCATE ( occ_mat (num_wann, num_wann, nkstot) )
  ALLOCATE (l_do_alpha(num_wann), group_alpha(num_wann) ) 
  !
  l_do_alpha = .TRUE.
  DO i = 1, num_wann; group_alpha(i)=i; ENDDO
  occ_mat = 0.D0
  !
  ! ... Set up the coulomb kernel. If l_vcut=.true. the Gygi Balderschi scheme is used.
  ! ... Otherwise the g=0 component is set to zero 
  !
  call setup_coulomb_exx()
  !
  DEALLOCATE ( nbnd_occ )  ! otherwise allocate_ph complains: FIXME
  !
  ! ... Read the q-point grid written by wann2kc
  !
  WRITE( stdout, '(/, 5X,"INFO: READING Wannier-orbital Densities ...")') 
  !
  iun_qlist = 127
  OPEN (iun_qlist, file = TRIM(tmp_dir)//TRIM(prefix)//'.qlist')
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
    tmp_dir_phq= TRIM (tmp_dir_ph) // TRIM(prefix) // '_q' &
                & // TRIM(int_to_char(iq))//'/'
    !
    DO i = 1, num_wann
      file_base=TRIM(tmp_dir_phq)//TRIM(prefix)//'.save/rhowann_iwann_'//TRIM(int_to_char(i))
      CALL read_rhowann( file_base, dffts, rhowann_aux )
      rhowann(:,i) = rhowann_aux(:)
    ENDDO
    !
    ! ... Save the rho_q on a direct access file
    !
    lrrho=num_wann*dffts%nnr
    CALL save_buffer (rhowann, lrrho, iurho_wann, iq)
    !
  ENDDO
  !
  ! ... Which orbital to compute
  !
  iorb_from_input = .FALSE. 
  !
  IF (i_orb == -1 ) THEN 
    !
    iorb_start = 1
    iorb_end = num_wann
    !
  ELSE IF  (i_orb .gt. num_wann ) THEN 
    !
    CALL errore('kc_setup','i_orb > num_wann' ,ABS(num_wann))
    !
  ELSE
    !
    iorb_start = i_orb
    iorb_end = i_orb
    iorb_from_input = .TRUE. 
    !
  ENDIF
  !
  IF (iorb_from_input) WRITE( stdout, '(5X,"INFO: Orbital to compute specified from input", i5)') i_orb
  WRITE( stdout, '(5X,"INFO: total number of wannier to compute", i5, " from ", i5, " to ", i5)')  &
  iorb_end-iorb_start+1, iorb_start, iorb_end
  !
  ! ... Group the orbitals as a function of the Self-Hartree 
  !
  IF( .NOT. iorb_from_input) CALL group_orbitals ( )
  !
  WRITE( stdout, '(/,5X,"INFO: PREPARING THE KCWANN CALCULATION ... DONE")')
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

  !
  ALLOCATE (done_iq(nqs) )
  done_iq = .FALSE.
  !
  CALL close_buffer  ( iuwfc, 'KEEP' )
  !
  CALL stop_clock ('kc_setup')
  !
  DEALLOCATE (rhowann, rhowann_aux)
  !
  RETURN
  !
END SUBROUTINE kc_setup_screen
